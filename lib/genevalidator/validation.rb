require 'genevalidator/blast'
require 'genevalidator/output'
require 'genevalidator/exceptions'
require 'genevalidator/tabular_parser'
require 'genevalidator/validation_length_cluster'
require 'genevalidator/validation_length_rank'
require 'genevalidator/validation_blast_reading_frame'
require 'genevalidator/validation_gene_merge'
require 'genevalidator/validation_duplication'
require 'genevalidator/validation_open_reading_frame'
require 'genevalidator/validation_alignment'
require 'genevalidator/validation_codon_bias'
require 'bio-blastxmlparser'
require 'net/http'
require 'open-uri'
require 'uri'
require 'io/console'
require 'yaml'
require 'thread'

Pair1 = Struct.new(:x, :y)

class Validation

  attr_reader :type
  attr_reader :fasta_filepath
  attr_reader :html_path
  attr_reader :yaml_path
  attr_reader :mafft_path
  attr_reader :filename
  attr_reader :raw_seq_file
  attr_reader :raw_seq_file_index
  attr_reader :raw_seq_file_load
  # current number of the querry processed
  attr_accessor :idx
  attr_reader :start_idx
  # array of indexes for the start offsets of each query in the fasta file
  attr_reader :query_offset_lst
  # list with all validation reports
  attr_reader :all_query_outputs

  attr_reader :vlist
  attr_reader :tabular_format
  attr_reader :overall_evaluation
  attr_reader :multithreading

  # global variables
  attr_reader :no_queries
  attr_reader :scores
  attr_reader :good_predictions
  attr_reader :bad_predictions
  attr_reader :nee
  attr_reader :no_mafft
  attr_reader :no_internet
  attr_reader :map_errors
  attr_reader :map_running_times 
  attr_reader :global_codon_bias
  attr_reader :local_codon_bias_vector 

  attr_reader :wq
  attr_reader :threads
  attr_reader :mutex
  attr_reader :mutex_yaml
  attr_reader :mutex_html
  attr_reader :mutex_array

  ##
  # Initilizes the object
  # Params:
  # +fasta_filepath+: fasta file with query sequences
  # +vlist+: list of validations
  # +tabular_format+: list of column names for parsing the tablar blast output
  # +xml_file+: name of the precalculated blast xml output (used in 'skip blast' case)
  # +mafft_path+: path of the MAFFT program installation
  # +start_idx+: number of the sequence from the file to start with
  # +overall_evaluation+: boolean variable for printing / not printing overall evaluation
  def initialize( fasta_filepath,
                  vlist = ["all"],
                  tabular_format = nil,
                  xml_file = nil,
                  raw_seq_file = nil,
                  mafft_path = nil,
                  start_idx = 1,
                  overall_evaluation = true,
                  multithreading = false)

    puts "\nDepending on your input and your computational "<<
           "resources, this may take a while. Please wait..."


    # start a worker thread
    @threads = []
    @mutex = Mutex.new
    @mutex_yaml = Mutex.new
    @mutex_html = Mutex.new
    @mutex_array = Mutex.new

    @fasta_filepath = fasta_filepath
   
    @xml_file = xml_file
    @vlist = vlist.map{|v| v.gsub(/^\s/,"").gsub(/\s\Z/,"").split(/\s/)}.flatten

    @idx = 0
    @start_idx = 1

    # global variables
    @no_queries = 0
    @scores = []
    @good_predictions = 0
    @bad_predictions = 0
    @nee = 0
    @no_mafft = 0
    @no_internet = 0
    @map_errors = Hash.new(0)
    @map_running_times = Hash.new(Pair1.new(0,0))
    @global_codon_bias = AllCodonsBias.new 
    @local_codon_bias_vector = Hash.new({})

    raise FileNotFoundException.new unless File.exists?(@fasta_filepath)
    raise FileNotFoundException.new unless File.file?(@fasta_filepath)

    fasta_content = IO.binread(@fasta_filepath);

    # the expected type for the sequences is the
    # type of the first query
       
    # autodetect the type of the sequence in the FASTA
    # also check if the fasta file contains a single type of queries
    @type = BlastUtils.type_of_sequences(fasta_content)

    # create a list of index of the queries in the FASTA
    @query_offset_lst = fasta_content.enum_for(:scan, /(>[^>]+)/).map{ Regexp.last_match.begin(0)}
    raise FileNotFoundException.new unless @query_offset_lst != []
    @query_offset_lst.push(fasta_content.length)
    fasta_content = nil # free memory for variable fasta_content
    GC.start
    @tabular_format = tabular_format

    if mafft_path == nil
      @mafft_path = which("mafft")
    else
      @mafft_path = mafft_path
    end

    begin

      # fasta file is not a file
      raise FileNotFoundException.new unless File.file?(@fasta_filepath)

      # index raw_sequence file
      if raw_seq_file != nil
        raise FileNotFoundException.new unless File.exists?(raw_seq_file)
        @raw_seq_file = raw_seq_file

        # leave only the identifiers in the fasta description
        content = File.open(raw_seq_file, "rb").read.gsub(/ .*/, "")
        File.open(raw_seq_file, 'w+') { |file| file.write(content)}

        #index the fasta file
        keys = content.scan(/>(.*)\n/).flatten
        values = content.enum_for(:scan, /(>[^>]+)/).map{ Regexp.last_match.begin(0)}

        # make an index hash
        index_hash = Hash.new
        keys.each_with_index do |k, i|
          start = values[i]
          if i == values.length - 1
            endf = content.length - 1
          else
            endf = values[i+1]
          end
          index_hash[k] = [start, endf]
        end

        # create FASTA index
        @raw_seq_file_index = "#{raw_seq_file}.idx"
        @raw_seq_file_load = index_hash

        File.open(@raw_seq_file_index, "w") do |f|
          YAML.dump(index_hash, f)
        end

        content = nil

        @multithreading = multithreading
        @overall_evaluation = overall_evaluation

      end
      rescue Exception => error
        $stderr.print "Error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
          "Possible cause: your file with raw sequences is not FASTA. Please use get_raw_sequences executable to create a correct one.\n"
    end
   
    # build the path of html folder output
    path = File.dirname(@fasta_filepath)
    @html_path = "#{fasta_filepath}.html"
    @yaml_path = path

    @filename = File.basename(@fasta_filepath)#.scan(/\/([^\/]+)$/)[0][0]
    @all_query_outputs = []

    # create 'html' directory
    if File.exists? @html_path
      $stderr.print "The output directory already exists for this fasta file. For a new validation please remove the following directory: #{@html_path}\n"
      exit
    end
    Dir.mkdir(@html_path)

    # copy auxiliar folders to the html folder
    FileUtils.cp_r(File.join(File.dirname(File.expand_path(__FILE__)), "../../aux/css"), @html_path)
    FileUtils.cp_r(File.join(File.dirname(File.expand_path(__FILE__)), "../../aux/js"), @html_path)
    FileUtils.cp_r(File.join(File.dirname(File.expand_path(__FILE__)), "../../aux/img"), @html_path)
    FileUtils.cp_r(File.join(File.dirname(File.expand_path(__FILE__)), "../../aux/font"), @html_path)
    FileUtils.cp_r(File.join(File.dirname(File.expand_path(__FILE__)), "../../aux/doc"), @html_path)

  rescue SequenceTypeError => error
    $stderr.print "Sequence Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
     "Possible cause: input file containes mixed sequence types.\n"
    exit
  rescue FileNotFoundException => error
    $stderr.print "File not found error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}."<<
     "Possible cause: input file does not exist.\n"
    exit
  end

  ##
  # Parse the blast output and run validations
  def validation
      if @xml_file == nil
        #file seek for each query
        @query_offset_lst[0..@query_offset_lst.length-2].each_with_index do |pos, i|
          if (i+1) >= @start_idx
            query = IO.binread(@fasta_filepath, @query_offset_lst[i+1] - @query_offset_lst[i], @query_offset_lst[i]);

            #call blast with the default parameters
            if type == :protein
              output = BlastUtils.call_blast_from_stdin("blastp", query, 11, 1)
            else
              output = BlastUtils.call_blast_from_stdin("blastx", query, 11, 1)
            end

            #parse output
            parse_output(output, :stream)
          else
            @idx = @idx + 1
          end
        end
      else

        #file = File.open(@xml_file, "rb").read
        #check the format of the input file
        parse_output(@xml_file)
      end

      # evaluate codon bias
      #@global_codon_bias = AllCodonsBias.new
      if @local_codon_bias_vector != []
        @local_codon_bias_vector.each do |key,codon_bias|
          @global_codon_bias.update_codon_bias(codon_bias)
        end

        @global_codon_bias.each{|codon| codon.print_bias}
      end

      if @overall_evaluation
        #Output.print_footer(@all_query_outputs, @html_path, @filename)
        Output.print_footer(@no_queries, @scores, @good_predictions, @bad_predictions, @nee, @no_mafft, @no_internet, @map_errors, @map_running_times, @html_path, @filename)
        #@map_running_times
      end

    rescue SystemCallError => error
      $stderr.print "Load error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: input file is not valid\n"
      exit
    rescue SequenceTypeError => error
      $stderr.print "Sequence Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: the blast output was not obtained against a protein database.\n"
      exit!
    rescue Exception => error
       $stderr.print "Error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}.\n"
       exit!
  end

  ##
  # Params:
  # +output+: filename or stream, according to the type
  # +type+: file or stream
  def parse_output(output, type=:file)

    # check the format of the blat input (xml or tabular)
    if type == :file
      begin
        iterator_xml = Bio::BlastXMLParser::XmlIterator.new(output).to_enum

        # the following instruction will raise exception if the file is not xml
        hits = BlastUtils.parse_next_query_xml(iterator_xml, @type)

        input_file_type = :xml
        iterator_xml = Bio::BlastXMLParser::XmlIterator.new(output).to_enum

      rescue SequenceTypeError => error
        $stderr.print "Sequence Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
          "Possible cause: the blast output was not obtained against a protein database.\n"
        exit!
      rescue Exception => error
        begin
          query = IO.binread(@fasta_filepath, @query_offset_lst[idx+1] - @query_offset_lst[idx], @query_offset_lst[idx])
          parse_query = query.scan(/>([^\n]*)\n([A-Za-z\n]*)/)[0]
          definition = parse_query[0].gsub("\n","")
          identifier = definition.gsub(/ .*/,"")

          iterator_tab = TabularParser.new(@xml_file, tabular_format, @type)

          # the following instruction will raise exception if the file is not tabular
          hits = iterator_tab.next(identifier)
 
          input_file_type = :tabular
          iterator_tab = TabularParser.new(@xml_file, tabular_format, @type)

          if @tabular_format == nil
            puts "Note: Please specify the --tabular argument if you used tabular format input with nonstandard columns.\n"
          end

        rescue SequenceTypeError => error
          $stderr.print "Sequence Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
            "Possible cause: the blast output was not obtained against a protein database.\n"
          exit!
        rescue Exception => error
          $stderr.print "Blast file error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
             "Possible cause: blast output file format is neihter xml nor tabular.\n"
          exit!
        end
      end
    else
      input_file_type = :xml
      iterator_xml = Bio::BlastXMLParser::NokogiriBlastXml.new(output).to_enum
    end

    begin

      if @idx+1 == @query_offset_lst.length
        break
      end

      # get info about the query
      # get the @idx-th sequence from the fasta file
      prediction = Sequence.new

      query = IO.binread(@fasta_filepath, @query_offset_lst[idx+1] - @query_offset_lst[idx], @query_offset_lst[idx])
      parse_query = query.scan(/>([^\n]*)\n([A-Za-z\n]*)/)[0]

      prediction.definition = parse_query[0].gsub("\n","")
      prediction.identifier = prediction.definition.gsub(/ .*/,"")
      prediction.type = @type
      prediction.raw_sequence = parse_query[1].gsub("\n","")
      prediction.length_protein = prediction.raw_sequence.length

      if @type == :nucleotide
        prediction.length_protein /= 3
      end

      @idx = @idx + 1
 
      if input_file_type == :tabular
        if @idx < @start_idx
          iterator_tab.jump_next
        else
          hits = iterator_tab.next(prediction.identifier)
        end
      else
        if input_file_type == :xml
          if @idx < @start_idx
            iter = iterator_xml.next
          else
            hits = BlastUtils.parse_next_query_xml(iterator_xml, @type)
          end
        end
      end

      if hits == nil
        @idx = @idx -1
        break
      end

      # the first validation should be treated separately
      if @idx == @start_idx
        validate(prediction, hits, idx)
      else
        if @multithreading
          @threads << Thread.new(prediction, hits, @idx){ |prediction_local, hits_local, idx_local|
            validate(prediction_local, hits_local, idx_local)
            GC.start
          }
        else
          validate(prediction, hits, idx)
        end
      end

      hits = nil # free memory
      GC.start # force garbage collector
      #puts `ps -o rss= -p #{Process.pid}`.to_i
    end while 1

    @threads.each {|t| t.join}

  end

  ##
  # Validate one query and create validation report
  # Params:
  # +prediction+: Sequence object
  # +hits+: Array of +Sequence+ objects
  # +idx+: the index number of the query
  def validate(prediction, hits, idx)

    query_output = do_validations(prediction, hits, idx)
    query_output.generate_html
    query_output.print_output_file_yaml
    query_output.print_output_console

    validations = query_output.validations

    no_mafft = 0
    no_internet = 0
    errors = []
    validations.each do |v| 
      if v.errors != nil
        no_mafft += v.errors.select{|e| e == NoMafftInstallationError}.length
        no_internet += v.errors.select{|e| e == NoInternetError}.length
      end
      if v.validation == :error
        errors.push(v.short_header)
      end 
    end

    no_evidence = validations.count{|v| v.result == :unapplicable or v.result == :warning} == validations.length 
    if no_evidence 
      nee = 1
    else
      nee = 0
    end

    #report = ValidationReportStat.new(query_output.overall_score, no_evidence, no_mafft, no_internet, errors) 

    good_predictions = 0
    bad_predictions = 0
    if query_output.overall_score >= 75
      good_predictions = 1
    else
      bad_predictions = 1
    end

    #@map_running_times = Hash.new(0)

    @mutex_array.synchronize {
#      @all_query_outputs.push(query_output)
      #@all_query_outputs.push(report)

      @no_queries += 1
      @scores.push(query_output.overall_score)
      @good_predictions += good_predictions 
      @bad_predictions += bad_predictions 
      @nee += nee
      @no_mafft += no_mafft
      @no_internet += no_internet
      errors.each{|err| @map_errors[err] += 1} 
      
      validations.each do |v|
        if v.running_time != 0 and 
          v.running_time != nil and 
          v.validation != :unapplicable and
          v.validation != :error
            p = Pair1.new(@map_running_times[v.short_header].x + v.running_time, @map_running_times[v.short_header].y + 1)
            @map_running_times[v.short_header] = p
        end
      end

      codon_bias_validation = validations.select {|v| v.class == CodonBiasValidationOutput }[0]

      if codon_bias_validation != []
        if codon_bias_validation.validation_result != :unapplicable
          codons = codon_bias_validation.codon_bias
          #codons.each{|codon| codon.print_bias}
          @local_codon_bias_vector[idx] = codons
        end
      end

    }
  
    return query_output

  end

  ##
  # Removes identical hits
  # Params:
  # +prediction+: Sequence object
  # +hits+: Array of +Sequence+ objects
  # Output:
  # new array of hit +Sequence+ objects
  def remove_identical_hits(prediction, hits)
    # remove the identical hits
    # identical hit means 100%coverage and >99% identity
    identical_hits = []
    hits.each do |hit|
      # check if all hsps have identity more than 99%
      low_identity = hit.hsp_list.select{|hsp| hsp.pidentity == nil or hsp.pidentity < 99}

      # check the coverage
      coverage = Array.new(prediction.length_protein,0)
      hit.hsp_list.each do |hsp|
         len = hsp.match_query_to - hsp.match_query_from + 1
         coverage[hsp.match_query_from-1..hsp.match_query_to-1] = Array.new(len, 1)
      end

      if low_identity.length == 0 and coverage.uniq.length == 1
        identical_hits.push(hit)
      end
    end

    identical_hits.each {|hit| hits.delete(hit)}
    return hits
  end
  
  ##
  # Runs all the validations and prints the outputs given the current
  # prediction query and the corresponding hits
  # Params:
  # +prediction+: Sequence object
  # +hits+: Array of +Sequence+ objects
  # +idx+: the index number of the query
  # Output:
  # +Output+ object
  def do_validations(prediction, hits, idx)

    begin
      hits = remove_identical_hits(prediction, hits)
      rescue Exception => error #NoPIdentError
    end

    # do validations
    query_output = Output.new(@mutex, @mutex_yaml, @mutex_html, @filename, @html_path, @yaml_path, idx, @start_idx)
    query_output.prediction_len = prediction.length_protein
    query_output.prediction_def = prediction.definition
    query_output.nr_hits = hits.length

    plot_path = "#{html_path}/#{filename}_#{idx}"

    validations = []
    validations.push LengthClusterValidation.new(@type, prediction, hits, plot_path)
    validations.push LengthRankValidation.new(@type, prediction, hits)
    validations.push BlastReadingFrameValidation.new(@type, prediction, hits)
    validations.push GeneMergeValidation.new(@type, prediction, hits, plot_path)
    validations.push DuplicationValidation.new(@type, prediction, hits, @mafft_path, @raw_seq_file, @raw_seq_file_index, @raw_seq_file_load)
    validations.push OpenReadingFrameValidation.new(@type, prediction, hits, plot_path, ["ATG"])
    validations.push AlignmentValidation.new(@type, prediction, hits, plot_path, @mafft_path, @raw_seq_file, @raw_seq_file_index, @raw_seq_file_load)
    validations.push CodonBiasValidation.new(@type, prediction, hits)

    # check the class type of the elements in the list
    validations.each do |v|
      raise ValidationClassError unless v.is_a? ValidationTest
    end

    # check alias duplication
    aliases = validations.map(&:cli_name)
    unless aliases.length == aliases.uniq.length
      raise AliasDuplicationError
    end

    if vlist.map{|v| v.strip.downcase}.include? "all"

      validations.map{|v| v.run}
      # check the class type of the validation reports

      validations.each do |v|
        raise ReportClassError unless v.validation_report.is_a? ValidationReport
      end

      query_output.validations = validations.map{|v| v.validation_report}

    else

      desired_validations = validations.select {|v| vlist.map{|vv| vv.strip.downcase}.include? v.cli_name.downcase }
      desired_validations.each do |v|
        v.run
        raise ReportClassError unless v.validation_report.is_a? ValidationReport
      end
      query_output.validations = desired_validations.map{|v| v.validation_report}
 
      if query_output.validations.length == 0
        raise NoValidationError
      end

    end
   
    # compute validation score
    validations = query_output.validations
    successes = validations.map{|v| v.result ==
      v.expected}.count(true)

    fails = validations.map{|v| v.validation != :unapplicable and
      v.validation != :error and
      v.result != v.expected}.count(true)

    lcv = validations.select{|v| v.class == LengthClusterValidationOutput}
    lrv = validations.select{|v| v.class == LengthRankValidationOutput}
    if lcv.length == 1 and lrv.length == 1
      score_lcv = (lcv[0].result == lcv[0].expected)
      score_lrv = (lrv[0].result == lrv[0].expected)
      # if both are true this should be counted as a single success
      if score_lcv == true and score_lrv == true
        successes = successes - 1
      else
      # if both are false this will be a fail
        if score_lcv == false and score_lrv == false
          fails = fails - 1
        else
          successes = successes - 0.5
          fails = fails - 0.5
        end
      end
    end

    query_output.successes = successes
    query_output.fails = fails
    query_output.overall_score = (successes*100/(successes + fails + 0.0)).round(0)

=begin
counts = Hash.new{ 0 }
ObjectSpace.each_object do |o|
  counts[o.class] += 1
end
puts counts.to_s
=end
    return query_output

  rescue ValidationClassError => error
     $stderr.print "Class Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: type of one of the validations is not ValidationTest\n"
    exit!
  rescue NoValidationError => error
     $stderr.print "Validation error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: your -v arguments are not valid aliases\n"
     exit!
  rescue ReportClassError => error
      $stderr.print "Class Type error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: type of one of the validation reports returned by the 'run' method is not ValidationReport\n"
      exit!
  rescue AliasDuplicationError => error
      $stderr.print "Alias Duplication error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}. "<<
        "Possible cause: At least two validations have the same CLI alias\n"
      exit!
  rescue Exception => error 
      puts error.backtrace
      $stderr.print "Error at #{error.backtrace[0].scan(/\/([^\/]+:\d+):.*/)[0][0]}.\n"
      exit!
  end

  ##
  # The ruby equivalent for 'which' command in unix
  def which(cmd)
    exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
    ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
      exts.each { |ext|
        exe = File.join(path, "#{cmd}#{ext}")
        return exe if File.executable? exe
      }
    end
    return nil
  end

end

