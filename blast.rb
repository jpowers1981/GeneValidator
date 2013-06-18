#!/usr/bin/env ruby

require './clustering'
require './sequences'
require 'bio-blastxmlparser'
require 'net/http'
require 'open-uri'
require 'uri'

class Blast

  # blast command: blastn or blastp
  attr_reader :command
  # result of executing command
  attr_reader :result
  #blast output
  attr_reader :xml_output
  #query sequence type
  attr_reader :type
  #query sequence fasta file
  attr_reader :fasta_file
  #predicted sequence
  attr_reader :predicted_seq
  #array list with the reference sequences
  attr_reader :ref_seq_list

  def initialize(fasta_file, type)
    @type = type
    @fasta_file = fasta_file
  end

  #call blast according to the type of the sequence
  def blast
    #call blast with the default parameters
    if type == 'protein'
      advanced_blast("blastp", @fasta_file, 11, 1)
    else
      advanced_blast("blastn", @fasta_file, 5, 2)
    end
  end

  #call blast with specific parameters
  def advanced_blast(command, filename, gapopen, gapextend)
    
    raise TypeError unless command.is_a? String and filename.is_a? String

    evalue = "1e-5"

    #blast output format:
    #0 = pairwise,
    #1 = query-anchored showing identities,
    #2 = query-anchored no identities,
    #3 = flat query-anchored, show identities,
    #4 = flat query-anchored, no identities,
    #5 = XML Blast output,
    #6 = tabular,
    #7 = tabular with comment lines,
    #8 = Text ASN.1,
    #9 = Binary ASN.1,
    #10 = Comma-separated values,
    #11 = BLAST archive format

    cmd = "#{command} -query #{filename} -db nr -remote -evalue #{evalue} -outfmt 5 -gapopen #{gapopen} -gapextend #{gapextend} "
    puts "Executing \"#{cmd}\"..."
    puts "This may take a while..."
    output = %x[#{cmd} 2>/dev/null]

    if output == ""
      $stderr.print "BLAST error. Did you add BLAST path to CLASSPATH?\n"
      exit
    end

    @xml_output = output
    output
  end

  #parse the xml blast output given as string parameter (optional parameter)
  def parse_output(output = nil)

    @predicted_seq = Sequence.new

    hits = Array.new

    if output == nil
      output = @xml_output
    end

    xml = Bio::BlastXMLParser::NokogiriBlastXml.new(output).to_enum #XmlIterator.new(filename).to_enum

    iter = xml.next
    @predicted_seq.xml_length = iter.field("Iteration_query-len")

    iter.each do | hit |
      hsp = hit.hsps.first
      hsp.field("Hsp_bit-score")
      seq = Sequence.new

      seq.object_type = "ref"
      seq.seq_type = @type
      seq.database = iter.field("BlastOutput_db")
      seq.id = hit.hit_id
      seq.definition = hit.hit_def
      seq.species = hit.hit_def.scan(/\[([^\]\[]+)\]$/)[0][0]
      seq.accession_no = hit.accession
      seq.e_value = hsp.evalue
      seq.xml_length = hit.len
      seq.hit_from = hsp.hit_from
      seq.hit_to = hsp.hit_to

      #get gene by accession number
      if @type == "protein"
        seq.raw_sequence = get_sequence_by_accession_no(seq.accession_no, "protein")
      else
        seq.raw_sequence = get_sequence_by_accession_no(seq.accession_no, "nucleotide")
      end
      seq.fasta_length = seq.raw_sequence.length

      align = Alignment.new
      align.query_seq = hsp.qseq
      align.hit_seq = hsp.hseq
      align.bit_score = hsp.bit_score
      align.score = hsp.score

      regex = align.hit_seq.gsub(/[+ -]/, '+' => '.', ' ' => '.', '-' => '')
      #puts "----\n#{regex}\n----"

      #seq.alignment_start_offset = seq.raw_sequence.index(/#{regex}/)
      seq.alignment = align

      hits.push(seq)
      #seq.print
      puts "getting sequence #{seq.accession_no}..."
    end     
    
    @ref_seq_list = hits	
    hits
  end

  #get gene by accession number from a givem database
  def get_sequence_by_accession_no(accno,db)

    uri = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=#{db}&retmax=1&usehistory=y&term=#{accno}/"
    result = Net::HTTP.get(URI.parse(uri))

    result2 = result
    queryKey = result2.scan(/<\bQueryKey\b>([\w\W\d]+)<\/\bQueryKey\b>/)[0][0]
    webEnv = result.scan(/<\bWebEnv\b>([\w\W\d]+)<\/\bWebEnv\b>/)[0][0]

    uri = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?rettype=fasta&retmode=text&retstart=0&retmax=1&db=#{db}&query_key=#{queryKey}&WebEnv=#{webEnv}"

    result = Net::HTTP.get(URI.parse(uri))

    #parse FASTA output
    rec=result
    nl = rec.index("\n")
    header = rec[0..nl-1]
    seq = rec[nl+1..-1]
    seq.gsub!(/\n/,'')

    seq

  end

  #clusterization by length from a list of sequences
  def clusterization_by_length(lst = nil)

    if lst == nil
      lst = @ref_seq_list
    end

    raise TypeError unless lst[0].is_a? Sequence

    contents = lst.map{ |x| x.xml_length.to_i }.sort{|a,b| a<=>b}

    clusters = hierarchical_clustering(contents)
    max_density = 0;
    max_density_cluster = 0;
    clusters.each_with_index do |item, i|
      if item.density > max_density
        max_density = item.density
        max_density_cluster = i;
      end
    end
      
    puts "Predicted sequence length: #{@predicted_seq.xml_length}"
    puts "Maximum sequence length: #{contents.max}"
    puts "Number of sequences: #{contents.length}"
    puts "\nMost dense cluster:"

    clusters[max_density_cluster].print_cluster

  end
end

#Main body
#Test certain methods of Blast class
=begin
b = Blast.new("ana","protein")
puts b.get_sequence_by_accession_no("EF100000","nucleotide")
file = File.open("/home/monique/GSoC2013/data/output_prot1_predicted.xml", "rb").read
b.parse_output(file)
=end

