require 'genevalidator/validation_report'
require 'genevalidator/codon_bias'
require 'genevalidator/validation_open_reading_frame'
require 'genevalidator/exceptions'

##
# Class that stores the validation output information
class CodonBiasValidationOutput < ValidationReport

  attr_reader :codon_bias

  def initialize (codon_bias, expected = :yes)

    @short_header = "CodonBias"
    @header = "Codon Bias"
    @description = "Check if the codon bias of the coding parts of the prediction"<<
    " fits the codon bias of the whole genome dataset."

    @codon_bias = codon_bias
    @result = validation
    @expected = expected
  end

  def print
    "pass"
  end

  def validation
    :yes
=begin
    if gaps < @threshold and extra_seq < @threshold and (1-consensus) < @threshold
      :yes
    else
      :no
    end
=end
  end

end

##
# This class contains the methods necessary for
# validations based on multiple alignment
class CodonBiasValidation < ValidationTest

  ##
  # Initilizes the object
  # Params:  
  # +type+: type of the predicted sequence (:nucleotide or :protein)
  # +prediction+: a +Sequence+ object representing the blast query
  # +hits+: a vector of +Sequence+ objects (usually representig the blast hits)
  def initialize(type, prediction, hits)
    super

    @short_header = "CodonBias"
    @header = "Codon Bias"
    @description = "Check if the codon bias of the coding parts of the prediction"<<
    " fits the codon bias of the whole genome dataset."
    @cli_name = "codon"

  end

  ##
  # Find gaps/extra regions based on the multiple alignment
  # of the first n hits
  # Output:
  # +AlignmentValidationOutput+ object
  def run

    begin
      raise Exception unless prediction.is_a? Sequence

      if type.to_s != "nucleotide"
        @validation_report = ValidationReport.new("", :unapplicable)
        return @validation_report
      end

    start = Time.now

    orfs = OpenReadingFrameValidation.get_orfs(100, prediction, ["ATG"], [])
    #get the reading frame
    if prediction.nucleotide_rf == nil
      prediction_len = prediction.raw_sequence.length
      max_ratio = orfs.map{|key, value| value.map{|orf| (orf[1]-orf[0])/(prediction_len + 0.0)}.max}.select{|e| e != nil}.max

      ratios = {}
      orfs.each do |key, value| 
        ratios[value.map{|orf| (orf[1]-orf[0])/(prediction_len + 0.0)}.max] = key
      end

      prediction.nucleotide_rf = ratios[max_ratio]

=begin      
      reading_frames = []
      ratios.each_with_index do |ratio, idx|
         if ratio != nil
           if ratio > 0.9
             reading_frames.push(idx)
           end
         end
      end

      puts reading_frames     

      if(reading_frames.length != 1)
        prediction.nucleotide_rf = 0 # muliple reading frames
      else
        prediction.nucleotide_rf = reading_frames[0]
      end
=end
    end  

    if prediction.nucleotide_rf == nil
        @validation_report = ValidationReport.new("", :unapplicable)
        return @validation_report
    end

    codons = AllCodonsBias.new
    map_codon_aa = codons.get_map_codon_aa

    # get the sequence corresponding to the open readin frames
    seq_orfs = orfs[prediction.nucleotide_rf]

    seq = ""
    seq_orfs.each do |orf|
      seq<<prediction.raw_sequence[orf[0]..orf[1]].upcase
    end

    if prediction.nucleotide_rf < 0
      seq = seq.reverse
    end

    seq_codons = seq.scan(/.{3}/)
    seq_codons.each do |codon|
      if map_codon_aa[codon] != nil
        map_codon_aa[codon].add_amino_acid(codon)
      end
    end

    codons = codons.select{|codon| codon.total_count!=0}

    @validation_report = CodonBiasValidationOutput.new(codons)
    @validation_report.running_time = Time.now - start

    return @validation_report

    # Exception is raised when blast founds no hits
    rescue  NotEnoughHitsError => error
      @validation_report = ValidationReport.new("Not enough evidence", :warning, @short_header, @header, @description)
      return @validation_report
    rescue ReadingFrameError => error
      puts error.backtrace
      @validation_report = ValidationReport.new("Multiple reading frames", :error, @short_header, @header, @description)
      return @validation_report
    rescue Exception => error
      puts error.backtrace
      @validation_report.errors.push "Unexpected Error"
      @validation_report = ValidationReport.new("Unexpected error", :error, @short_header, @header, @description)
      @validation_report.errors.push OtherError
      return @validation_report
    end
  end

end
