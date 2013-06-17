#!/usr/bin/env ruby

require './clustering'
require 'bio-blastxmlparser'
require 'net/http'
require 'open-uri'
require 'uri'

class Blast

  # blast command: blastn or blastp
  attr_reader :command

  # result of executing command
  attr_reader :result

  #BLAST output
  attr_reader :xml_output

  def initialize
  end

  def blast(command, filename, gapopen, gapextend)
    # we don't know what to do if the arguments ain't String
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
    output = %x[#{cmd} 2>/dev/null]

    if output == ""
      $stderr.print "BLAST error. Did you add BLAST path to CLASSPATH?\n"
      exit
    end

    @xml_output = output
    output
  end

  def parse_output
    output = @xml_output

  end

  def get_sequence_by_accession_no(accno,db)

    #url = URI.parse('http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&retmax=1&usehistory=y&term=EF100000')

    uri = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&retmax=1&usehistory=y&term=#{accno}/"

    #result = open(uri).read
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

  def clusterization
    output = @xml_output
    output2 = @xml_output

    contents = output.scan(/<\bHit_len\b>(\d+)<\/\bHit_len\b>/)
    contents = contents.map{ |x| x[0].to_i }.sort{|a,b| a<=>b}

    query_len = output2.scan(/<\bIteration_query-len\b>(\d+)<\/\bIteration_query-len\b>/)

    clusters = hierarchical_clustering(contents)
    max_density = 0;
    max_density_cluster = 0;
    clusters.each_with_index do |item, i|
      if item.density > max_density
        max_density = item.density
        max_density_cluster = i;
      end
    end
      
    puts "Predicted sequence length: #{query_len}"
    puts "Maximum sequence length: #{contents.max}"
    puts "Number of sequences: #{contents.length}"
    puts "\nMost dense cluster:"

    clusters[max_density_cluster].print_cluster
  end
end

b = Blast.new
out = b.get_sequence_by_accession_no("EF100000","nucleotide")
puts out
