require 'optparse'
require './clustering'
require './blast'
require './sequences'

# Argument validation

options = {}
opt_parser = OptionParser.new do |opt|
  opt.banner = "Usage: TYPE [SKIP_BLAST] FILE"
  opt.separator  ""
  opt.separator  "File: filename of the FASTA file containing the predicted sequences"
  opt.separator  ""
  opt.separator  "Options"

  opt.on("-t","--type TYPE","type of the predicted sequences: protein/mRNA") do |type|
    if type.to_s.downcase == 'protein' or type.to_s.downcase == 'mrna'
      options[:type] = type
    else 
      $stderr.print "Error: type may be protein or mRNA." + "\n"
      exit
    end
  end

  opt.on("-skip_blast","--skip_blast","skip blast-ing part and provide a blast xml output as input to this script") do
    options[:skip_blast] = true
  end

  opt.on("-h","--help","help") do
    puts opt_parser
  end
end

begin
  opt_parser.parse!(ARGV)

  unless options[:type]
    $stderr.puts "Error: you must specify --type option." + "\n"
    exit
  end 

  unless ARGV.length == 1
    $stderr.puts "Error: you must specify a single input file." + "\n"
    exit
  end 

  rescue OptionParser::ParseError
    $stderr.print "Error: " + $! + "\n"
    exit
end


# Main body

b = Blast.new(ARGV[0], options[:type].to_s.downcase)

unless options[:skip_blast]
  b.blast
  b.parse_output
  b.clusterization_by_length
  exit
end

# Skip the blast-ing part and provide a xml blast output file as argument to this ruby script
file = File.open(ARGV[0], "rb").read
b.parse_output(file)
b.plot_lengths
#b.clusterization_by_length
#b.plot_lengths
