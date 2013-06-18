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

unless options[:skip_blast]

  b = Blast.new(ARGV[0], options[:type].to_s.downcase)
  b.blast
  b.clusterization
  exit
end


# Skip the blast-ing part and provide a xml blast output file as argument to this ruby script

file = File.open(ARGV[0], "rb")

contents = file.read.scan(/<\bHit_len\b>(\d+)<\/\bHit_len\b>/)
contents = contents.map{ |x| x[0].to_i }.sort{|a,b| a<=>b}
 
clusters = hierarchical_clustering(contents)
max_density = 0;
max_density_cluster = 0;
clusters.each_with_index{|item, i|
        if item.density > max_density
                max_density = item.density
                max_density_cluster = i;
        end
}

file = File.open(ARGV[0], "rb")
query_len = file.read.scan(/<\bIteration_query-len\b>(\d+)<\/\bIteration_query-len\b>/)

puts "Predicted sequence length: #{query_len}"
puts "Maximum sequence length: #{contents.max}"
puts "Number of sequences: #{contents.length}"
puts "\nMost dense cluster:"
clusters[max_density_cluster].print_cluster



