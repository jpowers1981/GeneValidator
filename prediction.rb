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
else
  # Skip the blast-ing part and provide a xml blast output file as argument to this ruby script
  file = File.open(ARGV[0], "rb").read
  b.parse_output(file)
end


printf "\nQuery No | Definition | Predicted Length | Confidence Interval | Silhouette Score | Accepted (yes/no)\n"

idx = 0
begin
  idx = idx + 1
  seq = b.parse_next_query #return [hits, predicted_seq]

  if seq == nil
    break
  end

  hits = seq[0]
  predicted_seq = seq[1]

  rez = b.clusterization_by_length(hits, predicted_seq, false) #return [clusters, max_density_cluster_idx]
  clusters = rez[0]
  max_cluster_idx = rez[1]
  b.plot_histo_clusters(clusters, predicted_seq.xml_length, max_cluster_idx, "#{ARGV[0]}_#{idx}")
  b.plot_length(hits, predicted_seq, "#{ARGV[0]}_#{idx}")
  silhouette = b.sequence_silhouette(predicted_seq, max_cluster_idx, clusters)

  limits = clusters[max_cluster_idx].get_limits
  max_len = hits.map{|x| x.xml_length}.max
  predicted_len = predicted_seq.xml_length
  pval = clusters[0].t_test(clusters, predicted_len)

  if predicted_len <= limits[1] and predicted_len >= limits[0]
    printf "Query %3d: %-20s %6d [%6d:%6d]    %+.2f  yes  p-value = %f \n", idx, predicted_seq.definition[0, [predicted_seq.definition.length-1,20].min],
         predicted_len, limits[0], limits[1], silhouette, pval
  else
    printf "Query %3d: %-20s %6d [%6d:%6d]    %+.2f  no   p-value = %f\n", idx, predicted_seq.definition[0, [predicted_seq.definition.length-1,20].min],
         predicted_len, limits[0], limits[1], silhouette, pval


  end
end while 1

