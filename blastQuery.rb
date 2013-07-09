#/usr/bin/env ruby

require './clustering'
require './sequences'
require './blastQuery'
require 'rinruby.rb'

class QueryError < Exception
end

##
# This class stores all the data obtained from the blast query

class BlastQuery

  attr_reader :filename
  attr_reader :query_index
  attr_reader :hits
  attr_reader :prediction
  attr_reader :clusters
  attr_reader :max_density_cluster
  attr_reader :mean

##
# Initilizes the object
# Params:
# +hits+: a vector of +Sequence+ objects (usually representig the blast hits)
# +prediction+: a +Sequence+ object representing the blast query
# +filename+: name of the input file, used when generatig the plot files
# query_index: the number of the query in the blast output

  def initialize(hits, prediction, filename, query_index)
    begin
      raise QueryError unless hits[0].is_a? Sequence and prediction.is_a? Sequence and filename.is_a? String and query_index.is_a? Fixnum
      @hits = hits
      @prediction = prediction
      @filename = filename
      @query_index = query_index

    end
  end

##
# Calculates a precentage based on the rank of the predicion among the hit lengths
# Params:
# +hits+ (optional): a vector of +Sequence+ objects
# +prediction+ (optional): a +Sequence+ object
  def length_rank(debug = false, hits = @hits, prediction = @prediction)
    begin
      raise TypeError unless hits[0].is_a? Sequence and prediction.is_a? Sequence

      lengths = hits.map{ |x| x.xml_length.to_i }.sort{|a,b| a<=>b}
      len = lengths.length
      median = len % 2 == 1 ? lengths[len/2] : (lengths[len/2 - 1] + lengths[len/2]).to_f / 2

      predicted_len = prediction.xml_length.to_i
      if predicted_len < median
        rank = lengths.find_all{|x| x < predicted_len}.length
      else
        rank = lengths.find_all{|x| x > predicted_len}.length
      end

      percentage = rank / (len + 0.0)
      percentage.round(2)

    rescue TypeError
      $stderr.print "Type error. Possible cause: one of the arguments of 'length_rank' method has not the proper type.\n"
      exit
    end
  end

  def length_validation

      ret = clusterization_by_length  #returns [clusters, max_density_cluster_idx]

      @clusters = ret[0]
      @max_density_cluster = ret[1]
      @mean = @clusters[@max_density_cluster].mean      

      plot_histo_clusters(@filename)
      plot_length(@filename)

      silhouette = sequence_silhouette

      limits = @clusters[@max_density_cluster].get_limits
      max_len = @hits.map{|x| x.xml_length}.max
      predicted_len = @prediction.xml_length
      pval = @clusters[@max_density_cluster].t_test(@clusters, predicted_len)
      wval = @clusters[@max_density_cluster].wilcox_test(@clusters, predicted_len)
      deviation = @clusters[@max_density_cluster].deviation(@clusters, predicted_len)

      if predicted_len <= limits[1] and predicted_len >= limits[0]
        status = "YES"
      else
        status = "NO"
      end

      return [limits[0], limits[1]]
  end

  ## Test for reading frame inconsistency
  # Params:
  # +lst+: vector of +Sequence+ objects
  # Output:
  # yes/no answer
  def reading_frame_validation(lst = @hits)
    frames_histo = Hash[lst.group_by { |x| x.query_reading_frame }.map { |k, vs| [k, vs.length] }]
    rez = ""
    frames_histo.each do |x, y|
      rez << "#{x} #{y}; "
    end

    #if there are different reading frames of the same sign
    #count for positive reading frames
    count_p = 0
    count_n = 0
    frames_histo.each do |x, y|
      if x > 0
        count_p = count_p + 1
      else 
        if x < 0
          count_n = count_n + 1
        end
      end
    end

    if count_p > 1 or count_n > 1
      answ = "INVALID"
    else
      answ = "VALID"
    end

    return answ
  end

  def gene_merge_validation

    plot_matched_regions(@filename)

    # make a histogram with the middles of each matched subsequence from the predicted sequence
    middles = @hits.map { |hit| ((hit.match_query_from + hit.match_query_to)/2.0).to_i }

    #calculate the likelyhood to have a binomial distribution
    #split in two clusters
    #puts middles
    clusters = hierarchical_clustering(middles, 2, 1)

    max_freq = clusters.map{ |x| x.lengths.map{|y| y[1]}.max}.max

    R.eval "colors = c('red', 'blue', 'yellow', 'green', 'gray', 'orange')"
    R.eval "jpeg('#{@filename}_match_distr.jpg')"

    clusters.each_with_index do |cluster, i|
      cluster_values = cluster.lengths.sort{|a,b| a[0]<=>b[0]}.map{ |x| a = Array.new(x[1],x[0])}.flatten
      color = "colors[#{i%5+1}]"
         
      R.eval "hist(c#{cluster_values.to_s.gsub('[','(').gsub(']',')')},
                     breaks = 30,
                     xlim=c(#{middles.min-10},#{middles.max+10}),
                     ylim=c(0,#{max_freq}),
                     col=#{color},
                     border=#{color},
                     main='Predction match distribution (middle of the matches)', xlab='position idx', ylab='Frequency')"
      R.eval "par(new=T)"    
    end
    R.eval "dev.off()"
=begin
    R.eval "library(diptest)"
    R.eval "pval = dip.test(c#{middles.to_s.gsub('[','(').gsub(']',')')}, simulate.p.value = FALSE)$p.value"
    pval = R.pull("pval")
=end
    pval = 0
    
    wss1 = clusters[0].wss
    wss2 = clusters[1].wss   

    mean1 = clusters[0].mean
    mean2 = clusters[1].mean

    n1 = clusters[0].density
    n2 = clusters[1].density

    big_cluster = clusters[0].dup
    big_cluster.add(clusters[1])
    big_mean = big_cluster.mean
    big_wss = big_cluster.wss

    #puts "#{mean1} #{mean2} #{big_mean}"

    bss = n1 * (mean1 - big_mean) * (mean1 - big_mean)
    bss += n2 * (mean2 - big_mean) * (mean2 - big_mean)

    ratio = (wss1 + wss2) / (wss1 + wss2 + bss + 0.0)
    ratio_y = (wss1 + wss2) / (wss1 + wss2 + big_wss + 0.0)

    ratio.round(2)

  end

  ## 
  # Clusterization by length from a list of sequences
  # Params:
  # +lst+:: array of +Sequence+ objects
  # +predicted_seq+:: +Sequence+ objetc
  # +debug+ (optional):: true to display debug information, false by default (optional argument)
  # Output
  # output 1:: array of Cluster objects
  # output 2:: the index of the most dense cluster
  def clusterization_by_length(debug = false, lst = @hits, predicted_seq = @prediction)
    begin

      raise TypeError unless lst[0].is_a? Sequence and predicted_seq.is_a? Sequence

      contents = lst.map{ |x| x.xml_length.to_i }.sort{|a,b| a<=>b}

      clusters = hierarchical_clustering(contents,0)
      max_density = 0;
      max_density_cluster_idx = 0;
      clusters.each_with_index do |item, i|
        if item.density > max_density
          max_density = item.density
          max_density_cluster_idx = i;
        end
      end

      return [clusters, max_density_cluster_idx]

    rescue TypeError
      $stderr.print "Type error. Possible cause: one of the arguments of 'clusterization_by_length' method has not the proper type.\n"
      exit
    end
  end

  ######################################################################
  #calculates the silhouette score of the sequence
  #input1: Sequence object
  #input2: index of the target cluster within the clusters (see input 3)
  #input3: an array of Cluster objects
  #output: the silhouette of the sequence
  def sequence_silhouette (seq = @prediction, idx = @max_density_cluster, clusters = @clusters)
    seq_len = seq.xml_length    

    #the average dissimilarity of the sequence with other elements in idx cluster
    a = 0
    clusters[idx].lengths.each do |len, frecv|
      a = a + (len - seq_len).abs
    end
    a = a.to_f / clusters[idx].lengths.length
    
    b_vector = Array.new
    
    clusters.each_with_index do |cluster, i|
      #the average dissimilarity of the sequence with the members of cluster i
      if i != idx
        b = 0
        cluster.lengths.each do |len, frecv|
          b = b + (len - seq_len).abs
        end
        b = b.to_f / cluster.lengths.length
        unless b == 0
          b_vector.push(b)
        end
      end  
    end
    b = b_vector.min
    if b == nil
      b=0
    end
    silhouette = (b - a).to_f / [a,b].max
    return silhouette  
  end

  #########################################  
  #plots lines corresponding to each length
  #highlights the start and end hit offsets
  #input 1: lst = array of Sequence objects
  #input 2: predicted_seq = Sequence objetc
  def plot_length(output, lst = @hits, predicted_seq = @prediction)
    max_len = lst.map{|x| x.xml_length.to_i}.max
    lst = lst.sort{|a,b| a.xml_length<=>b.xml_length}

    max_plots = 120
    skip= lst.length/max_plots

    R.eval "jpeg('#{output}_len.jpg')"
    R.eval "plot(1:#{[lst.length-1,max_plots].min}, xlim=c(1,#{max_len}), xlab='Hit Length (black) vs part of the hit that matches the query (red)',ylab='Hit Number', col='white')"
    height = -1;
    lst.each_with_index do |seq,i|
      if skip == 0 or i%skip == 0
        height = height + 1
        R.eval "lines(c(1,#{seq.xml_length}), c(#{height}, #{height}), lwd=10)"
        R.eval "lines(c(#{seq.hit_from},#{seq.hit_to}), c(#{height}, #{height}), lwd=6, col='red')"
      end
    end
    R.eval "dev.off()"
  end

  #########################################  
  #plots lines corresponding to each hit
  #highlights the matched region of the prediction for each hit
  #input 1: lst = array of Sequence objects
  #input 2: predicted_seq = Sequence objetc
  def plot_matched_regions(output, lst = @hits, predicted_seq = @prediction)

    max_len = lst.map{|x| x.xml_length.to_i}.max

    max_plots = 120
    skip= lst.length/max_plots
    len = predicted_seq.xml_length

    R.eval "jpeg('#{output}_match.jpg')"
    R.eval "plot(1:#{[lst.length-1,max_plots].min}, xlim=c(1,#{len}), xlab='Prediction length (black) vs part of the prediction that matches hit x (yellow)',ylab='Hit Number', col='white')"
    height = -1;
    lst.each_with_index do |seq,i|
      if skip == 0 or i%skip == 0
        height = height + 1
        R.eval "lines(c(1,#{len}), c(#{height}, #{height}), lwd=10)"
        R.eval "lines(c(#{seq.match_query_from},#{seq.match_query_to}), c(#{height}, #{height}), lwd=6, col='yellow')"
      end
    end
    R.eval "dev.off()"
  end


  #############################################################################
  #plots a histogram of the length distribution given a list of Cluster objects
  #input 1: array of Cluster objects
  #input 2: length of the predicted sequence (number)
  #input 3: index from the clusters array of the most_dense_cluster_idx
  #input 4: name of the histogram file (optional argument, if missing the histogram will be displayed in a new window)
  def plot_histo_clusters(output, clusters = @clusters, predicted_length = @prediction.xml_length, 
                          most_dense_cluster_idx = @max_dense_cluster)
    begin
      raise TypeError unless clusters[0].is_a? Cluster and predicted_length.is_a? Fixnum

      lengths = clusters.map{ |c| c.lengths.sort{|a,b| a[0]<=>b[0]}.map{ |x| a = Array.new(x[1],x[0])}.flatten}.flatten
      lengths.push(predicted_length)

      max_freq = clusters.map{ |x| x.lengths.map{|y| y[1]}.max}.max
      #make the plot in a new process
        R.eval "colors = c('orange', 'blue', 'yellow', 'green', 'gray')"

        unless output == nil
          #puts "---- #{output}"
          #R.eval "dev.copy(png,'#{output}.png')"
          R.eval "jpeg('#{output}.jpg')"
        end

        clusters.each_with_index do |cluster, i|
          cluster_lengths = cluster.lengths.sort{|a,b| a[0]<=>b[0]}.map{ |x| a = Array.new(x[1],x[0])}.flatten

          if i == @max_density_cluster
            color = "'red'"
          else
            color = "colors[#{i%5+1}]"
          end

          R.eval "hist(c#{cluster_lengths.to_s.gsub('[','(').gsub(']',')')},
                     breaks = seq(#{lengths.min-10}, #{lengths.max+10}, 0.1),
                     xlim=c(#{lengths.min-10},#{lengths.max+10}),
                     ylim=c(0,#{max_freq}),
                     col=#{color},
                     border=#{color},
                     main='Histogram for length distribution', xlab='length\nblack = predicted sequence, red = most dense cluster')"
          R.eval "par(new=T)"
        end

        R.eval "abline(v=#{predicted_length})"

        unless output == nil
          R.eval "dev.off()"
        end

    rescue TypeError
      $stderr.print "Type error. Possible cause: one of the arguments of 'plot_histo_clusters' method has not the proper type.\n"
      exit
    end
  end

end


