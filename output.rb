
require 'yaml'

class Output

  attr_accessor :prediction_len
  attr_accessor :prediction_def
  attr_accessor :nr_hits

  attr_accessor :lv_cluster
  attr_accessor :length_cluster_limits  

  attr_accessor :length_rank_score
  attr_accessor :length_rank_msg

  attr_accessor :reading_frame_validation

  attr_accessor :merged_genes_score
  attr_accessor :duplication

  attr_accessor :filename
  attr_accessor :image_histo_len
  attr_accessor :image_plot_merge
  attr_accessor :image_histo_merge
  attr_accessor :idx

  def initialize(filename, idx)

    @prediction_len = 0
    @prediction_def = "no_definition"
    @nr_hits = 0

    @lv_cluster = "?"
    @length_cluster_limits = [0,0]
    
    @length_rank_score = 0
    @length_rank_msg = "?"
    
    @reading_frame_validation = "?"
    @merged_genes_score = 0
    @duplication = "?"

    @filename = filename
    @idx = idx

  end
  
  def print_output_console

    if @prediction_len >= @length_cluster_limits[0] and @prediction_len <= @length_cluster_limits[1]
      @lv_cluster = "YES"
    else
      @lv_cluster = "NO"
    end

    printf " %3d | %25s | %25s | %s(%s) | %15s| %15s | %10s |\n",              
              @idx,
              @prediction_def.scan(/([^ ]+)/)[0][0],
              "#{@prediction_len} #{@length_cluster_limits} #{@lv_cluster}", 
              @length_rank_msg, @length_rank_score, 
              @reading_frame_validation,
              @merged_genes_score, @duplication

  end

  def print_output_file_yaml
    File.open("#{@filename}.yaml", "w+") do |f|
      f.write({@prediction_def.scan(/([^ ]+)/)[0][0] => self}.to_yaml)
    end
  end

  def generate_html
    output = "<h1>#{@prediction_def}</h1>  \\
		<h3>Length Validation</h3>  \\
		<p><img alt=img src= #{@image_histo_len} style=\"width: 300px;\" /></p>		\\
		<ul>  \\
			<li>clustering : YES</li>  \\
			<li>ranking</li>  \\
		</ul> \\
		<h3> Reading Frame Validation:</h3> \\
		<h3> Gene Merge Validation:</h3> \\
		<p> \\
			<img alt=img src=#{@image_histo_len} style=\"width: 300px; height: 300px; float: left; margin-left: 20px; margin-right: 20px;\" /></p> \\
		<h3> \\
			<img alt=img src=#{@image_histo_len} style=\"width: 300px; height: 300px; float: left; margin-left: 20px; margin-right: 20px;\" /></h3>"
  end
end
