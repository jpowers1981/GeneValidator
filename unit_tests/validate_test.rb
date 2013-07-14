require "rubygems"
require "test/unit"
require "shoulda"
require 'yaml'
require './clustering'

class OutputValidationTestn < Test::Unit::TestCase

  context "Hierarchical clusterization" do

    #read reference test file YAML
    #log = File.open( "test_reference.yaml" )
    thing = YAML.load_file('test_reference.yml')
    puts thing.inspect
=begin
    yp = YAML::load_documents( log ) do |doc|
      puts "#{doc['valid_length']} #{doc['valid_rf']}"
    end
=end    
    should "make clusterization " do
      hc = HierarchicalClusterization.new(vec)
      assert_equal 2, hc.hierarchical_clustering(vec,2).length
    end
  end
end
