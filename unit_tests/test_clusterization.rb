require "rubygems"
require "test/unit"
require "shoulda"
require './clustering'

class TestHierarchicalClusterization < Test::Unit::TestCase

  context "Hierarchical clusterization" do

    vec = [4,5,8,11,11,14,15,15,15,15,15,16,17,17,20]   

    should "make clusterization " do
      hc = HierarchicalClusterization.new(vec)
      assert_equal 2, hc.hierarchical_clustering(vec,2).length
    end

    should "most dense cluster, method 1" do
      hc = HierarchicalClusterization.new(vec)
      hc.hierarchical_clustering(vec)
      result = {14=>1, 15=>5, 16=>1, 17=>2}
      assert_equal result , hc.most_dense_cluster.lengths
    end

    should "most dense cluster, method 2" do
      hc = HierarchicalClusterization.new(vec)
      hc.hierarchical_clustering(vec, 0, 1)
      result = {14=>1, 15=>5, 16=>1, 17=>2}
      assert_equal result , hc.most_dense_cluster.lengths
    end

    should "most dense cluster mean" do
      hc = HierarchicalClusterization.new(vec)
      hc.hierarchical_clustering(vec)
      assert_equal 15 , hc.most_dense_cluster.mean
    end

  end
end
