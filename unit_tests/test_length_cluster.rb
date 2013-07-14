require "rubygems"
require "calculator"
require "test/unit"
require "shoulda"

class TC_Calculator_Shoulda < Test::Unit::TestCase

  context "Gene prediction validation" do

    should "check if prediction lies in the densest cluster" do
      assert_equal :yes,Calculator.new(3,4).addition
    end

    should "subtraction of two numbers " do
     assert_equal 1,Calculator.new(4,3).subtraction
    end

    should "multiplication of two numbers" do
      assert_equal 12,Calculator.new(3,4).multiplication
    end

    should "division of two numbers" do
      assert_equal 4,Calculator.new(12,3).division
    end
  end

end
