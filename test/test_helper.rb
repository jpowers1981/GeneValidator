require 'codeclimate-test-reporter'
require 'simplecov'
require 'coveralls'

SimpleCov.formatter = SimpleCov::Formatter::MultiFormatter[
  CodeClimate::TestReporter::Formatter,
  Coveralls::SimpleCov::Formatter,
  SimpleCov::Formatter::HTMLFormatter
]
SimpleCov.start 'rails'
