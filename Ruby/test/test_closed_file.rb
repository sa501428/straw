require "minitest/autorun"
require_relative "../lib/hic_straw"

class ClosedFileTest < Minitest::Test
  def test_methods_reject_closed_native_handle
    path = File.expand_path("../../R/inst/extdata/test.hic", __dir__)
    file = HicStraw::File.new(path)
    assert_operator file.version, :>=, 6
    file.close
    file.close
    [-> { file.version }, -> { file.genome }, -> { file.chromosomes },
     -> { file.resolutions },
     -> { file.records("1", "2", resolution: 2_500_000) },
     -> { file.raw("1", "2", resolution: 2_500_000) }].each do |operation|
      assert_raises(IOError, &operation)
    end
  ensure
    file&.close
  end
end
