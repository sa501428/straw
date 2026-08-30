Gem::Specification.new do |spec|
  spec.name = "hic-straw"
  spec.version = "1.0.0"
  spec.authors = ["Aiden Lab"]
  spec.homepage = "https://github.com/aidenlab/straw"
  spec.summary = "Thin Ruby FFI bindings for libstraw"
  spec.license = "MIT"
  spec.files = Dir["lib/**/*.rb", "lib/hic_straw/native/**/*", "README.md"]
  spec.require_paths = ["lib"]
  spec.required_ruby_version = ">= 2.6"
  spec.add_dependency "ffi", ">= 1.16"
end
