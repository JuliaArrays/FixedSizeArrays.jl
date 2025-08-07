# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Version [v1.2.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.2.0) - 2025-08-07

### Added

* Pkgextension with [`Adapt.jl`](https://github.com/JuliaGPU/Adapt.jl), for converting other arrays to `FixedSizeArray`. ([#149](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/149))

## Version [v1.1.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.1.0) - 2025-08-05

### Added

* Support for broadcasting of view of `FixedSizeArray`s. ([#146](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/146), [#147](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/147))

## Version [v1.0.1](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.0.1) - 2025-06-08

### Fixed

* `rand!(A::FixedSizeArray)` now returns `A`, instead of its parent. ([#139](https://github.com/JuliaArrays/FixedSizeArrays.jl/issues/139), [#140](https://github.com/JuliaArrays/FixedSizeArrays.jl/pull/140))

## Version [v1.0.0](https://github.com/JuliaArrays/FixedSizeArrays.jl/releases/tag/v1.0.0) - 2025-06-04

Initial release.
