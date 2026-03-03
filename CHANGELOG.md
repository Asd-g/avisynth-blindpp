# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-03-03

### Added
- Support for 10, 12, 14, 16-bit integer pixel formats.
- Support for 32-bit floating point processing.
- SIMD code for SSE2, SSE4.1, AVX2, and AVX-512 instruction sets.
- Parameter `opt`.
- Support for Y and YUV444 planar formats.
- CMake building system.

### Changed
- Increased the maximum `quant` limit from 31 to 63 to allow for finer control over filtering strength.
- Replaced inline MMX assembly with C++ code and SIMD intrinsics.
