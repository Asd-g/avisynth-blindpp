## Description

An AviSynth+ plugin for MPEG-2 / MPEG-4 8x8 pixel blocks deblock and / or dering filter.

It is originally used as part of [DGDecode](https://www.rationalqm.us/dgmpgdec/DGDecodeManual.html#BlindPP).

### Requirements:

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

#### Usage:

```
BlindPP(clip clip, int "quant", int "cpu", string "cpu2", bool "iPP", int "moderate_h", int "moderate_v", int "opt")
```

#### Parameters:

##### ***`clip`***
A clip to process.<br>
It must be in a planar Y or YUV planar format. The frame dimensions must be a multiple of 16 (mod16).

##### ***`quant`***
The quantization level used to determine the strength of the deblocking and deringing filters.<br>
Must be between `0..63`.<br>
Default: `2`.

##### ***`cpu`***
Predefined post-processing levels. This parameter is ignored if `cpu2` is specified.<br>
* `0`: None.
* `1`: Luma horizontal deblocking.
* `2`: Luma horizontal and vertical deblocking.
* `3`: Luma deblocking + Chroma horizontal deblocking.
* `4`: Luma deblocking + Chroma horizontal and vertical deblocking.
* `5`: All deblocking + Luma deringing.
* `6`: All deblocking + Luma and Chroma deringing.
Default: `6`.

##### ***`cpu2`***
A custom configuration string of exactly 6 characters. Each character corresponds to a specific filter. Use `'x'` or `'X'` to enable a
filter, and any other character to disable it.<br>
* Position 1: Luma horizontal deblocking.
* Position 2: Luma vertical deblocking.
* Position 3: Chroma horizontal deblocking.
* Position 4: Chroma vertical deblocking.
* Position 5: Luma deringing.
* Position 6: Chroma deringing.
Example: `"xxxx--"` enables all deblocking but disables all deringing.<br>
Default: `""` (empty string, uses `cpu` instead).

##### ***`iPP`***
Enables interlaced post-processing. If true, the filter processes the clip as interlaced material by handling fields separately.<br>
Default: `false`.

##### ***`moderate_h`***
Sensitivity threshold for horizontal deblocking. Higher values make the filter more aggressive in detecting block edges.<br>
Must be between `0..255`.<br>
Default: `20`.

##### ***`moderate_v`***
Sensitivity threshold for vertical deblocking. Higher values make the filter more aggressive in detecting block edges.<br>
Must be between `0..255`.<br>
Default: `40`.

##### ***`opt`***
Sets the CPU optimization level.<br>
* `-1`: Auto-detect (uses the highest instruction set supported by your CPU).
* `0`: C++ (No SIMD optimization).
* `1`: SSE2.
* `2`: SSE4.1.
* `3`: AVX2.
* `4`: AVX-512.<br>
Default: `-1`.

### Building:

```
Requirements:
    - CMake
    - Ninja (optional)
    - C++20 compiler
```

```
Steps:
    git clone --depth 1 https://github.com/Asd-g/avisynth-blindpp
    cd avisynth-blindpp
    cmake -B build -G Ninja -DCMAKE_PREFIX_PATH="avisynth_headers"
    ninja -C build
```
