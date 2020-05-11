# WormFeatureDetector.jl
A collection of heuristics for locating various features of the worm.
### Worm curvature similarity heuristic
This heuristic (implemented by the `generate_elastix_difficulty_wormcurve` method) posits that two frames are similar to each other if the worm's curvature is similar, as this would result in a smaller amount of bending. It computes an estimate for the worm's centerline based on the images, and outputs its centerline fits as images which can be inspected for errors.

This heuristic requires data that has nuclear-localized fluorescent proteins in enough neurons to get an estimate of the worm shape, and has already been filtered (eg by the `GPUFilter.jl` package). It also requires you to have previously determined the worm's head location, such as in the `WormFeatureDetector.jl` package.

Example code:

```julia
generate_elastix_difficulty_wormcurve("/path/to/data", "MHD", "head_pos.txt", "img_prefix", 2, 1:100, "elastix_difficulty.txt", "worm_curves")
```

### HSN and nerve ring location heuristic
This heuristic (implemented by the `generate_elastix_difficulty_HSN_NR` method) tries to identify frames with similar HSN and nerve ring locations to be registered together. It only works on data taken with a non-nuclear-localized fluorescent protein expressed only in HSN, and it also requires you to have previously identified the HSN and nerve ring locations in each frame, which is implemented in the `WormFeatureDetector.jl` package.

Example code:

```julia
generate_elastix_difficulty_HSN_NR("/path/to/data", "hsn_locs.txt", "nr_locs.txt", 1:100, "elastix_difficulty.txt")
```