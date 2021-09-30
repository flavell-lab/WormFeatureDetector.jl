module WormFeatureDetector

using WormCurveFinder, MHDIO, Statistics, Plots, ProgressMeter, Images, ImageTransformations, HDF5, ImageSegmentation,
        LinearAlgebra, CoordinateTransformations, StaticArrays, ImageDataIO, FlavellBase, Interpolations, Rotations

include("worm_feature_detector.jl")
include("worm_curve_finder.jl")
include("heuristics.jl")
include("find_head.jl")

export
    find_gut_granules,
    find_head,
    find_head_unet,
    find_hsn,
    find_nerve_ring,
    in_conv_hull,
    curve_distance,
    elastix_difficulty_hsn_nr,
    elastix_difficulty_wormcurve!
end # module
