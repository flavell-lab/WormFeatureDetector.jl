module WormFeatureDetector

using WormCurveFinder, MHDIO, Statistics, Plots, ProgressMeter

include("worm_feature_detector.jl")
include("worm_curve_finder.jl")

export
    create_dir,
    output_gut_granules,
    find_hsn,
    find_nerve_ring,
    curve_distance,
    elastix_difficulty_hsn_nr,
    elastix_difficulty_wormcurve
end # module
