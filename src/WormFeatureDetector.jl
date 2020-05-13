module WormFeatureDetector

using WormCurveFinder, MHDIO, Statistics, Plots, ProgressMeter, Images, ImageTransformations

include("worm_feature_detector.jl")
include("worm_curve_finder.jl")
include("heuristics.jl")
include("find_head.jl")
include("crop_worm.jl")
include("centroids_io.jl")

export
    create_dir,
    output_gut_granules,
    find_hsn,
    find_nerve_ring,
    curve_distance,
    elastix_difficulty_hsn_nr,
    elastix_difficulty_wormcurve,
    find_head,
    write_centroids,
    read_centroids_transformix,
    read_centroids_roi,
    centroids_to_img,
    crop_rotate,
    crop_rotate_output,
    crop_rotate_images
end # module
