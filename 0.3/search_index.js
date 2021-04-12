var documenterSearchIndex = {"docs":
[{"location":"head/#Head-Detection-API","page":"Head Detection API","title":"Head Detection API","text":"","category":"section"},{"location":"head/","page":"Head Detection API","title":"Head Detection API","text":"find_head\nin_conv_hull","category":"page"},{"location":"head/#WormFeatureDetector.find_head","page":"Head Detection API","title":"WormFeatureDetector.find_head","text":"Finds the tip of the nose of the worm in each time point, and warns of bad time points. Uses a series of blob-approximations of the worm with different sensitivities, by using local convex hull. The convex hulls should be set up in increasing order (so the last convex hull is the most generous). The difference between the first two convex hulls is used to determine the direction of the worm's head. The third convex hull is used to find the tip of the worm's head.\n\nArguments\n\ncentroids: the locations of the neuron centroids\nimsize: the image size\n\nOptional keyword arguments\n\ntf (default [10,10,30]): threshold for required neuron density for convex hull i is (number of centroids) / tf[i]\nmax_d (default [30,50,50]): the maximum distance for a neuron to be counted as part of convex hull i is max_d[i]\nhd_threshold::Integer (default 100): if convex hulls 2 and 3 give head locations farther apart than this many pixels, set error flag.\nvc_threshold::Integer (default 300): if convex hulls 2 and 3 give tail locations farther apart than this many pixels, set error flag.\nnum_centroids_threshold::Integer (default 90): if there are fewer than this many centroids, set error flag.\nedge_threshold::Integer (default 5): if the boundary of the worm is closer than this to the edge of the frame, set error flag.\n\nOutputs a tuple (head_pos, q_flag, crop_x, crop_y, crop_z, theta, centroid)\n\nhead_pos: position of the worm's head.\nq_flag: array of error flags (empty means no errors were detected)\ncrop_x, crop_y, crop_z: cropping parameters\ntheta: amount by which to rotate the image to align it in the x-y plane\ncentroid: centroid of the worm\n\n\n\n\n\nFinds the tip of the nose of the worm in each time point, and warns of bad time points. Uses a series of blob-approximations of the worm with different sensitivities, by using local convex hull. The convex hulls should be set up in increasing order (so the last convex hull is the most generous). The difference between the first two convex hulls is used to determine the direction of the worm's head. The third convex hull is used to find the tip of the worm's head.\n\nArguments:\n\nparam_path::Dict: Dictionary of paths including the keys:\npath_head_pos: Path to head position output file\npath_dir_mhd_crop: Path to MHD input files\npath_dir_centroid: Path to centroid input files\nparam::Dict: Dictionary of parameter settings including the keys:\nhead_threshold: threshold for required neuron density for convex hull i is (number of centroids) / param[\"head_threshold\"][i]\nhead_max_distance: the maximum distance for a neuron to be counted as part of convex hull i is max_d[i]\nhead_err_threshold: if convex hulls 2 and 3 give head locations farther apart than this many pixels, set error flag.\nhead_vc_err_threshold: if convex hulls 2 and 3 give tail locations farther apart than this many pixels, set error flag.\nhead_edge_err_threshold: if the boundary of the worm is closer than this to the edge of the frame, set error flag.\nt_range: The time points to compute head location\nf_basename::Function: Function that takes as input a time point and a channel and gives the base name of the corresponding MHD file.\n\n\n\n\n\n","category":"function"},{"location":"head/#WormFeatureDetector.in_conv_hull","page":"Head Detection API","title":"WormFeatureDetector.in_conv_hull","text":"Computes whether a point is in the local convex hull of a collection of 2D points.\n\nArguments:\n\npoint: a point in 2D space.\ncentroids: a cloud of points in 2D space.\nmax_d::Real: the farthest a point in centroids can be from point and still count for the convex hull\n\n\n\n\n\n","category":"function"},{"location":"curve/#Worm-Curvature-API","page":"Worm Curvature API","title":"Worm Curvature API","text":"","category":"section"},{"location":"curve/","page":"Worm Curvature API","title":"Worm Curvature API","text":"elastix_difficulty_wormcurve!\ncurve_distance","category":"page"},{"location":"curve/#WormFeatureDetector.elastix_difficulty_wormcurve!","page":"Worm Curvature API","title":"WormFeatureDetector.elastix_difficulty_wormcurve!","text":"Computes registration difficulty between two time points based on the worm curvature heuristic. Requires that the data be filtered in some way (eg: total-variation filtering), and that the head position of the worm is known in each time point.\n\nArguments:\n\ndict_curve::Dict: dictionary of worm curves found so far. The method will attempt to find the worm's curvature in this dictionary,   and will compute and add it to the dictionary if not found.\nimg1::Array{<:AbstractFloat,3}: image 1 array (volume)\nimg2::Array{<:AbstractFloat,3}: image 2 array (volume)\nt1::Int: time point 1\nt2::Int: time point 2\nhead_pos_t1::Dict: head position dictionary at time point 1\nhead_pos_t2::Dict: head position dictionary at time point 2\n\nOther parameters (optional):\n\npath_dir_fig::Union{Nothing,String}: Path to save figures of worm curvature. If nothing, figures will not be generated.\ndownscale::Integer: log2(factor) by which to downscale the image before processing. Default 3 (ie: downscale by a factor of 8)\nnum_points::Integer: number of points (not including head) in generated curve. Default 9.\nheadpt::Integer: First position from head (in index of curves) to be aligned. Default 4.\ntailpt::Integer: Second position from head (in index of curves) to be aligned. Default 7.\n\n\n\n\n\nComputes registration difficulty between two time points based on the worm curvature heuristic. Requires that the data be filtered in some way (eg: total-variation filtering), and that the head position of the worm is known in each time point.\n\nArguments\n\ndict_curve::Dict: Dictionary of worm curves found so far. The method will attempt to find the worm's curvature in this dictionary,\n\nand will compute and add it to the dictionary if not found.\n\nparam::Dict: Dictionary of parameter settings, including:\nworm_curve_n_pts: number of points (not including head) in generated curve.\nworm_curve_head_idx: First position from head (in index of curves) to be aligned.\nworm_curve_tail_idx: Second position from head (in index of curves) to be aligned.\nworm_curve_downscale: log2(factor) by which to downscale the image before processing.\nparam_path_fixed::Dict: Dictionary of paths for the fixed (t1) time point, including:\npath_dir_mhd_filt: Path to filtered cropped MHD files\npath_head_pos: Path to head position\npath_dir_worm_curve: Path to save worm curve images\nget_basename: Function that takes as input a time point and a channel and gives the base name of the corresponding MHD file. \nparam_path_moving::Dict: Dictionary of paths for the moving (t2) time point, including the same keys as for the fixed dictionary.\nt1::Int: Fixed time point\nt2::Int: Moving time point\nch::Int: Channel\nsave_curve_fig::Bool (optional, default false): whether to save worm curve images\nmax_fixed_t::Union{Integer,Nothing} (optional, default nothing): If using two different data sets, the maximum time point in the fixed dataset.   Moving dataset time points will be incremented by this amount.\n\n\n\n\n\n","category":"function"},{"location":"curve/#WormFeatureDetector.curve_distance","page":"Worm Curvature API","title":"WormFeatureDetector.curve_distance","text":"Computes the difficulty of an elastix transform using the heuristic that more worm-unbending is harder.\n\nArguments:\n\nx1_c: Array of x-coordinates of first worm\ny1_c: Array of y-coordinates of first worm\nx2_c: Array of x-coordinates of second worm\ny2_c: Array of y-coordinates of second worm\nheadpt::Integer: First position from head (in index of curves) to be aligned. Default 4.\ntailpt::Integer: Second position from head (in index of curves) to be aligned. Default 7.\n\n\n\n\n\n","category":"function"},{"location":"#WormFeatureDetector.jl-Documentation","page":"WormFeatureDetector.jl Documentation","title":"WormFeatureDetector.jl Documentation","text":"","category":"section"},{"location":"","page":"WormFeatureDetector.jl Documentation","title":"WormFeatureDetector.jl Documentation","text":"Pages = [\"head.md\", \"curve.md\", \"hsn.md\"]","category":"page"},{"location":"hsn/#HSN-Feature-Detection-API","page":"HSN Feature Detection API","title":"HSN Feature Detection API","text":"","category":"section"},{"location":"hsn/","page":"HSN Feature Detection API","title":"HSN Feature Detection API","text":"find_gut_granules\nfind_hsn\nfind_nerve_ring\nelastix_difficulty_hsn_nr","category":"page"},{"location":"hsn/#WormFeatureDetector.find_gut_granules","page":"HSN Feature Detection API","title":"WormFeatureDetector.find_gut_granules","text":"Finds, returns, and outputs gut granule masks for a set of images. The masks filter out the gut granules, so 1 is not a gut granule, and 0 is a gut granule.\n\nArguments\n\npath::String: working directory path; all other directory inputs are relative to this\nnames::Array{String,1}: filenames of images to process\nthreshold::Real: pixel intensity brightness value. Pixels below this intensity are excluded\ndensity::Real: density of nearby pixels that must meet the threshold for the original pixel to be counted\nradius: distances in each dimension (in pixels) away from the original pixel that are counted as nearby.   For example, radius = [3, 2, 1] would allow a distance of three pixels in the x-direction, two pixels in the y-direction,   and one pixel in the z-direction.\n\nOptional keyword arguments\n\nmhd::String: path to MHD directory, where the image will be found. Default \"MHD\".\nout::String: path to output directory to store the masks. If left blank (default), images will not be written.\n\n\n\n\n\n","category":"function"},{"location":"hsn/#WormFeatureDetector.find_hsn","page":"HSN Feature Detection API","title":"WormFeatureDetector.find_hsn","text":"Finds location of the HSN soma in a frame. First threshold to remove densely-packed regions that might correspond to gut fluorescence or neuropil (by excluding too-dense regions), then threshold again to ensure high local density (by excluding not-dense regions). If multiple regions are still included, the larger region is chosen. Can optionally choose to output data to a file, for use with heuristics.\n\nArguments:\n\npath::String: working directory path; all other directory inputs are relative to this\nframes: frames of images to process\nmhd::String: path to MHD directory, where the image will be found.\nimg_prefix::String: image prefix not including the timestamp. It is assumed that each frame's filename   will be, eg, img_prefix_t0123_ch2.mhd for frame 123 with channel=2.\nchannel::Integer: channel being used.\nthreshold_outer::Real: pixel intensity brightness value. Pixels below this intensity are excluded\ndensity_outer::Real: density of nearby pixels that must meet the outer threshold for the original pixel to NOT be counted\nradius_outer: distances in each dimension (in pixels) away from the original pixel that are counted as nearby.\n\nFor example, radius = [3, 2, 1] would allow a distance of three pixels in the x-direction, two pixels in the y-direction, and one pixel in the z-direction.\n\nthreshold_inner::Real: pixel intensity brightness value. Pixels below this intensity are excluded\ndensity_inner::Real: density of nearby pixels that must meet the inner threshold for the original pixel to be counted\nradius_inner: distances in each dimension (in pixels) away from the original pixel that are counted as nearby.\n\nFor example, radius = [3, 2, 1] would allow a distance of three pixels in the x-direction, two pixels in the y-direction, and one pixel in the z-direction.\n\nradius_detection: If multiple locations remain as possible HSN locations after both thresholding steps,   the location with the most other such points within radius_detection of it is chosen. \n\nOptional keyword arguments\n\noutfile::String: path to HSN output file. If left blank (default), no output will be written.\n\n\n\n\n\n","category":"function"},{"location":"hsn/#WormFeatureDetector.find_nerve_ring","page":"HSN Feature Detection API","title":"WormFeatureDetector.find_nerve_ring","text":"Finds location of the nerve ring in a frame. Can optionally output nerve ring locations to a file for use with heuristics.\n\nArguments:\n\npath::String: working directory path; all other directory inputs are relative to this\nframes: frames of images to process\nmhd::String: path to MHD directory, where the image will be found.\nimg_prefix::String: image prefix not including the timestamp. It is assumed that each frame's filename   will be, eg, img_prefix_t0123_ch2.mhd for frame 123 with channel=2.\nchannel::Integer: channel being used.\nthreshold::Real: pixel intensity brightness value. Pixels below this intensity are excluded\nregion: region of the image that will be searched for the nerve ring. Generically, you should try to include the nerve ring   and exclude any other regions (such as gut granules, HSN soma, etc)\nradius: The location with the most other points that meet the threshold within radius of it is chosen   as the nerve ring location.\n\nOptional keyword arguments\n\noutfile::String: path to nerve ring output file. If left blank (default), no output will be written.\n\n\n\n\n\n","category":"function"},{"location":"hsn/#WormFeatureDetector.elastix_difficulty_hsn_nr","page":"HSN Feature Detection API","title":"WormFeatureDetector.elastix_difficulty_hsn_nr","text":"Computes registration difficulty between two frames based on the HSN and nerve ring locations\n\nArguments:\n\nrootpath::String: working directory path; all other directory inputs are relative to this\nframe1::Integer: first frame\nframe2::Integer: second frame\nhsn_location_file::String: path to file containing HSN locations\nnr_location_file::String: path to file containing nerve ring locations\n\nHeuristic parameters (optional):\n\nnr_weight::Real: weight of nerve ring location relative to HSN location. Default 1.\n\n\n\n\n\n","category":"function"}]
}