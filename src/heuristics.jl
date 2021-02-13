"""
Computes registration difficulty between two frames based on the HSN and nerve ring locations

# Arguments:
- `rootpath::String`: working directory path; all other directory inputs are relative to this
- `frame1::Integer`: first frame
- `frame2::Integer`: second frame
- `hsn_location_file::String`: path to file containing HSN locations
- `nr_location_file::String`: path to file containing nerve ring locations
# Heuristic parameters (optional):
- `nr_weight::Real`: weight of nerve ring location relative to HSN location. Default 1.
"""
function elastix_difficulty_hsn_nr(rootpath::String, frame1::Integer, frame2::Integer, 
        hsn_location_file::String, nr_location_file::String; nr_weight::Real=1)
    # load HSN locations
    hsn_locs = Dict()
    open(joinpath(rootpath, hsn_location_file)) do hsn
        for line in eachline(hsn)
            s = [parse(Int32, x) for x in split(line)]
            hsn_locs[s[1]] = s[2:4]
        end
    end
    # load nerve ring locations
    nr_locs = Dict()
    open(joinpath(rootpath, nr_location_file)) do nr
        for line in eachline(nr)
            s = [parse(Int32, x) for x in split(line)]
            nr_locs[s[1]] = s[2:4]
        end
    end

    return sqrt(sum((hsn_locs[frame1] .- hsn_locs[frame2]).^2))
                + nr_weight * sqrt(sum((nr_locs[frame1] .- nr_locs[frame2]).^2))

end

"""
Computes registration difficulty between two frames based on the worm curvature heuristic.
Requires that the data be filtered in some way (eg: total-variation filtering),
and that the head position of the worm is known in each frame.

# Arguments:
- `dict_curve::Dict`: dictionary of worm curves found so far. The method will attempt to find the worm's curvature in this dictionary,
    and will compute and add it to the dictionary if not found.
- `img1::Array{<:AbstractFloat,3}`: image 1 array (volume)
- `img2::Array{<:AbstractFloat,3}`: image 2 array (volume)
- `t1::Int`: index of frame 1
- `t2::Int`: index of frame 2
- `head_pos::Dict`: head position dictionaru

## Other parameters (optional):
- `path_dir_fig::Union{Nothing,String}`: Path to save figures of worm curvature. If `nothing`, figures will not be generated.
- `downscale::Integer`: log2(factor) by which to downscale the image before processing. Default 3 (ie: downscale by a factor of 8)
- `num_points::Integer`: number of points (not including head) in generated curve. Default 9.
- `headpt::Integer`: First position from head (in index of curves) to be aligned. Default 4.
- `tailpt::Integer`: Second position from head (in index of curves) to be aligned. Default 7.
"""
function elastix_difficulty_wormcurve!(dict_curve::Dict, img1::Array{<:AbstractFloat,3}, img2::Array{<:AbstractFloat,3},
        t1::Int, t2::Int, head_pos_t1::Dict, head_pos_t2::Dict; downscale::Int=3, num_points::Int=9, headpt::Int=4, tailpt::Int=7,
        path_dir_fig::Union{Nothing,String}=nothing)

    img1 = maxprj(img1, dims=3)
    img2 = maxprj(img2, dims=3)

    if haskey(dict_curve, t1)
        x1_c, y1_c = dict_curve[t1]
    else
        x1_c, y1_c = find_curve(img1, downscale, head_pos_t1[t1]./2^downscale, num_points)
        dict_curve[t1] = (x1_c, y1_c)
    end
    
    if haskey(dict_curve, t2)
        x2_c, y2_c = dict_curve[t2]
    else
        x2_c, y2_c = find_curve(img2, downscale, head_pos_t2[t2]./2^downscale, num_points)
        dict_curve[t2] = (x2_c, y2_c)
    end

    if !isnothing(path_dir_fig)
        create_dir(path_dir_fig)
        if !haskey(dict_curve, t1)
            fig = heatmap(transpose(img1), fillcolor=:grays, aspect_ratio=1, flip=false, showaxis=false, legend=false)
            scatter!(fig, x1_c.-1, y1_c.-1, color="red");
            scatter!(fig, [x1_c[1].-1], [y1_c[1].-1], color="cyan", markersize=5);
            savefig(fig, joinpath(path_dir_fig, "$(t1).png"));
        end
        if !haskey(dict_curve, t2)
            fig = heatmap(transpose(img2), fillcolor=:grays, aspect_ratio=1, flip=false, showaxis=false, legend=false)
            scatter!(fig, x2_c.-1, y2_c.-1, color="red");
            scatter!(fig, [x2_c[1].-1], [y2_c[1].-1], color="cyan", markersize=5);
            savefig(fig, joinpath(path_dir_fig, "$(t2).png"));
        end
    end

    return curve_distance(x1_c, y1_c, x2_c, y2_c, headpt=headpt, tailpt=tailpt)
end

# TODO: document
function elastix_difficulty_wormcurve!(dict_curve::Dict, param::Dict, param_path_fixed::Dict, param_path_moving::Dict, t1::Int, t2::Int, ch::Int;
        save_curve_fig=false, max_fixed_t::Union{Integer,Nothing}=nothing)

    worm_curve_n_pts = param["worm_curve_n_pts"] 
    worm_curve_tail_idx = param["worm_curve_tail_idx"]
    worm_curve_head_idx = param["worm_curve_head_idx"]
    worm_curve_downscale = param["worm_curve_downscale"]
    
    if isnothing(max_fixed_t)
        max_fixed_t = 0
    else
        if t1 > max_fixed_t || t2 <= max_fixed_t
            return Inf
        end
    end

    path_mhd_t1 = joinpath(param_path_fixed["path_dir_mhd_filt"], param_path_fixed["get_basename"](t1, ch) * ".mhd")
    path_mhd_t2 = joinpath(param_path_moving["path_dir_mhd_filt"], param_path_moving["get_basename"](t2 - max_fixed_t, ch) * ".mhd")

    img1 = Float64.(read_img(MHD(path_mhd_t1)))
    img2 = Float64.(read_img(MHD(path_mhd_t2)))
    
    head_pos_t1 = read_head_pos(param_path_fixed["path_head_pos"])
    head_pos_t2 = read_head_pos(param_path_moving["path_head_pos"])

    if !isnothing(max_fixed_t)
        head_pos_t2_shifted = Dict()
        for t in keys(head_pos_t2)
            head_pos_t2_shifted[t + max_fixed_t] = head_pos_t2[t]
        end
        head_pos_t2 = head_pos_t2_shifted
    end


    path_dir_worm_curve = save_curve_fig ? param_path_fixed["path_dir_worm_curve"] : nothing
    
    elastix_difficulty_wormcurve!(dict_curve, img1, img2, t1, t2, head_pos_t1, head_pos_t2,
        downscale=worm_curve_downscale, num_points=worm_curve_n_pts, headpt=param["worm_curve_head_idx"],
        tailpt=param["worm_curve_tail_idx"], path_dir_fig=path_dir_worm_curve)
end
