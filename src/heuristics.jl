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
Reads the worm head position from the file `head_path::String`.
Returns a dictionary mapping frame => head position of the worm at that frame.
"""
function read_head_pos(head_path::String)
    head_pos = Dict()
    open(head_path) do f
        for line in eachline(f)
            l = split(line)
            head_pos[parse(Int16, l[1])] = Tuple(map(x->parse(Float64, x), l[2:end]))
        end
    end
    return head_pos
end

"""
Computes registration difficulty between two frames based on the worm curvature heuristic.
Requires that the data be filtered in some way (eg: total-variation filtering),
and that the head position of the worm is known in each frame.

# Arguments:
- `rootpath::String`: working directory path; all other directory inputs are relative to this
- `frame1::Integer`: first frame
- `frame2::Integer`: second frame
- `mhd_path::String`: path to MHD image files
- `head_path::String`: path to a file containing positions of the worm's head at each frame.
- `img_prefix::String`: image prefix not including the timestamp. It is assumed that each frame's filename
    will be, eg, `img_prefix_t0123_ch2.mhd` for frame 123 with channel=2.
- `channel::Integer`: channel being used.

## Other parameters (optional):
- `figure_save_path::String`: Path to save figures of worm curvature. If left blank, figures will not be generated.
    If supplied, the figure for frame 1 will be generated.
- `downscale::Integer`: log2(factor) by which to downscale the image before processing. Default 3 (ie: downscale by a factor of 8)
- `num_points::Integer`: number of points (not including head) in generated curve. Default 9.
- `headpt::Integer`: First position from head (in index of curves) to be aligned. Default 4.
- `tailpt::Integer`: Second position from head (in index of curves) to be aligned. Default 7.
"""
function elastix_difficulty_wormcurve(rootpath::String, frame1::Integer, frame2::Integer, mhd_path::String, head_path::String,
        img_prefix::String, channel::Integer; figure_save_path::String="",
        downscale::Integer=3, num_points::Integer=9, headpt::Integer=4, tailpt::Integer=7)


    img1_path = joinpath(rootpath, mhd_path, img_prefix*"_t"*string(frame1, pad=4)*"_ch$(channel).mhd")
    img1 = Float64.(maxprj(read_img(MHD(img1_path)), dims=3))

    img2_path = joinpath(rootpath, mhd_path, img_prefix*"_t"*string(frame2, pad=4)*"_ch$(channel).mhd")
    img2 = Float64.(maxprj(read_img(MHD(img2_path)), dims=3))

    head_pos = read_head_pos(joinpath(rootpath, head_path))

    x1_c, y1_c = find_curve(img1, downscale, head_pos[frame1]./2^downscale, num_points)
    x2_c, y2_c = find_curve(img1, downscale, head_pos[frame2]./2^downscale, num_points)
    if figure_save_path != nothing
        create_dir(joinpath(rootpath, figure_save_path))
        fig = heatmap(transpose(img1), fillcolor=:grayscale, aspect_ratio=1, flip=false, showaxis=false, legend=false)
        scatter!(fig, x1_c.-1, y1_c.-1, color="red");
        scatter!(fig, [x1_c[1].-1], [y1_c[1].-1], color="cyan", markersize=5);
        savefig(fig, joinpath(rootpath, figure_save_path, "$(frames[i]).png"));
    end

    return curve_distance(x1_c, y1_c, x2_c, y2_c, headpt=headpt, tailpt=tailpt)
end
