"""
Rotates and then crops an image, along with its head and centroid locations.

# Arguments

`img`: image to transform (3D)
`crop_x`: crop amount in x-dimension
`crop_y`: crop amount in y-dimension
`crop_z`: crop amount in z-dimension
`theta`: rotation amount in xy-plane
`worm_centroid`: centroid of the worm (to rotate around). NOT the centroids of ROIs in the worm.
`head`: position of the worm's head.
`centroids`: the worm's ROI centroids

## Optional keyword argument
`fill`: what value to put in pixels that were rotated in from outside the original image.
    If kept at its default value "median", the median of the image will be used.
    Otherwise, it can be set to a numerical value.

Outputs a tuple `(new_img, new_head, new_centroids)` corresponding to transformed
versions of `img`, `head`, and `centroids`.
"""
function crop_rotate(img, crop_x, crop_y, crop_z, theta, worm_centroid, head, centroids; fill="median")
    new_img = zeros(Int16, (crop_x[2] - crop_x[1] + 1, crop_y[2] - crop_y[1] + 1, crop_z[2] - crop_z[1] + 1))
    if fill == "median"
        fill_val = Int16(round(median(img)))
    else
        fill_val = fill
    end
    
    imsize = size(img)
    # make sure we aren't trying to crop past image boundary
    cz = (max(crop_z[1], 1), min(crop_z[2], imsize[3]))
    cx = (crop_x[1], crop_x[2])
    cy = (crop_y[1], crop_y[2])
    tfm = recenter(RotMatrix(theta), worm_centroid)
    for z=cz[1]:cz[2]
        new_img_z = warp(img[:,:,z], tfm)
        # make sure we aren't trying to crop past image boundary - use harshest cropping
        cx = (max(cx[1], new_img_z.offsets[1]+1), min(cx[2], new_img_z.offsets[1] + size(new_img_z)[1]))
        cy = (max(cy[1], new_img_z.offsets[2]+1), min(cy[2], new_img_z.offsets[2] + size(new_img_z)[2]))
        for x=cx[1]:cx[2]
            for y=cy[1]:cy[2]
                val = new_img_z[x,y]
                # val is NaN
                if val != val || val == 0
                    val = fill_val
                end
                new_img[x-cx[1]+1,y-cy[1]+1,z-cz[1]+1] = Int16(round(val))
            end
        end
    end
    # find new head and worm_centroid positions
    h = SVector{2}(head)
    new_h = inv(tfm)(h)
    new_head = (Int32(round(new_h[1]))-cx[1]+1, Int32(round(new_h[2]))-cy[1]+1)
    new_centroids = []
    for cent in centroids
        c = SVector{2}(cent[1:2])
        new_c = inv(tfm)(c)
        push!(new_centroids, (round(new_c[1])-cx[1]+1, round(new_c[2])-cy[1]+1, round(cent[3]-cz[1]+1)))
    end
    return (new_img[1:cx[2]-cx[1]+1, 1:cy[2]-cy[1]+1, 1:cz[2]-cz[1]+1], new_head, new_centroids)
end

"""
Crops and rotates an image, and outputs the resulting image, along with its transformed centroids.

# Arguments
- `infile::String`: input MHD file
- `outdir::String`: output directory for transformed MHD file
- `centroid_out::String`: output file for transformed centroid file
`crop_x`: crop amount in x-dimension
`crop_y`: crop amount in y-dimension
`crop_z`: crop amount in z-dimension
`theta`: rotation amount in xy-plane
`worm_centroid`: centroid of the worm (to rotate around). NOT the centroids of ROIs in the worm.
`head`: position of the worm's head.
`centroids`: the worm's ROI centroids

## Optional keyword argument
`fill`: what value to put in pixels that were rotated in from outside the original image.
    If kept at its default value "median", the median of the image will be used.
    Otherwise, it can be set to a numerical value.
"""
function crop_rotate_output(infile::String, outdir::String, centroid_out::String, crop_x, crop_y, crop_z, theta, worm_centroid, head, centroids;
        fill="median", output_mhd=true, output_centroids=true)
    filename = split(split(infile, "/")[end], ".")[1]
    mhd = MHD(infile)
    img = read_img(mhd)
    spacing = split(mhd.mhd_spec_dict["ElementSpacing"], " ")
    new_img, new_head, new_centroids = crop_rotate(img, crop_x, crop_y, crop_z, theta, worm_centroid, head, centroids; fill=fill)
    if output_mhd
        write_raw(outdir*"/"*filename*".raw", new_img)
        write_MHD_spec(outdir*"/"*filename*".mhd", spacing[1], spacing[end], size(new_img)[1],
            size(new_img)[2], size(new_img)[3], filename*".raw")
    end
    if output_centroids
        write_centroids(new_centroids, centroid_out)
    end
    return new_head
end



"""
Finds head location of each worm frame, crops the frame, and outputs the new cropped image,
along with transformed centroids and head locations.
Returns a dictionary of error flags that arose during the computations.

# Arguments
- `rootpath::String`: working directory path; all other directory inputs are relative to this
- `frames`: frames to be cropped
- `MHD_in::String`: directory of the input (non-cropped) MHD files
- `MHD_out::String`: directory of the output cropped MHD files
- `img_prefix::String`: image prefix not including the timestamp. It is assumed that each frame's filename
    will be, eg, `img_prefix_t0123_ch2.mhd` for frame 123 with channel=2.
- `channel::Integer`: channel being used.
- `centroids_in::String`: directory of the input centroid files
- `centroids_out::String`: directory of the transformed centroid files
- `head_file::String`: name of the file for storing head locations. It will be generated within the `centroids_out` directory.

## Optional keyword arguments
- `tf` (default [10,10,30,30]): threshold for required neuron density for convex hull `i` is (number of centroids) / `tf[i]`
- `max_d` (default [30,50,50,100]): the maximum distance for a neuron to be counted as part of convex hull `i` is `max_d[i]`
- `hd_threshold::Integer` (default 100): if convex hulls 2 and 3 give head locations farther apart than this many pixels, set error flag.
- `vc_threshold::Integer` (default 300): if convex hulls 2 and 3 give tail locations farther apart than this many pixels, set error flag.
- `num_centroids_threshold::Integer` (default 90): if there are fewer than this many centroids, set error flag.
- `edge_threshold::Integer` (default 5): if the boundary of the worm is closer than this to the edge of the frame, set error flag.
- `crop_pad` (default [5,5,2]): pad the fourth convex hull by this much before cropping to it.
"""
function crop_rotate_images(rootpath::String, frames, MHD_in::String, MHD_out::String, img_prefix::String, channel::Integer,
        centroids_in::String, centroids_out::String, head_file::String; tf=[10,10,30,30], max_d=[30,50,50,100], hd_threshold::Integer=100, 
        vc_threshold::Integer=300, num_centroids_threshold::Integer=90, edge_threshold::Integer=5, crop_pad=[5,5,2])

    q_flags = Dict()
    create_dir(MHD_out)
    create_dir(centroids_out)

    # get image size
    img = read_img(MHD(joinpath(rootpath, MHD_in, img_prefix*"_t"*string(frames[1], pad=4)*"_ch$(channel).mhd")))
    imsize = size(img)

    n = length(frames)

    open(joinpath(rootpath, centroids_out, head_file), "w") do f
        @showprogress for i in 1:n
            frame = frames[i]
            try
                q_flags[i] = []
                centroids = read_centroids_roi(joinpath(rootpath, centroids_in, "$(i).txt"))
                head, q_flag, crop_x, crop_y, crop_z, theta, worm_centroid = find_head(centroids, imsize;
                        tf=tf, max_d=max_d, hd_threshold=hd_threshold, vc_threshold=vc_threshold,
                        num_centroids_threshold=num_centroids_threshold, edge_threshold=edge_threshold, crop_pad=crop_pad)
                append!(q_flags[i], q_flag)
                new_head = crop_output(joinpath(rootpath, MHD_in, img_prefix*"_t"*string(frame1, pad=4)*"_ch$(channel).mhd"),
                   joinpath(rootpath, MHD_out), joinpath(rootpath, centroids_out, "$(i).txt"), crop_x, crop_y, crop_z, theta, worm_centroid, head, centroids)
                write(f, string(i)*"    "*replace(string(new_head), r"\(|\,|\)" => "")*"\n")
            catch e
                append!(q_flags[i], "ERROR: $(e)")
            end
        end
    end
    return q_flags
end
