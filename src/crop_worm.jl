"""
Rotates and then crops an image, optionally along with its head and centroid locations.

# Arguments

- `img`: image to transform (3D)
- `crop_x`: crop amount in x-dimension
- `crop_y`: crop amount in y-dimension
- `crop_z`: crop amount in z-dimension
- `theta`: rotation amount in xy-plane
- `worm_centroid`: centroid of the worm (to rotate around). NOT the centroids of ROIs in the worm.

## Optional keyword argument
- `fill`: what value to put in pixels that were rotated in from outside the original image.
If kept at its default value "median", the median of the image will be used.
Otherwise, it can be set to a numerical value.
- `degree`: degree of the interpolation. Default `Linear()`; can set to `Constant()` for nearest-neighbors.
- `dtype`: type of data in resulting image

Outputs a new image that is the cropped and rotated version of `img`.
"""
function crop_rotate(img, crop_x, crop_y, crop_z, theta, worm_centroid; fill="median", degree=Linear(), dtype=Int16)
    new_img = zeros(dtype, (crop_x[2] - crop_x[1] + 1, crop_y[2] - crop_y[1] + 1, crop_z[2] - crop_z[1] + 1))
    if fill == "median"
        fill_val = dtype(round(median(img)))
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
        new_img_z = warp(img[:,:,z], tfm, degree)
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
                if dtype <: Integer
                    val = round(val)
                end
                new_img[x-cx[1]+1,y-cy[1]+1,z-cz[1]+1] = dtype(val)
            end
        end
    end

    return new_img[1:cx[2]-cx[1]+1, 1:cy[2]-cy[1]+1, 1:cz[2]-cz[1]+1]
end


"""
Generates cropping parameters from a frame by detecting the worm's location with thresholding and noise removal.

# Arguments
- `img`: Image to crop

# Optional keyword arguments
- `threshold`: Number of standard deviations above mean for a pixel to be considered part of the worm
- `size_threshold`: Number of adjacent pixels that must meet the threshold to be counted.
- `crop_pad`: Number of pixels to pad in each dimension.
"""
function get_cropping_parameters(img; threshold::Real=3, size_threshold=10, crop_pad=[3,3,3])
    # threshold image to detect worm
    thresh_img = consolidate_labeled_img(labels_map(fast_scanning(img .> mean(img) + threshold*std(frame), 0.2)), size_threshold);
    # extract worm points
    frame_worm_nonzero = map(x->collect(Tuple(x)), filter(x->frame_worm[x]!=0, CartesianIndices(frame_worm)))
    # get center of worm
    worm_centroid = reduce((x,y)->x.+y, frame_worm_nonzero) ./ length(frame_worm_nonzero)
    # get axis of worm
    deltas = map(x->collect(x.-worm_centroid), frame_worm_nonzero)
    cov_mat = cov(deltas)
    mat_eigvals = eigvals(cov_mat)
    mat_eigvecs = eigvecs(cov_mat)
    eigvals_order = sortperm(mat_eigvals, rev=true)
    # PC1 will be the long dimension of the worm
    # We only want to rotate in xy plane, so project to that plane
    long_axis = mat_eigvecs[:,eigvals_order[1]]
    long_axis[3] = 0
    long_axis = long_axis ./ sqrt(sum(long_axis .^ 2))
    short_axis = [long_axis[2], -long_axis[1], 0]
    theta = atan(long_axis[2]/long_axis[1])

    
    # get coordinates of points in worm axis-dimension
    distances = map(x->sum(x.*long_axis), deltas)
    distances_short = map(x->sum(x.*short_axis), deltas)
    distances_z = map(x->sum(x.*[0,0,1]), deltas)


    # get cropping parameters
    crop_x = (Int64(floor(minimum(distances) + worm_centroid[1])) - crop_pad[1], Int64(ceil(maximum(distances) + worm_centroid[1])) + crop_pad[1])
    crop_y = (Int64(floor(minimum(distances_short) + worm_centroid[2])) - crop_pad[2], Int64(ceil(maximum(distances_short) + worm_centroid[2])) + crop_pad[2])
    crop_z = (Int64(floor(minimum(distances_z) + worm_centroid[3])) - crop_pad[3], Int64(ceil(maximum(distances_z) + worm_centroid[3])) + crop_pad[3])

    return (crop_x, crop_y, crop_z, theta, worm_centroid)
end

