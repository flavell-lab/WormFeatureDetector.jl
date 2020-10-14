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
- `crop_pad`: amount to pad the cropping in each dimension. Default 3 pixels in each dimension.
- `min_crop_size`: minimum size of cropped image. Default [210,96,51], the UNet input size.

Outputs a new image that is the cropped and rotated version of `img`.
"""
function crop_rotate(img, crop_x, crop_y, crop_z, theta, worm_centroid; fill="median", degree=Linear(), dtype=Int16, crop_pad=[3,3,3], min_crop_size=[210,96,51])
    new_img = zeros(dtype, (crop_x[2] - crop_x[1] + 1, crop_y[2] - crop_y[1] + 1, crop_z[2] - crop_z[1] + 1))
    if fill == "median"
        fill_val = dtype(round(median(img)))
    else
        fill_val = fill
    end
    
    imsize = size(img)
    # make sure we aren't trying to crop past image boundary after adding crop padding
    cz = (max(crop_z[1]-crop_pad[3], 1), min(crop_z[2]+crop_pad[3], imsize[3]))
    cz = increase_crop_size(cz, imsize[3], min_crop_size[3])

    cx = nothing 
    cy = nothing 
    tfm = recenter(RotMatrix(theta), worm_centroid[1:2])

    for z=cz[1]:cz[2]
        new_img_z = warp(img[:,:,z], tfm, degree)
        # initialize x and y cropping parameters
        if cx == nothing
            cx = (crop_x[1]-crop_pad[1], crop_x[2]+crop_pad[1])
            cx = (max(cx[1], new_img_z.offsets[1]+1), min(cx[2], new_img_z.offsets[1] + size(new_img_z)[1]))
            cx = increase_crop_size(cx, imsize[1], min_crop_size[1])
        end
        
        if cy == nothing
            cy = (crop_y[1]-crop_pad[2], crop_y[2]+crop_pad[2])
            cy = (max(cy[1], new_img_z.offsets[2]+1), min(cy[2], new_img_z.offsets[2] + size(new_img_z)[2]))
            cy = increase_crop_size(cy, imsize[2], min_crop_size[2])
        end

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

""" Increases crop size of a `crop`. Requires the size of the image `img_size`, and the desired minimum crop size `crop_size`. """
function increase_crop_size!(crop, img_size, crop_size)
    if crop_size > img_size
        throw("Crop size cannot be larger than image size")
    end
    while crop[2] - crop[1] < crop_size
        d = argmax([crop[1] - 1, img_size - crop[2]])
        crop[d] += 2*d - 3
    end
    nothing
end

