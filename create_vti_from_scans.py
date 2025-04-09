import vtk
import cv2
# import imageio
import os
import glob
from vtkmodules.util.numpy_support import numpy_to_vtk  # Import numpy_to_vtk directly
import myVTKPythonLibrary as myvtk
import numpy as np 


def get_metada_PGM(
                    input_file                            :str              ,
                    metadata_fields                       :list             ,
                    header_size                           :int      = 32    ):

    print(input_file)
    metadata = []
    try:
        with open(input_file, "r") as img:
            i = 0
            for line in img:
                for field in metadata_fields:
                    if field in line:
                        for ch in "#@(){}[]!?:":                                                    # remove some specific characters from line
                            line = line.replace(ch, "")
                        line = line.strip()
                        line = line.split()
                        metadata.append(line[-1])
                i += 1
                if i == header_size: break    
        assert len(metadata) == len(metadata_fields), "Unable to read metadata from PGM image"      # Check that all metadata_fields have been found
                                                 # Converts the meta data in floats
    except:
        print("Unable to read metadata from PGM image")
    metadata = [float(i) for i in metadata]            
    return metadata


def pgm2array(
        input_files                                          : str              , 
        file_extension                                       : str      = '.pgm'):
        
    pgm_files               = glob.glob(input_files+"*"+file_extension)                                     
    pgm_files.sort()
    slices                  = [cv2.imread(file, cv2.IMREAD_GRAYSCALE) for file in pgm_files]

    image_array             = np.dstack(slices)
    image_array             = np.rot90(image_array, 1, axes=(1, 0))                                                   # rotation in plane XY for coherence with Collin's code
    image_shape             = slices[0].shape


    return pgm_files, image_array

def array2vti(
        image_shape,
        pixel_size,
        image_pos,
        input_array,
        field_name,
        output):

    output += '.vti'
    vtk_array=numpy_to_vtk(num_array=input_array, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    vtk_array.SetName(field_name)
    # Create image
    image_vti=vtk.vtkImageData()
    image_vti.SetDimensions(image_shape)
    image_vti.SetSpacing(pixel_size)
    image_vti.SetOrigin(image_pos)

    image_vti.GetPointData().SetScalars(vtk_array)
    # Write image
    print("Writing image:", output)
    myvtk.writeImage(image_vti, output)

    return


def get_z_metadata_flatten_image(   files_list          :list       = None      , 
                                    metadata            :list       = None      , 
                                    image_array         :array      = None      ,
                                    get_metadata        :bool       = True      ): 
    
    image_shape             = image_array.shape
    if get_metadata:
        z_abscisse_scd_image    = get_metada_PGM(
                                            input_file             = files_list[1],
                                            metadata_fields   = ["Slice_Location"]
                                        )

        metadata[0]             = "{0:.3f}".format(abs(z_abscisse_scd_image[0]-metadata[0]))
        metadata[0]             = float(metadata[0])


        pixel_size              = [metadata[1], metadata[1], metadata[0]]
        image_pos               = [ps/2 for ps in pixel_size]                                                   # center of the first voxel of the first slice
        return image_shape, pixel_size, image_pos, flatten_image_array
    else:
        flatten_image_array     = np.reshape(image_array, image_shape).flatten(order="F")
        return image_shape, flatten_image_array
    

