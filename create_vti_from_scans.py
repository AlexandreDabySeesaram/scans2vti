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


def sequence2array(
        input_files                                          : str              , 
        file_extension                                       : str      = '.pgm'):

    img_files               = glob.glob(input_files+"*"+file_extension)                                     
    img_files.sort()
    slices                  = [cv2.imread(file, cv2.IMREAD_GRAYSCALE) for file in img_files]

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
                                    image_array         :np.array      = None      ,
                                    get_metadata        :bool       = True      ): 
    
    image_shape             = image_array.shape
    flatten_image_array     = np.reshape(image_array, image_shape).flatten(order="F")
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
        return image_shape, flatten_image_array
    



def sign_masking_binary(
    input_name          : str   = None              , 
    suffix              : str   = "signed"          , 
    field_name          : str   = "pixel intensity" ,
    scalar2zero         : int   = 200               ,
    scalar_background   : int   = 0                 ,                                           # Initial background pixel intensity
    scalar_foreground   : int   = 100               ,                                           # Initial foreground pixel intensity
    target_value_bg     : float = 50                ,                                           # target background pixel intensity
    target_value_fg     : float = -50               ,                                           # target foreground pixel intensity
    target_type         : str   = "signed_char"     ,                                           # unsigned_char, signed_char, float
    ):                                         


    import os.path
    if os.path.isfile(input_name+"_"+suffix+".vti"):
        print(f'Masked VTI file {input_name+"_"+suffix+".vti"} already exists')
        return 
    input_file      = input_name+".vti"

    # Load the VTI file
    reader          = vtk.vtkXMLImageDataReader()
    reader.SetFileName(input_file)
    reader.Update()

    # Get the image data
    image_data      = reader.GetOutput()

    # Extract scalar data and convert to Numpy array
    scalars         = image_data.GetPointData().GetScalars()
    scalar_array    = vtk_to_numpy(scalars)

    scalar_array[(scalar_array == scalar2zero)]         = 0
    scalar_array[(scalar_array == scalar_background)]   = target_value_bg
    scalar_array[(scalar_array == scalar_foreground)]   = target_value_fg
    match target_type:
        case "unsigned_char":
            scalar_array = scalar_array.astype(np.uint8)
        case "signed_char":
            scalar_array = scalar_array.astype(np.int8)
        case "float":
            scalar_array = scalar_array.astype(np.float64)

    modified_scalars = numpy_to_vtk(scalar_array)
    modified_scalars.SetName(field_name)  # Set the name of the scalar field
    image_data.GetPointData().SetScalars(modified_scalars)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(input_name+"_"+suffix+".vti")
    writer.SetInputData(image_data)
    writer.Write()

    print("Done sign masking. "+input_name)


def PGM2vti(
            output_name         : int   = None,
            Raw_PGM_base_name   : int   = None,
            Bin_PGM_base_name   : int   = None,):
    
    import os.path
    if os.path.isfile(output_name+".vti"):
        print(f'VTI file {output_name+".vti"} already exists')
        return 

    if Bin_PGM_base_name == None:
        pgm_files_raw, image_array  = sequence2array(Raw_PGM_base_name)
    else:
        pgm_files, image_array      = sequence2array(Bin_PGM_base_name)
        pgm_files_raw, _            = sequence2array(Raw_PGM_base_name)
    metadata_fields             = ["Slice_Location", "Pixel_Size"]
    metadata                    = get_metada_PGM(
                                    input_file        = pgm_files_raw[0],
                                    metadata_fields   = metadata_fields
                                    )
    image_shape, pixel_size, image_pos, flatten_image_array = get_z_metadata_flatten_image(
                                                                    files_list  = pgm_files_raw, 
                                                                    metadata    = metadata, 
                                                                    image_array = image_array)
    array2vti(
            image_shape         = image_shape,
            pixel_size          = pixel_size,
            image_pos           = image_pos,
            input_array         = flatten_image_array,
            field_name          = 'pixel intensity',
            output              = output_name)

def TIFF2vti(
            output_name         : int   = None,
            images_base_name    : int   = None,
            pixel_size                  = 1):

    if not isinstance(pixel_size, list):
        pixel_size = [pixel_size]*3
    
    import os.path
    if os.path.isfile(output_name+".vti"):
        print(f'VTI file {output_name+".vti"} already exists')
        return 

    img_files, image_array      = sequence2array(images_base_name, ".tif")

    image_shape, flatten_image_array = get_z_metadata_flatten_image(
                                            image_array     = image_array, 
                                            get_metadata    = False)
    
    array2vti(
            image_shape         = image_shape,
            pixel_size          = pixel_size,
            image_pos           = [ps/2 for ps in pixel_size],
            input_array         = flatten_image_array,
            field_name          = 'pixel intensity',
            output              = output_name)




def copy_folder(source, destination):
    import os

    if not os.path.isdir(destination):
        os.system("mkdir -p "+destination)
        print("folder "+destination+" created")

    assert os.path.isdir(source), f"source folder {source} not found. Aborting"

    cp_comand = "rsync -azvp  "+source+"/ "+destination
    os.system(cp_comand)

def sequence2vti(
            output_name         : str   = None,
            Raw_base_name       : str   = None,
            Bin_base_name       : str   = None,
            image_ext                   = ".pgm"                                    ,
            pixel_size                  = 1, 
            ):

    match image_ext:
        case ".pgm":
            PGM2vti(
                output_name         = output_name       ,
                Raw_PGM_base_name   = Raw_base_name     )
        case ".tif":
            TIFF2vti(
                output_name         = output_name       ,
                images_base_name    = Raw_base_name     ,
                pixel_size          = pixel_size        )