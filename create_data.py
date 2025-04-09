import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from create_vti_from_scans import get_metada_PGM, pgm2array, array2vti, get_z_metadata_flatten_image


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
        pgm_files_raw, image_array  = pgm2array(Raw_PGM_base_name)
    else:
        pgm_files, image_array      = pgm2array(Bin_PGM_base_name)
        pgm_files_raw, _            = pgm2array(Raw_PGM_base_name)
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

    img_files, image_array      = pgm2array(images_base_name, ".tiff")

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
