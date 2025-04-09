import vtk
import nibabel as nib
import numpy as np
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import myVTKPythonLibrary as myvtk  # Replace with your library for writing VTI files


def nii2array(input_file):
    """
    Load a NIfTI file and return the image data as a NumPy array and the pixel dimensions.
    """
    # Load the nifti file
    nifti_image = nib.load(input_file)
    
    # Get the image data as a numpy array
    image_array = nifti_image.get_fdata()


    image_array = np.array([np.rot90(slice, 1, axes=(1, 0)) for slice in image_array])


    image_array = np.moveaxis(image_array, 0, 2)


    unique_segments = np.unique(image_array)
    print("Unique segment IDs:", unique_segments)

    # Check the nifti header for metadata
    header = nifti_image.header
    print("Header information:", header)


    # Get the pixel sizes from the header
    pixel_size = nifti_image.header.get_zooms()
    
    return image_array, pixel_size


def array2vti(image_shape, pixel_size, image_pos, input_array, field_name, output_file):
    """
    Convert a NumPy array to a VTK ImageData format and save as a VTI file.
    """
    # Flatten the NumPy array
    flattened_array = np.reshape(input_array, image_shape).flatten(order="F")
    
    # Convert NumPy array to VTK array
    vtk_array = numpy_to_vtk(num_array=flattened_array, deep=True, array_type=vtk.VTK_FLOAT)
    vtk_array.SetName(field_name)
    
    # Create VTK ImageData
    image_vti = vtk.vtkImageData()
    image_vti.SetDimensions(image_shape)
    image_vti.SetSpacing(pixel_size)
    image_vti.SetOrigin(image_pos)
    image_vti.GetPointData().SetScalars(vtk_array)
    
    # Write VTI file
    output_file += '.vti'
    print(f"Writing image to {output_file}")
    myvtk.writeImage(image_vti, output_file)



def nii2vti(input_name, output_name):
    import os.path
    if os.path.isfile(output_name+".vti"):
        print(f'Bin VTI file {output_name+".vti"} already exists')
        return 
    # Load NIfTI image data
    image_array, pixel_size = nii2array(input_name+".nii")

    pixel_size = [0.7*p for p in pixel_size]

    # Define metadata for VTI
    image_shape = image_array.shape
    image_pos = [ps / 2 for ps in pixel_size]  # Center of the first voxel

    # Convert and save to VTI
    array2vti(
        image_shape=image_shape,
        pixel_size=pixel_size,
        image_pos=image_pos,
        input_array=image_array,
        field_name="pixel intensity",
        output_file=output_name
    )




def coloured2bin(
    input_name          : str   = None              , 
    output_name         : str   = None              , 
    field_name          : str   = "pixel intensity" ,
    RL_value            : int   = 100               ,
    LL_value            : int   = 200               ,                                           
    integer_separation  : int   = 185               ,                                           
    trachea_value       : int   = 80                ,                                           
    ):                                         


    import os.path
    if os.path.isfile(output_name+".vti"):
        print(f'Bin VTI file {output_name+".vti"} already exists')
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

    scalar_array[(scalar_array == trachea_value) ]   = 0
    scalar_array[(scalar_array <= integer_separation) & (scalar_array > trachea_value)]   = RL_value
    scalar_array[(scalar_array >= integer_separation)]   = LL_value


    modified_scalars = numpy_to_vtk(scalar_array)
    modified_scalars.SetName(field_name)  # Set the name of the scalar field
    image_data.GetPointData().SetScalars(modified_scalars)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_name+".vti")
    writer.SetInputData(image_data)
    writer.Write()

    print("Done bin masking. "+input_name)




# # Input and Output Paths

import glob

input_files = glob.glob("/Users/daby/LargeFiles/masks_from_catalyn/nii/PAT*")
input_files = [file[:-4] for file in input_files]
output_files =  ["/Users/daby/LargeFiles/masks_from_catalyn/vti/"+file[-5:] for file in input_files]


for i in range(len(input_files)):
    nii2vti(input_files[i], output_files[i])


Alice_patients  = 9                                                                                     # Number of Alice patients
i = Alice_patients                                                                                                   

for file in output_files:
    i+=1
    # patient_number = file[-2:]
    output_name = "/Users/daby/LargeFiles/masks_from_catalyn/BIN/Pat"+str(i)+"_BIN"
    coloured2bin(
    input_name          = file              , 
    output_name         = output_name       , 
    field_name          = "pixel intensity" ,
    RL_value            = 100               ,
    LL_value            = 200               ,                                           
    integer_separation  = 185               ,                                           
    )


catalyn_patients_ID = list(range(Alice_patients+1, len(input_files)+Alice_patients+1))

for patient in catalyn_patients_ID:
    "./Images/PA"+str(patient)+"/BIN/Pat"+str(patient)+"_BIN",
    destination = "./Images/PA"+str(patient)+"/BIN"
    source = "/Users/daby/LargeFiles/masks_from_catalyn/BIN/Pat"+str(patient)+"_BIN"+".vti"
    if not os.path.isdir(destination):
        os.system("mkdir -p "+destination)
        print("folder "+destination+" created")
    rsync_cmd = "rsync -azvp  "+source+" "+destination
    import os
    os.system(rsync_cmd)
