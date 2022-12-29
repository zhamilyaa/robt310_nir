Change the name of the importing file to the proper from the images lists:
In the "samples" file: 
Infrared images are labeled as ..._nir.png

RGB images are labeled as ..._rgb.jpg

There are two source codes named robt310_project.m and nir_image_colorization.m. In robt310_project we implemented dehazing and BDPHE algorithms. In the second we have colorization code, region growing segmentation, median filtering and sharpening algorithms from scratch. 

Code is written in MATLAB, no additional libraries are needed. 

Colorization code runs approximately 5-10 mins, since it needs to segment and map the corresponding pixels form RGB to NIR image. 