#------------visualize_3d.py---------------------------------------------------#
#
# Purpose: This file will use blender to visualize the 3d data from GPUE and 
#          will output an image (or video, to be added later)
#
#   Notes: this function will rely on the gen_data.py function to create a 
#          .bvox file for blender to use.
#
#          To create an image with this, please use:
#              blender -b -P visualize_3d.py
#
#------------------------------------------------------------------------------#

# I think we'll probably need all the f(x)'s... 
# This requires toms trickery because we are using blender for this
# from gen_data.py import *

# library for blender python
import bpy

# other libraries
import numpy as np

# files for visualization
voxelfile = "test_wfc.bvox"
infile = open(voxelfile, "r")

#------------------------------------------------------------------------------#
# FUNCTIONS
#------------------------------------------------------------------------------#

# function to remove all (unnecessary) objects from the scene
def remove_obj(scene):
    for ob in scene.objects:
        if ob.name != "Camera":
            scene.objects.unlink(ob)

# Function to define the scene
def def_scene(box_length, res_stand, xres, yres, zres):
    # first, we need to relocate the camera
    x_cam = 3.0
    y_cam = 3.2
    z_cam = 2.0

    scene = bpy.context.scene

    # removing unnecessary obects (basically leaving only the camera)
    remove_obj(scene)

    # defining dummy location (empty) to point camera at.
    bpy.ops.object.add(type="EMPTY",
                       location=(box_length * xres * 0.5 / res_stand,
                                 box_length * yres * 0.5 / res_stand,
                                 box_length * zres * 0.5 / res_stand))

    context = bpy.context

    bpy.ops.object.select_pattern(pattern="Camera")
    bpy.context.scene.objects.active = bpy.context.scene.objects["Camera"]
    ob = bpy.data.objects["Camera"]
    bpy.ops.object.constraint_add(type = "TRACK_TO")
    target = bpy.data.objects.get("Empty", False)
    ob.constraints["Track To"].target=target

    ob.constraints["Track To"].track_axis = "TRACK_NEGATIVE_Z"
    ob.constraints["Track To"].up_axis = "UP_Y"

    scene.camera.location.x = box_length * x_cam * xres / res_stand
    scene.camera.location.y = box_length * y_cam * yres / res_stand
    scene.camera.location.z = box_length * z_cam * zres / res_stand

    # set's FOV
    scene.camera.data.angle = 50*(np.pi/180.0)
    bpy.data.cameras["Camera"].ortho_scale = 21.0

    # Set number of cores used
    scene.render.threads_mode = "FIXED"
    scene.render.threads = 8

    # Sets the BG to be black
    bpy.data.worlds["World"].horizon_color = (0,0,0)

    return scene

# function to create cube for data
def create_cube(box_length, res_stand, xres, yres, zres, step_size, 
                dens_scale, voxelfile, color_num):
    cube = bpy.ops.mesh.primitive_cube_add(
               location=((box_length * xres / (2*res_stand)),
                         (box_length * yres / (2*res_stand)),
                         (box_length * zres / (2*res_stand))),
               radius = box_length * 0.5)

    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.object.mode_set(mode="OBJECT")
    ob = bpy.context.object
    ob.scale = ((xres/res_stand, yres/res_stand, zres/res_stand))

    # setting voxel material
    me = ob.data
    mat = create_volume("MaterialVolume", xres, yres, zres, step_size,
                        dens_scale, voxelfile, color_num)
    me.materials.append(mat)

    return cube

# function to create voxel material for cube
def create_volume(passedName, xres, yres, zres, step_size, dens_scale, 
                  voxelfile, color_num):
    volMat = bpy.data.materials.new(passedName)
    volMat.type = "VOLUME"
    volMat.volume.density = 0.0
    volMat.volume.step_method = "CONSTANT"
    volMat.volume.depth_threshold = 0.01
    volMat.volume.density_scale = dens_scale
    matTex = volMat.texture_slots.add()
    voxTex = bpy.data.textures.new("VoxelData", type = "VOXEL_DATA")
    voxTex.voxel_data.file_format = "BLENDER_VOXEL"
    voxTex.use_color_ramp = True
    voxTex.color_ramp.color_mode = "RGB"
    ramp = voxTex.color_ramp

    values = [(0.0,(0,0,1,0)), (0.5,(1,0,1,0.3)), (1.0, (1,0,0,1))]

    for n,value in enumerate(values):
        ramp.elements.new((n+1)*0.2)
        elt = ramp.elements[n]
        (pos, color) = value
        elt.position = pos
        elt.color = color
    voxTex.voxel_data.filepath = voxelfile
    voxTex.voxel_data.resolution = (xres, yres, zres)
    matTex.texture = voxTex
    matTex.use_map_to_bounds = True
    matTex.texture_coords = 'ORCO'
    matTex.use_map_color_diffuse = True 
    matTex.use_map_emission = True 
    matTex.emission_factor = 1
    matTex.emission_color_factor = 1
    matTex.use_map_density = True 
    matTex.density_factor = 1
    return volMat

def createCage(passedName):
    cageMat = bpy.data.materials.new(passedName)
    cageMat.use_shadeless = True
    return cageMat

# Render Scene into image
def render_img(filename):
    bpy.data.scenes['Scene'].render.filepath = filename
    bpy.ops.render.render( write_still=True )

#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#

scene = def_scene(5, 64, 64, 64, 64)
create_cube(5, 64, 64, 64,64, 0.1, 0.5, voxelfile, 1)
render_img("image.png")
