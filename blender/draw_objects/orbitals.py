import bpy
import numpy as np
import os, ntpath, glob

from ..materials import get_orbital_material, get_wavy_material


def nbo_from_mcubes(name_tempate, material): # name_tempate = trino_48_plus.{type}iso
    objname = ntpath.basename(name_tempate).split('.')[0]
    vertices = list(np.load(name_tempate.format(type='vertices')))
    edges = list(np.load(name_tempate.format(type='edges')))
    triangles = list(np.load(name_tempate.format(type='triangles')))

    new_mesh = bpy.data.meshes.new(objname + "_mesh")
    new_mesh.from_pydata(vertices, edges, triangles)
    new_mesh.update()
    new_object = bpy.data.objects.new(objname, new_mesh)
    bpy.context.collection.objects.link(new_object)
    bpy.context.view_layer.objects.active = new_object
    bpy.ops.object.modifier_add(type='SMOOTH')
    bpy.data.objects[objname].modifiers["Smooth"].iterations = 1
    bpy.data.objects[objname].modifiers["Smooth"].factor = 1
    bpy.ops.object.shade_smooth()
    for p in new_object.data.polygons:
        p.use_smooth = True
    bpy.data.objects[objname].data.materials.append(bpy.data.materials[material])
    return objname


def get_mcubes_templates(nboname, nbodir):
    surf_types = []
    for file in glob.glob(os.path.join(nbodir, nboname + '*')):
        curtype = os.path.join(nbodir, ntpath.basename(file).split('.')[0]) + ".{type}iso"
        if curtype not in surf_types:
            surf_types.append(curtype)
    return surf_types


def plot_nbo(nboname, color="#377eb8", reverse=False, nbodir="./calcfiles"):
    files = get_mcubes_templates(nboname, nbodir)
    if len(files) > 2 or len(files) < 1:
        raise Exception("Unexpected number of %s_*.wrl files" % nboname)

    if isinstance(color, str):
        plusmat = get_orbital_material(color, 1)
        minusmat = get_orbital_material(color, -1)
    elif isinstance(color, list) and len(color) == 2:
        plusmat = get_wavy_material(color, 1)
        minusmat = get_wavy_material(color, -1)
    else:
        raise Exception("Unsupported color datatype")

    if len(files) == 1:
        if not reverse:
            a = nbo_from_mcubes(files[0], material=plusmat)
        else:
            a = nbo_from_mcubes(files[0], material=minusmat)
        return a
    elif len(files) == 2:
        if not reverse:
            a, b = nbo_from_mcubes(files[0], material=plusmat), nbo_from_mcubes(files[1], material=minusmat)
        else:
            a, b = nbo_from_mcubes(files[1], material=plusmat), nbo_from_mcubes(files[0], material=minusmat)
        return a, b


""" # Ancient code for drawing from Jmol files
def plot_wrl(file, material="orbital_template_plus", nbodir="../nbofiles/"):
    io_scene_x3d.import_x3d.load(bpy.context, nbodir + file)
    bpy.data.objects.remove(bpy.data.objects["Viewpoint"], do_unlink=True)
    surface_name = ntpath.basename(file).replace("_", "").replace(".wrl", "")
    bpy.data.objects["Shape_IndexedFaceSet"].name = surface_name
    bpy.context.view_layer.objects.active  = bpy.data.objects[surface_name]
    bpy.ops.object.modifier_add(type='SMOOTH')
    bpy.data.objects[surface_name].modifiers["Smooth"].iterations = 2
    bpy.data.objects[surface_name].modifiers["Smooth"].factor = 2
    bpy.ops.object.shade_smooth()
    bpy.data.objects[surface_name].data.materials.append(bpy.data.materials[material])
    return surface_name
def plot_nbo(nboname, color="#377eb8", reverse=False, nbodir="../nbofiles/"):
    files = glob.glob(nbodir + "%s_*.wrl" % nboname)
    if len(files) > 2 or len(files) < 1:
        raise Exception("Unexpected number of %s_*.wrl files" % nboname)
    group_name = "".join(ntpath.basename(files[0]).replace(".wrl", "").split("_")[:2])
    plusmat = get_material(REDCOLOR, 1)
    minusmat = get_material(BLUECOLOR, -1)
    if len(files) == 1:
        plot_wrl(files[0], material=plusmat)
        # bpy.data.collections[group_name].objects.link(bpy.data.objects[])
    elif len(files) == 2:
        if not reverse:
            a, b = plot_wrl(files[0], material=plusmat), plot_wrl(files[1], material=minusmat)
        else:
            a, b = plot_wrl(files[1], material=plusmat), plot_wrl(files[0], material=minusmat)
"""
