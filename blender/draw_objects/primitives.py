import bpy
import math
import numpy as np


def draw_point_xyz(point, size = 0.05, material=None):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=size, location=point)
    bpy.ops.object.shade_smooth()
    if material is not None:
        bpy.context.object.data.materials.append(bpy.data.materials[material])


def cylinder_between(start, end, r, material=None):
    rad_vec = end-start
    dist = np.linalg.norm(rad_vec)

    bpy.ops.mesh.primitive_cylinder_add(
        radius = r,
        vertices=258,
        depth = dist,
        location = rad_vec/2 + start
    )

    phi = math.atan2(rad_vec[1], rad_vec[0])
    theta = math.acos(rad_vec[2] / dist)

    bpy.context.object.rotation_euler[1] = theta
    bpy.context.object.rotation_euler[2] = phi
    bpy.ops.object.shade_smooth()
    if material is not None:
        bpy.context.object.data.materials.append(bpy.data.materials[material])


def cone_between(start, end, r, material=None):
    rad_vec = end-start
    dist = np.linalg.norm(rad_vec)

    bpy.ops.mesh.primitive_cone_add(
        radius1 = r,
        vertices=258,
        depth = dist,
        location = rad_vec/2 + start
    )

    phi = math.atan2(rad_vec[1], rad_vec[0])
    theta = math.acos(rad_vec[2]/dist)

    bpy.context.object.rotation_euler[1] = theta
    bpy.context.object.rotation_euler[2] = phi
    # bpy.ops.object.shade_smooth()
    if material is not None:
        bpy.context.object.data.materials.append(bpy.data.materials[material])
