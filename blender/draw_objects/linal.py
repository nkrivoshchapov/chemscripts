import bpy
import numpy as np
from numpy.linalg import norm
from math import sin, cos, acos

from .primitives import cylinder_between, cone_between, draw_point_xyz
from ..math import DEG2RAD


def draw_arrow(start, end, material=None, len_arrow = 0.4, rad_cyl = 0.07, cone_size = 1.4):
    rad_cone = rad_cyl * cone_size

    rad_vec = end.astype('float') - start.astype('float')
    rad_vec /= np.linalg.norm(rad_vec)
    cylinder_between(start, end - rad_vec * (len_arrow + 0.002), rad_cyl, material=material)
    cone_between(end - rad_vec*len_arrow, end, rad_cone, material=material)


def draw_curvedarrow(rA, rB, position=np.array([0, 0, 0]), leftgap=0, rightgap=0, material=None):
    leftgap *= DEG2RAD
    rightgap *= DEG2RAD
    print(repr(norm(rA)))
    beveldepth = 0.1
    points_per_deg = 10
    mynorm = norm(rA)
    z_axis = np.cross(rA, rB)
    z_axis /= norm(z_axis)
    x_axis = rA
    x_axis /= norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    basis_change = np.array([x_axis, y_axis, z_axis]).transpose()
    angle = acos(np.dot(rA / norm(rA), rB / norm(rB)))
    print("Angle = " + repr(angle))

    Npoints = int((angle - leftgap - rightgap) / DEG2RAD * points_per_deg)
    angle_step = 1 / points_per_deg

    # create the Curve Datablock
    curveData = bpy.data.curves.new('myCurve', type='CURVE')
    curveData.dimensions = '3D'
    curveData.resolution_u = 2

    # map coords to spline
    polyline = curveData.splines.new('POLY')
    polyline.points.add(Npoints - 1)
    curangle = leftgap
    sphere_pos = (basis_change @ np.array([cos(curangle), sin(curangle), 0])) * mynorm + position

    draw_point_xyz(sphere_pos, size=beveldepth, material=material)
    for i in range(Npoints):
        rot_vector = np.array([cos(curangle), sin(curangle), 0])
        myvec = (basis_change @ rot_vector) * mynorm
        polyline.points[i].co = (myvec[0], myvec[1], myvec[2], 1)
        curangle += angle_step * DEG2RAD

    conestart = (basis_change @ np.array([cos(curangle), sin(curangle), 0])) * mynorm + position
    conedirection = basis_change @ np.array([-sin(curangle), cos(curangle), 0])
    coneend = conestart + conedirection * 0.8
    cone_between(conestart, coneend, beveldepth * 2.5, material=material)

    # create Object
    curveOB = bpy.data.objects.new('myCurve', curveData)
    curveData.bevel_depth = beveldepth

    # attach to scene and validate context
    scn = bpy.context.scene
    scn.collection.objects.link(curveOB)
    curveOB.location = position
    if material is not None:
        curveOB.data.materials.append(bpy.data.materials[material])
