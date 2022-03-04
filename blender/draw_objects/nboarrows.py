import bpy
from copy import copy
import numpy as np
from numpy.linalg import norm, inv

from ..math import rotation_matrix_2axis, rotation_matrix_3axis
from ..materials import get_arrow_material
from .primitives import cone_between
from .linal import draw_arrow


def draw_aearrow(m, ael, start_x_shift=0, start_z_shift=0, start_ang_z=0, start_ang_y=0, end_shift=0, curv=1,
                 offplane_ang=0, offplane_st=0, offplane_acc=0, color=None):
    don_p = copy(m.atoms[ael[1]].location)
    donA = m.atoms[ael[0]].location - m.atoms[ael[1]].location
    donB = m.atoms[ael[2]].location - m.atoms[ael[1]].location
    don_dir = np.cross(donA, donB)
    don_dir /= norm(don_dir)
    acc_p = copy(m.atoms[ael[2]].location)
    acc_dir = m.atoms[ael[2]].location - m.atoms[ael[3]].location
    acc_dir /= norm(acc_dir)
    if np.dot(don_dir, acc_dir) < 0:
        don_dir[:] = -don_dir[:]

    if color is not None:
        newmat = get_arrow_material(color)
    else:
        newmat = None

    # draw_arrow(don_p, don_p + don_dir)
    # draw_arrow(acc_p, acc_p + acc_dir)

    startpos = don_p + start_x_shift * donB / norm(donB) + start_z_shift * don_dir
    av = copy(don_dir)
    bv = copy(donB) / norm(donB)
    cv = np.cross(av, bv)
    localcoord = np.array([av, bv, cv]).T
    startdir = localcoord @ rotation_matrix_2axis(start_ang_y) @ rotation_matrix_3axis(start_ang_z) @ np.array(
        [1, 0, 0])
    endpos = acc_p + end_shift * acc_dir

    startpos = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ (startpos - don_p) + don_p
    endpos = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ (endpos - don_p) + don_p
    startdir = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ startdir
    acc_dir = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ acc_dir

    startpos = localcoord @ rotation_matrix_2axis(offplane_st) @ inv(localcoord) @ (startpos - don_p) + don_p
    endpos = localcoord @ rotation_matrix_2axis(offplane_acc) @ inv(localcoord) @ (endpos - don_p) + don_p
    startdir = localcoord @ rotation_matrix_2axis(offplane_st) @ inv(localcoord) @ startdir
    acc_dir = localcoord @ rotation_matrix_2axis(offplane_acc) @ inv(localcoord) @ acc_dir

    bpy.ops.curve.primitive_bezier_curve_add(enter_editmode=False)
    curve = bpy.context.active_object
    curve.name = 'AE Arrow'

    bd = 0.1
    curve.data.use_fill_caps = True
    curve.data.bevel_mode = 'ROUND'
    curve.data.bevel_depth = bd

    bez_points = curve.data.splines[0].bezier_points
    for bez_point in bez_points:
        bez_point.handle_left_type = 'FREE'
        bez_point.handle_right_type = 'FREE'

    bez_points[0].co = startpos
    bez_points[0].handle_left = startpos + startdir * curv
    bez_points[0].handle_right = startpos + startdir * curv

    bez_points[1].co = endpos
    bez_points[1].handle_left = endpos + acc_dir * curv
    bez_points[1].handle_right = endpos + acc_dir * curv

    curve.data.materials.append(bpy.data.materials[newmat])
    cone_between(endpos - acc_dir * 0.005, endpos - acc_dir * 0.5, 2.0 * bd, material=newmat)

def draw_aearrow_pp(m, ael, start_x_shift=0, start_z_shift=0, start_ang_z=0, start_ang_y=0, end_shift=0, curv=1,
                 offplane_ang=0, offplane_st=0, offplane_acc=0, color=None):
    don_p = copy(m.atoms[ael[1]].location)
    donA = m.atoms[ael[0]].location - m.atoms[ael[1]].location
    donB = m.atoms[ael[2]].location - m.atoms[ael[1]].location
    don_dir = np.cross(donA, donB)
    don_dir /= norm(don_dir)

    acc_p = copy(m.atoms[ael[2]].location)
    accA = m.atoms[ael[1]].location - m.atoms[ael[2]].location
    accB = m.atoms[ael[3]].location - m.atoms[ael[2]].location
    acc_dir = np.cross(accA, accB)
    acc_dir /= norm(acc_dir)
    if np.dot(don_dir, acc_dir) < 0:
        acc_dir[:] = -acc_dir[:]

    if color is not None:
        newmat = get_arrow_material(color)
    else:
        newmat = None

    # draw_arrow(don_p, don_p + don_dir)
    # draw_arrow(acc_p, acc_p + acc_dir)

    startpos = don_p + start_x_shift * donB / norm(donB) + start_z_shift * don_dir
    av = copy(don_dir)
    bv = copy(donB) / norm(donB)
    cv = np.cross(av, bv)
    localcoord = np.array([av, bv, cv]).T
    startdir = localcoord @ rotation_matrix_2axis(start_ang_y) @ rotation_matrix_3axis(start_ang_z) @ np.array(
        [1, 0, 0])
    endpos = acc_p + end_shift * acc_dir

    startpos = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ (startpos - don_p) + don_p
    endpos = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ (endpos - don_p) + don_p
    startdir = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ startdir
    acc_dir = localcoord @ rotation_matrix_2axis(offplane_ang) @ inv(localcoord) @ acc_dir

    startpos = localcoord @ rotation_matrix_2axis(offplane_st) @ inv(localcoord) @ (startpos - don_p) + don_p
    endpos = localcoord @ rotation_matrix_2axis(offplane_acc) @ inv(localcoord) @ (endpos - don_p) + don_p
    startdir = localcoord @ rotation_matrix_2axis(offplane_st) @ inv(localcoord) @ startdir
    acc_dir = localcoord @ rotation_matrix_2axis(offplane_acc) @ rotation_matrix_3axis(-20) @ inv(localcoord) @ acc_dir

    bpy.ops.curve.primitive_bezier_curve_add(enter_editmode=False)
    curve = bpy.context.active_object
    curve.name = 'AE Arrow'

    endpos -= bv * 0.25

    bd = 0.1
    curve.data.use_fill_caps = True
    curve.data.bevel_mode = 'ROUND'
    curve.data.bevel_depth = bd

    bez_points = curve.data.splines[0].bezier_points
    for bez_point in bez_points:
        bez_point.handle_left_type = 'FREE'
        bez_point.handle_right_type = 'FREE'

    bez_points[0].co = startpos
    bez_points[0].handle_left = startpos + startdir * curv
    bez_points[0].handle_right = startpos + startdir * curv

    bez_points[1].co = endpos
    bez_points[1].handle_left = endpos + acc_dir * curv
    bez_points[1].handle_right = endpos + acc_dir * curv

    curve.data.materials.append(bpy.data.materials[newmat])
    cone_between(endpos - acc_dir * 0.005, endpos - acc_dir * 0.5, 2.0 * bd, material=newmat)
