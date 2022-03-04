import bpy
from .utils import hex_to_rgb


# MATERIALS FOR ORBITALS
def get_orbital_material(color, sign):
    if sign == 1:
        matname = color.replace("#", "") + "_orbital_plus"
    elif sign == -1:
        matname = color.replace("#", "") + "_orbital_minus"
    else:
        raise Exception("Forbidden sign value")
    exists = False
    for item in bpy.data.materials:
        if item.name == matname:
            exists = True
            break

    if not exists:
        if sign == 1:
            template_mat = bpy.data.materials["orbital_template_plus"]
        elif sign == -1:
            template_mat = bpy.data.materials["orbital_template_minus"]
        else:
            raise Exception("Forbidden sign value")
        copied_mat = template_mat.copy()
        copied_mat.name = matname
        copied_mat.node_tree.nodes["RGB"].outputs[0].default_value = hex_to_rgb(int(color.replace('#', '0x'), 16))
    return matname


def get_wavy_material(color, sign):
    if sign == 1:
        matname = "{firstcolor}_{secondcolor}_mixorbital_plus".format(
                                    firstcolor=color[0].replace("#", ""),
                                    secondcolor=color[1].replace("#", ""))
    elif sign == -1:
        matname = "{firstcolor}_{secondcolor}_mixorbital_minus".format(
                                    firstcolor=color[0].replace("#", ""),
                                    secondcolor=color[1].replace("#", ""))
    else:
        raise Exception("Forbidden sign value")
    exists = False
    for item in bpy.data.materials:
        if item.name == matname:
            exists = True
            break

    if not exists:
        if sign == 1:
            template_mat = bpy.data.materials["mixorbital_template_plus"]
        elif sign == -1:
            template_mat = bpy.data.materials["mixorbital_template_minus"]
        else:
            raise Exception("Forbidden sign value")
        copied_mat = template_mat.copy()
        copied_mat.name = matname
        copied_mat.node_tree.nodes["RGB"].outputs[0].default_value = hex_to_rgb(int(color[0].replace('#', '0x'), 16))
        copied_mat.node_tree.nodes["RGB.001"].outputs[0].default_value = hex_to_rgb(int(color[1].replace('#', '0x'), 16))
    return matname


# OTHER MATERIALS FOR NBO PLOTS
def get_arrow_material(color):
    matname = color.replace("#", "") + "_arrow"
    exists = False
    for item in bpy.data.materials:
        if item.name == matname:
            exists = True
            break

    if not exists:
        template_mat = bpy.data.materials["arrow_template"]
        copied_mat = template_mat.copy()
        copied_mat.name = matname
        copied_mat.node_tree.nodes["RGB"].outputs[0].default_value = hex_to_rgb(int(color.replace('#', '0x'), 16))
    return matname
