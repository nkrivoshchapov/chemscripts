import bpy
import fnmatch

# reserved_materials = ["orbital_template_plus", "orbital_template_minus", "mixorbital_template_plus", "mixorbital_template_minus", "arrow_template", "plane", "Enegative_arrow"]


def srgb_to_linearrgb(c):
    if   c < 0:       return 0
    elif c < 0.04045: return c/12.92
    else:             return ((c+0.055)/1.055)**2.4


def hex_to_rgb(h,alpha=1):
    r = (h & 0xff0000) >> 16
    g = (h & 0x00ff00) >> 8
    b = (h & 0x0000ff)
    return tuple([srgb_to_linearrgb(c/0xff) for c in (r,g,b)] + [alpha])


def cleanup(protected_items=("Camera*", "Light*", "Area"), protected_materials=("orbital_template_plus", "orbital_template_minus", "mixorbital_template_plus", "mixorbital_template_minus", "arrow_template", "plane", "Enegative_arrow")):
    bpy.context.scene.cursor.location = (0.0, 0.0, 0.0)

    for material in bpy.data.materials:
        if material.name not in protected_materials:
            bpy.data.materials.remove(material)

    for obj in bpy.data.objects:
        reserved = False
        for item in protected_items:
            if fnmatch.fnmatch(obj.name, item):
                reserved = True
                break
        if not reserved:
            bpy.data.objects.remove(obj, do_unlink=True)
    for item in bpy.data.collections:
        if item.name != "Collection":
            bpy.data.collections.remove(item)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.ops.object.select_by_type(type="MESH")
    bpy.ops.object.delete()
    bpy.ops.object.select_by_type(type="CURVE")
    bpy.ops.object.delete()
