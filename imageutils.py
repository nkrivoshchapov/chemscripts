import pandas as pd
import fnmatch, glob
import os
from PIL import Image, ImageEnhance, ImageChops
import imageio
import numpy as np
from plotnine import ggplot, geom_point, aes, geom_line, theme_bw, scale_color_manual, ggtitle


class TrimmingBox:
    def __init__(self):
        self.points = [None, None, None, None]

    @staticmethod
    def _get_trimming_box(im):
        bg = Image.new(im.mode, im.size, im.getpixel((0, 0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        return bbox

    @property
    def height(self):
        return self.points[3] - self.points[1]

    @property
    def width(self):
        return self.points[2] - self.points[0]

    @property
    def size(self):
        return (self.width, self.height)

    def resize(self, im):
        self.points = TrimmingBox._get_trimming_box(im)

    def extend(self, im):
        curbox = TrimmingBox._get_trimming_box(im)
        for i in [0, 1]:
            if self.points[i] is None or curbox[i] < self.points[i]:
                self.points[i] = curbox[i]
        for i in [2, 3]:
            if self.points[i] is None or curbox[i] > self.points[i]:
                self.points[i] = curbox[i]


class ImageLayout:
    # Insert types: 'rescale', 'top', 'middle', 'bottom'
    def __init__(self):
        self.frags = []

    @staticmethod
    def append_h(full_im, new_seg, pos, orient):
        if orient == 'top':
            full_im.paste(new_seg, (pos, 0))
        elif orient == 'middle':
            midpoint = full_im.height / 2
            startpos = midpoint - new_seg.height / 2
            full_im.paste(new_seg, (pos, int(startpos)))
        return full_im

    @staticmethod
    def append_v(full_im, new_seg, pos, orient):
        if orient == 'top':
            full_im.paste(new_seg, (0, pos))
        elif orient == 'middle':
            midpoint = full_im.width / 2
            startpos = midpoint - new_seg.width / 2
            full_im.paste(new_seg, (int(startpos), pos))
        return full_im

    @property
    def size(self):
        raise Exception(NotImplementedError)

    def insert(self, image, name=None, type='top'):
        self.frags.append({'name': name, 'image': image, 'type': type})

    def build(self):
        raise Exception(NotImplementedError)


class HImageLayout(ImageLayout):
    @property
    def size(self):
        mywidth = 0
        myheight = 0
        for frag in self.frags:
            mywidth += frag['image'].width
            if frag['image'].height > myheight:
                myheight = frag['image'].height
        return (mywidth, myheight)

    def build(self):
        res = Image.new('RGB', self.size, (255, 255, 255))
        curpos = 0
        for i in range(len(self.frags)):
            res = ImageLayout.append_h(res, self.frags[i]['image'], curpos, orient=self.frags[i]['type'])
            curpos += self.frags[i]['image'].width
        return res

class VImageLayout(ImageLayout):
    @property
    def size(self):
        mywidth = 0
        myheight = 0
        for frag in self.frags:
            myheight += frag['image'].height
            if frag['image'].width > mywidth:
                mywidth = frag['image'].width
        return (mywidth, myheight)

    def build(self):
        res = Image.new('RGB', self.size, (255, 255, 255))
        curpos = 0
        for i in range(len(self.frags)):
            res = ImageLayout.append_v(res, self.frags[i]['image'], curpos, orient=self.frags[i]['type'])
            curpos += self.frags[i]['image'].height
        return res