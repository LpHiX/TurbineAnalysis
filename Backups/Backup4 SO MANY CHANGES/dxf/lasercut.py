import pandas as pd
import numpy as np
import ezdxf
from ezdxf import colors
from ezdxf.enums import TextEntityAlignment
from ezdxf import recover
from ezdxf.addons.drawing import matplotlib

D_TSR = 5
R = 0.25 - 0.1
R0 = 0.025 + 0.02
maxcl = 1.4003
aoa = 7.5
B = 3


def twist(r):
    return 0

def chord(r):
    phi = np.arctan(2 / 3 / (5 * r / R))
    return (8 * np.pi * r * np.sin(phi)) / (3 * B * maxcl * (5 * r / R));

def transform(x, y, tw, ch):
    return (np.cos(tw) * x * ch + np.sin(tw) * y * ch, - np.sin(tw) * x * ch + np.cos(tw) * y * ch)


contour = pd.read_table("sg6043.dat",delim_whitespace=True,skiprows=[0],names=['x','y'],index_col=False)
doc = ezdxf.new(dxfversion="R2010")

msp = doc.modelspace()
x = contour["x"]
y = contour["y"]

allx = []
ally = []

N_sections = int(np.ceil((R - R0) / 0.0016))
r_list = np.linspace(R0, R, N_sections)
gridsize = int(np.ceil(np.sqrt(N_sections)))

# Add entities to a layout by factory methods: layout.add_...() 

for j in range(len(r_list)):
    r = r_list[j] 
    for i in range(1, len(contour)):
        (x1, y1) = transform(x[i], y[i], twist(r), chord(r))
        allx.append(x1)
        ally.append(y1)

minx = min(allx)
miny = min(ally)
maxx = max(allx)
maxy = max(ally)
ysize = maxy - miny
xsize = maxx - minx

N_x = 3;

for j in range(len(r_list)):
    r = r_list[j] 
    for i in range(1, len(contour) - 1):
        (x1, y1) = transform(x[i], y[i], twist(r), chord(r))
        (x2, y2) = transform(x[i+1], y[i+1], twist(r), chord(r))
        x1 = x1 + 1.1 * xsize * (j % N_x)
        x2 = x2 + 1.1 * xsize * (j % N_x)
        y1 = y1 + ysize * np.floor(j / N_x)
        y2 = y2 + ysize * np.floor(j / N_x)
        msp.add_line((x1, y1), (x2, y2), dxfattribs={"color": colors.YELLOW})

for j in range(len(r_list)):
    r = r_list[j] 
    for i in range(1, len(contour) - 1):
        (x1, y1) = transform(x[i], y[i], twist(r), chord(r))
        (x2, y2) = transform(x[i+1], y[i+1], twist(r), chord(r))
        x1 = xsize - x1 + 1.1 * xsize * (j % N_x)
        x2 = xsize - x2 + 1.1 * xsize * (j % N_x)
        y1 = ysize * np.floor(1.1 * N_sections / N_x) - (y1 + ysize * np.floor(j / N_x))
        y2 = ysize * np.floor(1.1 * N_sections / N_x) - (y2 + ysize * np.floor(j / N_x))
        msp.add_line((x1, y1), (x2, y2), dxfattribs={"color": colors.YELLOW})

# Save the DXF document.
doc.saveas("test.dxf")

doc, auditor = recover.readfile('test.dxf')
if not auditor.has_errors:
    matplotlib.qsave(doc.modelspace(), 'your.png')

print(xsize)