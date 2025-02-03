from PIL import Image, ImageFont, ImageDraw
import os
import numpy as np
import globals as const

font = ImageFont.truetype("/usr/share/fonts/truetype/asana-math/Asana-Math.otf",72)

p_gc = np.linspace(const.gc_start, const.gc_end, const.num_points) 
p_sc = np.linspace(const.sc_start, const.sc_end, const.num_points) 

p_gc = np.concatenate([[const.p_gc_i],p_gc])
p_sc = np.concatenate([[const.p_sc_i],p_sc])

prefix = 'top'
file_path = const.out_dir + "E1_" + str(const.E1_E3_gc) + "_E2_" + str(const.E2_gc) + const.outdirSuffix

i = 0
for i in range(0,len(p_gc)):
    Process.Mesh__System__Open('', file_path + '/GC_' + str(int(p_gc[i]*10)) + '.mdxm', 'No', 'No')
    Process.Mesh__Selection__Action__Clear_Selection('All')
    Process.Tools__System__Snapshot('tmp%03d.png'%i, 'No', '0', '0', '1.0', '95')
    
for i in range(0,len(p_gc)):
    in_file = 'tmp{0:03}'.format(i) + '.png'
        
    img = Image.open(in_file)
    overlay_img = Image.open("heatmap.png") 

    draw = ImageDraw.Draw(img)    

    overlay_size = overlay_img.size
    x = 100
    y = 600
    position = (x, y, x + overlay_size[0], y + overlay_size[1])
    
    draw.text((400,40),"GC Pressure = " + str(round(p_gc[i],2)) + " MPa", "white",font=font)
    draw.text((400,120),"SC  Pressure = " + str(round(p_sc[i],2)) + " MPa", "white",font=font)

    draw.text((10,1500),"Species: Barley", "white",font=font)
    draw.text((10,1580),"E1/E3 = 40 MPa", "white",font=font)
    draw.text((10,1660),"E2 = 75 MPa", "white",font=font)
    
    img.paste(overlay_img, position)

    out_file = prefix  + '{0:03}'.format(i) + '.png'
    img.save(out_file)
    
command = 'ffmpeg -i ' + prefix + '%03d.png ' + 'tmp.gif'
os.system(command)
command = 'convert -delay 32x100 ' + 'tmp.gif ' + prefix +'.gif'
os.system(command)
command = 'rm tmp*'
os.system(command)
command = 'mv *.png *.gif ' + file_path 
os.system(command)
command = 'mv ' + file_path + '/FemMembranes.png ./'
os.system(command)
command = 'mv ' + file_path + '/heatmap.png ./'
os.system(command)


