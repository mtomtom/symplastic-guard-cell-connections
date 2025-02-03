from PIL import Image, ImageFont, ImageDraw

font = ImageFont.truetype("/usr/share/fonts/truetype/asana-math/Asana-Math.otf",48)

pressure = [0.5, 1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]

prefix = 'cross' 
n = 0
for i in range(0,14):
    in_file = prefix + '{0:03}'.format(i) + '.png'
    
    img = Image.open(in_file)
    draw = ImageDraw.Draw(img)    
    draw.text((1000,0),"P_GC = " + str(pressure[i]) + " MPa", "white",font=font)
    out_file = prefix + '-time' + '{0:03}'.format(i) + '.png'
    img.save(out_file)




