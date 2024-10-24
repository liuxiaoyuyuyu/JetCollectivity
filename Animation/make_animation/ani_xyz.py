from PIL import Image 

image_files=[]
for i in range(0,500):
	image_files.append('./store_xyz/h_%d.png' % (i*1) )
	print(i)
images=[Image.open(x) for x in image_files]

images[0].save('xyz.gif', save_all=True, append_images=images[1:], optimize=False, duration=10, loop=0)