import plotly.graph_objects as go
import numpy as np
import plotly.io as pio


def plot(X,Y,Z,E,index):
	trace = go.Scatter3d(
    	x=X,
    	y=Y,
    	z=Z,
    	mode='markers',
    	marker=dict(
        	size=E,
        	color=E,
        	colorscale='hot',
        	opacity=0.8
    	),
    	showlegend=False
	)

	layout = go.Layout(
		width=800,
		height=800,

    	title='T=%dfm/c'%int(index/100),
    	scene=dict(
        	xaxis=dict(zeroline=False,range=[-5,5],mirror=True,linecolor='white',linewidth=5,showgrid=False, title='Y [fm]'),  
        	yaxis=dict(zeroline=False, range=[-5,5],mirror=True,linecolor='white',linewidth=5,showgrid=False, title='Z [fm]'),  
        	zaxis=dict(zeroline=False, range=[-5,5],mirror=True,linecolor='white',linewidth=5,showgrid=False, title='X [fm]')   
    	),
    	#paper_bgcolor='black',  
    	#plot_bgcolor='black'    
    	template="plotly_dark"
	)
	fig = go.Figure(data=[trace], layout=layout)
	fig.add_trace(go.Scatter3d(
        	x=[5,5],
        	y=[5,5],
        	z=[-5,5],
        	mode='lines', 
        	line=dict(width=5, color='white'),
        	showlegend=False  
	))

	fig.add_trace(go.Scatter3d(
        	x=[-5,5],
        	y=[-5,-5],
        	z=[-5,-5],
        	mode='lines',  
        	line=dict(width=2, color='white'),
        	showlegend=False  
	))

	fig.add_trace(go.Scatter3d(
        	x=[-5,-5],
        	y=[-5,5],
        	z=[-5,-5],
        	mode='lines',  
        	line=dict(width=2, color='white'),
        	showlegend=False  
	))

	fig.add_trace(go.Scatter3d(
        	x=[-5,-5],
        	y=[-5,-5],
        	z=[-5,5],
        	mode='lines',  
        	line=dict(width=2, color='white'),
        	showlegend=False  
	))

	fig.add_trace(go.Scatter3d(
        	x=[5,5],
        	y=[-5,5],
        	z=[5,5],
        	mode='lines', 
        	line=dict(width=5, color='white'),
        	showlegend=False  
	))
	fig.add_trace(go.Scatter3d(
        	x=[-5,5],
        	y=[5,5],
        	z=[5,5],
        	mode='lines', 
        	line=dict(width=5, color='white'),
        	showlegend=False 
	))	

	fig.update_layout(
    	scene=dict(
        	aspectmode='cube' 
    	)
	)


	fig.update_layout(
    	scene_camera=dict(
        	eye=dict(x=1.5, y=0.8, z=1.3)  
    	)

	)



	pio.write_image(fig, "store_xyz/h_%d.png"%index) 
	#fig.show()





file=open("animation_xyz.dat",'r')


for i in range(500):
	print(i)
	print(" ")
	line=file.readline()
	Npartons=int(line)
	X=[]
	Y=[]
	Z=[]
	E=[]
	#get position information
	for j in range(Npartons):
		line=file.readline()
		parts=line.split()
		if( float(parts[8] )>0.8 ):
			continue
		X.append( float(parts[0] ) )
		Y.append( float(parts[1] ) )
		Z.append( float(parts[2] ) )
		E.append( np.log( float(parts[6] ) )+4 )

	#end getting position information
	
	
	#plot
	#print(i)
	
	plot(Y,Z,X,E,i)

		
	


