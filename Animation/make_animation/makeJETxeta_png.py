import plotly.graph_objects as go
import numpy as np
import plotly.io as pio


def plot(X,Y,E,index):
	trace = go.Scatter(
    	x=X,
    	y=Y,
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
		height=400,

    	title='T=%dfm/c'%int(index/100),
    	
    	#paper_bgcolor='black',  
    	#plot_bgcolor='black'    
    	template="plotly_dark"
	)
	fig = go.Figure(data=[trace], layout=layout)
	
	
	fig.update_layout(
    	xaxis=dict(showgrid=False,zeroline=False,range=[-8,8],linecolor='white',linewidth=5, title=r"$\eta^{*}_{s} $",mirror=True), 
    	yaxis=dict(showgrid=False,zeroline=False,range=[-5,5],linecolor='white',linewidth=5, title='X* [fm]',mirror=True)   
	)



	pio.write_image(fig, "store_Jetxeta/h_%d.png"%index) 
	#fig.show()





file=open("animation_Jetxy.dat",'r')


for i in range(500):
	print(i)
	print(" ")
	line=file.readline()
	Npartons=int(line)
	X=[]
	ETA=[]
	
	E=[]
	#get position information
	for j in range(Npartons):
		line=file.readline()
		parts=line.split()

		X.append( float(parts[0] ) )
		ETA.append( float(parts[8] ) )
		
		E.append( np.log( float(parts[6] ) )+4 )

	#end getting position information
	
	
	#plot
	#print(i)
	
	plot(ETA,X,E,i)

		
	


