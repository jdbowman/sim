from mpl_toolkits.basemap import Basemap
import numpy as np

def Plotcat(ra_0, dec_0, fov_size=15., grid_size=2., allsky=False):
	pkscat = np.loadtxt('source_list_pks_140MHz',usecols=(0,1))
	pkscat[:,0] *= 15.
	
	sky = Basemap(projection='ortho', lon_0=-ra_0, lat_0=dec_0, celestial=True)
	if allsky==True:
		sky.drawmeridians(np.arange(0.,360.,grid_size))
		sky.drawparallels(np.arange(90,-90,-grid_size))
		sky.drawmapboundary(fill_color='White')
		catx, caty = sky(pkscat[:,0],pkscat[:,1])
		sky.scatter(catx,caty,3,marker='o',color='Black')
	else:
		cnreq = np.array([[ra_0+fov_size/2., dec_0-fov_size/2.], 
				[ra_0-fov_size/2., dec_0+fov_size/2.]])
			# = [[ll_ra, ll_dec], [ur_ra, ur_dec]]
		cnrxy = np.transpose(np.array(sky(cnreq[:,0],cnreq[:,1])))
			# = [[ll_x, ll_y], [ur_x, ur_y]]
		cenxy = np.array(sky(ra_0,dec_0))
		cnrmap = cnrxy-np.array([cenxy,cenxy])
		m = Basemap(projection='ortho', lon_0=-ra_0, lat_0=dec_0, 	
			celestial=True, llcrnrx=cnrxy[0,0]-cenxy[0], 
			llcrnry=cnrxy[0,1]-cenxy[1], urcrnrx=cnrxy[1,0]-cenxy[0], 
			urcrnry=cnrxy[1,1]-cenxy[1])
		m.drawmeridians(np.arange(cnreq[0,0], cnreq[1,0], -grid_size), 
						labels=[0,0,0,1],fmt=format_label)
		m.drawparallels(np.arange(cnreq[0,1], cnreq[1,1], grid_size), 
						labels=[1,0,0,0])
		m.drawmapboundary(fill_color='White')
		catx, caty = m(pkscat[:,0],pkscat[:,1])
		m.scatter(catx,caty,3,marker='o',color='Black')
		
def format_label(deg):
	return ('%02.0fh%02.0fm') % (np.floor(deg/15.), 
									np.remainder(deg,15.)*60./15.)
