#####################################################
#												   #
#   Dev. by Jonathan Biteau (biteau@in2p3.fr)	   #
#					   last edit: 2022-02-23	   #
#####################################################

# Data analysis
import numpy as np

# Utilities
import os.path

# Create spherical-map
import healpy as hp

# Manage & convert astro coordinates
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column

# Plotting
import matplotlib
import matplotlib.pyplot as plt

def MapToHealpyCoord(l, b):
	"""
	Transforms to healpy spherical coordinate system
	"""
	theta = np.pi/2 - b
	phi = l

	return phi, theta
	

def HealpyCoordToMap(phi, theta):
	"""
	Transforms to healpy spherical coordinate system
	"""	
	b = np.pi/2 - theta
	l = phi
	
	return l, b	

		
def DatasetLoad(fdir = "../Data",filename = "AugerApJS2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo"):
	"""
	Loads the Auger phase 1 high-energy data set. File is identified by its directory *fdir* and filename *filename*.
	"""
	fname = os.path.join(fdir,filename)
	if os.path.isfile(fname):
		table = Table.read(fname, format='ascii.basic', delimiter="\t", guess=False)
		#equatorial and Galactic coordinates
		sc_equatorial = SkyCoord(table["RA"], table["Dec"], unit=u.deg)
		table.add_column(sc_equatorial, name="sc_eq")
		table.add_column(sc_equatorial.galactic, name="sc_gal")
		table.add_column(sc_equatorial.supergalactic, name="sc_sgal")		
		#Spherical coordinates for healpy
		phi, theta = MapToHealpyCoord(table["sc_gal"].l.to(u.rad).value, table["sc_gal"].b.to(u.rad).value)
		table.add_column(Column(phi), name="phi")
		table.add_column(Column(theta), name="theta")		
		
		export = ["UTC","sc_eq","sc_gal","sc_sgal","phi", "theta", "E","Expo"]
		
		return table[export]
		
		
def LoadCountMap(dataset, nside, Eth):
	"""
	Loads the count map of the auger phase 1 high-energy data set. 
	"""	
	#Select events >= Eth
	sel = dataset["E"]>=Eth
	
	# Pixel index for each events of coordinates (theta, phi)
	index = hp.ang2pix(nside, np.array(dataset["theta"][sel]), np.array(dataset["phi"][sel]))

	# Count map of parameter nside
	npix = hp.nside2npix(nside)# npix corresponds to the number of pixel associated to a NSIDE healpy map
	count_map = np.histogram(index, bins=np.arange(npix + 1))[0]
	
	return count_map		
		
		
def DirectionalExposureLoad(expo_tot_km2yrsr, fdir = "../Data",filename = "exposure.fits"):
	"""
	Loads the directional exposure associated to the Auger phase 1 high-energy data set. File is identified by its directory *fdir* and filename *filename*.
	"""
	fname = os.path.join(fdir,filename)
	if os.path.isfile(fname):
		hpx, header = hp.read_map(fname, h=True, verbose=False)
		return expo_tot_km2yrsr*hpx/np.sum(hpx)


def PlotHPmap(hp_map, title, color_bar_title, save_filename=None, projection="hammer", cmap='afmhot', vmin=None, vmax=None, ax=None, src_list=None):
	"""
	Loads the directional exposure associated to the Auger phase 1 high-energy data set. File is identified by its directory *fdir* and filename *filename*.
	"""	
	# Transform healpix map into matplotlib map (grid_map)
	xsize = 2000 # grid size for matplotlib
	ysize = int(xsize/2.)
	theta = np.linspace(np.pi, 0, ysize)
	phi = np.linspace(-np.pi, np.pi, xsize)
	PHI, THETA = np.meshgrid(phi, theta)
	nside = hp.get_nside(hp_map)
	grid_pix = hp.ang2pix(nside, THETA, PHI)
	grid_map = hp_map[grid_pix]

	# Define the figure
	width = 12# width of the figure
	fig = plt.figure(figsize=(width,width/1.5))
	plt.subplots_adjust(left=0.06, right=0.94, top=0.83, bottom=0.02)
	fig.suptitle(title, size=30, y = 0.94)

	if ax is None:
		ax = fig.add_subplot(111, projection=projection)
	labelsize = 16
	ax.tick_params(axis='x', labelsize=labelsize)
	ax.tick_params(axis='y', labelsize=labelsize)
			
	# Set the size of the other fonts
	fontsize = 22
	font = {'size': fontsize}
	matplotlib.rc('font', **font)
	matplotlib.rcParams['figure.max_open_warning']=0
	
	# minimum and maximum values along the z-scale (colorbar)
	if vmax is None:
		vmax = np.max(hp_map)
	if vmin is None:
		vmin = np.min(hp_map)
	
	# Plot the map, reverse the longitude axis "[::-1]" (astronomical convention)
	longitude = np.radians(np.linspace(-180, 180, xsize))
	latitude = np.radians(np.linspace(-90, 90, ysize))
	image = ax.pcolormesh(longitude[::-1], latitude, grid_map, rasterized=True, 
						  cmap=plt.get_cmap(cmap), shading='auto', vmin=vmin, vmax=vmax)
	cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05)
	cb.set_label(color_bar_title, size=fontsize)

	#Supergalactic plane
	l = np.radians(np.linspace(-180, 180, xsize))
	l = (l + np.pi) % 2*np.pi - np.pi
	b = np.zeros_like(l)	
	c = SkyCoord(sgl=l, sgb=b, frame='supergalactic', unit="rad")
	l_sgp, b_sgp = -c.galactic.l.degree, c.galactic.b.degree
	l_sgp = (l_sgp + 180.) % 360. - 180.
	#breaks the array if angular distance is large
	i_br = np.where(np.abs(l_sgp[:-1]-l_sgp[1:])>30)[0]
	plt.plot(np.radians(l_sgp[i_br[0]+1:i_br[1]]),np.radians(b_sgp[i_br[0]+1:i_br[1]]), color="tab:gray",linewidth=1)	
	
	# Plot the labels considering if it is galactic or equatorial coordinates
	ax.set_xticklabels([r"150$\degree$", r"120$\degree$", r"90$\degree$", r"60$\degree$", r"30$\degree$", r"GC", r"330$\degree$", r"300$\degree$", r"270$\degree$", r"240$\degree$", r"210$\degree$"])
	ax.set_title("Galactic")
	ax.set_xlabel('longitude', size=labelsize)
	ax.set_ylabel('latitude', size=labelsize)
		
	ax.grid(True, alpha=0.25)
	
	if(src_list!=None):
		for src in src_list:
			name, l, b = src[0], src[1], src[2]
			l = -l
			l = (l + 180.) % 360. - 180.
			plt.annotate(name, (np.radians(l), np.radians(b)), color = "gray", fontsize='small')
			plt.scatter(np.radians(l), np.radians(b), s=20, facecolors='none', edgecolors='gray')
	
	if(save_filename!=None):
		plt.savefig(save_filename)
	
	
def LoadSmoothedMap(hp_map, radius_deg, top_hat = True):
	"""
	Return smoothed map with a top-hat beam of radius radius_deg
	"""	
	nside = hp.get_nside(expo_map)
	radius_rad = np.radians(radius_deg)	

	def top_hat_beam(radius_rad, nside):
		b = np.linspace(0.0, np.pi, 10000)
		bw = np.where(abs(b) <= radius_rad, 1, 0)
		return hp.sphtfunc.beam2bl(bw, b, lmax=nside*3) #beam in the spherical harmonics space

	def fisher_beam(radius_rad, nside):
		b = np.linspace(0.0, np.pi, 10000)
		kappa = 1./(2*(1.-np.cos(radius_rad)))
		bw = np.exp(kappa*np.cos(b))
		return hp.sphtfunc.beam2bl(bw, b, lmax=nside*3)#beam in the spherical harmonics space
	
	if(top_hat):
		solid_angle = 2.*np.pi*(1. - np.cos(radius_rad))
		return hp.smoothing(hp_map, beam_window=top_hat_beam(radius_rad, nside), verbose=False) / solid_angle
	else:
		return hp.smoothing(hp_map, beam_window=fisher_beam(radius_rad, nside), verbose=False) #note: this one is not normalized


def LoadModelMap(filename, alpha, theta_deg, expo_map):
	"""
	Load a model map from the file filename
	"""
	data = Table.read(filename, format='ascii.no_header', delimiter=' ', names=['name', 'ra', 'dec', 'w'])
	nside = hp.get_nside(expo_map)

	sc = SkyCoord(ra=data['ra'], dec=data['dec'], frame='icrs', unit="deg")
	phi, theta = MapToHealpyCoord(sc.galactic.l.to(u.rad).value, sc.galactic.b.to(u.rad).value)
	index = hp.ang2pix(hp.get_nside(expo_map), theta, phi)
	cumul_weights = np.histogram(index, weights = data['w'], bins=np.arange(hp.nside2npix(nside)+1))[0]

	src_map = LoadSmoothedMap(cumul_weights, theta_deg, top_hat = False)
	model_map = LoadSmoothedMap(src_map, radius_deg=25, top_hat = True)
	
	return model_map/np.max(model_map)


def	PlotFluxMap(expo_map, count_map, Emin_EeV, radius_deg, save_filename=None, display_OneOverOmega = False, index_list=None):
	"""
	Plot flux map based on exposure and count maps.
	"""	
	# Load count map above Emin
	count_map = LoadCountMap(table_data, nside = hp.get_nside(expo_map), Eth = Emin_EeV)
	print("Total counts >=",Emin_EeV,"EeV:",np.sum(count_map))		

	# Smooth exposure and count maps
	smoothed_expo_map = LoadSmoothedMap(expo_map, radius_deg)
	smoothed_count_map = LoadSmoothedMap(count_map, radius_deg)	
	
	# Flux map 
	flux_map = smoothed_count_map / smoothed_expo_map
	sel = expo_map>np.max(expo_map)/100#empty pixels where the exposure is zero
	flux_map[np.invert(sel)] = np.nan #uniform above maximum dec
	flux_min = max(0,np.min(flux_map[sel]))
	flux_max = np.max(flux_map[sel])	

	# Plot skymap
	title = r"$\Phi(E_{\rm Auger} \geq"+ f"{Emin_EeV}"+ r"\:{\rm EeV})$ - "+ f"$\Psi = {radius_deg}\\degree$"
	color_bar_title = r"Flux [$10^{-3}\:\rm km^{-2} \: sr^{-1} \: yr^{-1}$]"
	scaling = 1E3
	PlotHPmap(scaling*flux_map, title, color_bar_title, vmin=0., vmax=scaling*flux_max, save_filename=save_filename)


def PlotLiMa(expo_map, table_data, Emin_EeV = 41, radius_deg = 27, save_filename=None):
	"""
	Li & Ma significance map display
	"""

	def LiMaMap(nside, Non, Noff, alpha):
		"""
		Li & Ma significance map computation according to Eq. 17 in Li & Ma (1983)
		"""
		Non_log_term = (1. + alpha)*Non / (alpha*(Non + Noff))
		Noff_log_term = (1. + alpha)*Noff / (Non + Noff)
		
		sig2_ov_2 = np.zeros_like(Non)
		ind = np.where((Non > 0) & (alpha > 0))# ensures non negative log terms induced by smoothing
		sig2_ov_2[ind] += Non[ind]*np.log(Non_log_term[ind])
		ind = np.where(Noff > 0)# ensures non negative log terms induced by smoothing
		sig2_ov_2[ind] += Noff[ind]*np.log(Noff_log_term[ind])

		return np.sign(Non-alpha*Noff)*np.sqrt(np.abs(2*sig2_ov_2))
		
	# Load count map above Emin
	print(expo_map)
	nside = hp.get_nside(expo_map)
	count_map = LoadCountMap(table_data, nside, Eth = Emin_EeV)

	# Smooth exposure and count maps
	smoothed_expo_map = LoadSmoothedMap(expo_map, radius_deg)
	smoothed_count_map = LoadSmoothedMap(count_map, radius_deg)	
	
	# Normalize to the total number of events in the smoothing region
	solid_angle_per_pix = 4.*np.pi / hp.nside2npix(nside)
	solid_angle = 2.*np.pi * (1. - np.cos(np.radians(radius_deg)))
	npix_per_region = solid_angle/solid_angle_per_pix

	# alpha, Non, & Noff
	alpha = smoothed_expo_map / (np.sum(expo_map)/npix_per_region - smoothed_expo_map)
	Non = smoothed_count_map*npix_per_region
	Noff = np.sum(count_map) - Non

	# Significance map
	sig_map = LiMaMap(nside, Non, Noff, alpha) # Compute the significance map
	isig_max = np.argmax(sig_map)
	sel = expo_map>np.max(expo_map)/100#empty pixels where the exposure is zero
	sig_map[np.invert(sel)] = np.nan #uniform above maximum dec
	sig_min = np.min(sig_map[sel])
	sig_max = np.max(sig_map[sel])	
	sig_abs = max(np.abs(sig_max),np.abs(sig_min))

	# Plot skymap
	title = r"$\sigma(E_{\rm Auger} \geq"+ f"{Emin_EeV}"+ r"\:{\rm EeV})$ - "+ f"$\Psi = {radius_deg}\\degree$"
	color_bar_title = r"Li & Ma significance [$\sigma$]"
	PlotHPmap(sig_map, title, color_bar_title, vmin=-sig_abs, vmax=sig_abs, cmap='twilight_shifted', save_filename=save_filename)	


def PrintFluxRegions(coord_list, table_data, expo_map, Emin_EeV, radius_deg):
	"""
	Print the flux at the list of given cordinates as well as that along the Galactic and supergalactic plane. The 1/Omega approach is used to estimate the flux and associated uncertainty: see https://journals.jps.jp/doi/10.7566/JPSCP.19.01102
	"""
	
	def print_flux(text, flux, var_flux):
		rounded_flux = np.round(flux*1E4)/10
		rounded_rms = np.round(np.sqrt(var_flux)*1E4)/10
		print(text+": ("+str(rounded_flux)+" +/- "+str(rounded_rms)+") x 1E-3 km-2 sr-1 yr-1")
		
	nside = hp.get_nside(expo_map)
	index = hp.ang2pix(nside, table_data["theta"], table_data["phi"])
	expo_weights = 1./(expo_map[index])#exp in : km2 yr sr / pixel	
	table_data.add_column(Column(expo_weights), name="w")
	table_sel = table_data[table_data["E"]>=Emin_EeV]
		
	# Isotropic flux over the field of view
	sel_expo = expo_map>np.max(expo_map)/100
	npix_per_sky = len(expo_map[sel_expo])

	index = hp.ang2pix(nside, table_sel["theta"], table_sel["phi"])
	is_in_FoV = expo_map[index]>np.max(expo_map)/100
	
	flux_iso = np.sum(is_in_FoV*table_sel["w"])
	var_iso = np.sum(is_in_FoV*table_sel["w"]**2)
	print_flux("Isotropic flux", flux_iso/npix_per_sky, var_iso/npix_per_sky**2)
		
	# Planar flux
	solid_angle_per_pix = 4.*np.pi / hp.nside2npix(nside)
	solid_angle = 4.*np.pi * np.sin(np.radians(radius_deg))
	npix_per_plane = solid_angle/solid_angle_per_pix

	## Galactic plane
	sel = np.abs(table_sel["sc_gal"].b.to(u.deg).value)<radius_deg
	flux_plane = np.sum(is_in_FoV[sel]*table_sel["w"][sel])
	var_plane = np.sum(is_in_FoV[sel]*table_sel["w"][sel]**2)
	print_flux("Galactic plane flux", flux_plane/npix_per_plane, var_plane/npix_per_plane**2)

	## Supergalactic plane
	sel = np.abs(table_sel["sc_sgal"].sgb.to(u.deg).value)<radius_deg
	flux_plane = np.sum(is_in_FoV[sel]*table_sel["w"][sel])
	var_plane = np.sum(is_in_FoV[sel]*table_sel["w"][sel]**2)
	print_flux("Supergalactic plane flux", flux_plane/npix_per_plane, var_plane/npix_per_plane**2)	
			
	# Flux at test positions
	solid_angle = 2.*np.pi * (1. - np.cos(np.radians(radius_deg)))
	npix_per_region = solid_angle/solid_angle_per_pix
	for i_test, src in enumerate(coord_list):
		sc_testpos = SkyCoord(l=src[1]*u.degree, b=src[2]*u.degree, frame='galactic')
		sel = sc_testpos.separation(table_sel["sc_gal"])<radius_deg*u.deg
		flux_region = np.sum(is_in_FoV[sel]*table_sel["w"][sel])
		var_region = np.sum(is_in_FoV[sel]*table_sel["w"][sel]**2)
		print_flux(src[0], flux_region/npix_per_region, var_region/npix_per_region**2)		
		

### Main ##########################################################
if __name__ == "__main__":		

	# Load data set and exposure
	table_data = DatasetLoad()
	expo_map = DirectionalExposureLoad(expo_tot_km2yrsr = 122E3)
	
	# List of interesting targets: name, l[deg], b[deg]
	GC = 		["Galactic center", 0., 		0.]
	CenA = 		["Cen A",		 	309.516, 	19.417]
	NGC253 =	["NGC 253", 		97.364, 	-87.965]
	VirgoCl = 	["Virgo Cluster", 	279.676, 	74.460]
	NGC4945 = 	["NGC 4945", 		305.272, 	13.340]	
	NGC4151 = 	["NGC 4151", 		155.077, 	75.063]	
	NGC4736 = 	["NGC 4736", 		123.362, 	76.008]		
	Mrk421 = 	["Mrk 421", 		179.832,	65.031]	
	NGC1275 = 	["NGC 1275", 		150.576, 	-13.261]	
	M82 = 		["M 82", 			141.409,	40.567]	
	M81_82 = 	["M 81, M 82", 		141.8, 		40.7]	
	CenA_NGC4945 = ["Cen A, NGC 4945", 307.4, 	16.4]		

	# Save Li & Ma map
	Emin_EeV = 41
	radius_deg = 24
	PlotLiMa(expo_map, table_data, Emin_EeV = Emin_EeV, radius_deg = radius_deg, save_filename = f"auger_lima_map_{Emin_EeV}EeV_{radius_deg}deg.pdf")	

	# Save flux map
	radius_deg = 25	
	for Emin_EeV in range(32,81):
		save_filename = f"auger_flux_map_{Emin_EeV}EeV_{radius_deg}deg.pdf"
		PlotFluxMap(expo_map, table_data, Emin_EeV = Emin_EeV, radius_deg = radius_deg, save_filename = save_filename)	
		
	# Compute the flux
	Emin_EeV = 40
	radius_deg = 25
	
	coord_list = []#Target, l[deg], b[deg]
	coord_list.append(GC)
	coord_list.append(CenA)
	coord_list.append(NGC253)
	coord_list.append(VirgoCl)
	
	index_list = []#Target, ipix
	for i, coord in enumerate(coord_list):
		phi, theta = MapToHealpyCoord(np.radians(coord[1]), np.radians(coord[2]))
		index = hp.ang2pix(hp.get_nside(expo_map), theta, phi)
		index_list.append([coord[0], index])

	PrintFluxRegions(coord_list, table_data, expo_map, Emin_EeV, radius_deg)
	
	# Save model maps
	path = "../Catalog_Based/Catalogs/ModelsUHECR/"
	ext = f"_EPO1st_threshold40"
	ext_title = " - "+ r"$\Psi = 25\degree$"
	models = []#file, alpha[%], theta[deg] above 40 EeV
	source_list = []#name, l[deg], b[deg]
	cat_name = []
	
	models.append(["Starburst galaxies (radio)"+ext_title,	path+"sbg/sbg"+ext,							9.44, 15.3])
	source_list.append([M82, NGC4945, NGC253])
	cat_name.append("sbg")
	
	models.append([r"Jetted AGN ($\gamma$-rays)"+ext_title,	path+"fermi_lat_agn/fermi_lat_agn"+ext,		5.82, 14.1])
	source_list.append([Mrk421, CenA, NGC1275])	
	cat_name.append("gagn")
	
	models.append(["All AGN (hard X-rays)"+ext_title, 		path+"swift_bat_agn/swift_bat_agn"+ext,		7.43, 16.2])
	source_list.append([CenA_NGC4945, NGC4151])	
	cat_name.append("xagn")
	
	models.append(["Galaxies > 1 Mpc (IR)"+ext_title, 		path+"2mass_hyperleda/2mass_hyperleda"+ext,	15.7, 15.1])
	source_list.append([M81_82, NGC4736, VirgoCl, CenA_NGC4945, NGC253])		
	cat_name.append("2mrs")
	
	for i, m in enumerate(models):	
		model_map = LoadModelMap(filename=m[1], alpha=m[2]/100, theta_deg=m[3], expo_map=expo_map)
		title = m[0]
		color_bar_title = "Model flux, "+r"$\Phi(E_{\rm Auger} \geq 40\:{\rm EeV})$ [arb. unit]"
		save_filename = "fig10_"+cat_name[i]+".pdf"
		PlotHPmap(model_map, title, color_bar_title, vmin=0,vmax=1, src_list = source_list[i], save_filename=save_filename)
