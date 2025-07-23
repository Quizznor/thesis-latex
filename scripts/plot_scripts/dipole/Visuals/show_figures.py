import numpy as np
import ROOT
import ctypes
import os
##############################################################
#Make sure to run the codes ... before running this script 
###############################################################

isExist = os.path.exists("./Targeted")
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs("./Targeted")

isExist = os.path.exists("./Autocorrelation")
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs("./Autocorrelation")

isExist = os.path.exists("./FigureSets")
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs("./FigureSets")

##############################################################
#This loop creates the 4 plots in Figure 2 and the plot in figure 5
###############################################################
ROOT.gStyle.SetPalette(53)

filenames = ["Targeted/Results/Figure2SGP","Targeted/Results/Figure2GP","Targeted/Results/Figure2GC","Autocorrelation/Results/Figure2AC","Targeted/Results/Figure5"]
titles= ["Supergalactic plane","Galactic plane","Galactic center", "Autocorrelation","Centaurus region"]
for i in range(len(filenames)):

	inFileName = "../Targeted_Blind/"+filenames[i]+".root"
	inFile = ROOT.TFile.Open(inFileName, "READ")

	if( i<4):
		hP = inFile.Get("hProb")
	else:
		hP=  inFile.Get("hP")
	hP.SetStats(0);
	c1 = ROOT.TCanvas( 'c1', titles[i])
	c1.SetLeftMargin(0.15);
	c1.SetRightMargin(0.2);
	c1.SetTopMargin(0.1);
	c1.SetBottomMargin(0.15);
	c1.SetTickx()
	c1.SetTicky()
	c1.SetLogz()
	
	hP.SetTitle(titles[i])
	hP.SetTitleFont(11,"t")
	hP.GetYaxis().SetTitle("Threshold Energy, #it{E}_{th}  [EeV]")
	hP.GetXaxis().SetTitle("Search Angle, \Psi  [deg]")
	hP.GetZaxis().SetTitle("Local p-value")
	hP.GetXaxis().SetTitleOffset(1.2)
	hP.GetYaxis().SetTitleOffset(1.3)
	hP.GetZaxis().SetTitleOffset(1.3)
	hP.GetXaxis().CenterTitle()
	hP.GetYaxis().CenterTitle()
	hP.GetZaxis().CenterTitle()
	hP.GetYaxis().SetTitleFont(62)
	hP.GetYaxis().SetTitleSize(0.05)
	hP.GetXaxis().SetTitleFont(62)
	hP.GetXaxis().SetTitleSize(0.05)
	hP.GetZaxis().SetTitleFont(62)
	hP.GetZaxis().SetTitleSize(0.05)
	hP.GetYaxis().SetLabelFont(62)
	hP.GetYaxis().SetLabelSize(0.05)
	hP.GetXaxis().SetLabelFont(62)
	hP.GetXaxis().SetLabelSize(0.05)
	hP.GetZaxis().SetLabelFont(62)
	hP.GetZaxis().SetLabelSize(0.05)
	if(filenames[i]=="Targeted/Results/Figure5"):
		hP.GetZaxis().SetRangeUser(1E-7,1)
	else:
		hP.GetZaxis().SetRangeUser(1E-3,1)
	hP.SetTitleFont(hP.GetYaxis().GetTitleFont(),'t')

	hP.Draw("COLZ")

	binx, biny, binz = ctypes.c_int(),ctypes.c_int(),ctypes.c_int()
	hP.GetBinXYZ(hP.GetMinimumBin(),binx,biny,binz) or (binx, biny, binz)
	x,y,z=hP.GetXaxis().GetBinCenter(binx.value),hP.GetYaxis().GetBinCenter(biny.value),hP.GetZaxis().GetBinCenter(binz.value)
	MinMarker=ROOT.TMarker(x,y,34)
	MinMarker.SetMarkerColor(0)
	MinMarker.SetMarkerSize(2)
	MinMarker.Draw()

	c1.SaveAs(filenames[i].replace('/Results','')+".pdf")


##############################################################
#This part creates the plots in Figure 3 and 4
###############################################################

inFileName=r"../Catalog_Based/fig3.root"
inFile = ROOT.TFile.Open(inFileName ,"READ")
catnames=["Starburst galaxies (radio)","Jetted AGN (#gamma-rays)","All AGN (hard X-rays)","Galaxies > 1 Mpc (IR)"]
TS_canvas=inFile.Get("c_res")
leg=TS_canvas.FindObject("TPave")
leg.SetTextSize(0.05)
TS_canvas.SaveAs("fig3.pdf")

inFileName=r"../Catalog_Based/fig4.root"
inFile = ROOT.TFile.Open(inFileName ,"READ")
for i in range(4):
	hTS=TS_canvas.FindObject("G_TS"+str(i))
	xTS = hTS.GetX()
	yTS = hTS.GetY()
	binmax=np.argmax(yTS)
	hf=TS_canvas.FindObject("G_alpha"+str(i))
	yf = hf.GetY()
	hth=TS_canvas.FindObject("G_theta"+str(i))
	yth = hth.GetY()

	AbsMinMarker=ROOT.TMarker(yf[binmax]/100.,yth[binmax],20)
	AbsMinMarker.SetMarkerStyle(20)
	AbsMinMarker.SetMarkerColor(634)
	for Eth in range(32,80):
		cTS = inFile.Get("c2D_cat"+str(i)+"_Eth"+str(Eth))
		cTS.cd()
		if(Eth==xTS[binmax]):
			outname="./Figure4_"+catnames[i]+".pdf"
			cTS.SaveAs(outname)
		leg=cTS.FindObject("leg")
		AbsMinMarker.Draw()
		leg.AddEntry(AbsMinMarker,"Global best-fit","p")
		outname="./FigureSets/Figure4_"+str(i)+"_Eth"+str(Eth)+".pdf"
		cTS.SaveAs(outname)

