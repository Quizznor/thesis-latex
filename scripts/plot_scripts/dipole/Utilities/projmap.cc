#include <iostream>
#include <fstream>

#include "projmap.h"
#include "maptools.h"
#include "common.h"

#include "STClibrary.h"

// ROOT
#include "TLine.h"
#include "TPolyLine.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TPad.h"


static bool (*gCoordFcnAng2XY)(double, double, double&, double&, double, double, double);
static bool (*gCoordFcnXY2Ang)(double, double, double&, double&, double, double, double);
static char gObjName[1024];
static int gObjNumber = 0 ;
static char * GetObjName() {sprintf(gObjName,"SkyMap%d",gObjNumber++); return gObjName;}


// Background of the map
double gBackground = -999.;

// Color level of the palette
int gMaxPalette = 500;


TProjMap::TProjMap(const THealpixMap& map, int sizeX, int sizeY, double decLimit)
  : fMap(map), fSizeX(sizeX), fSizeY(sizeY)
{
  if( map.NPix() == 0 )
    {
      cout << "The map does not contain pixel. Set map size first. Returning." << endl;
      return;
    }
  fLongCenter = 0.; fLatCenter = 0.;
  fCosLatRef = 0.; fSinLatRef = 0.;
  fXMin = -180.; fXMax = 180; fYMin = -90.; fYMax = 90.;

	fXMarginLeft=0.12;
  fXMarginRight=0.08;
  fXMarginTop=0.16;
  fXMarginBottom=0.19;

  gCoordFcnXY2Ang = XYtoAngMollweide;
  gCoordFcnAng2XY = AngtoXYMollweide;

  fImage = map.Map2MollProj(sizeX, sizeY, gBackground, decLimit);
}


TProjMap::TProjMap(const THealpixMap& map, const THealpixMap& ptMap, int sizeX, int sizeY, double decLimit)
  : fMap(map), fSizeX(sizeX), fSizeY(sizeY)
{
  if( map.NPix() == 0 )
    {
      cout << "The map does not contain pixel. Set map size first. Returning." << endl;
      return;
    }
  fLongCenter = 0.; fLatCenter = 0.;
  fCosLatRef = 0.; fSinLatRef = 0.;
  fXMin = -180.; fXMax = 180; fYMin = -90.; fYMax = 90.;

	fXMarginLeft=0.12;
  fXMarginRight=0.08;
  fXMarginTop=0.16;
  fXMarginBottom=0.19;

  gCoordFcnXY2Ang = XYtoAngMollweide;
  gCoordFcnAng2XY = AngtoXYMollweide;

  fImage = map.Map2MollProj(sizeX, sizeY, gBackground, decLimit);
  fPtImage = ptMap.Map2MollProj(sizeX, sizeY, gBackground, decLimit);
}


TProjMap::TProjMap(const THealpixMap& map, int sizeX, int sizeY, double longCenter, double latCenter, double radius, double decLimit)
  : fMap(map), fSizeX(sizeX), fSizeY(sizeY), fLongCenter(longCenter), fLatCenter(latCenter)
{
  if( map.NPix() == 0 )
    {
      cout << "The map does not contain pixel. Set map size first. Returning." << endl;
      return;
    }
  fXMin = -radius; fXMax = radius; fYMin = -radius; fYMax = radius;
  fCosLatRef = cos(latCenter*DTOR); fSinLatRef = sin(latCenter*DTOR);

  gCoordFcnXY2Ang = XYtoAngLambertAzimuthal;
  gCoordFcnAng2XY = AngtoXYLambertAzimuthal;

  fImage = map.Map2LambertAzimuthalProj(sizeX, sizeY, longCenter, latCenter, radius, gBackground, decLimit);
}


TProjMap::~TProjMap()
{
  if(fCanvas) delete fCanvas;
  fImage.clear();
}


void TProjMap::SkyMap(string mapTitle)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = fMap.Max();
  double vMin = fMap.Min();
  for(int i = 1; i <= fSizeX; i++)
    {
      for(int j = 1; j <= fSizeY; j++)
        {
          double v = fImage[i-1][j-1];
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
        }
    }

  // Number of color levels
  mapProj->SetContour(gMaxPalette);
  mapProj->SetMaximum(vMax); mapProj->SetMinimum(vMin);

  // Display the color table only if the bins are not empty
  if( fMap.Max() != 0 )
    {
      mapProj->Draw("COLZ");
      fCanvas->Update();
      
      TPaletteAxis * colorTable = (TPaletteAxis *) mapProj->GetListOfFunctions()->FindObject("palette");
		colorTable->SetX1NDC(0.050); colorTable->SetX2NDC(0.075);
      colorTable->SetLabelFont(62);
    }
  // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
  else {mapProj->Draw("COL"); fCanvas->Update();}

	gPad->Update();
	TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
	title->SetTextAlign(11);
	title->SetTextFont(62);	
	title->SetTextSize(0.05);	
	title->SetBorderSize(0);
  title->SetX1NDC(0.01);
  title->SetY1NDC(0.9);
}

void TProjMap::SkyMap(string mapTitle, std::string sztitle, std::string sbeam_title)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  TH2D* mapPtProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = fMap.Max();
  double vMin = fMap.Min();
  for(int i = 1; i <= fSizeX; i++)
    {
      for(int j = 1; j <= fSizeY; j++)
        {
          double v = fImage[i-1][j-1];
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
					if(fPtImage.size()==fImage.size()) mapPtProj->SetBinContent(mapPtProj->GetBin(i,j),fPtImage[i-1][j-1]);
        }
    }

  // Number of color levels
  mapProj->SetContour(gMaxPalette);
  mapProj->SetMaximum(vMax); mapProj->SetMinimum(vMin);


	//Insert map
  TH2D* mapProj_insert = new TH2D(GetObjName(), mapTitle.c_str(), int(fSizeX/12.), -30., 30., int(fSizeY/6.), -30., 30.);
  mapProj_insert->GetXaxis()->SetLabelSize(0.); mapProj_insert->GetYaxis()->SetLabelSize(0.);
  mapProj_insert->GetXaxis()->SetAxisColor(0); mapProj_insert->GetYaxis()->SetAxisColor(0);
  mapProj_insert->SetTickLength(0,"X"); mapProj_insert->SetTickLength(0,"Y");
  mapProj_insert->SetStats(0);

	for(int i = 1; i<=mapProj_insert->GetNbinsX(); i++){
		for(int j = 1; j<=mapProj_insert->GetNbinsY(); j++){
			double coord_x = mapProj_insert->GetXaxis()->GetBinCenter(i);
			double coord_y = mapProj_insert->GetYaxis()->GetBinCenter(j);
			int pix0 = mapPtProj->FindBin(coord_x, coord_y);
			mapProj_insert->SetBinContent(mapProj_insert->GetBin(i,j),mapPtProj->GetBinContent(pix0));
		}
	}
	// Number of color levels
  mapProj_insert->SetContour(gMaxPalette);
  mapProj_insert->SetMaximum(vMax); mapProj_insert->SetMinimum(vMin);
  mapProj_insert->SetTitle("");

	double pad_width_X=0.135;
	double pad_width_Y=pad_width_X*fSizeX/(0.95*fSizeY) ;
  // Display the color table only if the bins are not empty
  if( fMap.Max() != 0 )
	{
		mapProj->Draw("COLZ");
		fCanvas->Update();

		//Draw the insert
		TPad *pad = new TPad("insert","", 0.99-pad_width_X,0.015,0.99,0.015+pad_width_Y);
		pad->SetLeftMargin(0.);
		pad->SetRightMargin(0.);
		pad->SetTopMargin(0.);
		pad->SetBottomMargin(0.);
		pad->Draw();
		pad->cd();
		mapProj_insert->Draw("COL");
		fCanvas->cd();
		fCanvas->Update();
      
		TPaletteAxis * colorTable = (TPaletteAxis *) mapProj->GetListOfFunctions()->FindObject("palette");
		colorTable->SetX1NDC(0.050); colorTable->SetX2NDC(0.075);
		colorTable->SetLabelFont(62);

		TLatex tlz;
		tlz.SetNDC();
		tlz.SetTextAlign(12);
		tlz.SetTextAngle(90.);
		tlz.SetTextSize(0.045);
		tlz.SetTextFont(62);
		tlz.DrawLatex(0.03,0.5-0.01*sztitle.size(),sztitle.c_str());

		TLatex tlbeam;
		tlbeam.SetNDC();
		tlbeam.SetTextAlign(12);
		tlbeam.SetTextSize(0.04);
		tlbeam.SetTextFont(62);

		tlbeam.DrawLatex(0.855,0.31,"Beam size");
		tlbeam.DrawLatex(0.855,0.265,sbeam_title.c_str());
	}
  // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
  else {mapProj->Draw("COL"); fCanvas->Update();}

	gPad->Update();
	TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
	title->SetTextAlign(11);
	title->SetTextFont(62);	
	title->SetTextSize(0.05);	
	title->SetBorderSize(0);
  title->SetX1NDC(0.01);
  title->SetY1NDC(0.9);
}



TGraph* GraphContour(TH2D* hist, double contour_level){

	TH2D* mapProj = (TH2D*)hist->Clone();
	double list_contours[1] = {contour_level};	
	mapProj->SetContour(1,list_contours);

	TCanvas can_tmp;
	can_tmp.cd();
	mapProj->Draw("CONT Z LIST");		
	can_tmp.Update();

	TObjArray *plah = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
	TList *list = (TList*)plah->At(0);
	TGraph* Gcont = (TGraph*)list->First();
	Gcont->SetLineWidth(2);

	return  (TGraph*)Gcont->Clone();
}

void TProjMap::SkyMapContours(string mapTitle, std::vector< double > vcontours, std::vector<Color_t> vcolor)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = fMap.Max();
  double vMin = fMap.Min();
  for(int i = 1; i <= fSizeX; i++)
    {
      for(int j = 1; j <= fSizeY; j++)
        {
          double v = fImage[i-1][j-1];
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
        }
    }

	std::vector< TGraph* > vGcont;
	for(unsigned int i=0; i<vcontours.size(); i++){
		vGcont.push_back(GraphContour(mapProj, vcontours[i]));
		vGcont.back()->SetLineColor(vcolor[i]);
	}
  // Number of color levels
  mapProj->SetContour(gMaxPalette);
  mapProj->SetMaximum(vMax); mapProj->SetMinimum(vMin);

  // Display the color table only if the bins are not empty
  if( fMap.Max() != 0 )
    {
      mapProj->Draw("COLZ");
      fCanvas->Update();
      
      TPaletteAxis * colorTable = (TPaletteAxis *) mapProj->GetListOfFunctions()->FindObject("palette");
		colorTable->SetX1NDC(0.050); colorTable->SetX2NDC(0.075);
      colorTable->SetLabelFont(62);
    }
  // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
  else {mapProj->Draw("COL"); fCanvas->Update();}

	gPad->Update();
	TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
	title->SetTextAlign(11);
	title->SetTextFont(62);	
	title->SetTextSize(0.05);	
	title->SetBorderSize(0);
  title->SetX1NDC(0.01);
  title->SetY1NDC(0.9);
}

void TProjMap::SkyMapContoursOnly(string mapTitle, std::vector< double > vcontours, std::vector<Color_t> vcolor, bool draw_same)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  TH2D* mapProj_plot = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj_plot->GetXaxis()->SetLabelSize(0.); mapProj_plot->GetYaxis()->SetLabelSize(0.);
  mapProj_plot->GetXaxis()->SetAxisColor(0); mapProj_plot->GetYaxis()->SetAxisColor(0);
  mapProj_plot->SetTickLength(0,"X"); mapProj_plot->SetTickLength(0,"Y");
  mapProj_plot->SetStats(0);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = fMap.Max();
  double vMin = fMap.Min();
  for(int i = 1; i <= fSizeX; i++){
		for(int j = 1; j <= fSizeY; j++)
		{
          double v = fImage[i-1][j-1];
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
          mapProj_plot->SetBinContent(mapProj->GetBin(i,j),0.);
        }
    }

	std::vector< TGraph* > vGcont;
	for(unsigned int i=0; i<vcontours.size(); i++){
		vGcont.push_back(GraphContour(mapProj, vcontours[i]));
		vGcont.back()->SetLineColor(vcolor[i]);
	}

	 // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
	if(!draw_same){
		fCanvas->cd();
		mapProj_plot->Draw("COL");
	}
	for(unsigned int i=0; i<vGcont.size(); i++) vGcont[i]->Draw("same l");
	if(!draw_same) fCanvas->Update();
}

void TProjMap::SkyMap(string mapTitle, double mini, double maxi)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = maxi;//max(fMap.Max(),maxi);
  double vMin = mini;//min(fMap.Min(),mini);
  for(int i = 1; i <= fSizeX; i++)
    {
      for(int j = 1; j <= fSizeY; j++)
        {
          double v = fImage[i-1][j-1];
          /*
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
           */ 
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
        }
    }

  // Number of color levels
  mapProj->SetContour(gMaxPalette);
  mapProj->SetMaximum(vMax); mapProj->SetMinimum(vMin);

  // Display the color table only if the bins are not empty
  if( fMap.Max() != 0 )
    {
      mapProj->Draw("COLZ");
      fCanvas->Update();
      
      TPaletteAxis * colorTable = (TPaletteAxis *) mapProj->GetListOfFunctions()->FindObject("palette");
		colorTable->SetX1NDC(0.050); colorTable->SetX2NDC(0.075);
      colorTable->SetLabelFont(62);
    }
  // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
  else {mapProj->Draw("COL"); fCanvas->Update();}

	gPad->Update();
	TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
	title->SetTextAlign(11);
	title->SetTextFont(62);	
	title->SetTextSize(0.05);	
	title->SetBorderSize(0);
  title->SetX1NDC(0.01);
  title->SetY1NDC(0.9);
}

void TProjMap::SkyMap(string mapTitle, double mini, double maxi, std::string sztitle, std::string sbeam_title)
{
  InitCanvas();

  // Two dimensional histogram containing the projected map
  TH2D* mapProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);
  mapProj->GetXaxis()->SetLabelSize(0.); mapProj->GetYaxis()->SetLabelSize(0.);
  mapProj->GetXaxis()->SetAxisColor(0); mapProj->GetYaxis()->SetAxisColor(0);
  mapProj->SetTickLength(0,"X"); mapProj->SetTickLength(0,"Y");
  mapProj->SetStats(0);

  // Two dimensional histogram containing the projected point map
  TH2D* mapPtProj = new TH2D(GetObjName(), mapTitle.c_str(), fSizeX, fXMin, fXMax, fSizeY, fYMin, fYMax);

  // Get the extremum values of the projected map
  // Fill the 2D histogram
  double vMax = max(fMap.Max(),maxi);
  double vMin = min(fMap.Min(),mini);
  for(int i = 1; i <= fSizeX; i++)
    {
      for(int j = 1; j <= fSizeY; j++)
        {
          double v = fImage[i-1][j-1];
          if( v != gBackground )
            {
              vMax = max(v,vMax);
              vMin = min(v,vMin);
            }
          mapProj->SetBinContent(mapProj->GetBin(i,j),fImage[i-1][j-1]);
					if(fPtImage.size()==fImage.size()) mapPtProj->SetBinContent(mapPtProj->GetBin(i,j),fPtImage[i-1][j-1]);
        }
    }

  // Number of color levels
  mapProj->SetContour(gMaxPalette);
  mapProj->SetMaximum(vMax); mapProj->SetMinimum(vMin);

	//Insert map
  TH2D* mapProj_insert = new TH2D(GetObjName(), mapTitle.c_str(), int(fSizeX/12.), -30., 30., int(fSizeY/6.), -30., 30.);
  mapProj_insert->GetXaxis()->SetLabelSize(0.); mapProj_insert->GetYaxis()->SetLabelSize(0.);
  mapProj_insert->GetXaxis()->SetAxisColor(0); mapProj_insert->GetYaxis()->SetAxisColor(0);
  mapProj_insert->SetTickLength(0,"X"); mapProj_insert->SetTickLength(0,"Y");
  mapProj_insert->SetStats(0);

	for(int i = 1; i<=mapProj_insert->GetNbinsX(); i++){
		for(int j = 1; j<=mapProj_insert->GetNbinsY(); j++){
			double coord_x = mapProj_insert->GetXaxis()->GetBinCenter(i);
			double coord_y = mapProj_insert->GetYaxis()->GetBinCenter(j);
			int pix0 = mapPtProj->FindBin(coord_x, coord_y);
			mapProj_insert->SetBinContent(mapProj_insert->GetBin(i,j),mapPtProj->GetBinContent(pix0));
		}
	}
	// Number of color levels
  mapProj_insert->SetContour(gMaxPalette);
  mapProj_insert->SetMaximum(vMax); mapProj_insert->SetMinimum(vMin);
  mapProj_insert->SetTitle("");

	double pad_width_X=0.135;
	double pad_width_Y=pad_width_X*fSizeX/(0.95*fSizeY) ;

  // Display the color table only if the bins are not empty
  if( fMap.Max() != 0 )
  {
		mapProj->Draw("COLZ");
		fCanvas->Update();

		//Draw the insert
		TPad *pad = new TPad("insert","", 0.99-pad_width_X,0.015,0.99,0.015+pad_width_Y);
		pad->SetLeftMargin(0.);
		pad->SetRightMargin(0.);
		pad->SetTopMargin(0.);
		pad->SetBottomMargin(0.);
		pad->Draw();
		pad->cd();
		mapProj_insert->Draw("COL");
		fCanvas->cd();
		fCanvas->Update();

		TPaletteAxis * colorTable = (TPaletteAxis *) mapProj->GetListOfFunctions()->FindObject("palette");
		colorTable->SetX1NDC(0.050); colorTable->SetX2NDC(0.075);
		colorTable->SetLabelFont(62);

		TLatex tlz;
		tlz.SetNDC();
		tlz.SetTextAlign(12);
		tlz.SetTextAngle(90.);
		tlz.SetTextSize(0.045);
		tlz.SetTextFont(62);
		tlz.DrawLatex(0.03,0.5-0.01*sztitle.size(),sztitle.c_str());

		TLatex tlbeam;
		tlbeam.SetNDC();
		tlbeam.SetTextAlign(12);
		tlbeam.SetTextSize(0.04);
		tlbeam.SetTextFont(62);

		tlbeam.DrawLatex(0.855,0.31,"Beam size");
		tlbeam.DrawLatex(0.855,0.265,sbeam_title.c_str());
	}
  // No color table displayed. Relevant to print only markers, which design for instance events or 
  // astrophysical objects.
  else {mapProj->Draw("COL"); fCanvas->Update();}

	gPad->Update();
	TPaveText *title = (TPaveText*)gPad->GetPrimitive("title");
	title->SetTextAlign(11);
	title->SetTextFont(62);	
	title->SetTextSize(0.05);	
	title->SetBorderSize(0);
  title->SetX1NDC(0.01);
  title->SetY1NDC(0.9);
}

void TProjMap::ShowGrid(double longStep, double latStep, int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  MakeLongitudesLatitudes(longStep, latStep, drawLine);
}

void TProjMap::ShowGrid_0_360(double longStep, double latStep, int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  MakeLongitudesLatitudes_0_360(longStep, latStep, drawLine);
}

void TProjMap::ShowFOV(double decLimit, int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  // Set up (l,b) for viewing range
  double raFOV = 0.;
  vector<double> lFOV(drawLine.fNPts),bFOV(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      raFOV = k*(360./drawLine.fNPts);
      radec2gal(raFOV/15.0, decLimit, &lFOV[k], &bFOV[k]); // 15.0 * convert from (24)hh (360)degrees
    }
  DrawRegion(lFOV, bFOV, drawLine);
}


void TProjMap::ShowFOV_in_Eq(double decLimit, int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  // Set up (l,b) for viewing range
  double raFOV = 0.;
  vector<double> lFOV(drawLine.fNPts),bFOV(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      raFOV = k*(360./drawLine.fNPts);
			lFOV[k] = raFOV;
			bFOV[k] = decLimit;
    }
  DrawRegion(lFOV, bFOV, drawLine);
}

void TProjMap::ShowFOV_in_SGC(double decLimit, int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  // Set up (l,b) for viewing range
  double raFOV = 0.;
  vector<double> lFOV(drawLine.fNPts),bFOV(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      raFOV = k*(360./drawLine.fNPts);
			double l, b;
      radec2gal(raFOV/15.0, decLimit, &l, &b);// 15.0 * convert from (24)hh (360)degrees
			gal2Sgal(l,b, lFOV[k], bFOV[k]);

    }
  DrawRegion(lFOV, bFOV, drawLine);
}

void TProjMap::ShowSGP(int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;

  double stepLongSGP = 360.0/drawLine.fNPts;
  // Set up (l,b) for SGP
  double sLongSGPmin = 0, sLongSGP = 0.0, sLatSGP = 0.0; 
  vector<double> lSGP(drawLine.fNPts),bSGP(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      sLongSGP = sLongSGPmin+k*stepLongSGP;
      Sgal2gal(sLongSGP,sLatSGP,lSGP[k],bSGP[k]);
    }
  DrawRegion(lSGP, bSGP, drawLine);
}

void TProjMap::ShowSGP_in_Eq_0_360(int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;

  double stepLongSGP = 360.0/drawLine.fNPts;
  // Set up (l,b) for SGP
  double sLongSGPmin = 0, sLongSGP = 0.0, sLatSGP = 0.0; 
  vector<double> lSGP(drawLine.fNPts),bSGP(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      sLongSGP = sLongSGPmin+k*stepLongSGP;
      double l, b;
      Sgal2gal(sLongSGP,sLatSGP,l,b);
      double ra, dec;
      gal2radec(l, b, &ra, &dec);
      ra*=kSTC::Hrtod;
      ra+=180.;
      lSGP[k] = ra;
      bSGP[k] = dec;
    }
  DrawRegion(lSGP, bSGP, drawLine);
}

void TProjMap::ShowGP_in_Eq_0_360(int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;

  double stepLongGP = 360.0/drawLine.fNPts;
  // Set up (l,b) for SGP
  double sLongGPmin = 0, sLongGP = 0.0, sLatGP = 0.0; 
  vector<double> lGP(drawLine.fNPts),bGP(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      sLongGP = sLongGPmin+k*stepLongGP;
      double ra, dec;
      gal2radec(sLongGP, sLatGP, &ra, &dec);
      ra*=kSTC::Hrtod;
      ra+=180.;
      lGP[k] = ra;
      bGP[k] = dec;
    }
  DrawRegion(lGP, bGP, drawLine);
}



void TProjMap::ShowGP_in_SGC(int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;

  double stepLongGP = 360.0/drawLine.fNPts;
  // Set up (l,b) for SGP
  double sLongGPmin = 0, sLongGP = 0.0, sLatGP = 0.0; 
  vector<double> lGP(drawLine.fNPts),bGP(drawLine.fNPts);
  for(unsigned int k = 0; k < drawLine.fNPts; k++)
    {
      sLongGP = sLongGPmin+k*stepLongGP;
      gal2Sgal(sLongGP,sLatGP,lGP[k],bGP[k]);
    }
  DrawRegion(lGP, bGP, drawLine);
}

void TProjMap::PlotSources(const vector<double>& lSources, const vector<double>& bSources, const vector<string> & nameSources, bool writeNameSources, int style, int color, double size)
{
  double xSources, ySources;
  char name[256];
  for(unsigned int i = 0; i < lSources.size(); i++)
    {
      gCoordFcnAng2XY(bSources[i]*DTOR, lSources[i]*DTOR, xSources, ySources, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
      xSources *= RTOD;
      ySources *= RTOD;
      if( IsInWindow(xSources, ySources) )
        {
          TGraph * plotSources = new TGraph(1,&xSources, &ySources);
          plotSources->SetMarkerStyle(style);
          plotSources->SetMarkerColor(color);
          plotSources->SetMarkerSize(size);
          plotSources->Draw("P");
          if( writeNameSources && nameSources.size() != 0 )
            {
              sprintf(name, "  %s", nameSources[i].c_str());
              TLatex * textPlotSources = new TLatex(xSources, ySources, name);
              textPlotSources->SetTextSize(0.04);
		          textPlotSources->SetTextFont(62);
              textPlotSources->SetLineWidth(1);
              textPlotSources->Draw();
            }
        }
    }
}


void TProjMap::PlotSource(double lSource, double bSource, string nameSource, bool writeNameSource, int style, int color, double size)
{
  double xSource, ySource;
  char name[256];
  gCoordFcnAng2XY(bSource*DTOR, lSource*DTOR, xSource, ySource, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
  xSource *= RTOD; ySource *= RTOD;
  if( IsInWindow(xSource, ySource) )
    {
      TGraph * plotSource = new TGraph(1,&xSource, &ySource);
      plotSource->SetMarkerStyle(style);
      plotSource->SetMarkerColor(color);
      plotSource->SetMarkerSize(size);
      plotSource->Draw("P");
      if( writeNameSource )
        {
          sprintf(name, "  %s", nameSource.c_str());
          TLatex * textPlotSource = new TLatex(xSource+0.45*size, ySource+0.45*size, name);
          textPlotSource->SetTextSize(0.04);
          textPlotSource->SetLineWidth(1);
          textPlotSource->SetTextFont(62);
          textPlotSource->SetTextColor(color);
          textPlotSource->Draw();
        }
    }
}



void TProjMap::PlotEvents(const vector<double>& lEvents, const vector<double>& bEvents, int style, int color, double size)
{
  double xEvents, yEvents;
  for(unsigned int i = 0; i < lEvents.size(); i++)
    {
      gCoordFcnAng2XY(bEvents[i]*DTOR, lEvents[i]*DTOR, xEvents, yEvents, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
      xEvents *= RTOD;
      yEvents *= RTOD;
      if( IsInWindow(xEvents, yEvents) )
        {
          TGraph * plotEvents = new TGraph(1,&xEvents, &yEvents);
          plotEvents->SetMarkerStyle(style);
          plotEvents->SetMarkerColor(color);
          plotEvents->SetMarkerSize(size);
          plotEvents->Draw("P");
        }
    }
}

void TProjMap::PlotEvents(const vector<double>& lEvents, const vector<double>& bEvents, const vector<double>& eEvents, int style, int color)
{
  double eMin, eMax;
  eMin = *min_element(eEvents.begin(),eEvents.end());
  eMax = *max_element(eEvents.begin(),eEvents.end());
  if( eMin == eMax ) eMin = eMax-1;
  double xEvents, yEvents;
  for(unsigned int i = 0; i < lEvents.size(); i++)
    {
      gCoordFcnAng2XY(bEvents[i]*DTOR, lEvents[i]*DTOR, xEvents, yEvents, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
      xEvents *= RTOD;
      yEvents *= RTOD;
      if( IsInWindow(xEvents, yEvents) )
        {
          TGraph * plotSources = new TGraph(1,&xEvents, &yEvents);
          plotSources->SetMarkerStyle(style);
          plotSources->SetMarkerColor(color);
          plotSources->SetMarkerSize(0.5+1.5*(eEvents[i]-eMin)/(eMax-eMin));
          plotSources->Draw("P");
        }
    }
}


void TProjMap::PlotCircle(const vector<double>& l, const vector<double>& b, double alpha,
			   int style, int color, int width)
{
  TDrawLine drawLine;
  drawLine.fNPts = 10000;
  drawLine.fStyle = style;
  drawLine.fColor = color;
  drawLine.fWidth = width;
  
  double cosAlpha = cos(alpha*DTOR), sinAlpha = sin(alpha*DTOR);
  double psi = 0.;
  vector<double> cosPsi(drawLine.fNPts), sinPsi(drawLine.fNPts);
  for(unsigned int i = 0; i < drawLine.fNPts; i++)
    {
      psi = TwoPi*i/(1.*drawLine.fNPts-1.);
      cosPsi[i] = cos(psi); 
      sinPsi[i] = sin(psi);
    }
  double theta = 0., phi = 0.;
  double x = 0., y = 0., z = 0.;
  vector<double> latPtsCercle(drawLine.fNPts), longPtsCercle(drawLine.fNPts);
  for(unsigned int i = 0; i < l.size(); i++)
    {
      theta = (90.-b[i])*DTOR; phi = l[i]*DTOR;
      for(unsigned int j = 0; j < drawLine.fNPts; j++)
        {
          x = cosAlpha*sin(theta)*cos(phi)+sinAlpha*(sin(phi)*sinPsi[j]-cos(theta)*cos(phi)*cosPsi[j]);
          y = cosAlpha*sin(theta)*sin(phi)-sinAlpha*(cos(phi)*sinPsi[j]+cos(theta)*sin(phi)*cosPsi[j]);
          z = cosAlpha*cos(theta)+sinAlpha*sin(theta)*cosPsi[j];
          xyz2thetaphi(x, y, z, latPtsCercle[j], longPtsCercle[j]);
          latPtsCercle[j] = (M_PI/2.-latPtsCercle[j])*RTOD;
          longPtsCercle[j] = longPtsCercle[j]*RTOD; 
        }
      DrawRegion(longPtsCercle,latPtsCercle,drawLine);
    }
}


void TProjMap::PlotMaxima(const vector<long>& ipAboveThreshold, int style, int color, double size)
{
  double xMaxima, yMaxima;
  vector<double> lMaxima, bMaxima;
  fMap.GiveLB(ipAboveThreshold, lMaxima, bMaxima);
  char number[256];
  for(unsigned int i = 0; i < ipAboveThreshold.size(); i++)
    {
      gCoordFcnAng2XY(bMaxima[i]*DTOR, lMaxima[i]*DTOR, xMaxima, yMaxima, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
      xMaxima *= RTOD;
      yMaxima *= RTOD;
      if( IsInWindow(xMaxima,yMaxima) )
        {
          TGraph* plotMaxima = new TGraph(1, &xMaxima, &yMaxima);
          plotMaxima->SetMarkerStyle(style);
          plotMaxima->SetMarkerColor(color);
          plotMaxima->SetMarkerSize(size);
          plotMaxima->Draw("P");
          sprintf(number,"   %i",i+1);

          TLatex * textMaxima = new TLatex(xMaxima, yMaxima-6., number);
          textMaxima->SetTextSize(0.02);
          textMaxima->SetLineWidth(1);
          textMaxima->Draw();
        }
    }
}


void TProjMap::Save(string fileName)
{
  fCanvas->Update();
  fCanvas->SaveAs(fileName.c_str());
}


void TProjMap::InitCanvas()
{
/*
  int * palette = GetPalette(gMaxPalette);
  gStyle->SetPalette(gMaxPalette,palette);
  delete [] palette;

  gStyle->SetPalette(53);
*/
  fCanvas = new TCanvas(GetObjName(), GetObjName(), (int)floor(1.3*fSizeX), (int)floor(1.3*fSizeY));
  fCanvas->SetFrameLineColor(0);
  fCanvas->SetLeftMargin(fXMarginLeft);
  fCanvas->SetRightMargin(fXMarginRight);
  fCanvas->SetTopMargin(fXMarginTop);
  fCanvas->SetBottomMargin(fXMarginBottom);
}

void TProjMap::MakeLongitudesLatitudes_0_360(double longStep, double latStep, const TDrawLine& drawLine)
{
  if( !fCanvas )
    {
      cout << "TProjMap::MakeLongitudesLatitudes. Call TProjMap::SkyMap first. Exiting." << endl;
      exit(0);
    }
  double latStart = -90, latStop = 90;
  double longStart = 0, longStop = 360;
  double dLat = latStop-latStart, dLong = longStop-longStart;
  int nLat = (int)floor(dLat/latStep)+1, nLong = (int)floor(dLong/longStep)+1;
  double sLat = dLat/(nLat-1), sLong = dLong/(nLong-1);
  double * longitudes = new double[nLong];
  double * latitudes = new double[nLat];
  double lambda0 = fLongCenter*DTOR;
  vector<string> latName(nLat);
  vector<string> longName(nLong+1); // +1 to solve the problem of the 180 to be drawn in each side of the map
  char tmp[128];
  double xProj, yProj;
  unsigned int count = 0;

  // Reference latitudes
  double * dLatRef = new double[nLat];
  int * dLatRefID = new int[nLat];
  bool isLatZero = false;
  for(int i = 0; i < nLat; i++)
    {
      latitudes[i] = (latStart+sLat*i)*DTOR;
      if( latitudes[i] == 0 ) isLatZero = true;
      sprintf(tmp,"%-10.0f",latitudes[i]*RTOD);
      latName[i] = tmp;
      dLatRef[i] = fabs(fLatCenter*DTOR-latitudes[i]);
    }
  TMath::Sort(nLat, dLatRef, dLatRefID, false);
  double latCenter = latitudes[dLatRefID[0]];
  delete [] dLatRef;
  delete [] dLatRefID;
  
  // Reference longitudes
  for(int i = 0; i < nLong; i++) longitudes[i] = (longStart+sLong*i)*DTOR;

  // Additionnal longitude to close the sky (at 180.001 degrees)
  double * longitudesComplete = new double[nLong+1];
  count = 0;
  for(int i = 0; i < nLong; i++)
    {
      longitudesComplete[count] = longitudes[i];
			double longi = longitudes[i];
			if(longi > M_PI) longi-=2.*M_PI;
      sprintf(tmp,"%-10.0f",longi*RTOD+180);
      longName[count] = tmp;
      count++;
      if( longitudes[i] == M_PI )
        {
          longitudesComplete[count] = 180.001*DTOR;
          longName[count] = "360";
          count++;
        }
    }
  if( latStart == -90 ) latitudes[0] = -89.9*DTOR;
  if( latStop == 90 ) latitudes[nLat-1] = 89.9*DTOR;
  
  // Draw isolatitude lines
  int align = 11;
  vector<double> longLine(drawLine.fNPts), latLine(drawLine.fNPts);
  for(unsigned int i = 0; i < drawLine.fNPts; i++) longLine[i] = longStart+dLong*i/(drawLine.fNPts-1);
  for(int i = 0; i < nLat; i++)
    {
      for(unsigned int j = 0; j < drawLine.fNPts; j++) latLine[j] = latitudes[i]/DTOR;
      DrawRegion(longLine, latLine, drawLine);
      gCoordFcnAng2XY(latitudes[i], lambda0, xProj, yProj, lambda0, fCosLatRef, fSinLatRef);
      xProj *= RTOD;
      yProj *= RTOD;
      TLatex * latText = new TLatex;
      if( IsInWindow(xProj, yProj) )
        {
          if( latitudes[i] == 0 && fMap.CoordSys() == 'G') latText->SetText(xProj, yProj, "GC");
          else latText->SetText(xProj, yProj, latName[i].c_str());
          latText->SetTextSize(0.04);
          latText->SetTextFont(62);
				  latText->SetTextColor(drawLine.fColor);
          latText->SetTextAlign(align);
          latText->Draw();
        }
    }
  
  // Draw galactic plane
  if( !isLatZero ) // Force galactic plane drawing
    {
      for(unsigned int i = 0; i < drawLine.fNPts; i++) longLine[i] = longStart+dLong*i/(drawLine.fNPts-1);
      for(unsigned int i = 0; i < drawLine.fNPts; i++) latLine[i] = 0; // GC
      DrawRegion(longLine, latLine, drawLine);
    }
  
  // Draw isolongitude lines
  longLine.resize(drawLine.fNPts+1);
  latLine.resize(drawLine.fNPts+1);
  for(unsigned int i = 0; i <= drawLine.fNPts; i++) latLine[i] = latStart+dLat*i/drawLine.fNPts;
  for(int i = 0; i <= nLong; i++)
    {
      for(unsigned int j = 0; j <= drawLine.fNPts; j++) longLine[j] = longitudesComplete[i]/DTOR;
      DrawRegion(longLine, latLine, drawLine);
      gCoordFcnAng2XY(latCenter, longitudesComplete[i], xProj, yProj, lambda0, fCosLatRef, fSinLatRef);
      xProj *= RTOD;
      yProj *= RTOD;
      TLatex * longText = new TLatex;
      if( IsInWindow(xProj,yProj) )
        {
          if( (int)floor(longitudes[i]*RTOD) != 0 )
            {
              longText->SetText(xProj, yProj, longName[i].c_str());
              longText->SetTextSize(0.04);
              longText->SetTextColor(drawLine.fColor);
              longText->SetTextFont(62);
              longText->SetTextAlign(align);
              longText->Draw();
            }
        }
    }
  delete [] longitudes;
  delete [] latitudes;
  delete [] longitudesComplete;
}

void TProjMap::MakeLongitudesLatitudes(double longStep, double latStep, const TDrawLine& drawLine)
{
  if( !fCanvas )
    {
      cout << "TProjMap::MakeLongitudesLatitudes. Call TProjMap::SkyMap first. Exiting." << endl;
      exit(0);
    }
  double latStart = -90, latStop = 90;
  double longStart = 0, longStop = 360;
  double dLat = latStop-latStart, dLong = longStop-longStart;
  int nLat = (int)floor(dLat/latStep)+1, nLong = (int)floor(dLong/longStep)+1;
  double sLat = dLat/(nLat-1), sLong = dLong/(nLong-1);
  double * longitudes = new double[nLong];
  double * latitudes = new double[nLat];
  double lambda0 = fLongCenter*DTOR;
  vector<string> latName(nLat);
  vector<string> longName(nLong+1); // +1 to solve the problem of the 180 to be drawn in each side of the map
  char tmp[128];
  double xProj, yProj;
  unsigned int count = 0;

  // Reference latitudes
  double * dLatRef = new double[nLat];
  int * dLatRefID = new int[nLat];
  bool isLatZero = false;
  for(int i = 0; i < nLat; i++)
    {
      latitudes[i] = (latStart+sLat*i)*DTOR;
      if( latitudes[i] == 0 ) isLatZero = true;
      sprintf(tmp,"%-10.0f",latitudes[i]*RTOD);
      latName[i] = tmp;
      dLatRef[i] = fabs(fLatCenter*DTOR-latitudes[i]);
    }
  TMath::Sort(nLat, dLatRef, dLatRefID, false);
  double latCenter = latitudes[dLatRefID[0]];
  delete [] dLatRef;
  delete [] dLatRefID;
  
  // Reference longitudes
  for(int i = 0; i < nLong; i++) longitudes[i] = (longStart+sLong*i)*DTOR;

  // Additionnal longitude to close the sky (at 180.001 degrees)
  double * longitudesComplete = new double[nLong+1];
  count = 0;
  for(int i = 0; i < nLong; i++)
    {
      longitudesComplete[count] = longitudes[i];
      sprintf(tmp,"%-10.0f",longitudes[i]*RTOD);
      longName[count] = tmp;
      count++;
      if( longitudes[i] == M_PI )
        {
          longitudesComplete[count] = 180.001*DTOR;
          longName[count] = "180";
          count++;
        }
    }
  if( latStart == -90 ) latitudes[0] = -89.9*DTOR;
  if( latStop == 90 ) latitudes[nLat-1] = 89.9*DTOR;
  
  // Draw isolatitude lines
  int align = 11;
  vector<double> longLine(drawLine.fNPts), latLine(drawLine.fNPts);
  for(unsigned int i = 0; i < drawLine.fNPts; i++) longLine[i] = longStart+dLong*i/(drawLine.fNPts-1);
  for(int i = 0; i < nLat; i++)
    {
      for(unsigned int j = 0; j < drawLine.fNPts; j++) latLine[j] = latitudes[i]/DTOR;
      DrawRegion(longLine, latLine, drawLine);
      gCoordFcnAng2XY(latitudes[i], lambda0, xProj, yProj, lambda0, fCosLatRef, fSinLatRef);
      xProj *= RTOD;
      yProj *= RTOD;
      TLatex * latText = new TLatex;
      if( IsInWindow(xProj, yProj) )
        {
          if( latitudes[i] == 0 && fMap.CoordSys() == 'G') latText->SetText(xProj, yProj, "GC");
          else latText->SetText(xProj, yProj, latName[i].c_str());
          latText->SetTextSize(0.04);
          latText->SetTextFont(62);
				  latText->SetTextColor(drawLine.fColor);
          latText->SetTextAlign(align);
          latText->Draw();
        }
    }
  
  // Draw galactic plane
  if( !isLatZero ) // Force galactic plane drawing
    {
      for(unsigned int i = 0; i < drawLine.fNPts; i++) longLine[i] = longStart+dLong*i/(drawLine.fNPts-1);
      for(unsigned int i = 0; i < drawLine.fNPts; i++) latLine[i] = 0; // GC
      DrawRegion(longLine, latLine, drawLine);
    }
  
  // Draw isolongitude lines
  longLine.resize(drawLine.fNPts+1);
  latLine.resize(drawLine.fNPts+1);
  for(unsigned int i = 0; i <= drawLine.fNPts; i++) latLine[i] = latStart+dLat*i/drawLine.fNPts;
  for(int i = 0; i <= nLong; i++)
    {
      for(unsigned int j = 0; j <= drawLine.fNPts; j++) longLine[j] = longitudesComplete[i]/DTOR;
      DrawRegion(longLine, latLine, drawLine);
      gCoordFcnAng2XY(latCenter, longitudesComplete[i], xProj, yProj, lambda0, fCosLatRef, fSinLatRef);
      xProj *= RTOD;
      yProj *= RTOD;
      TLatex * longText = new TLatex;
      if( IsInWindow(xProj,yProj) )
        {
          if( (int)floor(longitudes[i]*RTOD) != 0 )
            {
              longText->SetText(xProj, yProj, longName[i].c_str());
              longText->SetTextSize(0.04);
              longText->SetTextColor(drawLine.fColor);
              longText->SetTextFont(62);
              longText->SetTextAlign(align);
              longText->Draw();
            }
        }
    }
  delete [] longitudes;
  delete [] latitudes;
  delete [] longitudesComplete;
}


void TProjMap::DrawRegion(const vector<double>& longitudes, const vector<double>& latitudes, const TDrawLine& drawLine)
{
  // Line
  int nPts = int(longitudes.size());
  double ignore = -12345;
  double * xProj = new double[nPts+1];
  double * yProj = new double[nPts+1];
  xProj[nPts] = ignore; yProj[nPts] = ignore;
  int count(0);
  double xProjCurrent, yProjCurrent;
  double xProjPrevious(0), yProjPrevious(0);
  bool first = true;
  bool jump  = false;
  for(int k = 0; k < nPts; k++)
    {
      gCoordFcnAng2XY(latitudes[k]*DTOR,longitudes[k]*DTOR, xProjCurrent, yProjCurrent, fLongCenter*DTOR, fCosLatRef, fSinLatRef);
      xProjCurrent *= RTOD;
      yProjCurrent *= RTOD;
      if( first ) {first = false; jump = true;}
      else
        {
          jump = false;
          double dXProj = xProjCurrent-xProjPrevious, dYProj = yProjCurrent-yProjPrevious;
          double distance = sqrt(dXProj*dXProj+dYProj*dYProj);
          if( distance > 1 ) jump = true;
        }
      if( IsInWindow(xProjCurrent,yProjCurrent) && !jump )
        {
          xProj[k] = xProjCurrent;
          yProj[k] = yProjCurrent;
          count++;
        }
      else // mark discontinuity
        {
          xProj[k] = ignore;
          yProj[k] = ignore;
        }
      xProjPrevious = xProjCurrent;
      yProjPrevious = yProjCurrent;
    }
  if( count ) // some points are visible
    {
      if( count == nPts ) // all points are in the window
        {
          TPolyLine * region = new TPolyLine(nPts);
          region->SetPolyLine(nPts,xProj,yProj);
          region->SetLineWidth(drawLine.fWidth);
          region->SetLineColor(drawLine.fColor);
          region->SetLineStyle(drawLine.fStyle);
          region->Draw();
        }
      else // incomplete band
        {
          vector< vector<unsigned int> > bands = GetContiguousPoints(nPts,xProj,ignore);
          for(unsigned int i = 0; i < bands.size(); i++)
            {
              // draw lines one after the other
              TPolyLine * incompleteRegion = new TPolyLine;
              unsigned int size = bands[i].size();
              double * thex = new double[size];
              double * they = new double[size];
              for(unsigned int j = 0; j < size; j++)
                {
                  thex[j] = xProj[bands[i][j]];
                  they[j] = yProj[bands[i][j]];
                }
              incompleteRegion->SetPolyLine(size,thex,they);
              incompleteRegion->SetLineWidth(drawLine.fWidth);
              incompleteRegion->SetLineColor(drawLine.fColor);
              incompleteRegion->SetLineStyle(drawLine.fStyle);
              incompleteRegion->Draw();
              
              delete [] thex;
              delete [] they;
            }
        }
    }
}


bool TProjMap::IsInWindow(double xVal, double yVal)
{
  bool status = false;
  if( xVal <= fXMax && xVal >= fXMin && yVal <= fYMax && yVal >= fYMin ) status = true;
  return status;
}
