/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TNtuple.h>
#include <TObjString.h>
#include <TString.h>
#include <TStyle.h>

#include "AliFieldReader.h"
#include "AliMagFMaps.h"

ClassImp(AliFieldReader)



//_______________________________________________________________________
AliFieldReader::AliFieldReader():
    fField(0),
    fMap(0),
    fCatalogue(0),
    fHtmlMain(0),
    fStepSize(0.),
    fZStart(1383.),
    fDd(0.08),
    fDz(0.064),
    fPolarity(1.),
    fCatalogueName("goodfiles.list")
{
//
//  Constructor
//
}

AliFieldReader::AliFieldReader(const AliFieldReader& reader):
    TObject(reader),
    fField(0),
    fMap(0),
    fCatalogue(0),
    fHtmlMain(0),
    fStepSize(0.),
    fZStart(0.),
    fDd(0.),
    fDz(0.),
    fPolarity(0.),
    fCatalogueName(0)
{
//
//  Constructor
//
    reader.Copy(*this);
}

//_______________________________________________________________________
AliFieldReader::~AliFieldReader()
{
  //
  // Destructor
  //
}

//_______________________________________________________________________
void AliFieldReader::Init()
{   
//
// Initialize the reader
//
    // Calculated map
    fField = new AliMagFMaps("Maps","Maps", 2, 1., 10., 2);
    // Catalogue
    fCatalogue = fopen(fCatalogueName, "r");
    
    // HTML
    fHtmlMain = fopen("bmap.html", "w");
    MakeHtmlHeaderMain(fHtmlMain);
    
}


void AliFieldReader::ReadMap()
{   
//
// Read the measured dipole field map 
//    
    
    Float_t zA[450], bxzA[200], byzA[200], bzzA[200], bxzcA[200], byzcA[200], bzzcA[200];
    Float_t yA[450], bxyA[200], byyA[200], bzyA[200], bxycA[200], byycA[200], bzycA[200];
    Float_t xA[450], bxxA[200], byxA[200], bzxA[200], bxxcA[200], byxcA[200], bzxcA[200];
    
    Char_t sLine[255];
    Char_t fLine[255];
    
    Float_t xpos, ypos, zpos;
    Int_t iboxpos;
    
    Float_t h[3], x[3], b[3];
    Float_t temp;
    Float_t xnt[17];
    Int_t ires;
    
    Int_t iret, ipos;
    Int_t ic115 = 0;
    Float_t  dx = 0.;
    Float_t  dy = 0.;
    
    Init();
    // Read register
    ReadRegisterMap();
    // n-tuple
    fMap = new TNtuple("Field Map", "Map", 
		       "x:y:z:ix:iy:iz:bx:by:bz:bxc:byc:bzc:ifile:ireg:irot:temp:cal", 4000);
    
    //
    // Loop over files
    // 
    Int_t ifiles = 0;
    // File catalogue
    while ((fgets(fLine, 255, fCatalogue)) != NULL && ifiles <= 2000) {

	if (strncmp(fLine,"#HEADER", 7) == 0) {
	  iret = sscanf(&fLine[7],"%f", &fZStart);
	  continue;
	}	
	ifiles++;
	
//	if (ifiles != 87) continue;
	char fileName[32];
	iret = sscanf(fLine, "%s", fileName);
	printf("Reading File %s\n", fileName);
	
//  Get run name
	TString* tsFile          = new TString(fileName);
	TObjArray * tokens       =  tsFile->Tokenize(".");
	const char* runName      = (((TObjString*) tokens->At(0))->GetString()).Data();
	FILE* file               =  fopen(fileName, "r");
	
	
	Float_t bdl   = 0.;
	Float_t bdlc  = 0.;
	Int_t iz      = 0;
	Int_t izz     = 0;
	Int_t iy      = 0;
	Int_t ix      = 0;
	Int_t iheader = 0;
	Float_t stepsz = 0.;
    // Graphs go here
	TMultiGraph* bxzmg  = new TMultiGraph("bxzmg", "B_{x}");
	TMultiGraph* byzmg  = new TMultiGraph("byzmg", "B_{y}");
	TMultiGraph* bzzmg  = new TMultiGraph("bzzmg", "B_{z}");
	while ((fgets(sLine, 255, file)) != NULL) {
// read z-position
	    
	    if (strncmp(sLine," Current", 8) == 0) {
		Float_t current;
		iret = sscanf(&sLine[9],"%f", &current);
		printf("Current %f \n", current);
		
		fPolarity = current / 6000.;
	    }
	    if (strncmp(sLine," z", 2) == 0) {
		Float_t zmin, zmax;
		Int_t nsteps;
		iret = sscanf(&sLine[3],"%f %f %d", &zmin, &zmax, &nsteps);
		printf("zmin zmax %13.3f %13.3f %13.3f\n", 
		       zmin, zmax, TMath::Abs(zmax - zmin)/Float_t(nsteps));
		zmax = TMath::Max(zmin, zmax);
		stepsz  = TMath::Abs(zmax - zmin)/Float_t(nsteps);
	    }
	    
	    if (strncmp(sLine," X position", 11) == 0)
	    {
		TString string;
		TString* tsLine = new TString(sLine);
		TObjArray * tokens =  tsLine->Tokenize("=");	
		string = ((TObjString*) tokens->At(1))->GetString();
		iret = sscanf(string.Data(), "%f", &xpos);
		string = ((TObjString*) tokens->At(2))->GetString();
		iret = sscanf(string.Data(), "%f", &ypos);
		string = ((TObjString*) tokens->At(3))->GetString();
		iret = sscanf(string.Data(), "%d", &iboxpos);
		
		
		printf("This file is for  x = %13.3f y = %13.3f Box Position %1d\n", xpos, ypos, iboxpos);
		iheader++;
		
	    }
	    
	    if (strncmp(sLine," R	Z-pos",8) == 0)
	    {
		iret = sscanf(&sLine[8],"%e", &zpos);
		ipos  = 0;
		ic115 = 0;
		izz++;
	    }
	    
	    
	    if (strncmp(sLine,"C",1) == 0) {
		ipos++;
		Int_t ireg;
		iret = sscanf(&sLine[2],"%d %f %f %f %f %d", &ireg, &h[0], &h[1], &h[2], &temp, &ires);
		//
		// fix for address 115
		//
		if (ireg == 115) ic115++;
		if (ic115 == 2 && ireg == 115) ireg = 119;
		
		if (iheader == 1) {

		    Float_t bx = 0., by = 0., bz = 0.;
		    Int_t jx = fRegMap[ireg][0];
		    Int_t jy = fRegMap[ireg][1];
		    Int_t jz = fRegMap[ireg][2];

		    switch (iboxpos) {
		    case 0:
			bx =  h[1];
			by =  h[2];
			dx = -0.36 + jx * fDd;
			dy =  jy * fDd;
			break;
		    case 1:
			bx = -h[2];
			by =  h[1];
			dx =  0.36 + jy * fDd;
			dy = -(-0.36 + jx * fDd);
			break;
		    case 2:
			bx = -h[1];
			by = -h[2];
			dx =  0.36 - jx * fDd;
			dy = -jy * fDd;
			break;
		    case 3:
			bx =   h[2];
			by =  -h[1];
			dx =  -0.36 - jy * fDd;
			dy =   -(0.36 - jx * fDd);
			break;
		    } // switch

		    bz =  h[0];
		    if (jz == 0) {
			bz = -bz;
			if (iboxpos == 1 || iboxpos == 3) {
			    by = -by;
			} else {
			    bx = -bx;
			}
		    }
		    
		    
		    Float_t dz = (jz == 0)?  fDz/2. : - fDz/2.;
		    dz *= 100.;
		    

		    Float_t xc = xpos + dx;
		    Float_t yc = ypos + dy;

		    
		    x[0] =  - xc * 100.;
		    x[1] =  + yc * 100.;
		    x[2] =  - (-zpos * 100. + fZStart + dz);

		    fField->Field(x, b);
		    b[0] *= fPolarity;
		    b[1] *= fPolarity;
		    b[2] *= fPolarity;
		    
		    xnt[ 0] = xc;
		    xnt[ 1] = yc;
		    xnt[ 2] = -x[2] / 100.;
		    xnt[ 3] = Float_t (jx);
		    xnt[ 4] = Float_t (jy);
		    xnt[ 5] = Float_t (jz);
		    xnt[ 6] = bx;
		    xnt[ 7] = by;
		    xnt[ 8] = bz;
		    xnt[ 9] = b[0]/10.;
		    xnt[10] = b[1]/10.;
		    xnt[11] = b[2]/10.;		    
		    xnt[12] = Float_t (ifiles);
		    xnt[13] = Float_t (ireg);
		    xnt[14] = Float_t(iboxpos);
		    xnt[15] = temp;
		    xnt[16] = Float_t(ires);
		    
		    fMap->Fill(xnt);
		    

//
// Calculated field		

		    if (jy != -1 && jz == 1 && jy == 0 && izz == 40){
			x[1] = x[1] + 10.;
			fField->Field(x, b);
			yA  [jx]  = yc;
			bxyA[jx]  = bx;
			byyA[jx]  = by;
			bzyA[jx]  = bz;
			bxycA[jx] =   b[0] / 10.;
			byycA[jx] =   b[1] / 10.;
			bzycA[jx] =   b[2] / 10.;
			iy++;
		    } // if
		    
		    if (jy != -1 && jz == 1 && jy == 0 && izz == 1){
			x[1] = x[1] + 10.;
			fField->Field(x, b);
			xA  [jx]  = xc;
			bxxA[jx]  = bx;
			byxA[jx]  = by;
			bzxA[jx]  = bz;
			bxxcA[jx] =   b[0] / 10.;
			byxcA[jx] =   b[1] / 10.;
			bzxcA[jx] =   b[2] / 10.;
			ix++;
		    } // if


		    if (ireg == 181) {
			fField->Field(x, b);
//			printf("Field %f %f %f %f %f %f \n", x[0], x[1], x[2], b[0], b[1] , b[2]);
			
			bdl  +=  stepsz * bx;
			bdlc +=  stepsz * b[0];
			
			zA [iz]  = x[2];
			bxzA[iz]  = bx;
			byzA[iz]  = by;
			bzzA[iz]  = bz;
			bxzcA[iz] =   b[0] / 10.;
			byzcA[iz] =   b[1] / 10.;
			bzzcA[iz] =   b[2] / 10.;
			iz++;
		    }
		} // if 1st header 
	    } // if C
	} // next line
	gStyle->SetOptStat(0);
	char title[128];
	sprintf(title, "File#: %5d, X = %13.2f m , Y = %13.2f m, Box Orientation: %2d", ifiles, xpos, ypos, iboxpos);
	TCanvas* c2 = new TCanvas("c2", title, 1200, 800);
	c2->Divide(2,2);
	c2->cd(1);
	TGraph* bxg =  new TGraph(iz, zA, bxzA);
	TGraph* bxcg = new TGraph(iz, zA, bxzcA);    
	bxcg->SetLineColor(2);
	bxzmg->Add(bxg);
	bxzmg->Add(bxcg);
	bxzmg->Draw("lA");
	bxzmg->GetHistogram()->SetXTitle("z[m]");
	bxzmg->GetHistogram()->SetYTitle("B_{x} [T]");
	bxzmg->Draw("lA");
	
	c2->cd(2);
	TGraph* byg =  new TGraph(iz, zA, byzA);
	TGraph* bycg = new TGraph(iz, zA, byzcA);    
	bycg->SetLineColor(2);
	byzmg->Add(byg);
	byzmg->Add(bycg);
	byzmg->Draw("lA");
	byzmg->GetHistogram()->SetXTitle("z[m]");
	byzmg->GetHistogram()->SetYTitle("B_{y} [T]");
	byzmg->Draw("lA");
	
	c2->cd(3);
	TGraph* bzg =  new TGraph(iz, zA, bzzA);
	TGraph* bzcg = new TGraph(iz, zA, bzzcA);    
	bzcg->SetLineColor(2);
	bzzmg->Add(bzg);
	bzzmg->Add(bzcg);
	bzzmg->Draw("lA");
	bzzmg->GetHistogram()->SetXTitle("z[m]");
	bzzmg->GetHistogram()->SetYTitle("B_{z} [T]");
	bzzmg->Draw("lA");
	
	c2->Update();
	char pictFile[64];
	sprintf(pictFile, "%s.gif", runName);
	c2->SaveAs(pictFile);
	//
	// Html generation
	//
	char htmlFile[64];
	sprintf(htmlFile, "%s.html", runName);
	FILE* chtml = fopen(htmlFile, "w");
	MakeHtmlHeaderPict(chtml);
	MakeHtmlPict(chtml, pictFile);
	MakeHtmlTableEntry(fHtmlMain, fileName, htmlFile, xpos, ypos, iboxpos, bdl, ifiles);

	//
	//
	printf("Bdl [Tm] %f %f \n", 2. * bdl, 2 * bdlc / 10.);
    } // files
    MakeHtmlTrailor(fHtmlMain);
    TFile* out = new TFile("fmap.root", "recreate");
    fMap->Write();
    out->Close();
}

void AliFieldReader::ReadMapSolenoid(){
//
//  Read map for solenoid measurement
// 
    Float_t phiA[450], bzPhiA[200], brPhiA[200], btPhiA[200], bbPhiA[200];
    Float_t bzcPhiA[200], brcPhiA[200], btcPhiA[200], bbcPhiA[200];    
    Char_t sLine[255];
    Char_t fLine[255];
    
    Float_t zpos, phipos, skewing, temp;
    Int_t ical;
    Float_t zmin, zmax;
    Float_t h[3], x[3], b[3];
    Int_t iret;
    Int_t ipos = 0;
    
    Init();
    ReadRegisterMapSolenoid();
    fMap = new TNtuple("Field Map", "Map", 
		       "r:phi:z:br:bt:bz:brc:btc:bzc:ifile:ireg:temp:cal:arm", 4000);
    
    //
    // Loop over files
    // 
    Int_t ifiles = 0;
    // File catalogue
    while ((fgets(fLine, 255, fCatalogue)) != NULL && ifiles <= 2000) {

	if (strncmp(fLine,"#HEADER", 7) == 0) {
	  iret = sscanf(&fLine[7],"%f", &fZStart);
	  continue;
	}	
	ifiles++;
	
	char fileName[32];
	iret = sscanf(fLine, "%s", fileName);
	printf("Reading File %s\n", fileName);
	
//  Get run name
	TString* tsFile          = new TString(fileName);
	TObjArray * tokens       =  tsFile->Tokenize(".");
	Int_t n = tokens->GetEntries();
	char* runName = new  char[256]; 
	sprintf(runName, "%s", (((TObjString*) tokens->At(0))->GetString()).Data());
	if (n > 2) {
	    for (Int_t i = 1; i < n-1; i++)
	    {
		sprintf(runName, "%s.%s", 
			runName, (((TObjString*) tokens->At(i))->GetString()).Data());
	    }
	}
	FILE* file               =  fopen(fileName, "r");
	
	
	Float_t bdl  = 0.;
	Float_t bdlc = 0.;
	Int_t iA   = 0;
	Int_t izz  = 0;
	Int_t iphi = 0;
    // Graphs go here
	TMultiGraph* bxzmg  = new TMultiGraph("bxzmg", "B_{z}");
	TMultiGraph* byzmg  = new TMultiGraph("byzmg", "B_{r}");
	TMultiGraph* bzzmg  = new TMultiGraph("bzzmg", "B_{t}");
	TMultiGraph* bbmg   = new TMultiGraph("bbmg", "|B|");

	while ((fgets(sLine, 255, file)) != NULL) {
	    if (strncmp(sLine," z", 2) == 0) {
		Int_t nsteps;
		iret = sscanf(&sLine[3],"%f %f %d", &zmin, &zmax, &nsteps);
		printf("zmin zmax %13.3f %13.3f %13.3f\n", 
		       zmin, zmax, TMath::Abs(zmax - zmin)/Float_t(nsteps));
	    }
	    if (strncmp(sLine," R\tPOSITION NUMBER", 18) == 0)
	    {	    	      
		//
		// Current z-position
		TString string;
		TString* tsLine = new TString(sLine);
		TObjArray * tokens =  tsLine->Tokenize("=");	
		string = ((TObjString*) tokens->At(1))->GetString();
		iret = sscanf(string.Data(), "%f", &zpos);
		printf("POSITION NUMBER Z: %f\n", zpos);
		izz ++;
		delete tsLine;
	    }


	    if (strncmp(sLine," SKEWING ON Z:", 14) == 0)
	    {	    	      
		//
		// Skew in z
	      iret = sscanf(&sLine[14],"%f", &skewing);
	      printf("SKEWING ON Z: %f\n", skewing);
	    }


	    if (strncmp(sLine,"Phi", 3) == 0)
	    {
	      Float_t phiStart, phiStop;
	      
	      iret = sscanf(&sLine[3],"%f %f", &phiStart, &phiStop);

	      printf("phiStart phiStop %f %f\n", phiStart, phiStop);

	    }
	    

	    if (strncmp(sLine," R\tPhi-Angle",12) == 0)
	    {
		iret = sscanf(&sLine[12],"%e", &phipos);
		ipos = 0;
		iphi++;
	    }
	    

	    if (strncmp(sLine,"C",1) == 0) {
	
		Int_t ireg;
		iret = sscanf(&sLine[2],"%d %f %f %f %f %d", &ireg, &h[0], &h[1], &h[2], &temp, &ical);
		ipos++;

		Int_t ir = fRegMap[ireg][0];
		Int_t ia = fRegMap[ireg][1];
		Float_t rpos = 0.;

		
		if (ia == 0) {
		    rpos    = 0.2295 + ir * 0.16;
		} else {
		    if (ireg ==  81) rpos = 0.2295;
		    if (ireg ==  59) rpos = 1.0295;
		    if (ireg == 142) rpos = 2.1495;		    
		    if (ireg == 180) rpos = 3.1095;		    
		    if (ireg ==  69) rpos = 4.2295;		    

		    //		    if (ireg ==  55) rpos = 0.2295;
		    //if (ireg == 195) rpos = 1.0295;
		    //if (ireg == 129) rpos = 2.1495;		    
		    //if (ireg == 167) rpos = 3.1095;		    
		    //if (ireg == 142) rpos = 4.2295;		    
		}
		
		Float_t phi = phipos;
		if (ia == 1) {
		    phi += 180.;
		    if (phi > 360.) phi -= 360.;			    
		}
		
		phi = - phi * TMath::Pi() / 180.;
		Float_t xpos = rpos * TMath::Cos(phi);
		Float_t ypos = rpos * TMath::Sin(phi);
		x[0] = - xpos * 100.;
		x[1] = ypos * 100.;
		x[2] = -400. + zpos * 100.;
		
		fField->Field(x, b);
		Float_t phi0 = TMath::Pi() - phi;
		
		Float_t brc =   b[0] * TMath::Cos(phi0) +  b[1] * TMath::Sin(phi0);
		Float_t btc = - b[0] * TMath::Sin(phi0) +  b[1] * TMath::Cos(phi0);	
		Float_t bzc =   b[2];
		
		brc /= 10.;
		btc /= 10.;
		bzc /= 10.;
	
		fMap->Fill(rpos, -phi, -(615.5 - fZStart) / 100. + zpos, h[2], -h[1], h[0], brc, btc, bzc, ifiles, ireg, temp, Float_t(ical), Int_t(ia));
	
		if (ireg  == 174) {
		    printf("Field (Bx, By, Bz) at position %d: %13.3f %13.3f %13.3f %13.3f %13.3f %13.3f %5d\n", 
			   ipos, zpos, rpos, phipos, h[0], h[1], h[2], ia);	
		    if (izz == 1) {
			phiA[iA]   = phi;
			bzPhiA[iA] = -h[0];
			brPhiA[iA] =  h[2];
			btPhiA[iA] = -h[1];
			bbPhiA[iA] = TMath::Sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
			bzcPhiA[iA] = bzc;
			brcPhiA[iA] = brc;
			btcPhiA[iA] = btc;
			bbcPhiA[iA] = TMath::Sqrt(brc * brc + bzc * bzc + btc * btc);
			iA++;
		    }
		}
	    } // if R
	} // next line

	gStyle->SetOptStat(0);
	char title[128];
	sprintf(title, "Z = %13.2f m , Phi = %13.2f m", zpos, phipos);


	TCanvas* c2 = new TCanvas("c2", title, 1200, 800);
		
	c2->Divide(2,2);
	c2->cd(1);
	TGraph* bzg  =  new TGraph(iA, phiA, bzPhiA);
	TGraph* bzcg =  new TGraph(iA, phiA, bzcPhiA);
	bzcg->SetLineColor(2);
	bxzmg->Add(bzg);
	bxzmg->Add(bzcg);
	bxzmg->Draw("lA");
	bxzmg->GetHistogram()->SetXTitle("#phi[rad]");
	bxzmg->GetHistogram()->SetYTitle("B_{z} [T]");
	bxzmg->Draw("lA");

	c2->cd(2);
	TGraph* brg  =  new TGraph(iA, phiA, brPhiA);
	TGraph* brcg =  new TGraph(iA, phiA, brcPhiA);
	brcg->SetLineColor(2);
	byzmg->Add(brcg);
	byzmg->Add(brg);

	byzmg->SetMaximum(0.03);
	byzmg->SetMinimum(-0.03);	
	byzmg->Draw("lA");
	byzmg->GetHistogram()->SetXTitle("#phi[rad]");
	byzmg->GetHistogram()->SetYTitle("B_{r} [T]");
	byzmg->Draw("lA");

	c2->cd(3);
	TGraph* btg  =  new TGraph(iA, phiA, btPhiA);
	TGraph* btcg =  new TGraph(iA, phiA, btcPhiA);
	btcg->SetLineColor(2);
	bzzmg->Add(btcg);
	bzzmg->Add(btg);
	bzzmg->Draw("lA");
	bzzmg->SetMaximum(0.03);
	bzzmg->SetMinimum(-0.03);	
	bzzmg->GetHistogram()->SetXTitle("#phi[rad]");
	bzzmg->GetHistogram()->SetYTitle("B_{t} [T]");
	bzzmg->Draw("lA");
	

	c2->cd(4);
	TGraph* bg  =  new TGraph(iA, phiA, bbPhiA);
	TGraph* bcg =  new TGraph(iA, phiA, bbcPhiA);
	bcg->SetLineColor(2);
	bbmg->Add(bg);
	bbmg->Add(bcg);
	bbmg->Draw("lA");
	bbmg->GetHistogram()->SetXTitle("#phi[rad]");
	bbmg->GetHistogram()->SetYTitle("|B| [T]");
	bbmg->Draw("lA");
	



	char pictFile[64];
	sprintf(pictFile, "%s.gif", runName);
	c2->SaveAs(pictFile);

	//
	// Html generation
	//
	char htmlFile[64];
	sprintf(htmlFile, "%s.html", runName);
	FILE* chtml = fopen(htmlFile, "w");
	MakeHtmlHeaderPict(chtml);
	MakeHtmlPict(chtml, pictFile);
	MakeHtmlTableEntry(fHtmlMain, fileName, htmlFile, zmin, zmax, ifiles, 0., 0);

	//
	//
	printf("Bdl [Tm] %f %f \n", 2. * bdl, 2 * bdlc / 10.);
    } // files
    MakeHtmlTrailor(fHtmlMain);
    TFile* out = new TFile("fmap.root", "recreate");
    fMap->Write();
    out->Close();
}



void AliFieldReader::MakeHtmlHeaderMain(FILE* file)
{
//
//  Write the header of the heml output
//
    fprintf(file,"<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
    fprintf(file, "<html>\n");
    fprintf(file, "<head>\n");
    fprintf(file, "<meta http-equiv=\"content-type\"\n");
    fprintf(file, "content=\"text/html; charset=ISO-8859-1\"\n");
    fprintf(file, "<title>main.html</title>\n");
    fprintf(file, "<\\head>\n");
    fprintf(file, "<body>\n");
    fprintf(file, "<table cellpadding=\"1\" cellspacing=\"1\" border=\"1\"\n");
    fprintf(file, "style=\"text-align: left; width: 80;\">\n");
    fprintf(file, "<tbody>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">File#  <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">File Name  <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">X-Position <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">Y-Position <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">Box Orientation <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">B.dl <br></td>\n");
    fprintf(file, "<td style=\"vertical-align: top;\">Link to plots <br></td>\n");

}

void  AliFieldReader::MakeHtmlHeaderPict(FILE* file)
{
//
//  Write header for picture
//
    fprintf(file,"<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
    fprintf(file, "<html>\n");
    fprintf(file, "<head>\n");
    fprintf(file, "<meta http-equiv=\"content-type\"\n");
    fprintf(file, "content=\"text/html; charset=ISO-8859-1\"\n");
    fprintf(file, "<title>main.html</title>\n");
    fprintf(file, "</head>\n");
    fprintf(file, "<body>\n");
}

void AliFieldReader:: MakeHtmlPict(FILE* chtml, char* pictFile)
{
//
//  Write html for including picture
//
    fprintf(chtml, "<img src=\"./%s\" alt=\"%s\" style=\"width: 1196px; height: 772px;\">\n", 
	    pictFile, pictFile);
    
    fprintf(chtml, "<br> <br> <br>\n");	
    fprintf(chtml, "<a href=\"./bmap.html\">Back to main page</a><br>\n");
    
    fprintf(chtml, "</body>\n");
    fprintf(chtml, "</header>\n");	
    fclose(chtml);
}

void AliFieldReader::MakeHtmlTableEntry(FILE* htmlmain, char* fileName, char* htmlFile, Float_t x, Float_t y, Int_t i, Float_t bdl, Int_t ifile)
{	
    fprintf(htmlmain, "<tr>\n");
//
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%5d</td>\n", ifile);
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%s</td>\n", fileName);

//
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%13.2f</td>\n",x);
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%13.2f</td>\n",y);
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%3d   </td>\n",i);
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">%13.3f</td>\n",bdl);
//
    fprintf(htmlmain, "<td style=\"vertical-align: top;\">\n");
    fprintf(htmlmain, "<span style=\"text-decoration: underline;\"><a href=\"./%s\">gif</a></span></td>\n", htmlFile);
    fprintf(htmlmain, "</tr>\n");
}


void  AliFieldReader::MakeHtmlTrailor(FILE* htmlmain)
{
//
//  Write the html trailor
//
    fprintf(htmlmain, "</tbody>\n");
    fprintf(htmlmain, "</table>\n");	
    fprintf(htmlmain, "</body>\n");
    fprintf(htmlmain, "</header>\n");	
    fclose(htmlmain);
}

void AliFieldReader::ReadRegisterMap()
{
//
//  Read the register map
//
    FILE* regmap = fopen("register.map", "r");
    Int_t ireg;
    for (ireg = 0; ireg < 200; ireg++) {
	fRegMap[ireg][0] = -1;
	fRegMap[ireg][1] = -1;
	fRegMap[ireg][2] = -1;
    }
    for (Int_t iz = 0; iz < 2; iz++) {
	for (Int_t iy = 2; iy >= 0; iy--) {
	    for (Int_t ix = 0; ix < 10; ix++) {
		fscanf(regmap, "%d\n", &ireg);
		printf("Address %5d %5d %5d %5d \n", iz, iy, ix, ireg);
		    fRegMap[ireg][1] = iy;
		    fRegMap[ireg][2] = iz;
		if (iz == 1) {
		    fRegMap[ireg][0] = ix;
		} else {
		    fRegMap[ireg][0] = 9 - ix;
		}
	    } // ix
	} // iy
    } // iz
    fclose(regmap);
    printf("-> ReadRegisterMap()\n\n");
}


void AliFieldReader::ReadRegisterMapSolenoid()
{
//
//  Read the register map
//
    FILE* regmap = fopen("register.map", "r");
    Int_t ireg;
    
// Initialize
    for (ireg = 0; ireg < 200; ireg++) {
	fRegMap[ireg][0] = -1;
	fRegMap[ireg][1] = -1;
	fRegMap[ireg][2] = -1;
    }

// Main arm 
    for (Int_t ir = 0; ir < 33; ir++) {
	fscanf(regmap, "%d\n", &ireg);
	fRegMap[ireg][0] = ir;
	fRegMap[ireg][1] = 0;	
    }
// Opposite arm 
    for (Int_t ir = 0; ir < 5; ir++) {
	fscanf(regmap, "%d\n", &ireg);
	fRegMap[ireg][0] = ir;
	fRegMap[ireg][1] = 1;	
    }
    
    fclose(regmap);
    
}


AliFieldReader& AliFieldReader::operator=(const  AliFieldReader& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}


void AliFieldReader::Copy( TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}
