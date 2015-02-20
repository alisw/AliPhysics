#ifndef __CINT__

#include "TSystem.h"
#include "TGeoManager.h"
#include <AliRunLoader.h>
#include <AliGeomManager.h>
#include <AliITSgeom.h>
#include <AliTracker.h>
#include "TRandom.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <AliRun.h>
#include <TLegend.h>

#include <TGeoVolume.h>
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include <TGeoTube.h>
#include <AliITSUGeomTGeo.h>
#include <TPaveText.h>
#include <TText.h>

#endif

/*
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");

  .L GetMaterialBudget.C

  DrawMaterialBudget_SPLITTED(); 
  DrawMaterialBudget_FromTo();

  GetMaterialBudget_EtaVsPhi(0,0);
  GetMaterialBudget_inPhi(0,0);

*/

// Original version
// Chinorat Kobdaj <kobdaj@g.sut.ac.th> , Stefan Rossegger
// Revised and adapted to newer AliITSUv1 versions
// Mario Sitta <sitta@to.infn.it> - Nov 2014, Feb 2015


enum {
    kIBModelDummy=0,
    kIBModel0=1,
    kIBModel1=2, 
    kIBModel21=3,
    kIBModel22=4,
    kIBModel3=5,
    kIBModel4=10,
    kOBModelDummy=6,
    kOBModel0=7,
    kOBModel1=8, 
    kOBModel2=9 
};


//______________________________________________________________________
void PrintMaterialDefs(const char *geofile="geometry.root")
{
  //
  // Print material definition for all ITS materials/mediums
  // Rewritten - M.Sitta 15 Nov 2014
  //
  TGeoManager::Import(geofile);

  TGeoMedium *med;
  TGeoMaterial *mat;
  char mediaName[50], matName[50], shortName[50];

  Int_t nMedia = gGeoManager->GetListOfMedia()->GetEntries();

  printf("\n\n\n");
  printf("              ==== ITSU  Material  Properties ====\n\n");
  printf("    A      Z   d (g/cm3)  RadLen (cm)  IntLen (cm)\t Name\n");

  // Loop on media, select ITS materials, print their characteristics
  for (Int_t i = 0; i < nMedia; i++) {
    med = (TGeoMedium *)(gGeoManager->GetListOfMedia()->At(i));
    strcpy(mediaName, med->GetName());
    // The trick here: the name must begin with the string "ITS_"
    // so the pointer to the first occurrence of the substring as returned
    // by strstr() must match the beginning of the string
    if (strstr(mediaName,"ITS_") == mediaName) { // name begins with "ITS_"
      mat = med->GetMaterial();
      strcpy(matName, mat->GetName());
      strncpy(shortName, matName+4, strlen(matName)-5); // get rid of ITS_ , $
      shortName[strlen(matName)-5] = '\0';
      printf(" %5.1f %6.1f %8.3f %13.3f %11.3f\t %s\n", 
	     mat->GetA(), mat->GetZ(), mat->GetDensity(), mat->GetRadLen(),
	     mat->GetIntLen(), shortName);
    }
  }

  printf("\n");

//   printf("Values in Corrados estimate");
//   const Int_t nC=6;
//   TString medNC[nC]={"ITS_SI$","ITS_WATER$","ITS_CARBON$","ITS_KAPTON$","ITS_FLEXCABLE$","ITS_GLUE$"};
//   Double_t radl[nC]={9.36,36.1,25,28.6,13.3,44.37};
//   for (Int_t i=0; i<nC; i++) {
//     printf("%8.2lf (cm)  \t%s\n", radl[i],medNC[i].Data());
//   }
  
}


//______________________________________________________________________
Double_t ComputeMaterialBudget(Double_t rmin, Double_t rmax, Double_t etaPos,
			       Double_t phiMax,
			       TH1F* xOverX01d, TH2F* xOverX02d = 0)
{
  //
  // Ancillary function to compute material budget between rmin and rmax
  // and +/- etaPos in the 0--phiMax range, filling two histograms
  // Rewritten - M.Sitta 15 Nov 2014
  //
  Int_t n = 5e4; //testtracks

  TH1F *nParticles1d = new TH1F("", "", 300, 0, phiMax);
  TH2F *nParticles2d = new TH2F("", "", 300, -etaPos, etaPos, 300, 0, phiMax);

  Double_t mparam[7]={0.,0.,0.,0.,0.};
  Double_t point1[3],point2[3];

  for (Int_t it=0; it<n; it++) {

    // PHI VS ETA
    Double_t phi = gRandom->Rndm()*phiMax;
    Double_t eta = gRandom->Rndm()*2*etaPos - etaPos;
    Double_t theta = TMath::ATan(TMath::Exp(-eta))*2;
    point1[0] = rmin*TMath::Cos(phi);
    point1[1] = rmin*TMath::Sin(phi);
    point1[2] = rmin/TMath::Tan(theta);
    point2[0] = rmax*TMath::Cos(phi);
    point2[1] = rmax*TMath::Sin(phi);
    point2[2] = rmax/TMath::Tan(theta);
    
    AliTracker::MeanMaterialBudget(point1,point2,mparam);

    xOverX01d->Fill(phi, mparam[1]*100);
    nParticles1d->Fill(phi, 1.);

    if (xOverX02d != 0) {
      xOverX02d->Fill(eta, phi, mparam[1]*100);
      nParticles2d->Fill(eta, phi, 1.);
    }

    if (!(it%10000)) cout<<" : "<<mparam[1]*100<<"   phi:"<<phi<<endl;

  }

  // normalization to number of particles in case of phi vs eta
  Double_t theta = TMath::ATan(TMath::Exp(-etaPos/2))*2;
  printf("<eta>=%lf -> Sin(theta) %lf\n",etaPos/2,TMath::Sin(theta)); 

  for (Int_t ix = 1; ix<=nParticles1d->GetNbinsX(); ix++) {
    if (nParticles1d->GetBinContent(ix) > 0) 
      xOverX01d->SetBinContent(ix,xOverX01d->GetBinContent(ix)/nParticles1d->GetBinContent(ix)*TMath::Sin(theta));
    }

  if (xOverX02d) {
    for (Int_t ix = 1; ix<=nParticles2d->GetNbinsX(); ix++) {
      for (Int_t iy = 1; iy<=nParticles2d->GetNbinsY(); iy++) {
	if (nParticles2d->GetBinContent(ix,iy) > 0) 
	  xOverX02d->SetBinContent(ix,iy,xOverX02d->GetBinContent(ix,iy)/nParticles2d->GetBinContent(ix,iy));
      }
    }
  }

  Double_t mean = 0;
  for (Int_t ix = 1; ix<=xOverX01d->GetNbinsX(); ix++)
    mean+=xOverX01d->GetBinContent(ix);

  mean /= xOverX01d->GetNbinsX();

  return mean;
}


//______________________________________________________________________
void DrawMaterialBudget_Splitted(Int_t nLay = 0, Double_t rmin = 1.,
				 Double_t rmax = 5.)
{
  //
  // Function to factorize the percentage of each medium in the
  // total material budget between rmin and rmax for layer nLay
  // Cleaned up - M.Sitta 15 Nov 2014
  //
  TCanvas *c2 = new TCanvas("c2","c2");

  Double_t etaPos = 1.;

  Bool_t firstPlot = 1;

  // Check parameters
  if (nLay < 0 || nLay > 6) {
    printf("ERROR! Wrong layer number %d\n",nLay);
    return;
  }

  if (rmin < 0. || rmax < 0. || rmax < rmin) {
    printf("ERROR! Wrong radial values rmin = %f rmax = %f\n",rmin,rmax);
    return;
  }

  // In current version, build levels have different meaning for IB and OB
  const Int_t nScenariosIB = 6, nScenariosOB = 7;
  TString strDescripOB[nScenariosOB] = {"Copper", "Aluminum", "Glue", "Water",
				        "Kapton", "Carbon"  , "Silicon"      };
  //  Color codes
  //  Water                kBlack
  //  Kapton               kYellow
  //  Carbon               kRed
  //  Glue                 kBlue
  //  Aluminum             kCyan
  //  Si-Sensor            kGreen
  //  Copper               kGray
  Int_t colorsOB[nScenariosOB] = {kGray  , kCyan, kBlue, kBlack,
			          kYellow, kRed , kGreen       };

//   TString strDescripIB[nScenariosIB] = {"Carbon Structure", "Water",
// 				        "Cooling Pipe Walls and ColdPlate",
// 				        "Glue", "Flex Cable", "Pixel Chip"};
  TString strDescripIB[nScenariosIB] = {"Flex cable", "Glue",
				        "Carbon structure", "Water",
					"Cooling walls", "Pixel Chip"};
  //  Color codes
  //  Flex Cable           kBlack
  //  Glue                 kYellow
  //  Carbon Structure     kRed
  //  Water                kBlue
  //  Cooling Walls        kCyan 
  //  Si-Sensor            kGreen 
  Int_t colorsIB[nScenariosIB] = {kCyan, kBlue, kBlack, kYellow, kRed, kGreen};

  // Setup
  const Int_t maxNumberOfScenarios = nScenariosOB;

  Double_t meanX0[maxNumberOfScenarios];
  Double_t contribX0[maxNumberOfScenarios];
  TH1F *l[maxNumberOfScenarios+1];

  TString strDescrip[maxNumberOfScenarios];
  Int_t colors[maxNumberOfScenarios];

  // Choose which scenario based on Layer number and Model
  Int_t nScenarios;

  TGeoManager::Import(Form("geometry_0.root"));
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLad = gm->GetNStaves(nLay);

  Float_t phiMax = TMath::TwoPi()/nLad*2;

  char title[30];
  Int_t model, buildlevel;

  strcpy(title,
	 gGeoManager->GetVolume(Form("ITSULayer%d",nLay))->GetTitle());
  if (strlen(title) == 0) // No Model info -> old model schema
    model = -1;
  else  // New model schema -> extract model
    sscanf(title, "Model %d Build level %d", &model, &buildlevel);

  if (nLay < 3) { // IB has only 6 scenarios
    nScenarios = nScenariosIB;
    if (model == kIBModel4) // New model -> same name/colors as OB
      for (Int_t i=0; i<nScenarios; i++) {
	strDescrip[i] = strDescripOB[i+1];
	colors[i] = colorsOB[i+1];
      }
    else // Old model -> old names/colors
      for (Int_t i=0; i<nScenarios; i++) {
	strDescrip[i] = strDescripIB[i];
	colors[i] = colorsIB[i];
      }
  } else {
    if (model == kOBModel2) { // New OB
      nScenarios = nScenariosOB;
      for (Int_t i=0; i<nScenarios; i++) {
	strDescrip[i] = strDescripOB[i];
	colors[i] = colorsOB[i];
      }
    } else { // Old OB has only 6 scenarios like IB
      nScenarios = nScenariosIB;
      for (Int_t i=0; i<nScenarios; i++) {
	strDescrip[i] = strDescripOB[i+1];
	colors[i] = colorsOB[i+1];
      }
    }
  } // if (nLay < 3)

  delete gGeoManager;

  for (Int_t i=0; i<nScenarios; i++) {

    printf(" -> Loading geometry_%d.root .....\n",i);
    TGeoManager::Import(Form("geometry_%d.root",i)); 

    strcpy(title,
	   gGeoManager->GetVolume(Form("ITSULayer%d",nLay))->GetTitle());
    if (strlen(title) != 0) {
      sscanf(title, "Model %d Build level %d", &model, &buildlevel);
      if (i != buildlevel)
	printf("WARNING! Possible mismatch: file geometry_%d.root created with Build level %d\n",i,buildlevel);
    }

    TH1F *xOverX01d =  new TH1F("", "", 300, 0, phiMax);

    Double_t mean = ComputeMaterialBudget(rmin, rmax, etaPos, phiMax,
					  xOverX01d);
    meanX0[i] = mean;
    cout<<"Mean X/X0: " << meanX0[i] << " %" << endl;

    xOverX01d->GetXaxis()->SetTitle("#phi (rad)");
    xOverX01d->GetYaxis()->SetTitle(Form("X/X_{0} (%%) at #eta=0 "));
    xOverX01d->SetFillColor(colors[i]-7);
    if (i==0)  xOverX01d->SetFillColor(colors[i]);
    if ( (nScenarios == nScenariosIB && i==2) ||
	 (nScenarios == nScenariosOB && i==3) )
      xOverX01d->SetFillColor(colors[i]);
    xOverX01d->SetStats(0);

    l[i] = new TH1F("","",1,0,phiMax);
    l[i]->SetBinContent(1,mean);
    l[i]->SetStats(0);
    l[i]->SetLineColor(colors[i]-7);
    if (i==0) l[i]->SetLineColor(kGray);
    if ( (nScenarios == nScenariosIB && i==2) ||
	 (nScenarios == nScenariosOB && i==3) ) l[i]->SetLineColor(12);

    c2->cd();
    if (firstPlot) {
      xOverX01d->SetMinimum(0.);
      xOverX01d->DrawCopy(); 
      firstPlot=0;
    } else
      xOverX01d->DrawCopy("same"); 
     
    delete gGeoManager;
  }
  
  // Build meaningful legend
  TLegend *leg = new TLegend(0.76,0.77,0.99,0.99,"");
  leg->SetFillColor(0);

  for (Int_t i=0; i<nScenarios; i++) {
    // contribution in percent
    Double_t contr = 0;
    if (i == nScenarios-1) {
      contr = (meanX0[i])/meanX0[0]*100;
      contribX0[i] = meanX0[i];
    } else {
      contr = (meanX0[i]-meanX0[i+1])/meanX0[0]*100;
      contribX0[i] = meanX0[i]-meanX0[i+1];
    }
 
    strDescrip[i].Append(Form(" (%3.1lf%%)",contr));
    leg->AddEntry(l[i],strDescrip[i].Data(),"l");
  }

  TPaveText *pt = new TPaveText(0.76,0.70,0.99,0.77,"brNDC");
  pt->SetBorderSize(1); // no shadow
  pt->SetTextFont(12);
  pt->SetFillColor(0);
  pt->AddText(Form("Mean X/X0 = %4.3lf%%",meanX0[0]));
  pt->Draw();

  leg->Draw();

  l[nScenarios]=(TH1F*)l[0]->Clone();
  l[nScenarios]->SetLineColor(1);
  l[nScenarios]->SetLineWidth(2);
  l[nScenarios]->DrawCopy("same");
  
  printf("\n X/X0 contributions:\n");
  for (Int_t i=0; i<nScenarios; i++)
    printf("%s\t\t%f\n",strDescrip[i].Data(),contribX0[i]);

  c2->SaveAs("Material-details.pdf");
  
}


//______________________________________________________________________
void DrawMaterialBudget_FromTo(Int_t nLay = 0, Double_t rmin = 1.,
			       Double_t rmax = 5., Bool_t only2D = 1)
{
  //
  // Function to compute material budget between rmin and rmax
  // for layer nLay. If only2D is false a 1D histogram is also plotted
  // Cleaned up and simplified - M.Sitta 18 Nov 2014
  //

  Double_t etaPos = 0.25;

  // Check parameters
  if (nLay < 0 || nLay > 6) {
    printf("ERROR! Wrong layer number %d\n",nLay);
    return;
  }

  if (rmin < 0. || rmax < 0. || rmax < rmin) {
    printf("ERROR! Wrong radial values rmin = %f rmax = %f\n",rmin,rmax);
    return;
  }


  TGeoManager::Import("geometry.root"); 
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLad = gm->GetNStaves(nLay);

  Float_t phiMax = TMath::TwoPi()/nLad*2;

  TH2F *xOverX0   =  new TH2F("", "", 300, -1, 1., 300, 0, phiMax);
  TH1F *xOverX01d =  new TH1F("", "", 300,  0, phiMax);

  Double_t mean = ComputeMaterialBudget(rmin, rmax, -1, phiMax,
					xOverX01d, xOverX0);
  cout<<"Mean X/X0: " << mean << " %" << endl;

  TH1F *l = new TH1F("", "", 1, 0, phiMax);
  l->SetBinContent(1, mean);
  l->SetLineColor(2);
  l->SetStats(0);
  l->SetLineWidth(2);

  // Draw the histograms
  xOverX0->SetTitle(0);
  xOverX0->GetXaxis()->SetTitle("pseudorapidity #eta ");
  xOverX0->GetYaxis()->SetTitle("#phi (rad)");
  xOverX0->GetZaxis()->SetTitle("X/X_{0} (%)");
  xOverX0->SetStats(0);

  xOverX01d->SetTitle(Form("X/X_{0} (%%) within r=(%2.1lf-%2.1lf)cm average over #eta=(%1.1lf,%1.1lf)", rmin, rmax, -1., 1.));
  xOverX01d->GetXaxis()->SetTitle("#phi (rad)");
  xOverX01d->GetYaxis()->SetTitle("X/X_{0} (%) average over #eta=(-1,1)");
  xOverX01d->SetStats(0);

  if (!only2D) {
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->Divide(1,2);
    c1->cd(1); xOverX0->Draw("colz");
    c1->cd(2); xOverX01d->DrawCopy(); l->DrawCopy("same");
    c1->SaveAs("Material-1D.pdf");
  } else {
    TCanvas *c1 = new TCanvas("c1","c1");
    xOverX0->Draw("colz");
    c1->SaveAs("Material-2D.pdf");
  }

}


//______________________________________________________________________
void ComputeGeometricalParameters(const char *geofile, Double_t *rmin,
				  Double_t *rmax, Double_t *zPos,
				  Int_t *nlad, Bool_t fromGDML = 0)
{
  //
  // Ancillary function to retrieve various geometrical parameters
  // from a geometry file. The caller must ensure that all vectors
  // have proper dimensions: since this function is designed to be used
  // internally, no check is made on parameters!
  // Rewritten - M.Sitta 20 Nov 2014
  //

  if (fromGDML)  // from GDML
    TGeoManager::Import(geofile);
//    gGeoManager->ViewLeaves(true); 
//    gGeoManager->GetTopVolume()->Draw("ogl"); 
  else {         // from AliRoot simulation
        
//    gAlice=NULL;
//    AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
//    runLoader->LoadgAlice();
    
//    gAlice = runLoader->GetAliRun();
//    AliGeomManager::LoadGeometry(geofile);
    TGeoManager::Import(geofile);
    AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);

    TGeoVolume *pipeV = gGeoManager->GetVolume("IP_PIPE");
    if (!pipeV) {
      printf("Pipe volume %s is not in the geometry\n", "IP_PIPE");
      return;
    } else {
      printf("%s   ","IP_PIPE");
      TGeoTube *t = (TGeoTube*)(pipeV->GetShape());
      printf(" r = (%3.2lf %3.2lf) ", t->GetRmin(), t->GetRmax());
      printf(" z = %3.2lf\n", t->GetDz());

      rmin[0] = t->GetRmin();
      rmax[0] = t->GetRmax();
      nlad[0] = 0;
      zPos[0] = 0;
    }

    TGeoVolume *itsV = gGeoManager->GetVolume(gm->GetITSVolPattern());
    if (!itsV) {
      printf("ITS volume %s is not in the geometry\n", gm->GetITSVolPattern());
      return;
    }
    //
    // Loop on all ITSV nodes, count Layer volumes by checking names
    Int_t nNodes = itsV->GetNodes()->GetEntries();
    Int_t numberOfLayers = 0;
    for (Int_t j=0; j<nNodes; j++) {
      if ( strstr(itsV->GetNodes()->At(j)->GetName(),
		  gm->GetITSWrapVolPattern()) ) {

	TGeoNode *wrapN = (TGeoNode*)(itsV->GetNodes()->At(j));
	TGeoVolume *wrapV = wrapN->GetVolume();
 	Int_t nNodesWrap = wrapV->GetNodes()->GetEntries();

 	for (Int_t k=0; k<nNodesWrap; k++) {
 	  if ( strstr(wrapV->GetNodes()->At(k)->GetName(),
 		      gm->GetITSLayerPattern()) ) {

 	    Int_t l = numberOfLayers;
 	    numberOfLayers++;
	
 	    char laynam[30];
 	    snprintf(laynam, 30, "%s%d", gm->GetITSLayerPattern(), l);
 	    TGeoVolume* volLy = gGeoManager->GetVolume(laynam);
 	    printf("%s\t",volLy->GetName());

 	    // Layers are Assemblies not Tubes, so:
 	    // for rmax we assume the maximum extension of the assembly
 	    // for rmin: for OB (no Turbo) we take the stave height and
 	    // we subtract this from rmax; for IB (Turbo) we compute the
 	    // position of the most internal corner
 	    TGeoBBox *b = (TGeoBBox*)(volLy->GetShape());
 	    TGeoBBox *s = (TGeoBBox*)(gGeoManager->GetVolume(Form("%s%d",gm->GetITSStavePattern(),l))->GetShape());
 	    Double_t minr;

	    Double_t loc[3], mas[3];
 	    if (l < 3) {
 	      Double_t loc[3], mas[3];
  	      loc[0] = s->GetDX();
  	      loc[1] = s->GetDY();
 	      loc[2] = 0;
 	      volLy->GetNode(Form("%s%d_0",gm->GetITSStavePattern(),l))->GetMatrix()->LocalToMaster(loc,mas);
 	      minr = TMath::Sqrt(mas[0]*mas[0] + mas[1]*mas[1]);
	      zPos[l+1] = 0;
 	    } else {
 	      minr = b->GetDX() - 2*s->GetDX();
	      TGeoBBox *c = (TGeoBBox*)(gGeoManager->GetVolume(Form("%s%d",gm->GetITSChipPattern(),l))->GetShape());
	      zPos[l+1] = c->GetDZ();
 	    }
 	    rmin[l+1] = minr;
 	    rmax[l+1] = b->GetDX();
 	    nlad[l+1] = gm->GetNStaves(l);

 	    printf(" r = (%5.2lf , %5.2lf) ", rmin[l+1], rmax[l+1]);
 	    printf(" z = +/- %5.2lf ", b->GetDZ());
 	    printf(" #lad = %d \n", nlad[l+1]);

	  }
	}
      }
    }
  }

}


//______________________________________________________________________
void DrawMaterialBudget_inPhi(const char *geofile = "geometry.root",
			      Int_t lay = -1, Bool_t fromGDML = 0)
{
  //
  // Function to compute material budget as a function of phi
  // for layer lay (all layers and pipe if lay=-1). geofile is the
  // geometry file name. If fromGDML is true, use a GDML file, else
  // a standard geometry.root file.
  // Cleaned up and simplified - M.Sitta 18 Nov 2014
  //

  Double_t rmin[8], rmax[8], zPos[8];
  Int_t nlad[8];

  ComputeGeometricalParameters(geofile, rmin, rmax, zPos, nlad, fromGDML);

  Int_t n = 100000; //testtracks
  // Test Get material ?
  Double_t mparam[7]={0.,0.,0.,0.,0.};
  Double_t point1[3],point2[3];

  TH2F *xOvsPhi[8];
  Int_t ls = 0, le = 8;
  if (lay != -1) {
    ls = lay+1;
    le = lay+2;
  }

  for(Int_t i=ls; i<le; i++) {
    Double_t phiMax;
    if (i == 0) // PIPE doesn't have staves
      phiMax = TMath::TwoPi()/4.;
    else
      phiMax = TMath::TwoPi()/nlad[i]*5.; // show approx 5 ladders
    
//    xOvsPhi[i] = new TH2F("", "", 100, 0, phiMax, 100, 0.2, 0.8);
    if (i < 4)
      xOvsPhi[i] = new TH2F("", "", 100, 0, phiMax, 100, 0.0, 1.2);
    else
      xOvsPhi[i] = new TH2F("", "", 100, 0, phiMax, 100, 0.2, 2.2);

    for (Int_t it=0; it<n; it++) {
      Double_t phi = phiMax*gRandom->Rndm();
      Double_t z = zPos[i];
      point1[0] = rmin[i]*TMath::Cos(phi);
      point1[1] = rmin[i]*TMath::Sin(phi);
      point1[2] = z;
      point2[0] = rmax[i]*TMath::Cos(phi);
      point2[1] = rmax[i]*TMath::Sin(phi);
      point2[2] = z;
      AliTracker::MeanMaterialBudget(point1,point2,mparam);
      if (mparam[1] < 100)  // don't save fakes due to errors
	xOvsPhi[i]->Fill(phi,mparam[1]*100); //fxOverX0Layer
      if (!(it%10000)) cout << "layer" << i-1 << " : " << mparam[1]*100
			    << " r=(" <<rmin[i] << "," << rmax[i] << ")"
			    << " phi=" << phi <<endl;
    }

  }

  if (lay==-1) {
    TCanvas *c1 = new TCanvas("c1","Material Budget",700,800);
    c1->Divide(2,4);
    for(Int_t i=0; i<8; i++) {
      c1->cd(i+1);
      if (i>0) 
	xOvsPhi[i]->SetTitle(Form("Layer %d",i-1));
      else
	xOvsPhi[i]->SetTitle(Form("Beampipe"));
      xOvsPhi[i]->GetXaxis()->SetTitle("phi (rad)");
      xOvsPhi[i]->GetYaxis()->SetTitle("X/X_{0} (%)");
      xOvsPhi[i]->SetStats(0);
      xOvsPhi[i]->Draw("col");
    }
  } else {
    TCanvas *c1 = new TCanvas("c1","Material Budget");
    Int_t i = lay+1;
    xOvsPhi[i]->SetTitle(Form("Layer %d",lay));
    xOvsPhi[i]->GetXaxis()->SetTitle("phi (rad)");
    xOvsPhi[i]->GetYaxis()->SetTitle("X/X_{0} (%)");
    xOvsPhi[i]->SetStats(1111111);
    xOvsPhi[i]->Draw("col");
  }

}


//______________________________________________________________________
void DrawMaterialBudget_EtaVsPhi(const char *geofile = "geometry.root",
				 Int_t lay = -1, Bool_t fromGDML = 0)
{
  //
  // Function to compute material budget as a function of eta and phi
  // in eta +/- 1 for layer lay (all layers and pipe if lay=-1).
  // geofile is the geometry file name. If fromGDML is true, use a
  // GDML file, else a standard geometry.root file.
  // Cleaned up and simplified - M.Sitta 20 Nov 2014
  //

  Double_t rmin[10], rmax[10], zPos[10]; // zPos not used but needed for call
  Int_t nlad[10];

  ComputeGeometricalParameters(geofile, rmin, rmax, zPos, nlad, fromGDML);

  rmin[8] = rmin[1];
  rmax[8] = rmax[3];
  rmin[9] = rmin[4];
  rmax[9] = rmax[7];

  Int_t n = 100000; //testtracks
  // Test Get material ?
  Double_t mparam[7]={0.,0.,0.,0.,0.};
  Double_t point1[3],point2[3];

  TH2F *xOverX0[10];
  TH2F *nParticles[10];
  Int_t ls = 0, le = 10;
  if (lay != -1) {
    ls = lay+1;
    le = lay+2;
  }

  for(Int_t i=ls; i<le; i++) {
    Double_t phiMax;
    if (i == 0) // PIPE doesn't have staves
      phiMax = TMath::TwoPi()/4.;
    else if (i == 8) // Mean value of layer 0 to 2
      phiMax = TMath::TwoPi()/(0.5*(nlad[1]+nlad[3]))*5.;
    else if (i == 9) // Mean value of layer 3 to 6
      phiMax = TMath::TwoPi()/(0.5*(nlad[4]+nlad[7]))*5.;
    else
      phiMax = TMath::TwoPi()/nlad[i]*5.; // show approx 5 ladders

    xOverX0[i]    = new TH2F("", "", 100, -1., 1., 100, 0., phiMax);
    nParticles[i] = new TH2F("", "", 100, -1., 1., 100, 0., phiMax);

    for (Int_t it=0; it<n; it++) {
      Double_t phi = phiMax*gRandom->Rndm();
      Double_t eta = gRandom->Rndm()*2 - 1.; // +/- 1 eta
      Double_t theta = TMath::ATan(TMath::Exp(-eta))*2;
      point1[0] = rmin[i]*TMath::Cos(phi);
      point1[1] = rmin[i]*TMath::Sin(phi);
      point1[2] = rmin[i]/TMath::Tan(theta);
      point2[0] = rmax[i]*TMath::Cos(phi);
      point2[1] = rmax[i]*TMath::Sin(phi);
      point2[2] = rmax[i]/TMath::Tan(theta);

      AliTracker::MeanMaterialBudget(point1,point2,mparam);
      if (mparam[1] < 100) { // don't save fakes due to errors
	xOverX0[i]->Fill(eta,phi,mparam[1]*100); //fxOverX0Layer
	nParticles[i]->Fill(eta,phi,1.);
      }
      if (!(it%10000)) cout << "layer" << i-1 << " : " << mparam[1]*100
			    << " r=(" <<rmin[i] << "," << rmax[i] << ")"
			    << " phi=" << phi << " theta=" << theta <<endl;
    }

    // normalization to number of particles
    for (Int_t ix=1; ix<=nParticles[i]->GetNbinsX(); ix++) {
      for (Int_t iy=1; iy<=nParticles[i]->GetNbinsY(); iy++) {
	if (nParticles[i]->GetBinContent(ix,iy) > 0) 
	  xOverX0[i]->SetBinContent(ix,iy,xOverX0[i]->GetBinContent(ix,iy)/nParticles[i]->GetBinContent(ix,iy));
      }
    }
  }
  
  if (lay==-1) {
    TCanvas *c1 = new TCanvas("c1","Material Budget",750,900);
    c1->Divide(2,5);
    for(Int_t i=0; i<10; i++) {
      c1->cd(i+1);
      if (i==0) 
	xOverX0[i]->SetTitle(Form("Beampipe"));
      else if (i==8)
	xOverX0[i]->SetTitle(Form("Layer 0 to 2"));
      else if (i==9)
	xOverX0[i]->SetTitle(Form("Layer 3 to 6"));
      else
	xOverX0[i]->SetTitle(Form("Layer %d",i-1));
      
      
      xOverX0[i]->SetStats(0);
      xOverX0[i]->GetXaxis()->SetTitle("pseudorapidity #eta ");
      xOverX0[i]->GetYaxis()->SetTitle("#phi (rad)");
      xOverX0[i]->GetZaxis()->SetTitle("X/X_{0} (%)");
      xOverX0[i]->Draw("colz");
    }
  } else {
    TCanvas *c1 = new TCanvas("c1","Material Budget");
    Int_t i = lay+1;
    xOverX0[i]->SetTitle(Form("Layer %d",lay));
    xOverX0[i]->SetStats(0);
    xOverX0[i]->GetXaxis()->SetTitle("pseudorapidity #eta ");
    xOverX0[i]->GetYaxis()->SetTitle("#phi (rad)");
    xOverX0[i]->GetZaxis()->SetTitle("X/X_{0} (%)");
    xOverX0[i]->Draw("colz");
  }

}
