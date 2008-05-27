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

/// \ingroup macros
/// \file MUONSurveyCh8L.C
/// \brief Macro to process survey and photogrammetry data of chamber 8L
///  
///  Macro loads the survey data from .txt file using AliSurveyObj.
///  Macro MUONSurveyUtil.C is then loaded.
///
///  The transformations of the slats are obatained in 2 steps:
///  -  1. Fit a plane to the sticker targets -> psi, theta
///  -  2. Using above psi in theta obtain xc, yc, zc and phi by solving 
///        the equations from a local to global transformation of the
///        fixed button targets
///
///  Various histograms are filled and printed for monitoring.
///  MisAlignment object is then created.
/// 
/// \author Javier Castillo

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryMisAligner.h"

#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TROOT.h>
#include <TGeoManager.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>
#include <TF2.h>
#include <TH2.h>
#include <TGraph2DErrors.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPostScript.h>
#include <TPaveLabel.h>
#include <TStyle.h>

#include <fstream>

#endif

Bool_t MatrixToAngles(const Double_t *rot, Double_t *angles);
Double_t eqPlane(Double_t *x, Double_t *par);
Double_t xpCenter(Double_t *x, Double_t *par);
Double_t xnCenter(Double_t *x, Double_t *par);
Double_t ypCenter(Double_t *x, Double_t *par);
Double_t ynCenter(Double_t *x, Double_t *par);
Double_t zpCenter(Double_t *x, Double_t *par);
Double_t znCenter(Double_t *x, Double_t *par);
Double_t phixpp(Double_t *x, Double_t *par);
Double_t phixpn(Double_t *x, Double_t *par);
Double_t phixnp(Double_t *x, Double_t *par);
Double_t phixnn(Double_t *x, Double_t *par);
Double_t phiypp(Double_t *x, Double_t *par);
Double_t phiypn(Double_t *x, Double_t *par);
Double_t phiynp(Double_t *x, Double_t *par);
Double_t phiynn(Double_t *x, Double_t *par);
AliMUONGeometryTransformer *ReAlign(const AliMUONGeometryTransformer * transformer, 
				    int rMod, TGeoCombiTrans deltaDetElemTransf[], Bool_t verbose);
void MUONSurveyCh8L() {
  
  char str[100];
  char filename[100];
  
  int saveps = 1;
  const int cWidth = (int)(700*(29./21.));
  const int cHeight = 700;
//   const int   lineColor = 1;
//   const int   fillColor = 4;
//   const int   filetype  = 111; // portrait  
  const int   filetype  = 112; // landscape  
  
  sprintf(filename,"surveyChamber8L.ps");
  
  Int_t nSlats = 13;

  gROOT->LoadMacro("MUONSurveyUtil.C+");

  AliSurveyObj *so = new AliSurveyObj();
  
  Int_t size = so->GetEntries();
  printf("-> %d\n", size);
  
  so->FillFromLocalFile("Alice_MuonSystem_Chamber8LCavern_3561b.txt");
  size = so->GetEntries();
  printf("--> %d\n", size);

  Printf("Title: \"%s\"", so->GetReportTitle().Data());
  Printf("Date: \"%s\"", so->GetReportDate().Data());
  Printf("Detector: \"%s\"", so->GetDetector().Data());
  Printf("URL: \"%s\"", so->GetURL().Data());
  Printf("Number: \"%d\"", so->GetReportNumber());
  Printf("Version: \"%d\"", so->GetReportVersion());
  Printf("Observations: \"%s\"", so->GetObservations().Data());
  Printf("Coordinate System: \"%s\"", so->GetCoordSys().Data());
  Printf("Measurement Units: \"%s\"", so->GetUnits().Data());
  Printf("Nr Columns: \"%d\"", so->GetNrColumns());

  TObjArray *colNames = so->GetColumnNames();
  for (Int_t i = 0; i < colNames->GetEntries(); ++i)
    Printf("  Column %d --> \"%s\"", i, ((TObjString *) colNames->At(i))->GetString().Data());

  // Get Array of surveyed points
  Printf("Points:");
  TObjArray *points = so->GetData();
  
  for (Int_t i = 0; i < points->GetEntries(); ++i)
    Printf("  Point %d --> \"%s\"  %s ", i, ((AliSurveyPoint *) points->At(i))->GetPointName().Data(), points->At(i)->GetName());


  // Slats (#1 - #13) button targets local coordinates
  Double_t lSBTLoc6[13][2][3] = {{{  -412.50,  0.0, -(11.75+ 8.20+20.00)},
				  {   412.50,  0.0, -(11.75+16.50+20.00)}},
				 {{  -612.50,  0.0,  (11.75+ 8.20+20.00)}, 
				  {   612.50,  0.0,  (11.75+16.50+20.00)}},
				 {{ - 812.50,  0.0, -(11.75+ 8.20+20.00)}, 
				  {   812.50,  0.0, -(11.75+16.50+20.00)}},
				 {{ -1012.50,  0.0,  (11.75+ 8.20+20.00)}, 
				  {  1012.50,  0.0,  (11.75+16.50+20.00)}},
				 {{ -1012.50,  0.0, -(11.75+ 8.20+20.00)}, 
				  {  1012.50,  0.0, -(11.75+16.50+20.00)}},
				 {{ -1212.50,  5.0,  (11.75+ 8.20+20.00)}, 
				  {  1212.50,  0.0,  (11.75+16.50+20.00)}},
				 {{ -1012.50,  0.0, -(11.75+ 8.20+20.00)}, 
				  {  1012.50,  0.0, -(11.75+16.50+20.00)}},
				 {{ -1212.50,  5.0,  (11.75+ 8.20+20.00)}, 
				  {  1212.50,  0.0,  (11.75+16.50+20.00)}},
				 {{ -1012.50,  0.0, -(11.75+ 8.20+20.00)}, 
				  {  1012.50,  0.0, -(11.75+16.50+20.00)}},
				 {{ -1012.50,  0.0,  (11.75+ 8.20+20.00)}, 
				  {  1012.50,  0.0,  (11.75+16.50+20.00)}},
				 {{  -812.50,  0.0, -(11.75+ 8.20+20.00)}, 
				  {   812.50,  0.0, -(11.75+16.50+20.00)}},
				 {{  -612.50,  0.0,  (11.75+ 8.20+20.00)}, 
				  {   612.50,  0.0,  (11.75+16.50+20.00)}},
				 {{  -412.50,  0.0, -(11.75+ 8.20+20.00)},
				  {   412.50,  0.0, -(11.75+16.50+20.00)}}};
						    
  
  AliSurveyPoint *pointSST = 0;
  AliSurveyPoint *pointCST = 0;
  AliSurveyPoint *pointCPST = 0;
  AliSurveyPoint **pointSBT = new AliSurveyPoint*[2];
  
  char sPointName[10] = "5000"; 
  
  // Print length of slats 
  cout << "Slat lengths:" << endl;
  TVector3 vTemp(0., 0., 0.);
  TVector3 vSBT(0., 0., 0.);

  for (Int_t iSlat=0; iSlat<nSlats; iSlat++){
    // Get button targets survey points
    vTemp.SetXYZ(0., 0., 0.);
    if (iSlat+1<10) {
      sprintf(sPointName,"60%d%d",iSlat+1,1);
      pointSBT[0] = (AliSurveyPoint *)points->FindObject(sPointName);
      if(!pointSBT[0]) {
	cout << "Error! No button targets ... " << endl;
	break;
      }
      vSBT.SetXYZ(pointSBT[0]->GetX(),pointSBT[0]->GetY(),pointSBT[0]->GetZ());
      vTemp+=vSBT;
      sprintf(sPointName,"60%d%d",iSlat+1,2);
      pointSBT[1] = (AliSurveyPoint *)points->FindObject(sPointName);
      if(!pointSBT[1]) {
	cout << "Error! No button targets ... " << endl;
	break;
      }
      vSBT.SetXYZ(pointSBT[1]->GetX(),pointSBT[1]->GetY(),pointSBT[1]->GetZ());
      vTemp-=vSBT;
    }
    else {
      sprintf(sPointName,"6%d%d",iSlat+1,1);
      pointSBT[0] = (AliSurveyPoint *)points->FindObject(sPointName);
      if(!pointSBT[0]) {
	cout << "Error! No button targets ... " << endl;
	break;
      }
      vSBT.SetXYZ(pointSBT[0]->GetX(),pointSBT[0]->GetY(),pointSBT[0]->GetZ());
      vTemp+=vSBT;
      sprintf(sPointName,"6%d%d",iSlat+1,2);
      pointSBT[1] = (AliSurveyPoint *)points->FindObject(sPointName);
      if(!pointSBT[1]) {
	cout << "Error! No button targets ... " << endl;
	break;
      }
      vSBT.SetXYZ(pointSBT[1]->GetX(),pointSBT[1]->GetY(),pointSBT[1]->GetZ());
      vTemp-=vSBT;
    }
    cout << "Slat " << iSlat+1 << ": " << vTemp.Mag() << endl;
  }
 
  // Histograms for monitoring
  TH2F *hCPSTry = new TH2F("hCPSTry","hCPSTry",28,-600,2200,52,0,5200);
  TH2F *hCPSTly = new TH2F("hCPSTly","hCPSTly",28,-600,2200,52,0,5200);
  TH2F *hSSTry = new TH2F("hSSTry","hSSTry",70,-600,2200,130,0,5200);
  TH2F *hSSTly = new TH2F("hSSTly","hSSTly",70,-600,2200,130,0,5200);
  TH2F *hCSTy = new TH2F("hCSTy","hCSTy",70,-600,2200,130,0,5200);
  TH2F *hSSTrpy = new TH2F("hSSTrpy","hSSTrpy",70,-600,2200,130,0,5200); 
  TH2F *hSSTlpy = new TH2F("hSSTlpy","hSSTlpy",70,-600,2200,130,0,5200); 

  // Chamber plane sticker targets
  for (int iPoint=0; iPoint<9; iPoint++) {
    sprintf(sPointName,"700%d",iPoint+1);
    pointCPST = (AliSurveyPoint *)points->FindObject(sPointName);
    if(!pointCPST) {
      printf("Point %s is missing ...\n",sPointName);
      break;       
    }
    hCPSTry->Fill(pointCPST->GetX(),pointCPST->GetZ(),pointCPST->GetY());
  }
  for (int iPoint=9; iPoint<18; iPoint++) {
    sprintf(sPointName,"70%d",iPoint+1);
    pointCPST = (AliSurveyPoint *)points->FindObject(sPointName);
    if(!pointCPST) {
      printf("Point %s is missing ...\n",sPointName);
      break;       
    }   
    hCPSTly->Fill(pointCPST->GetX(),pointCPST->GetZ(),pointCPST->GetY());
  }

  // Chamber Side Targets
  for (int iPoint=0; iPoint<25; iPoint++) {
    if (iPoint+1<10) {   
      sprintf(sPointName,"800%d",iPoint+1);
    }
    else {
      sprintf(sPointName,"80%d",iPoint+1);
    }
    pointCST = (AliSurveyPoint *)points->FindObject(sPointName);
    if(!pointCPST) {
      printf("Point %s is missing ...\n",sPointName);
      break;
    }
    hCSTy->Fill(pointCST->GetX(),pointCST->GetZ(),pointCST->GetY());
  }

  
  // Graphs used for plane fitting
  TGraph2DErrors **gSST5 = new TGraph2DErrors*[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    gSST5[iSlat] = new TGraph2DErrors();
  }
  
  // Keep the id of slat sticker targets next to slat button targets
  Int_t iSSBT[13][2] = {{0}};
  
  // Fill graph with sticker target positions
  for (int iSlat=0; iSlat<nSlats; iSlat++){
    for (int iPoint=0; iPoint<9; iPoint++) {
      // Second sticker target is next to first button target
      // Previous to last sticker target is next to second button target
      iSSBT[iSlat][0] = 2;
      iSSBT[iSlat][1] = 8;
      if (iSlat+1<10) {
	sprintf(sPointName,"50%d%d",iSlat+1,iPoint+1);
      }
      else {
	sprintf(sPointName,"5%d%d",iSlat+1,iPoint+1);
      }
      pointSST = (AliSurveyPoint *)points->FindObject(sPointName);
      if(!pointSST) {
	printf("%s\n",sPointName);
	cout << iSlat << " " << iPoint  << " " << pointSST << endl;
	// Previous to last sticker target is next to second button target
	iSSBT[iSlat][1] = iPoint+1-2;
	break;
      }
     
      gSST5[iSlat]->SetPoint(iPoint,-pointSST->GetX(),pointSST->GetZ(),pointSST->GetY());
      gSST5[iSlat]->SetPointError(iPoint,pointSST->GetPrecisionX(),pointSST->GetPrecisionZ(),pointSST->GetPrecisionY());
      
      // Fill histograms of sticker targets. For monitoring purposes.
      if((iSlat+1)%2==0){
	hSSTly->Fill(pointSST->GetX(),pointSST->GetZ(),pointSST->GetY());
      }
      else {
	hSSTry->Fill(pointSST->GetX(),pointSST->GetZ(),pointSST->GetY());
      }
    }
  }

  Float_t xMin = -2200.;
  Float_t xMax =  2200.;
  Float_t yMin = -5200.;
  Float_t yMax =  5200.;
  Float_t zMin = -200.;
  Float_t zMax =  200.;

  Double_t xMinSlat, xMaxSlat;
  Double_t yMinSlat, yMaxSlat;  

  // Slat plane function
  char fsName[100] = "fSlat00"; 
  TF2 **fSlat = new TF2*[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    sprintf(fsName,"fSlat%d",iSlat+1);
    fSlat[iSlat] = new TF2(fsName,eqPlane,xMin,xMax,yMin,yMax,3);
  }

  // Xcenter functions
  char fxcName[100] = "fXcnSlat00"; 
  TF2 ***fXcSlat = new TF2**[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    fXcSlat[iSlat] = new TF2*[2];
    sprintf(fxcName,"fXcnSlat%d",iSlat+1);
    fXcSlat[iSlat][0] = new TF2(fxcName,xnCenter,xMin,xMax,yMin,yMax,7);
    sprintf(fxcName,"fXcpSlat%d",iSlat+1);
    fXcSlat[iSlat][1] = new TF2(fxcName,xpCenter,xMin,xMax,yMin,yMax,7);   
  }

  // Ycenter functions
  char fycName[100] = "fYcnSlat00"; 
  TF2 ***fYcSlat = new TF2**[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    fYcSlat[iSlat] = new TF2*[2];
    sprintf(fycName,"fYcnSlat%d",iSlat+1);
    fYcSlat[iSlat][0] = new TF2(fycName,ynCenter,yMin,yMax,yMin,yMax,8);
    sprintf(fycName,"fYcpSlat%d",iSlat+1);
    fYcSlat[iSlat][1] = new TF2(fycName,ypCenter,yMin,yMax,yMin,yMax,8);   
  }

  // Zcenter functions
  char fzcName[100] = "fZcnSlat00"; 
  TF2 ***fZcSlat = new TF2**[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    fZcSlat[iSlat] = new TF2*[2];
    sprintf(fzcName,"fZcnSlat%d",iSlat+1);
    fZcSlat[iSlat][0] = new TF2(fzcName,znCenter,zMin,zMax,zMin,zMax,8);
    sprintf(fzcName,"fZcpSlat%d",iSlat+1);
    fZcSlat[iSlat][1] = new TF2(fzcName,zpCenter,zMin,zMax,zMin,zMax,8);   
  }

  // Phi rotation using xglobal coords functions
  char fphixName[100] = "fPhiXnnSlat00"; 
  TF2 ****fPhiXSlat = new TF2***[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    fPhiXSlat[iSlat] = new TF2**[2];
    for (Int_t iX =0; iX<2; iX++)
      fPhiXSlat[iSlat][iX] = new TF2*[2];
    sprintf(fphixName,"fPhiXnnSlat%d",iSlat+1);
    fPhiXSlat[iSlat][0][0] = new TF2(fphixName,phixnn,xMin,xMax,xMin,xMax,7);
    sprintf(fphixName,"fPhixnpSlat%d",iSlat+1);
    fPhiXSlat[iSlat][0][1] = new TF2(fphixName,phixnp,xMin,xMax,xMin,xMax,7);   
    sprintf(fphixName,"fPhiXpnSlat%d",iSlat+1);
    fPhiXSlat[iSlat][1][0] = new TF2(fphixName,phixpn,xMin,xMax,xMin,xMax,7);
    sprintf(fphixName,"fPhixppSlat%d",iSlat+1);
    fPhiXSlat[iSlat][1][1] = new TF2(fphixName,phixpp,xMin,xMax,xMin,xMax,7);   
  }

  // Phi rotation using yglobal coords functions
  char fphiyName[100] = "fPhiYnnSlat00"; 
  TF2 ****fPhiYSlat = new TF2***[nSlats];
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    fPhiYSlat[iSlat] = new TF2**[2];
    for (Int_t iY =0; iY<2; iY++)
      fPhiYSlat[iSlat][iY] = new TF2*[2];
    sprintf(fphiyName,"fPhiYnnSlat%d",iSlat+1);
    fPhiYSlat[iSlat][0][0] = new TF2(fphiyName,phiynn,yMin,yMax,yMin,yMax,7);
    sprintf(fphiyName,"fPhiYnpSlat%d",iSlat+1);
    fPhiYSlat[iSlat][0][1] = new TF2(fphiyName,phiynp,yMin,yMax,yMin,yMax,7);   
    sprintf(fphiyName,"fPhiYpnSlat%d",iSlat+1);
    fPhiYSlat[iSlat][1][0] = new TF2(fphiyName,phiypn,yMin,yMax,yMin,yMax,7);
    sprintf(fphiyName,"fPhiYppSlat%d",iSlat+1);
    fPhiYSlat[iSlat][1][1] = new TF2(fphiyName,phiypp,yMin,yMax,yMin,yMax,7);   
  }

//   Double_t *xce = new Double_t[nSlats]; 
//   Double_t *yce = new Double_t[nSlats]; 
//   Double_t *zce = new Double_t[nSlats]; 
  Double_t *psi = new Double_t[nSlats]; 
  Double_t *tht = new Double_t[nSlats]; 
//   Double_t *phi = new Double_t[nSlats]; 

  Double_t **lCenSlat = new Double_t*[nSlats];
  Double_t **lRotSlat = new Double_t*[nSlats];
  Double_t **lDiffCenSlat0 = new Double_t*[nSlats];
  Double_t **lDiffRotSlat0 = new Double_t*[nSlats];
  Double_t **lDiffThCenSlat0 = new Double_t*[nSlats];
  Double_t **lDiffThRotSlat0 = new Double_t*[nSlats];
  Double_t **lDeltaDiffCenSlat0 = new Double_t*[nSlats];
  Double_t **lDeltaDiffRotSlat0 = new Double_t*[nSlats];

  for (int iSlat=0; iSlat<nSlats; iSlat++) {
    lCenSlat[iSlat] = new Double_t[3];
    lRotSlat[iSlat] = new Double_t[3];

    lDiffCenSlat0[iSlat] = new Double_t[3];
    lDiffRotSlat0[iSlat] = new Double_t[3];
    lDiffThCenSlat0[iSlat] = new Double_t[3];
    lDiffThRotSlat0[iSlat] = new Double_t[3];
    lDeltaDiffCenSlat0[iSlat] = new Double_t[3];
    lDeltaDiffRotSlat0[iSlat] = new Double_t[3];
  }


  TGeoTranslation transSlat[nSlats];
  TGeoRotation rotSlat[nSlats];
  TGeoCombiTrans trfSlat[nSlats];

  TGeoTranslation dtransSlat[nSlats];
  TGeoRotation drotSlat[nSlats];
  TGeoCombiTrans dtrfSlat[nSlats];

  TGeoTranslation transTemp;
  TGeoRotation rotTemp;
  TGeoCombiTrans trfTemp;


  Double_t lCenTemp[3];
  Double_t lRotTemp[3];

  Double_t lDiffTemp[9];
  Double_t lDiffMin[9];

  double tempDiff = 0.;
  double tempDiff1 = 0.;
  double tempDiff2 = 0.;

  AliSurveyPoint **pointSSBT = new AliSurveyPoint*[2];

  //
  // Get Slat transformation. 
  // Psi and Theta are obtained by fitting a plane to the sticker targets.
  // Then Xc, Yc, Zc and Phi are obtained by solving the equations to the ref.
  // syst. transformation of the button targets
  //
  for (int iSlat=0; iSlat<nSlats; iSlat++) {
    sprintf(fsName,"fSlat%d",iSlat+1);
    cout << "Fitting Slat" << iSlat+1 << " ..." << endl;

    // Fit a plane to the sticker targets
    gSST5[iSlat]->Fit(fsName,"","same");

    psi[iSlat] = TMath::ATan(fSlat[iSlat]->GetParameter(1));
    tht[iSlat] = TMath::ATan(fSlat[iSlat]->GetParameter(0));
    if (iSlat==5)
      psi[iSlat] += TMath::Pi(); // Rotated slat

    lRotSlat[iSlat][0] = psi[iSlat];
    lRotSlat[iSlat][1] = tht[iSlat];

    for(Int_t iS=0; iS<2; iS++){
      fXcSlat[iSlat][iS]->SetParameters(lSBTLoc6[iSlat][0][0],lSBTLoc6[iSlat][0][1],lSBTLoc6[iSlat][0][2],lSBTLoc6[iSlat][1][0],lSBTLoc6[iSlat][1][1],lSBTLoc6[iSlat][1][2],tht[iSlat]);
      fYcSlat[iSlat][iS]->SetParameters(lSBTLoc6[iSlat][0][0],lSBTLoc6[iSlat][0][1],lSBTLoc6[iSlat][0][2],lSBTLoc6[iSlat][1][0],lSBTLoc6[iSlat][1][1],lSBTLoc6[iSlat][1][2],psi[iSlat],tht[iSlat]);
      fZcSlat[iSlat][iS]->SetParameters(lSBTLoc6[iSlat][0][0],lSBTLoc6[iSlat][0][1],lSBTLoc6[iSlat][0][2],lSBTLoc6[iSlat][1][0],lSBTLoc6[iSlat][1][1],lSBTLoc6[iSlat][1][2],psi[iSlat],tht[iSlat]);
      for(Int_t jS=0; jS<2; jS++){
	fPhiXSlat[iSlat][iS][jS]->SetParameters(lSBTLoc6[iSlat][0][0],lSBTLoc6[iSlat][0][1],lSBTLoc6[iSlat][0][2],lSBTLoc6[iSlat][1][0],lSBTLoc6[iSlat][1][1],lSBTLoc6[iSlat][1][2],tht[iSlat]);
	fPhiYSlat[iSlat][iS][jS]->SetParameters(lSBTLoc6[iSlat][0][0],lSBTLoc6[iSlat][0][1],lSBTLoc6[iSlat][0][2],lSBTLoc6[iSlat][1][0],lSBTLoc6[iSlat][1][1],lSBTLoc6[iSlat][1][2],psi[iSlat],tht[iSlat]);
      }
    }

    //
    // Calculate Slat Center from button targets
    //

    // Get button targets survey points
    for (Int_t iPoint=0; iPoint<2; iPoint++) {
      if (iSlat+1<10) {
	sprintf(sPointName,"60%d%d",iSlat+1,iPoint+1);
	pointSBT[iPoint] = (AliSurveyPoint *)points->FindObject(sPointName);
	if(!pointSBT[iPoint]) {
	  cout << "Error! No button targets ... " << endl;
	  break;       
	}
	sprintf(sPointName,"50%d%d",iSlat+1,iSSBT[iSlat][iPoint]);
	pointSSBT[iPoint] = (AliSurveyPoint *)points->FindObject(sPointName);
	if(!pointSSBT[iPoint]) {
	  cout << "Error! No sticker target ... " << sPointName << endl;
	  break;       
	}
      }
      else {
	sprintf(sPointName,"6%d%d",iSlat+1,iPoint+1);
	pointSBT[iPoint] = (AliSurveyPoint *)points->FindObject(sPointName);
	if(!pointSBT[iPoint]) {
	  cout << "Error! No button targets ... " << endl;
	  break;
	}
	sprintf(sPointName,"5%d%d",iSlat+1,iSSBT[iSlat][iPoint]);
	pointSSBT[iPoint] = (AliSurveyPoint *)points->FindObject(sPointName);
	if(!pointSSBT[iPoint]) {
	  cout << "Error! No sticker targets ... " << sPointName << endl;
	  break;
	}
      }
    }

    tempDiff += TMath::Power(-1,iSlat)*((pointSBT[1]->GetY() - pointSSBT[1]->GetY())-(pointSBT[0]->GetY() - pointSSBT[0]->GetY()));
    tempDiff1 += TMath::Abs(pointSBT[0]->GetY() - pointSSBT[0]->GetY())-20;
    tempDiff2 += TMath::Abs(pointSBT[1]->GetY() - pointSSBT[1]->GetY())-20;
    cout << "BSdiff: " << TMath::Abs(pointSBT[0]->GetY() - pointSSBT[0]->GetY()) << " " << TMath::Abs(pointSBT[1]->GetY() - pointSSBT[1]->GetY()) << " " << tempDiff1/(iSlat+1) << " " << tempDiff2/(iSlat+1) << " " << tempDiff/(iSlat+1) << endl;


    //    Double_t p0l[3] = {0};
    Double_t p1l[3] = {0};
    Double_t p2l[3] = {0};
    //    Double_t p0g[3] = {0};
    Double_t p1g[3] = {0};
    Double_t p2g[3] = {0};

    //    p0l[2] = lSBTLoc6[iSlat][0][2];
    // Button targets local coordinates
    for(Int_t iCor=0; iCor<3; iCor++){
      p1l[iCor]= lSBTLoc6[iSlat][0][iCor];
      p2l[iCor]= lSBTLoc6[iSlat][1][iCor];
    }
    for(Int_t i=0; i<9; i++){
      lDiffMin[i]=1000000.;
    }
    // Trying 2x*2y*2z*2phi possibilities
    for(Int_t iX=0; iX<2; iX++){
      for(Int_t iY=0; iY<2; iY++){
	for(Int_t iZ=0; iZ<2; iZ++){
	  lCenTemp[0] = fXcSlat[iSlat][iX]->Eval(-pointSBT[0]->GetX(),-pointSBT[1]->GetX());
	  lCenTemp[1] = fYcSlat[iSlat][iY]->Eval(pointSBT[0]->GetZ(),pointSBT[1]->GetZ());
	  lCenTemp[2] = fZcSlat[iSlat][iZ]->Eval(pointSBT[0]->GetY(),pointSBT[1]->GetY());
	  //	 lCenTemp[2] = fZcSlat[iSlat][iZ]->Eval(pointSSBT[0]->GetY(),pointSSBT[1]->GetY());
	  lRotTemp[0] = psi[iSlat];
	  lRotTemp[1] = tht[iSlat];
	  for(Int_t iP=0; iP<2; iP++){
	    lRotTemp[2] = fPhiXSlat[iSlat][iX][iP]->Eval(-pointSBT[0]->GetX(),-pointSBT[1]->GetX());

	    trfTemp.SetTranslation(transTemp);
	    trfTemp.SetRotation(rotTemp);
	    trfTemp.Clear();
	    trfTemp.RotateZ(TMath::RadToDeg()*lRotTemp[2]);
	    trfTemp.RotateY(TMath::RadToDeg()*lRotTemp[1]);
	    trfTemp.RotateX(TMath::RadToDeg()*lRotTemp[0]);
	    trfTemp.SetTranslation(lCenTemp);
	   
	    // 	   trfTemp.LocalToMaster(p0l, p0g);
	    // 	   lCenTemp[2]= fSlat[iSlat]->Eval(p0g[0],p0g[1]) - trfTemp.GetRotationMatrix()[8]*p0l[2];
	    // 	   trfTemp.SetTranslation(lCenTemp);

	    trfTemp.LocalToMaster(p1l, p1g);
	    trfTemp.LocalToMaster(p2l, p2g);
	   
	    lDiffTemp[0] = (-pointSBT[0]->GetX()-p1g[0]);
	    lDiffTemp[1] = (pointSBT[0]->GetZ()-p1g[1]);
	    lDiffTemp[2] = (pointSBT[0]->GetY()-p1g[2]);
	    //	   lDiffTemp[2] = (pointSSBT[0]->GetY()-p1g[2]);
	    lDiffTemp[3] = TMath::Sqrt(lDiffTemp[0]*lDiffTemp[0]+lDiffTemp[1]*lDiffTemp[1]+lDiffTemp[2]*lDiffTemp[2]);	   
	    lDiffTemp[4] = (-pointSBT[1]->GetX()-p2g[0]);
	    lDiffTemp[5] = (pointSBT[1]->GetZ()-p2g[1]);
	    lDiffTemp[6] = (pointSBT[1]->GetY()-p2g[2]);
	    //	   lDiffTemp[6] = (pointSSBT[1]->GetY()-p2g[2]);
	    lDiffTemp[7] = TMath::Sqrt(lDiffTemp[4]*lDiffTemp[4]+lDiffTemp[5]*lDiffTemp[5]+lDiffTemp[6]*lDiffTemp[6]);
	    lDiffTemp[8] = TMath::Sqrt(lDiffTemp[3]*lDiffTemp[3]+lDiffTemp[7]*lDiffTemp[7]);
	   
	    if(lDiffTemp[8]<lDiffMin[8]){
	      cout << "Diffs" ;
	      for(Int_t i=0; i<9; i++){
		cout << " " << lDiffTemp[i]; 
	      }
	      cout << endl;
	      cout << "Slat" << iSlat+1 << " : mycenX" << iX << "Y" << iY << "Z" << iZ << "(" << lCenTemp[0] << "," << lCenTemp[1] << "," << lCenTemp[2] << "); rotx" << iP << "(" << lRotTemp[0] << "," << lRotTemp[1] << "," << lRotTemp[2] << ")"  << endl;  	   
	      cout << p1g[0] << " " << p1g[1] << " " << p1g[2] << " " << p2g[0] << " " << p2g[1] << " " << p2g[2] << endl;
	      cout << "Transformation improved ..." << endl;
	      lCenSlat[iSlat][0] = lCenTemp[0]; lCenSlat[iSlat][1] = lCenTemp[1]; lCenSlat[iSlat][2] = lCenTemp[2]; 
	      lRotSlat[iSlat][2] = lRotTemp[2];
	      for(Int_t i=0; i<9; i++){
		lDiffMin[i]=lDiffTemp[i];
	      }
	      if((lDiffMin[3]*lDiffMin[3]<0.1*0.1+0.1*0.1+0.1*0.1)&&
		 (lDiffMin[7]*lDiffMin[7]<0.1*0.1+0.1*0.1+0.1*0.1)){
		cout << "Correct Transformation found X " << iX  << " Y " << iY << " Z " << iZ << " xP " << iP << endl;
		lCenSlat[iSlat][0] = lCenTemp[0]; lCenSlat[iSlat][1] = lCenTemp[1]; lCenSlat[iSlat][2] = lCenTemp[2]; 
		lRotSlat[iSlat][2] = lRotTemp[2];
	      }
	    }
	  }
	  for(Int_t iP=0; iP<2; iP++){
	    lRotTemp[2] = fPhiYSlat[iSlat][iY][iP]->Eval(pointSBT[0]->GetZ(),pointSBT[1]->GetZ());
	    Double_t lPhi = TMath::ATan2((pointSBT[1]->GetZ()-pointSBT[0]->GetZ()),-(pointSBT[1]->GetX()-pointSBT[0]->GetX()));
	   
	    trfTemp.Clear();
	    trfTemp.RotateZ(TMath::RadToDeg()*lRotTemp[2]);
	    trfTemp.RotateY(TMath::RadToDeg()*lRotTemp[1]);
	    trfTemp.RotateX(TMath::RadToDeg()*lRotTemp[0]);
	    trfTemp.SetTranslation(lCenTemp);

	    // 	   trfTemp.LocalToMaster(p0l, p0g);
	    // 	   lCenTemp[2]= fSlat[iSlat]->Eval(p0g[0],p0g[1]) - trfTemp.GetRotationMatrix()[8]*p0l[2];
	    // 	   trfTemp.SetTranslation(lCenTemp);
	   
	    trfTemp.LocalToMaster(p1l, p1g);
	    trfTemp.LocalToMaster(p2l, p2g);

	    lDiffTemp[0] = (-pointSBT[0]->GetX()-p1g[0]);
	    lDiffTemp[1] = (pointSBT[0]->GetZ()-p1g[1]);	
	    lDiffTemp[2] = (pointSBT[0]->GetY()-p1g[2]);
	    //	   lDiffTemp[2] = (pointSSBT[0]->GetY()-p1g[2]);
	    lDiffTemp[3] = TMath::Sqrt(lDiffTemp[0]*lDiffTemp[0]+lDiffTemp[1]*lDiffTemp[1]+lDiffTemp[2]*lDiffTemp[2]);	   
	    lDiffTemp[4] = (-pointSBT[1]->GetX()-p2g[0]);
	    lDiffTemp[5] = (pointSBT[1]->GetZ()-p2g[1]);	
	    lDiffTemp[6] = (pointSBT[1]->GetY()-p2g[2]);
	    //	   lDiffTemp[6] = (pointSSBT[1]->GetY()-p2g[2]);
	    lDiffTemp[7] = TMath::Sqrt(lDiffTemp[4]*lDiffTemp[4]+lDiffTemp[5]*lDiffTemp[5]+lDiffTemp[6]*lDiffTemp[6]);
	    lDiffTemp[8] = TMath::Sqrt(lDiffTemp[3]*lDiffTemp[3]+lDiffTemp[7]*lDiffTemp[7]);
	   
	    if(lDiffTemp[8]<lDiffMin[8]){
	      cout << "Diffs" ;
	      for(Int_t i=0; i<9; i++){
		cout << " " << lDiffTemp[i]; 
	      }
	      cout << endl;	     
	      cout << "Slat" << iSlat+1 << " : mycenX" << iX << "Y" << iY << "Z" << iZ << "(" << lCenTemp[0] << "," << lCenTemp[1] << "," << lCenTemp[2] << "); roty" << iP << "(" << lRotTemp[0] << "," << lRotTemp[1] << "," << lRotTemp[2] << "(" << lPhi << "))"  << endl;
	      cout << p1g[0] << " " << p1g[1] << " " << p1g[2] << " " << p2g[0] << " " << p2g[1] << " " << p2g[2] << endl;	     
	      cout << "Transformation improved ..." << endl;
	      lCenSlat[iSlat][0] = lCenTemp[0]; lCenSlat[iSlat][1] = lCenTemp[1]; lCenSlat[iSlat][2] = lCenTemp[2]; 
	      lRotSlat[iSlat][2] = lRotTemp[2];

	      for(Int_t i=0; i<9; i++){
		lDiffMin[i]=lDiffTemp[i];
	      }
	      if((lDiffMin[3]*lDiffMin[3]<0.1*0.1+0.1*0.1+0.1*0.1)&&
		 (lDiffMin[7]*lDiffMin[7]<0.1*0.1+0.1*0.1+0.1*0.1)){
		cout << "Correct Transformation found X " << iX  << " Y " << iY << " Z " << iZ << " yP " << iP << endl;
		lCenSlat[iSlat][0] = lCenTemp[0]; lCenSlat[iSlat][1] = lCenTemp[1]; lCenSlat[iSlat][2] = lCenTemp[2]; 
		lRotSlat[iSlat][2] = lRotTemp[2];
	      }
	    }
	  } 
	} 
      } 
    }

    // Fill slat plane for fit monitor.
    xMinSlat = hSSTrpy->GetXaxis()->FindBin(pointSBT[0]->GetX()); 
    xMaxSlat = hSSTrpy->GetXaxis()->FindBin(pointSBT[1]->GetX()); 
    yMinSlat = hSSTrpy->GetYaxis()->FindBin(pointSBT[0]->GetZ()-200.); 
    yMaxSlat = hSSTrpy->GetYaxis()->FindBin(pointSBT[0]->GetZ()+200.); 
   
    for (int i=(int)xMinSlat; i<=(int)xMaxSlat; i++) {
      for (int j=(int)yMinSlat; j<=(int)yMaxSlat; j++) {
	Double_t zSlat = fSlat[iSlat]->Eval(-hSSTrpy->GetXaxis()->GetBinCenter(i),
					    hSSTrpy->GetYaxis()->GetBinCenter(j));
	if((iSlat+1)%2==0){
	  hSSTlpy->SetBinContent(i,j,zSlat);
	}
	else {
	  hSSTrpy->SetBinContent(i,j,zSlat);
	}
      }
    }   
  }

  //
  // Compare transformations to expected ones 
  //
  Int_t iSlatToPos[13] = {0, 11, 2, 9, 4, 7, 6, 5, 8, 3, 10, 1, 12};
 
  // Theoretical differences with respect to Slat 513 which is here Slat07
  lDiffThCenSlat0[0][0]  = (2-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[1][0]  = (3-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[2][0]  = (4-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[3][0]  = (5-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[4][0]  = (5-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[5][0]  = (6-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[6][0]  = (5-5)*400./2. +12.5 -375. -25. +375. +25. -12.5;
  lDiffThCenSlat0[7][0]  = (6-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[8][0]  = (5-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[9][0]  = (5-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[10][0] = (4-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[11][0] = (3-5)*400./2. +12.5 -375. -25.;
  lDiffThCenSlat0[12][0] = (2-5)*400./2. +12.5 -375. -25.;

  lDiffThCenSlat0[12][1] = 382.0 +378.5 +375.5 +294.0 +370.0 +286.0; 
  lDiffThCenSlat0[1][1]  = 382.0 +378.5 +375.5 +294.0 +370.0; 
  lDiffThCenSlat0[10][1] = 382.0 +378.5 +375.5 +294.0; 
  lDiffThCenSlat0[3][1]  = 382.0 +378.5 +375.5; 
  lDiffThCenSlat0[8][1]  = 382.0 +378.5; 
  lDiffThCenSlat0[5][1]  = 382.0; 
  lDiffThCenSlat0[6][1]  = 0.0; 
  lDiffThCenSlat0[7][1]  = -382.0; 
  lDiffThCenSlat0[4][1]  = -382.0 -378.5; 
  lDiffThCenSlat0[9][1]  = -382.0 -378.5 -375.5; 
  lDiffThCenSlat0[2][1]  = -382.0 -378.5 -375.5 -294.0; 
  lDiffThCenSlat0[11][1] = -382.0 -378.5 -375.5 -294.0 -370.0; 
  lDiffThCenSlat0[0][1]  = -382.0 -378.5 -375.5 -294.0 -370.0 -286.0; 

  lDiffThCenSlat0[12][2] =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[1][2]  = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[10][2] =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[3][2]  = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[8][2]  =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[5][2]  = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[6][2]  =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[7][2]  = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[4][2]  =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[9][2]  = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[2][2]  =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[11][2] = -(42.5)-(42.5); //-(42.5-1.175)-(42.5-1.175);
  lDiffThCenSlat0[0][2]  =  (42.5)-(42.5); // (42.5-1.175)-(42.5-1.175);

  TGraph *gDeltaDiffCenXSlat0 = new TGraph(nSlats);
  TGraph *gDeltaDiffCenYSlat0 = new TGraph(nSlats);
  TGraph *gDeltaDiffCenZSlat0 = new TGraph(nSlats);
  TGraph *gDeltaDiffPsiSlat0 = new TGraph(nSlats);
  TGraph *gDeltaDiffThtSlat0 = new TGraph(nSlats);
  TGraph *gDeltaDiffPhiSlat0 = new TGraph(nSlats);

  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    trfSlat[iSlat].SetTranslation(transSlat[iSlat]);
    trfSlat[iSlat].SetRotation(rotSlat[iSlat]);
    trfSlat[iSlat].Clear();
    trfSlat[iSlat].RotateZ(TMath::RadToDeg()*lRotSlat[iSlat][2]);
    trfSlat[iSlat].RotateY(TMath::RadToDeg()*lRotSlat[iSlat][1]);
    trfSlat[iSlat].RotateX(TMath::RadToDeg()*lRotSlat[iSlat][0]);
    trfSlat[iSlat].SetTranslation(lCenSlat[iSlat]); 
  }

  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    dtrfSlat[iSlat].SetTranslation(dtransSlat[iSlat]);
    dtrfSlat[iSlat].SetRotation(drotSlat[iSlat]);
    dtrfSlat[iSlat].Clear();
    dtrfSlat[iSlat] = trfSlat[6].Inverse()*trfSlat[iSlat]; 
    dtrfSlat[iSlat].Print();
    lDiffCenSlat0[iSlat] = (Double_t*)dtrfSlat[iSlat].GetTranslation();
    MatrixToAngles(dtrfSlat[iSlat].GetRotationMatrix(),lDiffRotSlat0[iSlat]);

    lDeltaDiffCenSlat0[iSlat][0] = lDiffCenSlat0[iSlat][0]-lDiffThCenSlat0[iSlat][0];
    lDeltaDiffCenSlat0[iSlat][1] =  -lDiffCenSlat0[iSlat][1]-lDiffThCenSlat0[iSlat][1];
    lDeltaDiffCenSlat0[iSlat][2] =  -lDiffCenSlat0[iSlat][2]-lDiffThCenSlat0[iSlat][2];
    lDeltaDiffRotSlat0[iSlat][0] = lDiffRotSlat0[iSlat][0];
    lDeltaDiffRotSlat0[iSlat][1] = lDiffRotSlat0[iSlat][1];
    lDeltaDiffRotSlat0[iSlat][2] = lDiffRotSlat0[iSlat][2];
    gDeltaDiffCenXSlat0->SetPoint(iSlat,lDeltaDiffCenSlat0[iSlat][0],iSlatToPos[iSlat]+1);
    gDeltaDiffCenYSlat0->SetPoint(iSlat,lDeltaDiffCenSlat0[iSlat][1],iSlatToPos[iSlat]+1);
    gDeltaDiffCenZSlat0->SetPoint(iSlat,lDeltaDiffCenSlat0[iSlat][2],iSlatToPos[iSlat]+1);
    gDeltaDiffPsiSlat0->SetPoint(iSlat,1e3*lDeltaDiffRotSlat0[iSlat][0],iSlatToPos[iSlat]+1);
    gDeltaDiffThtSlat0->SetPoint(iSlat,1e3*lDeltaDiffRotSlat0[iSlat][1],iSlatToPos[iSlat]+1);
    gDeltaDiffPhiSlat0->SetPoint(iSlat,1e3*lDeltaDiffRotSlat0[iSlat][2],iSlatToPos[iSlat]+1);
  }

  // Import TGeo geometry 
  char* geoFilename = "geometry.root";
  cout << "geometry imported" << endl;
  if ( ! AliGeomManager::GetGeometry() ) {
    AliGeomManager::LoadGeometry(geoFilename);
    if (! AliGeomManager::GetGeometry() ) {
      printf("MUONSurveyCh8L: getting geometry from file %s failed\n", geoFilename);
      return;
    }
  }

  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer();
//   transform->ReadGeometryData("volpath.dat", gGeoManager);
  transform->LoadGeometryData();
  cout << "geometry data read" << endl;
  AliMUONGeometryTransformer *newTransform = ReAlign(transform,11,dtrfSlat,true); 
  newTransform->WriteTransformations("transform2ReAlign.dat");
  cout << "newtransform read" << endl;
  // Generate realigned data in local cdb
  const TClonesArray* array = newTransform->GetMisAlignmentData();
   
  // CDB manager
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://ReAlignCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, 0); 
  cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);

  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    cout << lDeltaDiffCenSlat0[iSlat][0] << " " << lDiffCenSlat0[iSlat][0] << " " << lDiffThCenSlat0[iSlat][0] << endl;
  }
  cout << endl;
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    //    lDiffCenSlat0[iSlat][1] = lCenSlat[iSlat][1]-lCenSlat[6][1];
    //    lDeltaDiffCenSlat0[iSlat][1] = lDiffCenSlat0[iSlat][1]-lDiffThCenSlat0[iSlat][1];
    cout << lDeltaDiffCenSlat0[iSlat][1] << " " << lDiffCenSlat0[iSlat][1] << " " << lDiffThCenSlat0[iSlat][1] << endl;
  }
  cout << endl;
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    cout << lDeltaDiffCenSlat0[iSlat][2] << " " << lDiffCenSlat0[iSlat][2] << " " << lDiffThCenSlat0[iSlat][2] << endl;
  }
  cout << endl;
  cout << endl;
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    cout << lDeltaDiffRotSlat0[iSlat][0] << " " << lDiffRotSlat0[iSlat][0] << " " << lDiffThRotSlat0[iSlat][0] << endl;
  }
  cout << endl;
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    cout << lDeltaDiffRotSlat0[iSlat][1] << " " << lDiffRotSlat0[iSlat][1] << " " << lDiffThRotSlat0[iSlat][1] << endl;
  }
  cout << endl;
  for(Int_t iSlat=0; iSlat<nSlats; iSlat++){
    cout << lDeltaDiffRotSlat0[iSlat][2] << " " << lDiffRotSlat0[iSlat][2] << " " << lDiffThRotSlat0[iSlat][2] << endl;
  }

  TH1F *mySlatDeltaDiffCenX = new TH1F("mySlatDeltaDiffCenX","mySlatDeltaDiffCenX",100,-10,10);
  mySlatDeltaDiffCenX->SetMaximum(15);
  mySlatDeltaDiffCenX->SetMinimum(0);
  TH1F *mySlatDeltaDiffCenY = new TH1F("mySlatDeltaDiffCenY","mySlatDeltaDiffCenY",100,-10,10);
  mySlatDeltaDiffCenY->SetMaximum(15);
  mySlatDeltaDiffCenY->SetMinimum(0);
  TH1F *mySlatDeltaDiffCenZ = new TH1F("mySlatDeltaDiffCenZ","mySlatDeltaDiffCenZ",100,-20,20);
  mySlatDeltaDiffCenZ->SetMaximum(15);
  mySlatDeltaDiffCenZ->SetMinimum(0);

  TH1F *mySlatDeltaDiffRotX = new TH1F("mySlatDeltaDiffRotX","mySlatDeltaDiffRotX",100,-10,10);
  mySlatDeltaDiffRotX->SetMaximum(15);
  mySlatDeltaDiffRotX->SetMinimum(0);
  TH1F *mySlatDeltaDiffRotY = new TH1F("mySlatDeltaDiffRotY","mySlatDeltaDiffRotY",100,-10,10);
  mySlatDeltaDiffRotY->SetMaximum(15);
  mySlatDeltaDiffRotY->SetMinimum(0);
  TH1F *mySlatDeltaDiffRotZ = new TH1F("mySlatDeltaDiffRotZ","mySlatDeltaDiffRotZ",100,-5,5);
  mySlatDeltaDiffRotZ->SetMaximum(15);
  mySlatDeltaDiffRotZ->SetMinimum(0);
  //
  // ******** Starting plots 
  //
  TCanvas *canvas;
  TPad *pad;
  TPaveLabel *theTitle;

  TPostScript *ps = 0;

  if( saveps ){
    ps = new TPostScript(filename,filetype); 
    ps->NewPage();
  }
 
  // Inv Mass, Multiplicity
  sprintf(str,"Chamber 8L");
  TCanvas *cvn0 = new TCanvas("cvn0",str,cWidth,cHeight);
  canvas = cvn0;
  canvas->Range(0,0,21,29);
 
  TPaveLabel *theTitle0 = new TPaveLabel(3,27.0,18,28.5," Deformations of chamber 8L ","br");
  theTitle = theTitle0;
  theTitle->SetFillColor(18);
  theTitle->SetTextFont(32);
  theTitle->SetTextSize(0.4);
  theTitle->SetTextColor(1);
  theTitle->Draw();
 
  TPad *pad0 = new TPad("pad0","pad0",0.01,0.01,0.98,0.91,0);
  pad = pad0;
  pad->Draw();
  pad->Divide(2,2);

  pad->cd(1);
  gStyle->SetPalette(1);
  hCPSTry->SetMinimum(100);
  hCPSTry->SetMaximum(120);
  hCPSTry->Draw("lego2z");

  pad->cd(2);
  gStyle->SetPalette(1);
  hSSTry->SetMinimum(60);
  hSSTry->SetMaximum(80);
  hSSTry->Draw("lego2z");

  pad->cd(3);
  gStyle->SetPalette(1);
  hCSTy->SetMinimum(110);
  hCSTy->SetMaximum(130);
  hCSTy->Draw("lego2z");

  pad->cd(4);
  gStyle->SetPalette(1);
  hSSTly->SetMinimum(165);
  hSSTly->SetMaximum(185);
  hSSTly->Draw("lego2z");

  if(saveps){
    ps->NewPage();
  }

  // Inv Mass, Multiplicity
  sprintf(str,"Chamber 8L");
  TCanvas *cvn1 = new TCanvas("cvn1",str,cWidth,cHeight);
  canvas = cvn1;
  canvas->Range(0,0,21,29);
  
  TPaveLabel *theTitle1 = new TPaveLabel(3,27.0,18,28.5," Deformations of chamber 8L ","br");
  theTitle = theTitle1;
  theTitle->SetFillColor(18);
  theTitle->SetTextFont(32);
  theTitle->SetTextSize(0.4);
  theTitle->SetTextColor(1);
  theTitle->Draw();
 
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.98,0.91,0);
  pad = pad1;
  pad->Draw();
  pad->Divide(2,2);

  pad->cd(1);
  gStyle->SetPalette(1);
  hCPSTly->SetMinimum(120);
  hCPSTly->SetMaximum(140);
  hCPSTly->Draw("lego2z");

  pad->cd(2);
  gStyle->SetPalette(1);
  hSSTrpy->SetMinimum(60);
  hSSTrpy->SetMaximum(80);
  hSSTrpy->Draw("surf2z");

  pad->cd(3);
  gStyle->SetPalette(1);
  hCSTy->SetMinimum(110);
  hCSTy->SetMaximum(130);
  hCSTy->Draw("lego2z");

  pad->cd(4);
  gStyle->SetPalette(1);
  hSSTlpy->SetMinimum(165);
  hSSTlpy->SetMaximum(185);
  hSSTlpy->Draw("surf2z");

  // Inv Mass, Multiplicity
  sprintf(str,"Chamber 8L");
  TCanvas *cvn2 = new TCanvas("cvn2",str,cWidth,cHeight);
  canvas = cvn2;
  canvas->Range(0,0,21,29);
  
  TPaveLabel *theTitle2 = new TPaveLabel(3,27.0,18,28.5," Deformations of chamber 8L ","br");
  theTitle = theTitle2;
  theTitle->SetFillColor(18);
  theTitle->SetTextFont(32);
  theTitle->SetTextSize(0.4);
  theTitle->SetTextColor(1);
  theTitle->Draw();
 
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.98,0.91,0);
  pad = pad2;
  pad->Draw();
  pad->Divide(3,2);

  pad->cd(1);
  mySlatDeltaDiffCenX->Draw();
  mySlatDeltaDiffCenX->SetXTitle("#Delta[(xc_{i}^{m}-xc_{0}^{m})-(xc_{i}^{th}-xc_{0}^{th})] (mm)");
  mySlatDeltaDiffCenX->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffCenXSlat0->SetMarkerStyle(20);
  gDeltaDiffCenXSlat0->Draw("P");

  pad->cd(2);
  mySlatDeltaDiffCenY->Draw();
  mySlatDeltaDiffCenY->SetXTitle("#Delta[(yc_{i}^{m}-yc_{0}^{m})-(yc_{i}^{th}-yc_{0}^{th})] (mm)");
  mySlatDeltaDiffCenY->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffCenYSlat0->SetMarkerStyle(20);
  gDeltaDiffCenYSlat0->Draw("P");

  pad->cd(3);
  mySlatDeltaDiffCenZ->Draw();
  mySlatDeltaDiffCenZ->SetXTitle("#Delta[(zc_{i}^{m}-zc_{0}^{m})-(zc_{i}^{th}-zc_{0}^{th})] (mm)");
  mySlatDeltaDiffCenZ->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffCenZSlat0->SetMarkerStyle(20);
  gDeltaDiffCenZSlat0->Draw("P");

  pad->cd(4);
  mySlatDeltaDiffRotX->Draw();
  mySlatDeltaDiffRotX->SetXTitle("#Delta[(#psi_{i}^{m}-#psi_{0}^{m})-(#psi_{i}^{th}-#psi_{0}^{th})] (mrad)");
  mySlatDeltaDiffRotX->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffPsiSlat0->SetMarkerStyle(20);
  gDeltaDiffPsiSlat0->Draw("P");

  pad->cd(5);
  mySlatDeltaDiffRotY->Draw();
  mySlatDeltaDiffRotY->SetXTitle("#Delta[(#theta_{i}^{m}-#theta_{0}^{m})-(#theta_{i}^{th}-#theta_{0}^{th})] (mrad)");
  mySlatDeltaDiffRotY->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffThtSlat0->SetMarkerStyle(20);
  gDeltaDiffThtSlat0->Draw("P");

  pad->cd(6);
  mySlatDeltaDiffRotZ->Draw();
  mySlatDeltaDiffRotZ->SetXTitle("#Delta[(#phi_{i}^{m}-#phi_{0}^{m})-(#phi_{i}^{th}-#phi_{0}^{th})] (mrad)");
  mySlatDeltaDiffRotZ->SetYTitle("Slat ascending vertical ordering");
  gDeltaDiffPhiSlat0->SetMarkerStyle(20);
  gDeltaDiffPhiSlat0->Draw("P");

  if( saveps ){
    ps->Close();
  }
}
