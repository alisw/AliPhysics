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

// $Id$

/// \ingroup macros
/// \file MUONChamberMaterialBudget.C
/// \brief Utility macro to check tracking chamber material properties. 
///
/// \author Philippe Pillot, Subatech, Oct. 2009  

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TH2F.h>
#include <TFile.h>
#include <TStyle.h>
#include <TBox.h>
#include <TString.h>
#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoMaterial.h>
#include <TGeoShape.h>

#endif

void MUONChamberMaterialBudget(Bool_t warn = kFALSE)
{
  /// draw the local chamber thickness over x0 (x/x0) used in the computation of Multiple Coulomb Scattering effets.
  /// Compute <x> and <x/x0> in a limited area (displayed on the histograms) to avoid edge effets.
  
  Double_t OneOverX0MeanCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t totalLengthCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t OneOverX0Mean = 0.;
  Double_t totalLength = 0.;
  
  // Import TGeo geometry
  if (!gGeoManager) {
    TGeoManager::Import("geometry.root");
    if (!gGeoManager) {
      cout<<"getting geometry from file geometry.root failed"<<endl;
      return;
    }
  }
  
  // z intervals where to find the stations
  Double_t zIn[5] = {-510., -600., -850., -1250., -1350.};
  Double_t zOut[5] = {-600., -750., -1150., -1350., -1470.};
  
  // transverse area where to compute locally x and x/x0
  Double_t xIn0[5] = {0., 0., 0., 0., 0.};
  Double_t yIn0[5] = {0., 0., 0., 0., 0.};
  Int_t ixMax[5] = {90, 120, 180, 240, 270};
  Int_t iyMax[5] = {90, 120, 180, 240, 270};
  
  // transverse area where to compute <x> and <x/x0> for each chamber
  Double_t xIn0ForMean[5] = { 5.,  0.,  40.,  45.,  50.};
  Double_t yIn0ForMean[5] = {21., 30.,   0.,   0.,   0.};
  Int_t ixMaxForMean[5] = {  50,  70,   80,  115,  150 };
  Int_t iyMaxForMean[5] = {  60,  70,  110,  150,  150 };
  
  // define output histograms
  gStyle->SetPalette(1);
  TFile *f = TFile::Open("MaterialBudget.root","RECREATE");
  TH2F* hXOverX0[10];
  TBox* bXOverX0[10];
  for (Int_t i=0; i<10; i++) {
    Int_t st = i/2;
    hXOverX0[i] = new TH2F(Form("hXOverX0_%d",i+1), Form("x/x0 on ch %d (%%)",i+1),
			   ixMax[st], xIn0[st], xIn0[st]+ixMax[st],
			   iyMax[st], yIn0[st], yIn0[st]+iyMax[st]);
    hXOverX0[i]->SetOption("COLZ");
    bXOverX0[i] = new TBox(xIn0ForMean[st], yIn0ForMean[st],
			   xIn0ForMean[st]+ixMaxForMean[st], yIn0ForMean[st]+iyMaxForMean[st]);
    bXOverX0[i]->SetLineStyle(2);
    bXOverX0[i]->SetLineWidth(2);
    bXOverX0[i]->SetFillStyle(0);
    hXOverX0[i]->GetListOfFunctions()->Add(bXOverX0[i]);
  }
  
  // loop over stations
  for (Int_t ist=0; ist<5; ist++) {
    
    Int_t nPoints = 0;
    
    // loop over position in non bending direction (by step of 1cm)
    for (Int_t ix=0; ix<ixMax[ist]; ix++) {
      
      Double_t xIn = xIn0[ist] + 1.*ix;
      
      // loop over position in bending direction (by step of 1cm)
      for (Int_t iy=0; iy<iyMax[ist]; iy++) {
	
	Double_t yIn = yIn0[ist] + 1.*iy;
	
	// Initialize starting point and direction
	Double_t trackXYZIn[3] = {xIn, yIn, zIn[ist]};
	Double_t trackXYZOut[3] = {1.000001*xIn, 1.000001*yIn, zOut[ist]};
	//  Double_t trackXYZIn[3] = {35., -45., -510.};
	//  Double_t trackXYZOut[3] = {100., -138., -1500.};
	Double_t pathLength = TMath::Sqrt((trackXYZOut[0] - trackXYZIn[0])*(trackXYZOut[0] - trackXYZIn[0])+
					  (trackXYZOut[1] - trackXYZIn[1])*(trackXYZOut[1] - trackXYZIn[1])+
					  (trackXYZOut[2] - trackXYZIn[2])*(trackXYZOut[2] - trackXYZIn[2]));
	Double_t b[3];
	b[0] = (trackXYZOut[0] - trackXYZIn[0]) / pathLength;
	b[1] = (trackXYZOut[1] - trackXYZIn[1]) / pathLength;
	b[2] = (trackXYZOut[2] - trackXYZIn[2]) / pathLength;
	TGeoNode *currentnode = gGeoManager->InitTrack(trackXYZIn, b);
	if (!currentnode) {
	  if (warn) cout<<"start point out of geometry"<<endl;
	  break;
	}
	
	Bool_t OK = kTRUE;
	Double_t x0 = 0.;  // radiation-length (cm-1)
	Double_t localOneOverX0MeanCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t localTotalLengthCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t localPathLength = 0.;
	Double_t remainingPathLength = pathLength;
	Int_t chId = 2*ist;
	Bool_t enterCh = kFALSE;
	do {
	  // Get material properties
	  TGeoMaterial *material = currentnode->GetVolume()->GetMedium()->GetMaterial();
	  TString mName(material->GetName());
	  mName.ToUpper();
	  x0 = material->GetRadLen();
	  
	  // Get path length within this material
	  gGeoManager->FindNextBoundary(remainingPathLength);
	  localPathLength = gGeoManager->GetStep() + 1.e-6;
	  // Check if boundary within remaining path length. If so, make sure to cross the boundary to prepare the next step
	  if (localPathLength >= remainingPathLength) localPathLength = remainingPathLength;
	  else {
	    currentnode = gGeoManager->Step();
	    if (!currentnode) {
	      if (warn) cout<<"navigation failed(1)"<<endl;
	      OK = kFALSE;
	      break;
	    }
	    if (!gGeoManager->IsEntering()) {
	      // make another small step to try to enter in new absorber slice
	      gGeoManager->SetStep(0.001);
	      currentnode = gGeoManager->Step();
	      if (!gGeoManager->IsEntering() || !currentnode) {
		if (warn) cout<<"navigation failed(2)"<<endl;
		OK = kFALSE;
		break;
	      }
	      localPathLength += 0.001;
	    }
	  }
	  
	  // check if entering a new chamber
	  if (localPathLength > 10.) {
	    if(!mName.Contains("AIR")) {
	      if (warn) {
		material->Print();
	        cout<<"changing chamber while crossing material"<<endl;
	      }
	      OK = kFALSE;
	      break;
	    }
	    if (enterCh) chId++;
	    //if (chId > 2*ist+1) break;
	    enterCh = kFALSE;
	    remainingPathLength -= localPathLength;
	    continue;
	  } else {
	    enterCh = kTRUE;
	  }
	  
	  // print parameters
	  if(!mName.Contains("AIR")) {
	    if(mName.Contains("DIPO")) {
	      if (warn) cout<<"the track cross the dipole"<<endl;
	      OK = kFALSE;
	      break;
	    }
	    //material->Print();
	    //cout<<"localPathLength = "<<localPathLength<<endl;
	    if (chId <= 2*ist+1) {
	      localOneOverX0MeanCh[chId] += localPathLength / x0;
	      localTotalLengthCh[chId] += localPathLength;
	    }
	  }
	  
	  // prepare next step
	  remainingPathLength -= localPathLength;
	} while (remainingPathLength > TGeoShape::Tolerance());
	
	// account for the local material characteristic if computed successfully
	if (OK && (chId == 2*ist+2)) {
	  
	  // fill histograms in the full space
	  hXOverX0[2*ist]->Fill(xIn,yIn,100.*localOneOverX0MeanCh[2*ist]);
	  hXOverX0[2*ist+1]->Fill(xIn,yIn,100.*localOneOverX0MeanCh[2*ist+1]);
	  
	  // limit the chamber region for the computation of <x> and <x/x0>
	  if (xIn > xIn0ForMean[ist]-0.5 && xIn < xIn0ForMean[ist]+1.*ixMaxForMean[ist]-0.5 &&
	      yIn > yIn0ForMean[ist]-0.5 && yIn < yIn0ForMean[ist]+1.*iyMaxForMean[ist]-0.5) {
	    nPoints++;
	    OneOverX0MeanCh[2*ist] += localOneOverX0MeanCh[2*ist];
	    OneOverX0MeanCh[2*ist+1] += localOneOverX0MeanCh[2*ist+1];
	    totalLengthCh[2*ist] += localTotalLengthCh[2*ist];
	    totalLengthCh[2*ist+1] += localTotalLengthCh[2*ist+1];
	  }
	  
	}
	
      }
      
    }
    
    // normalize <x> and <x/x0> to the number of data points
    for (Int_t i=2*ist; i<=2*ist+1; i++) {
      OneOverX0MeanCh[i] /= (Double_t) nPoints;
      totalLengthCh[i] /= (Double_t) nPoints;
    }
    
  }
  
  // print results
  cout<<endl<<endl;
  cout<<"chamber   thickness (cm)    x/x0 (%)"<<endl;
  cout<<"------------------------------------"<<endl;
  for (Int_t i=0; i<10; i++) {
    printf("  %2d          %4.2f            %4.2f\n",i+1,totalLengthCh[i],100.*OneOverX0MeanCh[i]);
    totalLength += totalLengthCh[i];
    OneOverX0Mean += OneOverX0MeanCh[i];
  }
  cout<<"------------------------------------"<<endl;
  printf("  tot         %4.1f            %4.1f\n",totalLength,100.*OneOverX0Mean);
  cout<<endl;
  
  // save histograms
  f->Write();
  f->Close();
  
}

