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

void MUONChamberMaterialBudget(const char* geoFilename = "geometry.root", Int_t segmentationLevel = 1)
{
  /// Draw the local chamber thickness over x0 (x/x0) used in the computation of Multiple Coulomb Scattering effets.
  /// Compute <x> and <x/x0> in a limited area (displayed on the histograms) to avoid edge effets.
  /// The resolution can be changed by changing the sementation level: resolution = 1 cm / segmentationLevel.
  
  const char* chamberName[10] = {"SC01", "SC02", "SC03", "SC04", "SC05", "SC06", "SC07", "SC08", "SC09", "SC10"};
  Double_t OneOverX0MeanCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t totalLengthCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t OneOverX0Mean = 0.;
  Double_t totalLength = 0.;
  
  // Import TGeo geometry
  if (!gGeoManager) {
    TGeoManager::Import(geoFilename);
    if (!gGeoManager) {
      cout<<"getting geometry from file geometry.root failed"<<endl;
      return;
    }
  }
  
  // z intervals where to find the stations
  Double_t zIn[5] =  {-510., -600.,  -800., -1150., -1350.};
  Double_t zOut[5] = {-600., -800., -1150., -1350., -1470.};
  
  // transverse area where to compute locally x and x/x0
  Double_t xIn0[5] = {0., 0., 0., 0., 0.};
  Double_t yIn0[5] = {0., 0., 0., 0., 0.};
  Int_t ixMax[5] = {90, 120, 165, 250, 260};
  Int_t iyMax[5] = {90, 120, 180, 250, 270};
  
  // transverse area where to compute <x> and <x/x0> for each chamber
  Double_t xIn0ForMean[5] = { 5.,  5.,  35.,  40.,  40.};
  Double_t yIn0ForMean[5] = {20., 25.,   0.,   0.,   0.};
  Int_t ixMaxForMean[5] = {  50,  65,   85,  120,  160 };
  Int_t iyMaxForMean[5] = {  60,  70,  110,  150,  150 };
  
  // define output histograms
  gStyle->SetPalette(1);
  TFile *f = TFile::Open("MaterialBudget.root","RECREATE");
  TH2F* hXOverX0[10];
  TBox* bXOverX0[10];
  for (Int_t i=0; i<10; i++) {
    Int_t st = i/2;
    hXOverX0[i] = new TH2F(Form("hXOverX0_%d",i+1), Form("x/x0 on ch %d (%%)",i+1),
			   segmentationLevel*ixMax[st], xIn0[st], xIn0[st]+ixMax[st],
			   segmentationLevel*iyMax[st], yIn0[st], yIn0[st]+iyMax[st]);
    hXOverX0[i]->SetOption("COLZ");
    hXOverX0[i]->SetStats(kFALSE);
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
    for (Int_t ix=0; ix<segmentationLevel*ixMax[ist]; ix++) {
      
      Double_t xIn = xIn0[ist] + ((Double_t)ix+0.5) / ((Double_t)segmentationLevel);
      
      // loop over position in bending direction (by step of 1cm)
      for (Int_t iy=0; iy<segmentationLevel*iyMax[ist]; iy++) {
	Int_t permilDone = 1000 * (ix * segmentationLevel*iyMax[ist] + iy + 1) /
	                   (segmentationLevel*segmentationLevel*ixMax[ist]*iyMax[ist]);
	if (permilDone%10 == 0) cout<<"\rStation "<<ist+1<<": processing... "<<permilDone/10<<"%"<<flush;
	
	Double_t yIn = yIn0[ist] + ((Double_t)iy+0.5) / ((Double_t)segmentationLevel);
	
	// Initialize starting point and direction
	Double_t trackXYZIn[3] = {xIn, yIn, zIn[ist]};
	Double_t trackXYZOut[3] = {1.000001*xIn, 1.000001*yIn, zOut[ist]};
	Double_t pathLength = TMath::Sqrt((trackXYZOut[0] - trackXYZIn[0])*(trackXYZOut[0] - trackXYZIn[0])+
					  (trackXYZOut[1] - trackXYZIn[1])*(trackXYZOut[1] - trackXYZIn[1])+
					  (trackXYZOut[2] - trackXYZIn[2])*(trackXYZOut[2] - trackXYZIn[2]));
	Double_t b[3];
	b[0] = (trackXYZOut[0] - trackXYZIn[0]) / pathLength;
	b[1] = (trackXYZOut[1] - trackXYZIn[1]) / pathLength;
	b[2] = (trackXYZOut[2] - trackXYZIn[2]) / pathLength;
	TGeoNode *currentnode = gGeoManager->InitTrack(trackXYZIn, b);
	if (!currentnode) break;
	
	Bool_t OK = kTRUE;
	Double_t x0 = 0.;  // radiation-length (cm-1)
	Double_t localOneOverX0MeanCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t localTotalLengthCh[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Double_t localPathLength = 0.;
	Double_t remainingPathLength = pathLength;
	do {
	  // Get material properties
	  TGeoMaterial *material = currentnode->GetVolume()->GetMedium()->GetMaterial();
	  TString currentNodePath(gGeoManager->GetPath());
	  x0 = material->GetRadLen();
	  
	  // Get path length within this material
	  gGeoManager->FindNextBoundary(remainingPathLength);
	  localPathLength = gGeoManager->GetStep() + 1.e-6;
	  
	  // Check if boundary within remaining path length.
	  // If so, make sure to cross the boundary to prepare the next step
	  if (localPathLength >= remainingPathLength) localPathLength = remainingPathLength;
	  else {
	    currentnode = gGeoManager->Step();
	    if (!currentnode) { OK = kFALSE; break; }
	    if (!gGeoManager->IsEntering()) {
	      // make another small step to try to enter in new slice
	      gGeoManager->SetStep(0.001);
	      currentnode = gGeoManager->Step();
	      if (!gGeoManager->IsEntering() || !currentnode) { OK = kFALSE; break; }
	      localPathLength += 0.001;
	    }
	  }
	  remainingPathLength -= localPathLength;

	  // check if entering a chamber of the current station or go to next step
	  Int_t chId;
	  if (currentNodePath.Contains(chamberName[2*ist])) chId = 2*ist;
	  else if (currentNodePath.Contains(chamberName[2*ist+1])) chId = 2*ist+1;
	  else continue;
	  
	  // add current material budget
	  localOneOverX0MeanCh[chId] += localPathLength / x0;
	  localTotalLengthCh[chId] += localPathLength;
	  
	} while (remainingPathLength > TGeoShape::Tolerance());
	
	// account for the local material characteristic if computed successfully
	if (OK) {
	  
	  // fill histograms in the full space
	  hXOverX0[2*ist]->Fill(xIn,yIn,100.*localOneOverX0MeanCh[2*ist]);
	  hXOverX0[2*ist+1]->Fill(xIn,yIn,100.*localOneOverX0MeanCh[2*ist+1]);
	  
	  // computation of <x> and <x/x0> in a limited chamber region
	  if (xIn > xIn0ForMean[ist] && xIn < xIn0ForMean[ist]+1.*ixMaxForMean[ist] &&
	      yIn > yIn0ForMean[ist] && yIn < yIn0ForMean[ist]+1.*iyMaxForMean[ist]) {
	    nPoints++;
	    OneOverX0MeanCh[2*ist] += localOneOverX0MeanCh[2*ist];
	    OneOverX0MeanCh[2*ist+1] += localOneOverX0MeanCh[2*ist+1];
	    totalLengthCh[2*ist] += localTotalLengthCh[2*ist];
	    totalLengthCh[2*ist+1] += localTotalLengthCh[2*ist+1];
	  }
	  
	}
	
      }
      
    }
    cout<<"\rStation "<<ist+1<<": processing... 100%"<<endl;
    
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

