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

// Read beam halo background particles from a boundary source
// Boundary source is in the LHCb format http://www.hep.manchester.ac.uk/u/robert/LHC_backgrounds/Note-MIBStudies.pdf
// and has been provided by Robert Appleby
// Author: andreas.morsch@cern.ch


#include <stdlib.h>

#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPDGCode.h>
#include <TSystem.h>

#include "AliGenHalo.h"
#include "AliRun.h"
#include "AliLog.h"

ClassImp(AliGenHalo)

AliGenHalo::AliGenHalo()
    :AliGenerator(-1), 
     fFile(0),
     fFileName(0),
     fSide(1),
     fRunPeriod(kY3D90),
     fTimePerEvent(1.e-4),
     fNskip(0),
     fZ1(0),
     fZ2(0),
     fG1(0),
     fG2(0),
     fGPASize(0),
     fLossID(0),   
     fLossA(0),   
     fPdg(0),     
     fLossT0(0),
     fLossZ(0), 
     fLossW(0), 
     fXS(0),    
     fYS(0),    
     fZS(0),    
     fDX(0),    
     fDY(0),    
     fEkin(0),  
     fTS(0),    
     fWS(0)     
{
// Constructor
    
    fName  = "Halo";
    fTitle = "Halo from LHC Beam";
//
//  Read all particles
    fNpart = -1;
    SetAnalog(0);
}

AliGenHalo::AliGenHalo(Int_t npart)
    :AliGenerator(npart),
     fFile(0),
     fFileName(0),
     fSide(1),
     fRunPeriod(kY3D90),
     fTimePerEvent(1.e-4),
     fNskip(0),
     fZ1(0),
     fZ2(0),
     fG1(0),
     fG2(0),
     fGPASize(0),
     fLossID(0),   
     fLossA(0),   
     fPdg(0),     
     fLossT0(0),
     fLossZ(0), 
     fLossW(0), 
     fXS(0),    
     fYS(0),    
     fZS(0),    
     fDX(0),    
     fDY(0),    
     fEkin(0),  
     fTS(0),    
     fWS(0)     
{
// Constructor
    fName = "Halo";
    fTitle= "Halo from LHC Beam";
//
    fNpart   = npart;
//
    SetAnalog(0);
}

//____________________________________________________________
AliGenHalo::~AliGenHalo()
{
// Destructor
}

//____________________________________________________________
void AliGenHalo::Init() 
{
// Initialisation
    fFile = fopen(fFileName,"r");
    Int_t ir = 0;
    
    
    if (fFile) {
	printf("\n File %s opened for reading, %p ! \n ",  fFileName.Data(), (void*)fFile);
    } else {
	printf("\n Opening of file %s failed,  %p ! \n ",  fFileName.Data(), (void*)fFile);
	return;
    }

    if (fNskip > 0) {
      // Skip the first fNskip events
      SkipEvents();
    }
//
//
//
//    Read file with gas pressure values
    char *name = 0;
    if (fRunPeriod < 5) {
	name = gSystem->ExpandPathName("$(ALICE_ROOT)/LHC/gasPressure.dat" );
	fGPASize = 21;
	fG1 = new Float_t[fGPASize];
	fG2 = new Float_t[fGPASize];
	fZ1 = new Float_t[fGPASize];
	fZ2 = new Float_t[fGPASize];
    } else if (fRunPeriod == 5) {
	name = gSystem->ExpandPathName("$(ALICE_ROOT)/LHC/pressure_2003_startup.dat");
	fGPASize = 18853;
	fG1 = new Float_t[fGPASize];
	fZ1 = new Float_t[fGPASize];
    } else if (fRunPeriod ==6) {
	name = gSystem->ExpandPathName("$(ALICE_ROOT)/LHC/pressure_2003_conditioned.dat");
	fGPASize = 12719;
	fG1 = new Float_t[fGPASize];
	fZ1 = new Float_t[fGPASize];
    } else {
	Fatal("Init()", "No gas pressure file for given run period !");
    }

    
    FILE* file = 0;
    if (name) file = fopen(name, "r");
    if (!file) {
	AliError("No gas pressure file");
	return;
    }
    
    Float_t z;
    Int_t i;
    Float_t p[5];    

    const Float_t kCrossSection = 0.094e-28;      // m^2
    const Float_t kFlux         = 1.e11 / 25.e-9; // protons/s
    Float_t pFlux[5] = {0.2, 0.2, 0.3, 0.3, 1.0};

    if (fRunPeriod < 5) {
//
//  Ring 1   
// 

	for (i = 0; i < fGPASize; i++)
	{
	    ir = fscanf(file, "%f %f %f %f %f %f", &z, &p[0], &p[1], &p[2] , &p[3], &p[4]);
	    if (ir == 0) break;
	    
	    fG1[i] = p[fRunPeriod];
	    
	    if (i > 0) {
		fZ1[i] = fZ1[i-1] + z;
	    } else {
		fZ1[i] = 20.;
	    }
	}
//
// Ring 2
//
	for (i = 0; i < fGPASize; i++)
	{
	    ir = fscanf(file, "%f %f %f %f %f %f", &z, &p[0], &p[1], &p[2] , &p[3], &p[4]);
	    if (ir == 0) break;

	    fG2[i] = p[fRunPeriod];
	    if (i > 0) {
		fZ2[i] = fZ2[i-1] + z;
	    } else {
		fZ2[i] = 20.;
	    }
	}
//
// Interaction rates
//
	for (i = 0; i < fGPASize; i++)  
	{
	    fG1[i] = fG1[i] * kCrossSection * pFlux[fRunPeriod] * kFlux; // 1/m/s 
	    fG2[i] = fG2[i] * kCrossSection * pFlux[fRunPeriod] * kFlux; // 1/m/s
	}

    } else {
	for (i = 0; i < fGPASize; i++) 
	{
	    ir = fscanf(file, "%f %e %e %e %e %e", &z, &p[0], &p[1], &p[2], &p[3], &p[4]);
	    if (ir == 0) break;
	    z /= 1000.;
	    fG1[i] = p[4] * kCrossSection * kFlux;             // 1/m/s
	    // 1/3 of nominal intensity at startup
	    if (fRunPeriod ==  kLHCPR674Startup) fG1[i] /= 3.;
	    fZ1[i] = z;
	}
    }
    


    
//
//  Transform into interaction rates
//


    

    Float_t sum1 = 0.;
    Float_t sum2 = 0.;
    
    for (Int_t iz = 0; iz < 300; iz++) {
	Float_t zpos = 20. + iz * 1.;
	zpos *= 100;
	Float_t wgt1 = GasPressureWeight( zpos);
	Float_t wgt2 = GasPressureWeight(-zpos);
	sum1 += wgt1;
	sum2 += wgt2;
    }
    sum1/=250.;
    sum2/=250.;
    printf("\n %f %f \n \n", sum1, sum2);
    delete file;
}

//____________________________________________________________
void AliGenHalo::Generate()
{
// Generate by reading particles from input file
 
  Float_t polar[3]= {0., 0., 0.};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t tz, txy;
  Float_t mass;
  //
  Int_t nt;
  static Bool_t first = kTRUE;
  static Int_t  oldID = -1;
//

  if (first && (fNskip == 0)) ReadNextParticle();
  first = kFALSE;
  oldID = fLossID;
  
  while(1) {
      // Push particle to stack
      mass = TDatabasePDG::Instance()->GetParticle(fPdg)->Mass();
      p0  = TMath::Sqrt(fEkin * fEkin + 2. * mass * fEkin);
      txy = TMath::Sqrt(fDX * fDX + fDY * fDY);
      if (txy > 1.) {
	  tz = 0.;
      } else {
	  tz = - TMath::Sqrt(1. - txy);
      }
 
      p[0] =  p0 * fDX;
      p[1] =  p0 * fDY;
      p[2] =  p0 * tz;

      origin[0] = fXS;
      origin[1] = fYS;
      origin[2] = 1950.;

      PushTrack(fTrackIt , -1, fPdg , p, origin, polar, fTS - 1950./2.9979e10, kPNoProcess, nt, fWS);
      
      Int_t nc = ReadNextParticle();
      
      if (fLossID != oldID || nc == 0) {
	  oldID = fLossID;
	  break;
      }
  }
  SetHighWaterMark(nt);
}
 

Float_t AliGenHalo::GasPressureWeight(Float_t zPrimary)
{
//
// Return z-dependent gasspressure weight = interaction rate [1/m/s].
//
    Float_t weight = 0.;
    zPrimary /= 100.;        // m
    if (fRunPeriod < 5) {
	Float_t zAbs = TMath::Abs(zPrimary);
	if (zPrimary > 0.) 
	{
	    if (zAbs > fZ1[20]) {
		weight = 2.e4;
	    } else {
		for (Int_t i = 1; i < 21; i++) {
		    if (zAbs < fZ1[i]) {
			weight = fG1[i];
			break;
		    }
		}
	    }
	} else {
	    if (zAbs > fZ2[20]) {
		weight = 2.e4;
	    } else {
		for (Int_t i = 1; i < 21; i++) {
		    if (zAbs < fZ2[i]) {
		    weight = fG2[i];
		    break;
		    }
		}
	    }
	}
    } else {
	Int_t index = TMath::BinarySearch(fGPASize, fZ1, zPrimary);
	weight = fG1[index];
    }
    return weight;
}

void AliGenHalo::Draw(Option_t *)
{
// Draws the gas pressure distribution
    Float_t z[400];
    Float_t p[400];
    
    for (Int_t i = 0; i < 400; i++)
    {
	z[i] = -20000. + Float_t(i) * 100;
	p[i] = GasPressureWeight(z[i]);
    }
    
    TGraph*  gr = new TGraph(400, z, p);   
    TCanvas* c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    c1->cd();
    gr->Draw("AL");
}

Int_t AliGenHalo::ReadNextParticle()
{
    // Read the next particle from the file
    Int_t ncols = fscanf(fFile,"%d %f %f %d %f %d %f %f %f %f %f %f %f %f",
		   &fLossID, &fLossT0, &fLossZ, &fLossA, &fLossW, &fPdg, &fXS, &fYS, &fZS, &fDX, &fDY, &fEkin, &fTS, &fWS);
    fLossZ /= 10.;
    fXS    /= 10.;
    fYS    /= 10.; 
    fZS    /= 10.;   
    fTS    *= 1.e-9;
    return (ncols);
}

void AliGenHalo::SkipEvents()
{
  //
  // Skip the first fNskip events
  Int_t skip = fNskip;
  Int_t oldID = -1;

  while (skip >= 0)
    {
      ReadNextParticle();
      if (oldID != fLossID) {
	oldID = fLossID;
	skip--;
      }
    } 
}
void AliGenHalo::CountEvents()
{
    // Count total number of events
    Int_t nev = 0;
    Int_t oldID = -1;
    Int_t nc = 1;
    while (nc != -1)
    {
	nc = ReadNextParticle();
	if (oldID != fLossID) {
	    oldID = fLossID;
	    nev++;
	    printf("Number of events %10d %10d \n", nev, oldID);
	}
    }


    rewind(fFile);
}


