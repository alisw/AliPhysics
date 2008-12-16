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

//------------------------------------------------------------------------
// Magnetic field composed by 4 maps: the L3 magnet (inside and outside measured region), 
// extended region, and dipole magnet.
// Used in the configuration macros (macros/Config.C, etc.)
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//------------------------------------------------------------------------

#include <TClass.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliFieldMap.h"
#include "AliMagFMapsV1.h"

ClassImp(AliMagFMapsV1)
    

//_______________________________________________________________________
AliMagFMapsV1::AliMagFMapsV1():
  AliMagFMaps(),
  fMeasuredMap(0) 
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliMagFMapsV1::AliMagFMapsV1(const char *name, const char *title, Int_t integ, 
			     Float_t factor, Float_t fmax, Int_t map, 
			     Int_t l3):
    AliMagFMaps(name, title, integ, factor, fmax, map, l3),
    fMeasuredMap(0) 
{
    //
    // Constructor
    //
    char* fname;
    fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/mfcheb.root");
    TFile* file = new TFile(fname);
    if (fMap == k2kG) {
	fMeasuredMap = dynamic_cast<AliMagFCheb*>(file->Get("Sol12_Dip6_Hole"));
	fSolenoid = 0.2; // T
    } else if (fMap == k5kG) {
	fMeasuredMap = dynamic_cast<AliMagFCheb*>(file->Get("Sol30_Dip6_Hole"));
	fSolenoid = 0.5; // T
    } else if (fMap == k4kG){
	fMeasuredMap = 0;
	fSolenoid = 0.4; // T
    }
    
    
    
    file->Close();
    delete file;
}


//_______________________________________________________________________
AliMagFMapsV1::AliMagFMapsV1(const AliMagFMapsV1 &magf):
  AliMagFMaps(magf),
  fMeasuredMap(0)
{
  //
  // Copy constructor
  //
  magf.Copy(*this);
}

//_______________________________________________________________________
AliMagFMapsV1::~AliMagFMapsV1()
{
    // Destructor
    delete fMeasuredMap;
}

//_______________________________________________________________________
void AliMagFMapsV1::Field(float *x, float *b) const
{
  //
  // Method to calculate the magnetic field at position x
  //
    const float kRmax2 = 500. * 500.;
    const float kZmax  = 550.; 
    const float kTeslaTokG = 10.;
    const float kScale = 0.98838; // matching factor
    
    // Check if position inside measured map
    float r2 = x[0] * x[0] + x[1] * x[1];
    if (fMeasuredMap              &&
	r2 < kRmax2               && 
	TMath::Abs(x[2]) < kZmax
	) 
    {
	fMeasuredMap->Field(x, b);
	b[0] *= kTeslaTokG;
	b[1] *= kTeslaTokG;
	b[2] *= kTeslaTokG;
    } else {
	AliMagFMaps::Field(x, b);
	// Match to measure map
	b[0] = - b[0] * kScale;
	b[2] = - b[2] * kScale;
	b[1] = - b[1] * kScale;
    }
}

//_______________________________________________________________________
void AliMagFMapsV1::Field(double *x, double *b) const
{
  //
  // Method to calculate the magnetic field at position x
  //
    const double kRmax2 = 500. * 500.;
    const double kZmax  = 550.; 
    const double kTeslaTokG = 10.;
    const double kScale = 0.98838; // matching factor
    
    // Check if position inside measured map
    double r2 = x[0] * x[0] + x[1] * x[1];
    if (fMeasuredMap              &&
	r2 < kRmax2               && 
	TMath::Abs(x[2]) < kZmax
	) 
    {
	fMeasuredMap->Field(x, b);
	b[0] *= kTeslaTokG;
	b[1] *= kTeslaTokG;
	b[2] *= kTeslaTokG;
    } else {
	AliMagFMaps::Field(x, b);
	// Match to measure map
	b[0] = - b[0] * kScale;
	b[2] = - b[2] * kScale;
	b[1] = - b[1] * kScale;
    }
}


Float_t AliMagFMapsV1::SolenoidField() const
{
  //
  // Returns max. L3 (solenoid) field strength 
  // according to field map setting 
  //
	return fSolenoid;
}

void AliMagFMapsV1::Copy(TObject &fm) const
{
  // dummy since the base class constructor is dummy
  AliMagFMaps::Copy(fm);
  AliMagFMapsV1* fmp = (AliMagFMapsV1*) &fm;
  fmp->SetMeasuredMap(fMeasuredMap);
}
