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

/*
$Log$
Revision 1.2  2002/02/19 16:14:35  morsch
Reading of 0.2 T solenoid field map enabled.

Revision 1.1  2002/02/14 11:41:28  morsch
Magnetic field map for ALICE for L3+muon spectrometer stored in 3 seperate
root files.

*/

//
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//

#include <TFile.h>
#include <TSystem.h>
#include "AliFieldMap.h"
#include "AliMagFMaps.h"


ClassImp(AliMagFMaps)

//________________________________________
AliMagFMaps::AliMagFMaps(const char *name, const char *title, const Int_t integ, 
		     const Float_t factor, const Float_t fmax, const Int_t map)
  : AliMagF(name,title,integ,factor,fmax)
{
  //
  // Standard constructor
  //
  fType         = kConMesh;
  fFieldMap[0]  = 0;
  char* fname;
  
  fMap = map;
  TFile* file = 0;
  if (fMap == k2kG) {
      if (fL3Option) {
	  fFieldMap[0] = new AliFieldMap();
	  fFieldMap[0]->SetLimits(-800., 800., -800., 800., -700., 700.);
      } else {
	  fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B02.root");
	  file = new TFile(fname);
	  fFieldMap[0] = (AliFieldMap*) file->Get("L3B02");
	  file->Close();
	  delete file;
	  
	  fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB02.root");
	  file = new TFile(fname);
	  fFieldMap[1] = (AliFieldMap*) file->Get("DipB02");
	  file->Close();
	  delete file;;
	  
	  fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB02.root");
	  file = new TFile(fname);
	  fFieldMap[2] = (AliFieldMap*) file->Get("ExtB02");
	  file->Close();
	  delete file;
      }
	fSolenoid = 2.;
    } else if (fMap == k4kG) {
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B04.root");
	file = new TFile(fname);
	fFieldMap[0] = (AliFieldMap*) file->Get("L3B04");
	file->Close();
	delete file;

	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB04.root");
	file = new TFile(fname);
	fFieldMap[1] = (AliFieldMap*) file->Get("DipB04");
	file->Close();
	delete file;;

	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB04.root");
	file = new TFile(fname);
	fFieldMap[2] = (AliFieldMap*) file->Get("ExtB04");
	file->Close();
	delete file;
	fSolenoid = 4.;
    } else if (fMap == k5kG) {
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B05.root");
	file = new TFile(fname);
	fFieldMap[0] = (AliFieldMap*) file->Get("L3B05");
	file->Close();
	delete file;

	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB05.root");
	file = new TFile(fname);
	fFieldMap[1] = (AliFieldMap*) file->Get("DipB05");
	file->Close();
	delete file;;

	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB05.root");
	file = new TFile(fname);
	fFieldMap[2] = (AliFieldMap*) file->Get("ExtB05");
	file->Close();
	delete file;

	fSolenoid = 5.;
    }
    SetL3ConstField(0);
}

//________________________________________
AliMagFMaps::AliMagFMaps(const AliMagFMaps &magf)
{
  //
  // Copy constructor
  //
  magf.Copy(*this);
}

AliMagFMaps::~AliMagFMaps()
{
//
//  Destructor
//
    delete fFieldMap[0];
    delete fFieldMap[1];
    delete fFieldMap[2];    
}


Float_t AliMagFMaps::SolenoidField() const
{
//
// Returns max. L3 (solenoid) field strength 
// according to field map setting

    return fSolenoid;
}

    

//________________________________________
void AliMagFMaps::Field(Float_t *x, Float_t *b)
{
  //
  // Method to calculate the magnetic field
  //
  const Double_t kone=1;
  // --- find the position in the grid ---
  
  b[0]=b[1]=b[2]=0;
  AliFieldMap* map = 0;

  if (fFieldMap[0]->Inside(x[0], x[1], x[2])) {
      map = fFieldMap[0];
      if (fL3Option) {
//
//     Constant L3 field, if this option was selected
//
	  b[2] = fSolenoid;
	  return;
      }
  } else if (fFieldMap[1]->Inside(x[0], x[1], x[2])) {
      map = fFieldMap[1];
  } else if (fFieldMap[2]->Inside(x[0], x[1], x[2])) {
      map = fFieldMap[2];
  }
  
  if(map){
      map->Field(x,b);
  } else {
//This is the ZDC part
      Float_t rad2=x[0]*x[0]+x[1]*x[1];
      if(x[2]>kCORBEG2 && x[2]<kCOREND2){
	  if(rad2<kCOR2RA2){
	      b[0] = kFCORN2;
	  }
      }
      else if(x[2]>kZ1BEG && x[2]<kZ1END){  
	  if(rad2<kZ1RA2){
	      b[0] = -kG1*x[1];
	      b[1] = -kG1*x[0];
	  }
      }
      else if(x[2]>kZ2BEG && x[2]<kZ2END){  
	  if(rad2<kZ2RA2){
	      b[0] = kG1*x[1];
	      b[1] = kG1*x[0];
	  }
      }
      else if(x[2]>kZ3BEG && x[2]<kZ3END){  
	  if(rad2<kZ3RA2){
	      b[0] = kG1*x[1];
	      b[1] = kG1*x[0];
	  }
      }
      else if(x[2]>kZ4BEG && x[2]<kZ4END){  
	  if(rad2<kZ4RA2){
	      b[0] = -kG1*x[1];
	      b[1] = -kG1*x[0];
	  }
      }
      else if(x[2]>kD1BEG && x[2]<kD1END){ 
	  if(rad2<kD1RA2){
	      b[1] = -kFDIP;
	  }
      }
      else if(x[2]>kD2BEG && x[2]<kD2END){
	  if(((x[0]-kXCEN1D2)*(x[0]-kXCEN1D2)+(x[1]-kYCEN1D2)*(x[1]-kYCEN1D2))<kD2RA2
	     || ((x[0]-kXCEN2D2)*(x[0]-kXCEN2D2)+(x[1]-kYCEN2D2)*(x[1]-kYCEN2D2))<kD2RA2){
	      b[1] = kFDIP;
	  }
      }
  }
  if(fFactor!=1) {
      b[0]*=fFactor;
      b[1]*=fFactor;
      b[2]*=fFactor;
  }
}

//________________________________________
void AliMagFMaps::Copy(AliMagFMaps & /* magf */) const
{
  //
  // Copy *this onto magf -- Not implemented
  //
  Fatal("Copy","Not implemented!\n");
}

//________________________________________
AliMagFMaps & AliMagFMaps::operator =(const AliMagFMaps &magf)
{
  magf.Copy(*this);
  return *this;
}
