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
Revision 1.8  2002/10/14 14:57:32  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.5.6.2  2002/08/01 15:35:30  alibrary
Updating VirtualMC

Revision 1.7  2002/08/01 12:34:48  morsch
Initialize maps in first call to Field() after streaming in.

Revision 1.6  2002/04/11 11:17:48  morsch
ReadField() used in Constructor.

Revision 1.5  2002/02/26 09:48:14  morsch
Extra argument in constructor for l3 map choice.

Revision 1.4  2002/02/22 14:00:20  morsch
Protection against replication of fieldmap data in gAlice.

Revision 1.3  2002/02/21 09:23:41  morsch
Create dummy field map for L3 in case no detailed map is needed.

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

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps():
  fSolenoid(0),
  fL3Option(0),
  fFieldRead(0)
{
  //
  // Default constructor
  //
  fFieldMap[0] = fFieldMap[1] = fFieldMap[2] = 0;
}

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps(const char *name, const char *title, const Int_t integ, 
                         const Float_t factor, const Float_t fmax, const Int_t map, 
                         const Int_t l3):
  AliMagF(name,title,integ,factor,fmax),
  fSolenoid(0),
  fL3Option(l3),
  fFieldRead(0)
{
  //
  // Standard constructor
  //
  fType         = kConMesh;
  fFieldMap[0]  = 0;
  fMap          = map;
  fL3Option     = l3;

  ReadField();
  fFieldRead = 1;
  //
  // Don't replicate field information in gAlice
  for (Int_t i = 0; i < 3; i++)  fFieldMap[i]->SetWriteEnable(0);
  //
}

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps(const AliMagFMaps &magf):
  AliMagF(magf),
  fSolenoid(0),
  fL3Option(0),
  fFieldRead(0)
{
  //
  // Copy constructor
  //
  magf.Copy(*this);
}

//_______________________________________________________________________
AliMagFMaps::~AliMagFMaps()
{
  //
  //  Destructor
  //
  delete fFieldMap[0];
  delete fFieldMap[1];
  delete fFieldMap[2];    
}

//_______________________________________________________________________
void AliMagFMaps::ReadField()
{
  //  Read Field Map from file
  //
  //  don't read twice
  //
  if (fFieldRead) return;
  fFieldRead = 1;
  //    
  char* fname;
  TFile* file = 0;
  if (fMap == k2kG) {
	if (fL3Option) {
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B02.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B02"));
      file->Close();
      delete file;
	}
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB02.root");
	file = new TFile(fname);
	fFieldMap[1] = dynamic_cast<AliFieldMap*>(file->Get("DipB02"));
	file->Close();
	delete file;;
	
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB02.root");
	file = new TFile(fname);
	fFieldMap[2] = dynamic_cast<AliFieldMap*>(file->Get("ExtB02"));
	file->Close();
	delete file;
	fSolenoid = 2.;
  } else if (fMap == k4kG) {
	if (fL3Option) {
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B04.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B04"));
      file->Close();
      delete file;
	}
	
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB04.root");
	file = new TFile(fname);
	fFieldMap[1] = dynamic_cast<AliFieldMap*>(file->Get("DipB04"));
	file->Close();
	delete file;
	
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB04.root");
	file = new TFile(fname);
	fFieldMap[2] = dynamic_cast<AliFieldMap*>(file->Get("ExtB04"));
	file->Close();
	delete file;
	fSolenoid = 4.;
  } else if (fMap == k5kG) {
	if (fL3Option) {
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B05.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B05"));
      file->Close();
      delete file;
	}
	
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/DipB05.root");
	file = new TFile(fname);
	fFieldMap[1] = dynamic_cast<AliFieldMap*>(file->Get("DipB05"));
	file->Close();
	delete file;
	
	fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/ExtB05.root");
	file = new TFile(fname);
	fFieldMap[2] = dynamic_cast<AliFieldMap*>(file->Get("ExtB05"));
	file->Close();
	delete file;
	
	fSolenoid = 5.;
  }
  
  if (!fL3Option) {
    //
    // Dummy L3 map
	fFieldMap[0] = new AliFieldMap();
	fFieldMap[0] -> SetLimits(-800., 800., -800., 800., -700., 700.);
  }
}

//_______________________________________________________________________
Float_t AliMagFMaps::SolenoidField() const
{
  //
  // Returns max. L3 (solenoid) field strength 
  // according to field map setting 
  //
  return fSolenoid;
}

//_______________________________________________________________________
void AliMagFMaps::Field(Float_t *x, Float_t *b)
{
  //
  // Method to calculate the magnetic field
  //
  // --- find the position in the grid ---
  
  if (!fFieldRead) ReadField();
  
  b[0]=b[1]=b[2]=0;
  AliFieldMap* map = 0;
  if (fFieldMap[0]->Inside(x[0], x[1], x[2])) {
    map = fFieldMap[0];
    if (!fL3Option) {
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

//_______________________________________________________________________
void AliMagFMaps::Copy(AliMagFMaps & /* magf */) const
{
  //
  // Copy *this onto magf -- Not implemented
  //
  Fatal("Copy","Not implemented!\n");
}

//_______________________________________________________________________
void AliMagFMaps::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliMagFMaps.
  if (R__b.IsReading()) {
    AliMagFMaps::Class()->ReadBuffer(R__b, this);
    fFieldRead = 0;
  } else {
    AliMagFMaps::Class()->WriteBuffer(R__b, this);
  }
}
