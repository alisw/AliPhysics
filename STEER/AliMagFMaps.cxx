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
// Magnetic field composed by 3 maps: the L3 magnet, extended region, and
// dipole magnet
// Used in the configuration macros (macros/Config.C, etc.)
// Author: Andreas Morsch <andreas.morsch@cern.ch>
//------------------------------------------------------------------------

#include <TClass.h>
#include <TFile.h>
#include <TMath.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliFieldMap.h"
#include "AliMagFMaps.h"

ClassImp(AliMagFMaps)
    

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps():
  fSolenoid(0),
  fSolenoidUser(0.),
  fL3Option(0)
{
  //
  // Default constructor
  //
  //
  // Don't replicate field information in gAlice
  fFieldMap[0] = fFieldMap[1] = fFieldMap[2] = 0;
}

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps(const char *name, const char *title, Int_t integ, 
                         Float_t factor, Float_t fmax, Int_t map, 
                         Int_t l3):
  AliMagFC(name,title,integ,factor,fmax),
  fSolenoid(0),
  fSolenoidUser(0),
  fL3Option(l3)
{
  //
  // Standard constructor
  //
  fType         = kConMesh;
  fFieldMap[0]  = 0;
  fMap          = map;
  fL3Option     = l3;
  ReadField();
  //
  // Don't replicate field information in gAlice
  for (Int_t i = 0; i < 3; i++)  fFieldMap[i]->SetWriteEnable(0);
  //
}

//_______________________________________________________________________
AliMagFMaps::AliMagFMaps(const AliMagFMaps &magf):
  AliMagFC(magf),
  fSolenoid(0),
  fSolenoidUser(0),
  fL3Option(0)
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
  fgReadField = kTRUE;
}

//_______________________________________________________________________
void AliMagFMaps::ReadField()
{
  //  Read Field Map from file
  //
  //  don't read twice
  //
    if (!fgReadField) return;
    fgReadField = 0;
  //    
  char* fname;
  TFile* file = 0;
  if (fMap == k2kG) {
      fSolenoid = 2.;
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B02.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B02"));
      file->Close();
      delete file;
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
  } else if (fMap == k4kG) {
      fSolenoid = 4.;
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B04.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B04"));
      file->Close();
      delete file;
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
  } else if (fMap == k5kG) {
      fSolenoid = 5.;
      fname = gSystem->ExpandPathName("$(ALICE_ROOT)/data/maps/L3B05.root");
      file = new TFile(fname);
      fFieldMap[0] = dynamic_cast<AliFieldMap*>(file->Get("L3B05"));
      file->Close();
      delete file;
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
  }

  if (!fL3Option) fSolenoidUser = fSolenoid;

}

//_______________________________________________________________________
Float_t AliMagFMaps::SolenoidField() const
{
  //
  // Returns max. L3 (solenoid) field strength 
  // according to field map setting 
  //
  return -Factor()*fSolenoid;
}

//_______________________________________________________________________
void AliMagFMaps::Field(Float_t *x, Float_t *b) const
{
  //
  // Method to calculate the magnetic field
  //
  // --- find the position in the grid ---
  

    
 //    if (!fFieldRead) ReadField();

  //
  // Field Maps have been calculated for the coordinate system in which 
  // the Muon Spectrometer is placed at z > 0
  // Transform coordinates corresponingly 
  //
  
  b[0]=b[1]=b[2]=0;
  Float_t xm[3];
  xm[0] = - x[0];
  xm[1] =   x[1];
  xm[2] = - x[2];
  
  AliFieldMap* map = 0;
  if (fFieldMap[0]->Inside(xm[0], xm[1], xm[2])) {
      map = fFieldMap[0];
      Float_t r = TMath::Sqrt(xm[0] * xm[0] + xm[1] * xm[1]);
      
      if (!fL3Option && TMath::Abs(xm[2]) < 370. && r < 550.) {
      //
      //     Constant L3 field , if this option was selected
      //
	b[2] = (- fSolenoid)*fFactor;
	  return;
      } 
  } else if (fFieldMap[1]->Inside(xm[0], xm[1], xm[2])) {
    map = fFieldMap[1];
  } else if (fFieldMap[2]->Inside(xm[0], xm[1], xm[2])) {
    map = fFieldMap[2];
  }
  
  if(map){
    map->Field(xm,b);
    b[0] = - b[0];
    b[2] = - b[2];

  } else {
      //This is the ZDC part
      ZDCField(x, b);
  }

  
  if(fFactor!=1) {
      b[0]*=fFactor;
      b[1]*=fFactor;
      b[2]*=fFactor;
  }
}

//_______________________________________________________________________
void AliMagFMaps::Copy(TObject & /* magf */) const
{
  //
  // Copy *this onto magf -- Not implemented
  //
  AliFatal("Not implemented!");
}

//_______________________________________________________________________
void AliMagFMaps::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliMagFMaps.
  if (R__b.IsReading()) {
    AliMagFMaps::Class()->ReadBuffer(R__b, this);
    ReadField();
  } else {
    AliMagFMaps::Class()->WriteBuffer(R__b, this);
  }
}
