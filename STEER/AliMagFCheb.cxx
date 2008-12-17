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


#include <TClass.h>
#include <TFile.h>
#include <TSystem.h>

#include "AliMagFCheb.h"
#include "AliMagWrapCheb.h"
#include "AliLog.h"

ClassImp(AliMagFCheb)
    

//_______________________________________________________________________
AliMagFCheb::AliMagFCheb():
  AliMagFC(),
  fMeasuredMap(0),
  fSolenoid(5.)
{
  // Default constructor
  //
}

//_______________________________________________________________________
AliMagFCheb::AliMagFCheb(const char *name, const char *title, Int_t integ, 
			 Float_t factor, Float_t fmax, Int_t map, 
			 Bool_t dipoleON,const char* path):
  AliMagFC(name, title, integ, factor, fmax),
  fMeasuredMap(0),
  fSolenoid(5.)
{
  //
  fMap = map;
  char* fname = gSystem->ExpandPathName(path);
  TFile* file = TFile::Open(fname);
  if (!file) {
    AliError(Form("Failed to open magnetic field data file %s\n",fname)); 
    return;
  }
  const char* parname = 0;
  if        (fMap == k2kG) {
    fSolenoid = 2.;
    parname = dipoleON ? "Sol12_Dip6_Hole":"Sol12_Dip0_Hole";
  } else if (fMap == k5kG) {
    fSolenoid = 5.;
    parname = dipoleON ? "Sol30_Dip6_Hole":"Sol30_Dip0_Hole";
  } else {
    AliError(Form("Unknown field identifier %d is requested\n",fMap)); 
    return;
  }
  //
  fMeasuredMap = dynamic_cast<AliMagWrapCheb*>(file->Get(parname));
  if (!fMeasuredMap) {
    AliError(Form("Did not find field %s in %s\n",parname,fname)); 
    return;
  }
  file->Close();
  delete file;
}


//_______________________________________________________________________
AliMagFCheb::AliMagFCheb(const AliMagFCheb &src):
  AliMagFC(src),
  fMeasuredMap(0),
  fSolenoid(src.fSolenoid)
{
  if (src.fMeasuredMap) fMeasuredMap = new AliMagWrapCheb(*src.fMeasuredMap);
}

//_______________________________________________________________________
AliMagFCheb::~AliMagFCheb()
{
  delete fMeasuredMap;
}

//_______________________________________________________________________
void AliMagFCheb::GetTPCInt(const Float_t *xyz, Float_t *b) const
{
  // Method to calculate the integral of magnetic integral from xyz to nearest cathode plane
  //
  b[0]=b[1]=b[2]=0.0;
  if (fMeasuredMap) fMeasuredMap->GetTPCInt(xyz,b);
  for (int i=3;i--;) b[i] *= fFactor;
}

//_______________________________________________________________________
void AliMagFCheb::GetTPCIntCyl(const Float_t *rphiz, Float_t *b) const
{
  // Method to calculate the integral of magnetic integral from point to nearest cathode plane
  // in cylindrical coordiates ( -pi<phi<pi convention )
  b[0]=b[1]=b[2]=0.0;
  if (fMeasuredMap) fMeasuredMap->GetTPCIntCyl(rphiz,b);
  for (int i=3;i--;) b[i] *= fFactor;
}

//_______________________________________________________________________
void AliMagFCheb::Field(const Float_t *xyz, Float_t *b) const
{
  // Method to calculate the field at point  xyz
  //
  b[0]=b[1]=b[2]=0.0;
  if (xyz[2] > 919. || xyz[2] < -1972.) {
    ZDCField(xyz, b);
  } else {
    if (fMeasuredMap && fFactor !=0.) {
      fMeasuredMap->Field(xyz,b);
      for (int i=3;i--;) b[i] *= fFactor;
    }
  }
}

//_______________________________________________________________________
void AliMagFCheb::Field(const Double_t *xyz, Double_t *b) const
{
  // Method to calculate the field at point  xyz
  //
  b[0]=b[1]=b[2]=0.0;
  if (xyz[2] > 919. || xyz[2] < -1972.) {
    ZDCField(xyz, b);
  } 
  else {
    if (fMeasuredMap && fFactor !=0.) {
      fMeasuredMap->Field(xyz,b);
      for (int i=3;i--;) b[i] *= fFactor;
    }
  }
}

//_______________________________________________________________________
AliMagFCheb& AliMagFCheb::operator=(const AliMagFCheb& maps)
{
  fSolenoid=maps.fSolenoid;
  if (this != &maps && maps.fMeasuredMap) { 
    if (fMeasuredMap) delete fMeasuredMap;
    fMeasuredMap = new AliMagWrapCheb(*maps.fMeasuredMap);
  }
  return *this;
}

//_______________________________________________________________________
void AliMagFCheb::SetMeasuredMap(AliMagWrapCheb* parm) 
{
  if (fMeasuredMap) delete fMeasuredMap; 
  fMeasuredMap = parm;
}
