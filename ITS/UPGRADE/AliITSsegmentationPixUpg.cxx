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

/* $Id: AliITSsegmentationPixUpg.cxx 47180 2011-02-08 09:42:29Z masera $ */
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TObjArray.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include "AliITSsegmentationPixUpg.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Segmentation class for pixels                                                                          //
// Questions to solve: are guardrings needed and do they belong to the sensor or to the module in TGeo    //
//                     At the moment assume that the local coord syst. is located at bottom left corner   //
//                     of the ACTIVE matrix. If the guardring to be accounted in the local coords, in     //
//                     the Z and X conversions one needs to first subtract the  fGuardLft and fGuardBot   //
//                     from the local Z,X coordinates                                                     //
//                                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(AliITSsegmentationPixUpg)

const char* AliITSsegmentationPixUpg::fgkSegmListName = "ITSUpgradeSegmentations";

//_____________________________________________________________________________RS
AliITSsegmentationPixUpg::AliITSsegmentationPixUpg(int nchips,int ncol,int nrow,
						   double pitchX,double pitchZ,
						   double thickness,
						   double pitchLftC,double pitchRgtC,
						   double edgL,double edgR,double edgT,double edgB)
: AliITSsegmentation()
  ,fGuardLft(edgL)
  ,fGuardRgt(edgR)
  ,fGuardTop(edgT)
  ,fGuardBot(edgB)
  ,fPitchX(pitchX)
  ,fPitchZ(pitchZ)
  ,fPitchZLftCol(pitchLftC<0 ? pitchZ:pitchLftC)
  ,fPitchZRgtCol(pitchRgtC<0 ? pitchZ:pitchRgtC)
  ,fChipDZ(0)
  ,fNChips(nchips)
  ,fNColPerChip(nchips>0 ? ncol/nchips:0)
  ,fNRow(nrow)
  ,fNCol(ncol)
{
  // Default constructor, sizes in microns
  fChipDZ = (fNColPerChip-2)*fPitchZ + fPitchZLftCol + fPitchZRgtCol;
  SetDetSize( fNRow*fPitchX /*+fGuardTop+fGuardBot*/,
	      fNChips*fChipDZ /*+fGuardLft+fGuardRgt*/,
	      thickness);
  //
}

//_____________________________________________________________________________RS
void AliITSsegmentationPixUpg::GetPadIxz(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const 
{
  //  Returns pixel coordinates (ix,iz) for given real local coordinates (x,z)
  //  expects x, z in microns
  ix = int(x/fPitchX) + 1;     
  iz = int(Z2Col(z)) + 1;
  //  
  if (iz >  fNCol) iz= fNCol;
  if (ix >  fNRow) ix= fNRow;
  //
}

//_____________________________________________________________________________RS
void AliITSsegmentationPixUpg::GetPadTxz(Float_t &x,Float_t &z) const
{
  //  local transformation of real local coordinates (x,z)
  //  expects x, z in microns
  x /= fPitchX;
  z = Z2Col(z);
  //
}

//_____________________________________________________________________________RS
void AliITSsegmentationPixUpg::GetPadCxz(Int_t ix,Int_t iz,Float_t &x,Float_t&z) const
{
  // Transform from pixel to real local coordinates
  // returns x, z in microns
  x = (ix>0) ? Float_t((ix-0.5)*fPitchX) : Float_t((ix+0.5)*fPitchX);
  z = Col2Z(iz);
  //
}

//_____________________________________________________________________________RS
Float_t AliITSsegmentationPixUpg::Z2Col(Float_t z) const 
{
  // get column number (from 0) from local Z
  int chip = z/fChipDZ;
  float col = chip*fNColPerChip;
  z -= chip*fChipDZ;
  if (z<fPitchZLftCol) col += z/fPitchZLftCol;
  else {
    z = fPitchZLftCol;
    col += 1;
    if (z<(fChipDZ-fPitchZRgtCol)) col += 1+z/fPitchZ;
    else col += 1+(z - (fNColPerChip-2)*fPitchZ)/fPitchZRgtCol;
  }
  return col;
}

//_____________________________________________________________________________RS
Float_t AliITSsegmentationPixUpg::Col2Z(Int_t col) const 
{
  // convert column number (from 0) to Z coordinate
  int nchip = col/fNColPerChip;
  col %= fNColPerChip;
  float z = nchip*fChipDZ;
  if (!col) z -= fPitchZRgtCol/2;
  else if (col==1) z += fPitchZLftCol/2;
  else if (col==fNColPerChip-1) z += fChipDZ - fPitchZRgtCol/2;
  else    z += fPitchZLftCol + (col-0.5)*fChipDZ;
  return z;
  //
}

//______________________________________________________________________RS
AliITSsegmentationPixUpg& AliITSsegmentationPixUpg::operator=(const AliITSsegmentationPixUpg &src)
{
  // = operator
  if(this==&src) return *this;
  AliITSsegmentation::operator=(src);
  fNCol  = src.fNCol;
  fNRow  = src.fNRow;
  fNColPerChip  = src.fNColPerChip;
  fNChips = src.fNChips;
  fChipDZ = src.fChipDZ;
  fPitchZRgtCol = src.fPitchZRgtCol;
  fPitchZLftCol = src.fPitchZLftCol;
  fPitchZ = src.fPitchZ;
  fPitchX = src.fPitchX;
  //
  fGuardBot = src.fGuardBot;
  fGuardTop = src.fGuardTop;
  fGuardRgt = src.fGuardRgt;
  fGuardLft = src.fGuardLft;
  //
  return *this;
}

//____________________________________________________________________________RS
AliITSsegmentationPixUpg::AliITSsegmentationPixUpg(const AliITSsegmentationPixUpg &src) :
  AliITSsegmentation(src)
  ,fGuardLft(src.fGuardLft)
  ,fGuardRgt(src.fGuardRgt)
  ,fGuardTop(src.fGuardTop)
  ,fGuardBot(src.fGuardBot)
  ,fPitchX(src.fPitchX)
  ,fPitchZ(src.fPitchZ)
  ,fPitchZLftCol(src.fPitchZLftCol)
  ,fPitchZRgtCol(src.fPitchZRgtCol)
  ,fChipDZ(src.fChipDZ)
  ,fNChips(src.fNChips)
  ,fNColPerChip(src.fNColPerChip)
  ,fNRow(src.fNRow)
  ,fNCol(src.fNCol)  
{
}

//____________________________________________________________________________RS
Float_t AliITSsegmentationPixUpg::Dpx(Int_t ) const 
{
  //returs x pixel pitch for a give pixel
  return fPitchX;
}

//____________________________________________________________________________RS
Float_t AliITSsegmentationPixUpg::Dpz(Int_t col) const 
{
  // returns z pixel pitch for a given pixel (cols starts from 0)
  col %= fNColPerChip;
  if (!col) return fPitchZLftCol;
  if (col==fNColPerChip-1) return fPitchZRgtCol;
  return fPitchZ;
  //
}

//------------------------------
void AliITSsegmentationPixUpg::Neighbours(Int_t iX, Int_t iZ, Int_t* nlist, Int_t xlist[8], Int_t zlist[8]) const 
{
  // returns the neighbouring pixels for use in Cluster Finders and the like.
  //
  *nlist=8;
  xlist[0]=xlist[1]=iX;
  xlist[2]=iX-1;
  xlist[3]=iX+1;
  zlist[0]=iZ-1;
  zlist[1]=iZ+1;
  zlist[2]=zlist[3]=iZ;
  // Diagonal elements
  xlist[4]=iX+1;
  zlist[4]=iZ+1;
  //  
  xlist[5]=iX-1;
  zlist[5]=iZ-1;
  //
  xlist[6]=iX-1;
  zlist[6]=iZ+1;
  //
  xlist[7]=iX+1;
  zlist[7]=iZ-1;
  //
}

//______________________________________________________________________
Bool_t AliITSsegmentationPixUpg::LocalToDet(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const 
{
  // Transformation from Geant detector centered local coordinates (cm) to
  // Pixel cell numbers ix and iz.
  // Input:
  //    Float_t   x        detector local coordinate x in cm with respect to
  //                       the center of the sensitive volume.
  //    Float_t   z        detector local coordinate z in cm with respect to
  //                       the center of the sensitive volulme.
  // Output:
  //    Int_t    ix        detector x cell coordinate. Has the range 
  //                       0<=ix<fNRow.
  //    Int_t    iz        detector z cell coordinate. Has the range 
  //                       0<=iz<fNCol.
  // Return:
  //   kTRUE if point x,z is inside sensitive volume, kFALSE otherwise.
  //   A value of -1 for ix or iz indecates that this point is outside of the
  //   detector segmentation as defined.
  x = x*kCM2MC+0.5*Dx();
  z = z*kCM2MC+0.5*Dz();
  ix = iz = -1;
  if(x<0 || x>Dx()) return kFALSE; // outside x range.
  if(z<0 || z>Dz()) return kFALSE; // outside z range.
  ix = int(x/fPitchX);
  iz = Z2Col(z);
  return kTRUE; // Found ix and iz, return.
}

//______________________________________________________________________
void AliITSsegmentationPixUpg::DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const
{
// Transformation from Detector cell coordiantes to Geant detector centerd 
// local coordinates (cm).
// Input:
// Int_t    ix        detector x cell coordinate. Has the range 0<=ix<fNRow.
// Int_t    iz        detector z cell coordinate. Has the range 0<=iz<fNCol.
// Output:
// Float_t   x        detector local coordinate x in cm with respect to the
//                    center of the sensitive volume.
// Float_t   z        detector local coordinate z in cm with respect to the
//                    center of the sensitive volulme.
// If ix and or iz is outside of the segmentation range a value of -0.5*Dx()
// or -0.5*Dz() is returned.
  //
  x = -0.5/kCM2MC*Dx(); // default value.
  z = -0.5/kCM2MC*Dz(); // default value.
  // RS: to check: why we don't do strict check for [0:n)
  if(ix<0 || ix>=fNRow) return; // outside of detector 
  if(iz<0 || iz>=fNCol) return; // outside of detctor
  x += (ix+0.5)*fPitchX*(1./kCM2MC);       // RS: we go to the center of the pad, i.e. + pitch/2, not to the boundary as in SPD
  z += Col2Z(iz)*(1./kCM2MC); 
  return; // Found x and z, return.
}

//______________________________________________________________________
void AliITSsegmentationPixUpg::CellBoundries(Int_t ix,Int_t iz,Double_t &xl,Double_t &xu,Double_t &zl,Double_t &zu) const
{
  // Transformation from Detector cell coordiantes to Geant detector centerd 
  // local coordinates (cm).
  // Input:
  // Int_t    ix        detector x cell coordinate. Has the range 0<=ix<fNRow.
  // Int_t    iz        detector z cell coordinate. Has the range 0<=iz<fNCol.
  // Output:
  // Double_t   xl       detector local coordinate cell lower bounds x in cm
  //                    with respect to the center of the sensitive volume.
  // Double_t   xu       detector local coordinate cell upper bounds x in cm 
  //                    with respect to the center of the sensitive volume.
  // Double_t   zl       detector local coordinate lower bounds z in cm with
  //                    respect to the center of the sensitive volulme.
  // Double_t   zu       detector local coordinate upper bounds z in cm with 
  //                    respect to the center of the sensitive volulme.
  // If ix and or iz is outside of the segmentation range a value of -0.5*Dx()
  // and -0.5*Dx() or -0.5*Dz() and -0.5*Dz() are returned.
  Float_t x,z;
  DetToLocal(ix,iz,x,z);
  //
  if( ix<0 || ix>=fNRow || iz<0 || iz>=fNCol) {
    xl = xu = -0.5/kCM2MC*Dx(); // default value.
    zl = zu = -0.5/kCM2MC*Dz(); // default value.
    return; // outside of detctor
  }
  float zpitchH = Dpz(iz)/2./kCM2MC;
  float xpitchH = fPitchX/2./kCM2MC;
  xl -= xpitchH;
  xu += xpitchH;
  zl -= zpitchH;
  zu += zpitchH;
  return; // Found x and z, return.
}

//______________________________________________________________________
Int_t AliITSsegmentationPixUpg::GetChipFromChannel(Int_t, Int_t iz) const 
{
  // returns chip number (in range 0-4) starting from channel number
  if(iz>=fNCol  || iz<0 ){
    AliWarning("Bad cell number");
    return -1;
  }
  return iz/fNColPerChip;
}

//______________________________________________________________________
Int_t AliITSsegmentationPixUpg::GetChipFromLocal(Float_t, Float_t zloc) const 
{
  // returns chip number (in range 0-4) starting from local coordinates
  Int_t ix0,iz;
  if (!LocalToDet(0,zloc,ix0,iz)) {
    AliWarning("Bad local coordinate");
    return -1;
  } 
  return GetChipFromChannel(ix0,iz);
}

//______________________________________________________________________
Int_t AliITSsegmentationPixUpg::GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, Float_t, Float_t) const 
{
  // returns the number of chips containing a road defined by given local coordinate limits
  //
  const Float_t kconv = 1./kCM2MC; // converts microns to cm.
  //
  if (zmin>zmax) {
    AliWarning("Bad coordinate limits: zmin>zmax!");
    return -1;
  } 
  //
  Int_t nChipInW = 0;
  //
  Float_t zminDet = -0.5*kconv*Dz();
  Float_t zmaxDet =  0.5*kconv*Dz();
  if(zmin<zminDet) zmin=zminDet;
  if(zmax>zmaxDet) zmax=zmaxDet;

  Int_t n1 = GetChipFromLocal(0,zmin);
  array[nChipInW] = n1;
  nChipInW++;

  Int_t n2 = GetChipFromLocal(0,zmax);

  if(n2!=n1){
    Int_t imin=TMath::Min(n1,n2);
    Int_t imax=TMath::Max(n1,n2);
    for(Int_t ichip=imin; ichip<=imax; ichip++){
      if(ichip==n1) continue;
      array[nChipInW]=ichip;
      nChipInW++;
    }
  }
  //
  return nChipInW;
}

//______________________________________________________________________
void AliITSsegmentationPixUpg::Init()
{
  // init settings
}

//______________________________________________________________________
Bool_t AliITSsegmentationPixUpg::StoreWithID(UInt_t id, const char* outf)
{
  // store in the special list under given ID
  TString fns = outf;
  gSystem->ExpandPathName(fns);
  if (fns.IsNull()) {AliFatal("No file name provided"); return kFALSE;}
  TFile* fout = TFile::Open(fns.Data(),"update");
  if (!fout) {AliFatal(Form("Failed to open output file %s",outf)); return kFALSE;}
  TObjArray* arr = (TObjArray*)fout->Get(fgkSegmListName);
  if (!arr) arr = new TObjArray();
  else {
    int nent = arr->GetEntriesFast();
    for (int i=nent;i--;) {
      AliITSsegmentationPixUpg* segm = dynamic_cast<AliITSsegmentationPixUpg*>(arr->At(i));
      if (segm && segm->GetUniqueID()==id) {
	AliFatal(Form("Segmenation %d already exists in file %s",id,outf)); 
	return kFALSE;
      }
    }
  }
  //
  this->SetUniqueID(id);  
  arr->AddLast(this);
  arr->SetOwner(kTRUE);
  fout->WriteObject(arr,fgkSegmListName,"kSingleKey");
  fout->Close();
  delete fout;
  arr->Remove(this);
  delete arr;
  AliInfo(Form("Stored segmentation %d in %s",id,outf));
  return kTRUE;
  //
}

//______________________________________________________________________
AliITSsegmentationPixUpg* AliITSsegmentationPixUpg::LoadWithID(UInt_t id, const char* inpf)
{
  // store in the special list under given ID
  TString fns = inpf;
  gSystem->ExpandPathName(fns);
  if (fns.IsNull()) {AliFatalGeneral("LoadWithID","No file name provided"); return 0;}
  TFile* finp = TFile::Open(fns.Data());
  if (!finp) {AliFatalGeneral("LoadWithID",Form("Failed to open file %s",inpf)); return 0;}
  TObjArray* arr = (TObjArray*)finp->Get(fgkSegmListName);
  if (!arr) {
    AliFatalGeneral("LoadWithID",Form("Failed to find segmenation array %s in %s",fgkSegmListName,inpf)); 
    return 0;
  }
  AliITSsegmentationPixUpg* segm = 0;
  int nent = arr->GetEntriesFast();
  for (int i=nent;i--;) {
    segm = dynamic_cast<AliITSsegmentationPixUpg*>(arr->At(i));
    if (segm && segm->GetUniqueID()==id) {arr->RemoveAt(i); break;}
    segm = 0;
  }
  //
  if (!segm) {AliFatalGeneral("LoadWithID",Form("Failed to find segmenation %d in %s",id,inpf)); return 0;}
  //
  arr->SetOwner(kTRUE); // to not leave in memory other segmenations
  finp->Close();
  delete finp;
  delete arr;
  //
  return segm;
}
