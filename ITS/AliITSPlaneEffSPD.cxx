/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
///////////////////////////////////////////////////////////////////////////
//  Plane Efficiency class for ITS                      
//  It is used for chip by chip efficiency of the SPD,        
//  evaluated by tracks
//  (Inherits from AliITSPlaneEff)
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TMath.h>
#include "AliITSPlaneEffSPD.h"
#include "AliLog.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
//#include "AliCDBRunRange.h"
#include "AliITSCalibrationSPD.h"

ClassImp(AliITSPlaneEffSPD)	
//______________________________________________________________________
AliITSPlaneEffSPD::AliITSPlaneEffSPD():
  AliITSPlaneEff(){
//  for (UInt_t im=0; im<kNModule; im++){
//  for (UInt_t ic=0; ic<kNChip; ic++){
//    fFound[im][ic]=0;
//    fTried[im][ic]=0;
//  }}
  for (UInt_t i=0; i<kNModule*kNChip; i++){
    fFound[i]=0;
    fTried[i]=0;
  }
  // default constructor
  AliDebug(1,Form("Calling default constructor"));
}
//______________________________________________________________________
AliITSPlaneEffSPD::~AliITSPlaneEffSPD(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.
}
//______________________________________________________________________
AliITSPlaneEffSPD::AliITSPlaneEffSPD(const AliITSPlaneEffSPD &s) : AliITSPlaneEff(s) //,
//fHis(s.fHis),
{
    //     Copy Constructor
    // Inputs:
    //    AliITSPlaneEffSPD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

}
//_________________________________________________________________________
AliITSPlaneEffSPD& AliITSPlaneEffSPD::operator+=(const AliITSPlaneEffSPD &add){
    //    Add-to-me operator
    // Inputs:
    //    const AliITSPlaneEffSPD &add  simulation class to be added
    // Outputs:
    //    none.
    // Return:
    //    none
    for (UInt_t i=0; i<kNModule*kNChip; i++){
      fFound[i] += add.fFound[i];
      fTried[i] += add.fTried[i];
    }
    return *this;
}
//______________________________________________________________________
AliITSPlaneEffSPD&  AliITSPlaneEffSPD::operator=(const
                                           AliITSPlaneEffSPD &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSPD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:
 
    if(this==&s) return *this;
    s.Copy(*this);
//    if(&s == this) return *this;
//    for (UInt_t i=0; i<kNModule*kNChip; i++){
//      this->fFound[i] = s.fFound[i];
//      this->fTried[i] = s.fTried[i];
//    }
    return *this;
}
//______________________________________________________________________
void AliITSPlaneEffSPD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSPlaneEff::Copy(obj);
  //((AliITSPlaneEffSPD& ) obj).fNpx  = fNpx;
  for(Int_t i=0;i<kNModule*kNChip;i++) {
      ((AliITSPlaneEffSPD& ) obj).fFound[i] = fFound[i];
      ((AliITSPlaneEffSPD& ) obj).fTried[i] = fTried[i];
  }
}
//______________________________________________________________________
AliITSPlaneEff&  AliITSPlaneEffSPD::operator=(const
                                           AliITSPlaneEff &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSPD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

    if(&s == this) return *this;
    Error("AliITSPlaneEffSPD","Not allowed to make a = with "
          "AliITSPlaneEffSPD","Using default creater instead");

    return *this;
}
//_______________________________________________________________________
Int_t AliITSPlaneEffSPD::GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr,
          UInt_t im, UInt_t ic) const {
   
  //   Estimate the number of tracks still to be collected to attain a 
  //   given efficiency eff, with relative error RelErr
  //   Inputs:
  //         eff    -> Expected efficiency (e.g. those from actual estimate)
  //         RelErr -> tollerance [0,1] 
  //         im     -> module number [0,249]
  //         ic     -> chip number [0,4]
  //   Outputs: none
  //   Return: the estimated n. of tracks 
  //
if (im>=kNModule || ic>=kNChip) 
 {Error("AliITSPlaneEffSPD","you asked for a non existing chip");
 return -1;}
else return GetNTracksForGivenEff(eff,RelErr)-fTried[GetKey(im,ic)];
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSPD::PlaneEff(const UInt_t im,const UInt_t ic) const {
// Compute the efficiency for a basic block, 
// Inputs:
//        im     -> module number [0,249]
//        ic     -> chip number [0,4] 
if (im>=kNModule || ic>=kNChip) 
 {Error("AliITSPlaneEffSPD","you asked for a non existing chip"); return -1.;}
 Int_t nf=fFound[GetKey(im,ic)];
 Int_t nt=fTried[GetKey(im,ic)];
return AliITSPlaneEff::PlaneEff(nf,nt);
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSPD::ErrPlaneEff(const UInt_t im,const UInt_t ic) const {
    // Compute the statistical error on efficiency for a basic block,
    // using binomial statistics 
    // Inputs:
    //        im     -> module number [0,249]
    //        ic     -> chip number [0,4] 
if (im>=kNModule || ic>=kNChip) 
 {Error("AliITSPlaneEffSPD","you asked for a non existing chip"); return -1.;}
Int_t nf=fFound[GetKey(im,ic)];
Int_t nt=fTried[GetKey(im,ic)];
return AliITSPlaneEff::ErrPlaneEff(nf,nt);
} 
//_________________________________________________________________________
Bool_t AliITSPlaneEffSPD::UpDatePlaneEff(const Bool_t Kfound,
                                         const UInt_t im, const UInt_t ic) {
  // Update efficiency for a basic block
if (im>=kNModule || ic>=kNChip) 
 {Error("AliITSPlaneEffSPD","you asked for a non existing chip"); return kFALSE;}
 fTried[GetKey(im,ic)]++;
 if(Kfound) fFound[GetKey(im,ic)]++;
 return kTRUE;
}
//_________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetChipFromCol(const UInt_t col) const {
  // get chip given the column
if(col>=kNCol*kNChip) 
 {Error("AliITSPlaneEffSPD","you asked for a non existing column"); return 10;}
return col/kNCol;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetKey(const UInt_t mod, const UInt_t chip) const {
  // get key given a basic block
if(mod>=kNModule || chip>=kNChip)
  {Error("AliITSPlaneEffSPD::GetKey","you asked for a non existing block"); return 99999;}
return mod*kNChip+chip;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetModFromKey(const UInt_t key) const {
  // get mod. from key
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetModFromKey","you asked for a non existing key"); return 9999;}
return key/kNChip;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetChipFromKey(const UInt_t key) const {
  // retrieves chip from key
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetChipFromKey","you asked for a non existing key"); return 999;}
return (key%(kNModule*kNChip))%kNChip;
}
//__________________________________________________________________________
void AliITSPlaneEffSPD::GetModAndChipFromKey(const UInt_t key,UInt_t& mod,UInt_t& chip) const {
  // get module and chip from a key
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetModAndChipFromKey","you asked for a non existing key"); 
  mod=9999;
  chip=999;
  return;}
mod=key/kNChip;
chip=(key%(kNModule*kNChip))%kNChip;
return;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSPD::LivePlaneEff(UInt_t key) const {
  // returns plane efficieny after adding the fraction of sensor which is bad
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::LivePlaneEff","you asked for a non existing key");
   return -1.;}
Double_t leff=AliITSPlaneEff::LivePlaneEff(0); // this just for the Warning
leff=PlaneEff(key)+GetFracBad(key);
return leff>1?1:leff;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSPD::ErrLivePlaneEff(UInt_t key) const {
  // returns error on live plane efficiency
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::ErrLivePlaneEff","you asked for a non existing key");
   return -1.;}
Int_t nf=fFound[key];
Double_t triedInLive=GetFracLive(key)*fTried[key];
Int_t nt=TMath::Max(nf,TMath::Nint(triedInLive));
return AliITSPlaneEff::ErrPlaneEff(nf,nt); // for the time being: to be checked
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSPD::GetFracLive(const UInt_t key) const {
  // returns the fraction of the sensor which is OK
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetFracLive","you asked for a non existing key");
   return -1.;}
    // Compute the fraction of bad (dead+noisy) detector 
UInt_t dead=0,noisy=0;
GetDeadAndNoisyInChip(key,dead,noisy);
Double_t live=dead+noisy;
live/=(kNRow*kNCol);
return 1.-live;
}
//_____________________________________________________________________________
void AliITSPlaneEffSPD::GetDeadAndNoisyInChip(const UInt_t key,
      UInt_t& nrDeadInChip, UInt_t& nrNoisyInChip) const {
  // returns the number of dead and noisy pixels
nrDeadInChip=0;
nrNoisyInChip=0;
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip","you asked for a non existing key");
   return;}
    // Compute the number of bad (dead+noisy) pixel in a chip
//
if(!fInitCDBCalled) 
  {Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip","CDB not inizialized: call InitCDB first");
   return;};
AliCDBManager* man = AliCDBManager::Instance();
// retrieve map of dead Pixel 
AliCDBEntry *cdbSPDDead = man->Get("ITS/Calib/SPDDead", fRunNumber);
TObjArray* spdDead;
if(cdbSPDDead) {
  spdDead = (TObjArray*)cdbSPDDead->GetObject();
  if(!spdDead) 
  {Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip"," SPDDead not found in CDB");
   return;}
} else {
  Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip","Did not find Calib/SPDDead.");
  return;
}
// retrieve map of noisy Pixel 
AliCDBEntry *cdbSPDNoisy = man->Get("ITS/Calib/SPDNoisy", fRunNumber);
TObjArray* spdNoisy;
if(cdbSPDNoisy) {
  spdNoisy = (TObjArray*)cdbSPDNoisy->GetObject();
  if(!spdNoisy) 
  {Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip"," SPDNoisy not found in CDB");
   return;}
} else {
  Error("AliITSPlaneEffSPD::GetDeadAndNoisyInChip","Did not find Calib/SPDNoisy.");
  return;
}
//
UInt_t mod=GetModFromKey(key);
UInt_t chip=GetChipFromKey(key);
// count number of dead
AliITSCalibrationSPD* calibSPD=(AliITSCalibrationSPD*) spdDead->At(mod);
UInt_t nrDead = calibSPD->GetNrBad();
for (UInt_t index=0; index<nrDead; index++) {
  if(GetChipFromCol(calibSPD->GetBadColAt(index))==chip) nrDeadInChip++;
}
calibSPD=(AliITSCalibrationSPD*) spdNoisy->At(mod);
UInt_t nrNoisy = calibSPD->GetNrBad();
for (UInt_t index=0; index<nrNoisy; index++) {
  if(GetChipFromCol(calibSPD->GetBadColAt(index))==chip) nrNoisyInChip++;
}
return;
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSPD::GetFracBad(const UInt_t key) const {
  // returns 1-fractional live
if(key>=kNModule*kNChip)
  {Error("AliITSPlaneEffSPD::GetFracBad","you asked for a non existing key");
   return -1.;}
return 1.-GetFracLive(key);
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSPD::WriteIntoCDB() const {
// write onto CDB
if(!fInitCDBCalled)
  {Error("AliITSPlaneEffSPD::WriteIntoCDB","CDB not inizialized: call InitCDB first");
   return kFALSE;}
// to be written properly: now only for debugging 
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetObjectClassName("AliITSPlaneEff");
  md->SetResponsible("Giuseppe Eugenio Bruno");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion("head 19/11/07"); //root version
  AliCDBId id("ITS/PlaneEff/PlaneEffSPD",0,AliCDBRunRange::Infinity()); 
  AliITSPlaneEffSPD eff; 
  eff=*this;
  Bool_t r=AliCDBManager::Instance()->GetDefaultStorage()->Put(&eff,id,md);
  delete md;
  return r;
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSPD::ReadFromCDB() {
// read from CDB
if(!fInitCDBCalled)
  {Error("AliITSPlaneEffSPD::ReadFromCDB","CDB not inizialized: call InitCDB first");
   return kFALSE;}
//if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
//    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
//  }
AliCDBEntry *cdbEntry = AliCDBManager::Instance()->Get("ITS/PlaneEff/PlaneEffSPD",fRunNumber);
AliITSPlaneEffSPD* eff= (AliITSPlaneEffSPD*)cdbEntry->GetObject();
if(this==eff) return kFALSE;
eff->Copy(*this);
return kTRUE;
}
