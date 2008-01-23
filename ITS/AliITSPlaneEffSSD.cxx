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
//  It is used for module by module efficiency of the SSD,        
//  evaluated by tracks
//  (Inherits from AliITSPlaneEff)
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////

/*  $Id$ */

#include <TMath.h>
#include "AliITSPlaneEffSSD.h"
#include "AliLog.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
//#include "AliCDBRunRange.h"
#include "AliITSCalibrationSSD.h"

ClassImp(AliITSPlaneEffSSD)	
//______________________________________________________________________
AliITSPlaneEffSSD::AliITSPlaneEffSSD():
  AliITSPlaneEff(){
  for (UInt_t i=0; i<kNModule; i++){
    fFound[i]=0;
    fTried[i]=0;
  }
  // default constructor
  AliDebug(1,Form("Calling default constructor"));
}
//______________________________________________________________________
AliITSPlaneEffSSD::~AliITSPlaneEffSSD(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.
}
//______________________________________________________________________
AliITSPlaneEffSSD::AliITSPlaneEffSSD(const AliITSPlaneEffSSD &s) : AliITSPlaneEff(s) //,
//fHis(s.fHis),
{
    //     Copy Constructor
    // Inputs:
    //    AliITSPlaneEffSSD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:
}
//_________________________________________________________________________
AliITSPlaneEffSSD& AliITSPlaneEffSSD::operator+=(const AliITSPlaneEffSSD &add){
    //    Add-to-me operator
    // Inputs:
    //    const AliITSPlaneEffSSD &add  simulation class to be added
    // Outputs:
    //    none.
    // Return:
    //    none
    for (UInt_t i=0; i<kNModule; i++){
      fFound[i] += add.fFound[i];
      fTried[i] += add.fTried[i];
    }
    return *this;
}
//______________________________________________________________________
AliITSPlaneEffSSD&  AliITSPlaneEffSSD::operator=(const
                                           AliITSPlaneEffSSD &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSSD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:
 
    if(this==&s) return *this;
    s.Copy(*this);
    return *this;
}
//______________________________________________________________________
void AliITSPlaneEffSSD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSPlaneEff::Copy(obj);
  for(Int_t i=0;i<kNModule;i++) {
      ((AliITSPlaneEffSSD& ) obj).fFound[i] = fFound[i];
      ((AliITSPlaneEffSSD& ) obj).fTried[i] = fTried[i];
  }
}
//______________________________________________________________________
AliITSPlaneEff&  AliITSPlaneEffSSD::operator=(const
                                           AliITSPlaneEff &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSSD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

    if(&s == this) return *this;
    AliError("operator=: Not allowed to make a =, use default creater instead");
    return *this;
}
//_______________________________________________________________________
Int_t AliITSPlaneEffSSD::GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr,
          UInt_t im) const {
   
  //   Estimate the number of tracks still to be collected to attain a 
  //   given efficiency eff, with relative error RelErr
  //   Inputs:
  //         eff    -> Expected efficiency (e.g. those from actual estimate)
  //         RelErr -> tollerance [0,1] 
  //         im     -> module number [0,1697]
  //   Outputs: none
  //   Return: the estimated n. of tracks 
  //
if (im>=kNModule) 
 {AliError("GetMissingTracksForGivenEff: you asked for a non existing module");
 return -1;}
else return GetNTracksForGivenEff(eff,RelErr)-fTried[GetKey(im)];
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSSD::PlaneEff(const UInt_t im) const {
// Compute the efficiency for a basic block, 
// Inputs:
//        im     -> module number [0,1697]
if (im>=kNModule) 
 {AliError("PlaneEff(UInt_t): you asked for a non existing module"); return -1.;}
 Int_t nf=fFound[GetKey(im)];
 Int_t nt=fTried[GetKey(im)];
return AliITSPlaneEff::PlaneEff(nf,nt);
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSSD::ErrPlaneEff(const UInt_t im) const {
    // Compute the statistical error on efficiency for a basic block,
    // using binomial statistics 
    // Inputs:
    //        im     -> module number [0,1697]
if (im>=kNModule) 
 {AliError("ErrPlaneEff(UInt_t): you asked for a non existing module"); return -1.;}
Int_t nf=fFound[GetKey(im)];
Int_t nt=fTried[GetKey(im)];
return AliITSPlaneEff::ErrPlaneEff(nf,nt);
} 
//_________________________________________________________________________
Bool_t AliITSPlaneEffSSD::UpDatePlaneEff(const Bool_t Kfound, const UInt_t im) {
  // Update efficiency for a basic block
if (im>=kNModule) 
 {AliError("UpDatePlaneEff: you asked for a non existing module"); return kFALSE;}
 fTried[GetKey(im)]++;
 if(Kfound) fFound[GetKey(im)]++;
 return kTRUE;
}
//_________________________________________________________________________
UInt_t AliITSPlaneEffSSD::GetKey(const UInt_t mod) const {
  // get key given a basic block
if(mod>=kNModule)
  {AliError("GetKey: you asked for a non existing block"); return 99999;}
return mod;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSSD::GetModFromKey(const UInt_t key) const {
  // get mod. from key
if(key>=kNModule)
  {AliError("GetModFromKey: you asked for a non existing key"); return 9999;}
return key;
}
//__________________________________________________________________________
Double_t AliITSPlaneEffSSD::LivePlaneEff(UInt_t key) const {
  // returns plane efficieny after adding the fraction of sensor which is bad
if(key>=kNModule)
  {AliError("LivePlaneEff: you asked for a non existing key");
   return -1.;}
Double_t leff=AliITSPlaneEff::LivePlaneEff(0); // this just for the Warning
leff=PlaneEff(key)+GetFracBad(key);
return leff>1?1:leff;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSSD::ErrLivePlaneEff(UInt_t key) const {
  // returns error on live plane efficiency
if(key>=kNModule)
  {AliError("ErrLivePlaneEff: you asked for a non existing key");
   return -1.;}
Int_t nf=fFound[key];
Double_t triedInLive=GetFracLive(key)*fTried[key];
Int_t nt=TMath::Max(nf,TMath::Nint(triedInLive));
return AliITSPlaneEff::ErrPlaneEff(nf,nt); // for the time being: to be checked
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSSD::GetFracLive(const UInt_t key) const {
  // returns the fraction of the sensor which is OK
if(key>=kNModule)
  {AliError("GetFracLive: you asked for a non existing key");
   return -1.;}
    // Compute the fraction of bad (dead+noisy) detector 
UInt_t bad=0;
GetBadInModule(key,bad);
Double_t live=bad;
live/=(kNChip*kNSide*kNStrip);
return 1.-live;
}
//_____________________________________________________________________________
void AliITSPlaneEffSSD::GetBadInModule(const UInt_t key, UInt_t& nrBadInMod) const {
  // returns the number of dead and noisy pixels
nrBadInMod=0;
if(key>=kNModule)
  {AliError("GetBadInModule: you asked for a non existing key");
   return;}
    // Compute the number of bad (dead+noisy) pixel in a module
//
if(!fInitCDBCalled) 
  {AliError("GetBadInModule: CDB not inizialized: call InitCDB first");
   return;};
AliCDBManager* man = AliCDBManager::Instance();
// retrieve map of dead Pixel 
AliCDBEntry *cdbSSD = man->Get("ITS/Calib/BadChannelsSSD", fRunNumber);
TObjArray* ssdEntry;
if(cdbSSD) {
  ssdEntry = (TObjArray*)cdbSSD->GetObject();
  if(!ssdEntry) 
  {AliError("GetBadInChip: SSDEntry not found in CDB");
   return;}
} else {
  AliError("GetBadInChip: did not find Calib/BadChannelsSSD");
  return;
}
//
UInt_t mod=GetModFromKey(key);
//
AliITSBadChannelsSSD* badchannels=(AliITSBadChannelsSSD*) ssdEntry->At(mod);
// count the  number of bad channels on the p side
nrBadInMod += (badchannels->GetBadPChannelsList()).GetSize();
// add the  number of bad channels on the s side
nrBadInMod += (badchannels->GetBadNChannelsList()).GetSize();
return;
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSSD::GetFracBad(const UInt_t key) const {
  // returns 1-fractional live
if(key>=kNModule)
  {AliError("GetFracBad: you asked for a non existing key");
   return -1.;}
return 1.-GetFracLive(key);
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSSD::WriteIntoCDB() const {
// write onto CDB
if(!fInitCDBCalled)
  {AliError("WriteIntoCDB: CDB not inizialized. Call InitCDB first");
   return kFALSE;}
// to be written properly: now only for debugging 
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetObjectClassName("AliITSPlaneEff");
  md->SetResponsible("Giuseppe Eugenio Bruno");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion("head 02/01/08"); //root version
  AliCDBId id("ITS/PlaneEff/PlaneEffSSD",0,AliCDBRunRange::Infinity()); 
  AliITSPlaneEffSSD eff; 
  eff=*this;
  Bool_t r=AliCDBManager::Instance()->GetDefaultStorage()->Put(&eff,id,md);
  delete md;
  return r;
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSSD::ReadFromCDB() {
// read from CDB
if(!fInitCDBCalled)
  {AliError("ReadFromCDB: CDB not inizialized. Call InitCDB first");
   return kFALSE;}
//if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
//    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
//  }
AliCDBEntry *cdbEntry = AliCDBManager::Instance()->Get("ITS/PlaneEff/PlaneEffSSD",fRunNumber);
AliITSPlaneEffSSD* eff= (AliITSPlaneEffSSD*)cdbEntry->GetObject();
if(this==eff) return kFALSE;
eff->Copy(*this);
return kTRUE;
}
//_____________________________________________________________________________
UInt_t AliITSPlaneEffSSD::GetKeyFromDetLocCoord(Int_t ilay, Int_t idet,
                                                Float_t, Float_t) const {
// method to locate a basic block from Detector Local coordinate (to be used in tracking)
UInt_t key=999999;
if(ilay<4 || ilay>5)
  {AliError("GetKeyFromDetLocCoord: you asked for a non existing layer");
   return key;}
if(ilay==4 && (idet<0 || idet>747))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}
if(ilay==5 && (idet<0 || idet>949))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}

UInt_t mod=idet;
if(ilay==1) mod+=748;
key=GetKey(mod);
return key;
}
