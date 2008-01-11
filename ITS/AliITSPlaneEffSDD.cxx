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
//  It is used for chip by chip efficiency (eventually with subwing division 
//  along the drift direction ) of the SDD,        
//  evaluated by tracks
//  (Inherits from AliITSPlaneEff)
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TMath.h>
#include "AliITSPlaneEffSDD.h"
#include "AliLog.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
//#include "AliCDBRunRange.h"
#include "AliITSCalibrationSDD.h"

ClassImp(AliITSPlaneEffSDD)	
//______________________________________________________________________
AliITSPlaneEffSDD::AliITSPlaneEffSDD():
  AliITSPlaneEff(){
  // Default constructor
  for (UInt_t i=0; i<kNModule*kNChip*kNWing*kNSubWing; i++){
    fFound[i]=0;
    fTried[i]=0;
  }
  // default constructor
  AliDebug(1,Form("Calling default constructor"));
}
//______________________________________________________________________
AliITSPlaneEffSDD::~AliITSPlaneEffSDD(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.
}
//______________________________________________________________________
AliITSPlaneEffSDD::AliITSPlaneEffSDD(const AliITSPlaneEffSDD &s) : AliITSPlaneEff(s) //,
//fHis(s.fHis),
{
    //     Copy Constructor
    // Inputs:
    //    AliITSPlaneEffSDD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

}
//_________________________________________________________________________
AliITSPlaneEffSDD& AliITSPlaneEffSDD::operator+=(const AliITSPlaneEffSDD &add){
    //    Add-to-me operator
    // Inputs:
    //    const AliITSPlaneEffSDD &add  simulation class to be added
    // Outputs:
    //    none.
    // Return:
    //    none
    for (UInt_t i=0; i<kNModule*kNChip*kNWing*kNSubWing; i++){
      fFound[i] += add.fFound[i];
      fTried[i] += add.fTried[i];
    }
    return *this;
}
//______________________________________________________________________
AliITSPlaneEffSDD&  AliITSPlaneEffSDD::operator=(const
                                           AliITSPlaneEffSDD &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSDD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:
 
    if(this==&s) return *this;
    s.Copy(*this);
//    if(&s == this) return *this;
//    for (UInt_t i=0; i<kNModule*kNChip*kNWing*kNSubWing; i++){
//      this->fFound[i] = s.fFound[i];
//      this->fTried[i] = s.fTried[i];
//    }
    return *this;
}
//______________________________________________________________________
void AliITSPlaneEffSDD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSPlaneEff::Copy(obj);
  for(Int_t i=0;i<kNModule*kNChip*kNWing*kNSubWing;i++) {
      ((AliITSPlaneEffSDD& ) obj).fFound[i] = fFound[i];
      ((AliITSPlaneEffSDD& ) obj).fTried[i] = fTried[i];
  }
}
//______________________________________________________________________
AliITSPlaneEff&  AliITSPlaneEffSDD::operator=(const
                                           AliITSPlaneEff &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSDD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

    if(&s == this) return *this;
    Error("AliITSPlaneEffSDD","Not allowed to make a = with "
          "AliITSPlaneEffSDD","Using default creater instead");

    return *this;
}
//_______________________________________________________________________
Int_t AliITSPlaneEffSDD::GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr,
          UInt_t im, UInt_t ic, UInt_t iw, UInt_t is) const {
   
  //   Estimate the number of tracks still to be collected to attain a 
  //   given efficiency eff, with relative error RelErr
  //   Inputs:
  //         eff    -> Expected efficiency (e.g. those from actual estimate)
  //         RelErr -> tollerance [0,1] 
  //         im     -> module number [0,259]
  //         ic     -> chip number [0,3]
  //         iw     -> wing number [0,1]
  //         is     -> chip number [0,kNSubWing-1]
  //   Outputs: none
  //   Return: the estimated n. of tracks 
  //
if (im>=kNModule || ic>=kNChip || iw>=kNWing || is>=kNSubWing) 
 {Error("AliITSPlaneEffSDD","you asked for a non existing block");
 return -1;}
else return GetNTracksForGivenEff(eff,RelErr)-fTried[GetKey(im,ic,iw,is)];
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSDD::PlaneEff(const UInt_t im,const UInt_t ic,
                                      const UInt_t iw,const UInt_t is) const {
// Compute the efficiency for a basic block, 
// Inputs:
//        im     -> module number [0,259]
//        ic     -> chip number [0,3]
//        iw     -> wing number [0,1]
//        is     -> chip number [0,kNSubWing-1]
if (im>=kNModule || ic>=kNChip || iw>=kNWing || is>=kNSubWing) 
 {Error("AliITSPlaneEffSDD","you asked for a non existing block"); return -1.;}
 Int_t nf=fFound[GetKey(im,ic,iw,is)];
 Int_t nt=fTried[GetKey(im,ic,iw,is)];
return AliITSPlaneEff::PlaneEff(nf,nt);
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSDD::ErrPlaneEff(const UInt_t im,const UInt_t ic,
                                         const UInt_t iw,const UInt_t is) const {
    // Compute the statistical error on efficiency for a basic block,
    // using binomial statistics 
    // Inputs:
    //        im     -> module number [0,259]
    //        ic     -> chip number [0,3]
    //        iw     -> wing number [0,1]
    //        is     -> chip number [0,kNSubWing-1]
if (im>=kNModule || ic>=kNChip || iw>=kNWing || is>=kNSubWing) 
 {Error("AliITSPlaneEffSDD","you asked for a non existing block"); return -1.;}
Int_t nf=fFound[GetKey(im,ic,iw,is)];
Int_t nt=fTried[GetKey(im,ic,iw,is)];
return AliITSPlaneEff::ErrPlaneEff(nf,nt);
} 
//_________________________________________________________________________
Bool_t AliITSPlaneEffSDD::UpDatePlaneEff(const Bool_t Kfound,
                                         const UInt_t im, const UInt_t ic,
                                         const UInt_t iw, const UInt_t is) {
  // Update efficiency for a basic block
if (im>=kNModule || ic>=kNChip || iw>=kNWing || is>=kNSubWing) 
 {Error("AliITSPlaneEffSDD","you asked for a non existing block"); return kFALSE;}
 fTried[GetKey(im,ic,iw,is)]++;
 if(Kfound) fFound[GetKey(im,ic,iw,is)]++;
 return kTRUE;
}
//_________________________________________________________________________
void AliITSPlaneEffSDD::ChipAndWingFromAnode(const UInt_t anode, UInt_t& chip,
                                      UInt_t& wing) const {
  // Retun the chip number [0,3] and the wing number [0,1] given the anode number
  // input: anode number [0,511]
if(anode>=kNAnode*kNChip*kNWing)
 {Error("AliITSPlaneEffSDD::ChipAndWingFromAnode","you asked for a non existing anode"); 
  chip=999;
  wing=99;
  return;}
wing=0;
chip=anode/kNAnode;
if(anode>=kNChip*kNAnode) wing=1;
if(wing==1) chip-=kNChip;
return;
}
//_________________________________________________________________________
UInt_t AliITSPlaneEffSDD::ChipFromAnode(const UInt_t anode) const {
  // Retun the chip number [0,3] given the anode number 
  // input: anode number [0,511]
if(anode>=kNAnode*kNChip*kNWing) 
 {Error("AliITSPlaneEffSDD::ChipFromAnode","you asked for a non existing anode"); return 999;}
Int_t wing=0;
Int_t chip=anode/kNAnode;
if(anode>=kNChip*kNAnode) wing=1;
if(wing==1)chip-=kNChip;
return chip;
}
//_________________________________________________________________________
UInt_t AliITSPlaneEffSDD::WingFromAnode(const UInt_t anode) const {
  // return the wing number [0,1] given the anode number 
  // input: anode number [0,511]
if(anode>=kNAnode*kNChip*kNWing) 
 {Error("AliITSPlaneEffSDD::GetWingFromAnode","you asked for a non existing anode"); return 99;}
Int_t wing=0;
if(anode>=kNChip*kNAnode) wing=1;
return wing;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetKey(const UInt_t mod, const UInt_t chip,
                                 const UInt_t wing, const UInt_t subw) const {
  // get key given a basic block
if(mod>=kNModule || chip>=kNChip || wing>= kNWing || subw>=kNSubWing)
  {Error("AliITSPlaneEffSDD::GetKey","you asked for a non existing block"); return 99999;}
return mod*kNChip*kNWing*kNSubWing+chip*kNWing*kNSubWing+wing*kNSubWing+subw;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetModFromKey(const UInt_t key) const {
  // get mod. from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetModFromKey","you asked for a non existing key"); return 9999;}
return key/(kNChip*kNWing*kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetChipFromKey(const UInt_t key) const {
  // retrieves chip from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetChipFromKey","you asked for a non existing key"); return 999;}
return (key%(kNChip*kNWing*kNSubWing))/(kNWing*kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetWingFromKey(const UInt_t key) const {
  // retrieves wing from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetWingFromKey","you asked for a non existing key"); return 99;}
return ((key%(kNChip*kNWing*kNSubWing))%(kNWing*kNSubWing))/(kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetSubWingFromKey(const UInt_t key) const {
  // retrieves sub-wing from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetSubWingFromKey","you asked for a non existing key"); return 9;}
return ((key%(kNChip*kNWing*kNSubWing))%(kNWing*kNSubWing))%(kNSubWing);
}
//__________________________________________________________________________
void AliITSPlaneEffSDD::GetAllFromKey(const UInt_t key,UInt_t& mod,UInt_t& chip,
                                      UInt_t& wing,UInt_t& subw) const {
  // get module, chip, wing and subwing from a key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetAllFromKey","you asked for a non existing key"); 
  mod=9999;
  chip=999;
  wing=99;
  subw=9;
  return;}
mod=GetModFromKey(key);
chip=GetChipFromKey(key);
wing=GetWingFromKey(key);
subw=GetSubWingFromKey(key);
return;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSDD::LivePlaneEff(UInt_t key) const {
  // returns plane efficieny after adding the fraction of sensor which is bad
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::LivePlaneEff","you asked for a non existing key");
   return -1.;}
Double_t leff=AliITSPlaneEff::LivePlaneEff(0); // this just for the Warning
leff=PlaneEff(key)+GetFracBad(key);
return leff>1?1:leff;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSDD::ErrLivePlaneEff(UInt_t key) const {
  // returns error on live plane efficiency
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::ErrLivePlaneEff","you asked for a non existing key");
   return -1.;}
Int_t nf=fFound[key];
Double_t triedInLive=GetFracLive(key)*fTried[key];
Int_t nt=TMath::Max(nf,TMath::Nint(triedInLive));
return AliITSPlaneEff::ErrPlaneEff(nf,nt); // for the time being: to be checked
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSDD::GetFracLive(const UInt_t key) const {
  // returns the fraction of the sensor which is OK
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetFracLive","you asked for a non existing key");
   return -1.;}
    // Compute the fraction of bad (dead+noisy) detector 
UInt_t bad=0;
GetBadInBlock(key,bad);
Double_t live=bad;
live/=(kNAnode);
return 1.-live;
}
//_____________________________________________________________________________
void AliITSPlaneEffSDD::GetBadInBlock(const UInt_t key, UInt_t& nrBadInBlock) const {
  // Compute the number of bad (dead+noisy) anodes inside a block
  // (it depends on the chip, not on the sub-wing)
nrBadInBlock=0;
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetBadInBlock","you asked for a non existing key");
   return;}
//
if(!fInitCDBCalled) 
  {Error("AliITSPlaneEffSDD::GetBadInBlock","CDB not inizialized: call InitCDB first");
   return;};
AliCDBManager* man = AliCDBManager::Instance();
// retrieve map of dead Pixel 
AliCDBEntry *cdbSDD = man->Get("ITS/Calib/CalibSDD", fRunNumber);
TObjArray* sddEntry;
if(cdbSDD) {
  sddEntry = (TObjArray*)cdbSDD->GetObject();
  if(!sddEntry) 
  {Error("AliITSPlaneEffSDD::GetBadInBlock"," SDDEntry not found in CDB");
   return;}
} else {
  Error("AliITSPlaneEffSDD::GetBadInBlock","Did not find Calib/CalibSDD");
  return;
}
//
UInt_t mod=GetModFromKey(key);
UInt_t chip=GetChipFromKey(key);
UInt_t wing=GetWingFromKey(key);
// count number of dead
AliITSCalibrationSDD* calibSDD=(AliITSCalibrationSDD*) sddEntry->At(mod);
UInt_t nrBad = calibSDD-> GetDeadChannels();
for (UInt_t index=0; index<nrBad; index++) {
  if(ChipFromAnode(calibSDD->GetBadChannel(index))==chip &&
     WingFromAnode(calibSDD->GetBadChannel(index))==wing    ) nrBadInBlock++;
}
return;
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSDD::GetFracBad(const UInt_t key) const {
  // returns 1-fractional live
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {Error("AliITSPlaneEffSDD::GetFracBad","you asked for a non existing key");
   return -1.;}
return 1.-GetFracLive(key);
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSDD::WriteIntoCDB() const {
// write onto CDB
if(!fInitCDBCalled)
  {Error("AliITSPlaneEffSDD::WriteIntoCDB","CDB not inizialized: call InitCDB first");
   return kFALSE;}
// to be written properly: now only for debugging 
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetObjectClassName("AliITSPlaneEff");
  md->SetResponsible("Giuseppe Eugenio Bruno");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion("head 02/01/08"); //root version
  AliCDBId id("ITS/PlaneEff/PlaneEffSDD",0,AliCDBRunRange::Infinity()); 
  AliITSPlaneEffSDD eff; 
  eff=*this;
  Bool_t r=AliCDBManager::Instance()->GetDefaultStorage()->Put(&eff,id,md);
  delete md;
  return r;
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSDD::ReadFromCDB() {
// read from CDB
if(!fInitCDBCalled)
  {Error("AliITSPlaneEffSDD::ReadFromCDB","CDB not inizialized: call InitCDB first");
   return kFALSE;}
//if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
//    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
//  }
AliCDBEntry *cdbEntry = AliCDBManager::Instance()->Get("ITS/PlaneEff/PlaneEffSDD",fRunNumber);
AliITSPlaneEffSDD* eff= (AliITSPlaneEffSDD*)cdbEntry->GetObject();
if(this==eff) return kFALSE;
eff->Copy(*this);
return kTRUE;
}
