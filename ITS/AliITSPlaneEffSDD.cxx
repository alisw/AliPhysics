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
//  It is used for chip by chip efficiency (eventually with sui-bwing division 
//  along the drift direction) of the SDD,        
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
#include "AliITSgeom.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSsegmentationSDD.h"

ClassImp(AliITSPlaneEffSDD)	
//______________________________________________________________________
AliITSPlaneEffSDD::AliITSPlaneEffSDD():
  AliITSPlaneEff(){
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
    AliError("operator=: Not allowed to make a =, use default creater instead");
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
 {AliError("GetMissingTracksForGivenEff: you asked for a non existing block");
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
 {AliError("PlaneEff(UInt_t,UInt_t,UInt_t,UInt_t): you asked for a non existing block"); return -1.;}
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
 {AliError("ErrPlaneEff(UInt_t,UInt_t,UInt_t,UInt_t): you asked for a non existing block"); return -1.;}
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
 {AliError("UpDatePlaneEff: you asked for a non existing block"); return kFALSE;}
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
 {AliError("ChipAndWingFromAnode: you asked for a non existing anode"); 
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
 {AliError("ChipFromAnode: you asked for a non existing anode"); return 999;}
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
 {AliError("WingFromAnode: you asked for a non existing anode"); return 99;}
Int_t wing=0;
if(anode>=kNChip*kNAnode) wing=1;
return wing;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetKey(const UInt_t mod, const UInt_t chip,
                                 const UInt_t wing, const UInt_t subw) const {
  // get key given a basic block
if(mod>=kNModule || chip>=kNChip || wing>= kNWing || subw>=kNSubWing)
  {AliError("GetKey: you asked for a non existing block"); return 99999;}
return mod*kNChip*kNWing*kNSubWing+chip*kNWing*kNSubWing+wing*kNSubWing+subw;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetModFromKey(const UInt_t key) const {
  // get mod. from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("GetModFromKey: you asked for a non existing key"); return 9999;}
return key/(kNChip*kNWing*kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetChipFromKey(const UInt_t key) const {
  // retrieves chip from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("GetChipFromKey: you asked for a non existing key"); return 999;}
return (key%(kNChip*kNWing*kNSubWing))/(kNWing*kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetWingFromKey(const UInt_t key) const {
  // retrieves wing from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("GetWingFromKey: you asked for a non existing key"); return 99;}
return ((key%(kNChip*kNWing*kNSubWing))%(kNWing*kNSubWing))/(kNSubWing);
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetSubWingFromKey(const UInt_t key) const {
  // retrieves sub-wing from key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("GetSubWingFromKey: you asked for a non existing key"); return 9;}
return ((key%(kNChip*kNWing*kNSubWing))%(kNWing*kNSubWing))%(kNSubWing);
}
//__________________________________________________________________________
void AliITSPlaneEffSDD::GetAllFromKey(const UInt_t key,UInt_t& mod,UInt_t& chip,
                                      UInt_t& wing,UInt_t& subw) const {
  // get module, chip, wing and subwing from a key
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("GetAllFromKey: you asked for a non existing key"); 
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
  {AliError("LivePlaneEff: you asked for a non existing key");
   return -1.;}
Double_t leff=AliITSPlaneEff::LivePlaneEff(0); // this just for the Warning
leff=PlaneEff(key)+GetFracBad(key);
return leff>1?1:leff;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSDD::ErrLivePlaneEff(UInt_t key) const {
  // returns error on live plane efficiency
if(key>=kNModule*kNChip*kNWing*kNSubWing)
  {AliError("ErrLivePlaneEff: you asked for a non existing key");
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
  {AliError("GetFracLive: you asked for a non existing key");
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
  {AliError("GetBadInBlock: you asked for a non existing key");
   return;}
//
if(!fInitCDBCalled) 
  {AliError("GetBadInBlock: CDB not inizialized: call InitCDB first");
   return;};
AliCDBManager* man = AliCDBManager::Instance();
// retrieve map of dead Pixel 
AliCDBEntry *cdbSDD = man->Get("ITS/Calib/CalibSDD", fRunNumber);
TObjArray* sddEntry;
if(cdbSDD) {
  sddEntry = (TObjArray*)cdbSDD->GetObject();
  if(!sddEntry) 
  {AliError("GetBadInBlock: SDDEntry not found in CDB");
   return;}
} else {
  AliError("GetBadInBlock: Did not find Calib/CalibSDD");
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
  {AliError("GetFracBad: you asked for a non existing key");
   return -1.;}
return 1.-GetFracLive(key);
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSDD::WriteIntoCDB() const {
// write onto CDB
if(!fInitCDBCalled)
  {AliError("WriteIntoCDB: CDB not inizialized: call InitCDB first");
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
  {AliError("ReadFromCDB: CDB not inizialized: call InitCDB first");
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
//_____________________________________________________________________________
UInt_t AliITSPlaneEffSDD::GetKeyFromDetLocCoord(Int_t ilay, Int_t idet,
                                                Float_t locx, Float_t locz) const {
// method to locate a basic block from Detector Local coordinate (to be used in tracking)
//
// If kNSubWing = 1, i.e. no sub-wing subdivision, then the numbering scheme of the 
// unique key is the following, e.g. for the first detector (idet=0,  ilayer=2)
//
//			      ^x_loc (cm)
//			      |
//   _________________________|__________________________ 3.5085
//  |            |            |            |            |
//  |            |            |            |            |
//  |            |            |            |            |
//  |   key=1    |   key=3    |   key=5    |   key=7    |
//  |            |            |            |            |
//  |____________|____________|____________|____________|_0_____\  local z (cm)
//  |            |            |            |            |       /
//  |            |            |            |            |
//  |   key=0    |   key=2    |   key=4    |   key=6    |
//  |            |            |            |            |
//  |            |            |            |            |
//  |____________|____________|____________|____________| -3.5085
//-3.7632     -1.8816         0          1.1186      3.7632 
//
// for the second detector (idet=2, ilay=2), first key is 8 (bottom-left), 
// last one is 15 (upper-right), and so on. 
//
// If subwing division has been applied, then you count in each wing, starting from 
// the one with negative local ,  from the anode side (outer part) towards the 
// cathod strip (center). 
//                    first column: 
//                      bottom wing   (from below): 0,1,..,kNSubWing-1, 
//                      upper wing    (from up):    kNSubWing, ... , 2*kNSubWing-1
//                    2nd column  
//                      bottom wing   (from below): 2*kNSubWing, .. , 3*kNSubWing-1
//                      upper wing    (from up):    3*kNSubWing, ... ,4*kNSubWing-1
//                      ... 
//                    4nd (last) column :   
//                      bottom wing   (from below): 6*kNSubWing, .. , 7*kNSubWing-1
//                      upper wing    (from up):    7*kNSubWing, ... ,8*kNSubWing-1
//
// E.g. kNSubWing=2.
//
//			      ^x_loc (cm)
//		              |   
//   _________________________|__________________________ 3.5085
//  |            |            |            |            |
//  |   key=2    |   key=6    |   key=10   |   key=14   |
//  |____________|____________|____________|____________|
//  |            |            |            |            |
//  |   key=3    |   key=7    |   key=11   |   key=15   |
//  |____________|____________|____________|____________|_0_____\  local z (cm)
//  |            |            |            |            |       /
//  |   key=1    |   key=5    |   key=9    |   key=13   |
//  |____________|____________|____________|____________|
//  |            |            |            |            |
//  |   key=0    |   key=4    |   key=8    |   key=12   |
//  |____________|____________|____________|____________| -3.5085
//-3.7632     -1.8816         0          1.1186      3.7632 
//
//___________________________________________________________________________
//
UInt_t key=999999;
if(ilay<2 || ilay>3)
  {AliError("GetKeyFromDetLocCoord: you asked for a non existing layer");
   return key;}
if(ilay==2 && (idet<0 || idet>83))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}
if(ilay==3 && (idet<0 || idet>175))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}
UInt_t mod=idet;
if(ilay==3) mod+=84;
UInt_t chip=0,wing=0,subw=0;
ChipAndWingAndSubWingFromLocCoor(locx,locz,chip,wing,subw);
key=GetKey(mod,chip,wing,subw);
return key;
}
//_____________________________________________________________________________
void AliITSPlaneEffSDD::ChipAndWingAndSubWingFromLocCoor(Float_t xloc, Float_t zloc, 
                                UInt_t& chip, UInt_t& wing, UInt_t& subw) const {
AliITSgeom* geom=NULL;
//AliITSsegmentationSDD* sdd=new AliITSsegmentationSDD(geom);
AliITSsegmentationSDD sdd=AliITSsegmentationSDD(geom);
sdd.SetDriftSpeed(sdd.GetDriftSpeed()); // this only for setting fSetDriftSpeed=kTRUE !!!
Int_t ix,iz;
Int_t ntb;
UInt_t anode=0;
if(sdd.LocalToDet(xloc,zloc,ix,iz)) {  
  anode+=iz; 
  ChipAndWingFromAnode(anode,chip,wing);
  if(sdd.LocalToDet(0.,0.,ntb,iz)) {  // in this way the sub-division along time coordinate
    subw=SubWingFromTimeBin(ix,ntb); }  // is purely geometrical one and it does not 
  else {				// depen on the drift-velocity. 
    AliError("ChipAndWingAndSubWingFromLocCoor: cannot calculate n. of time bins for SubWing.");
    subw=9;
  }
} else {
  AliError("ChipAndWingAndSubWingFromLocCoor: cannot calculate anode number and time bin."); 
  chip=999; 
  wing=99; 
  subw=9;
}
delete geom;
}
//__________________________________________________________________________________
UInt_t AliITSPlaneEffSDD::SubWingFromTimeBin(const Int_t tb, const Int_t ntb) const {
if(tb<0 || tb>ntb || ntb<0) {
 AliError(Form("SubWingFromTimeBin: you asked time bin = %d with %d n. of bins",tb,ntb));
 return 9;
}
//AliDebug(Form("tb = %d, ntb= %d , NSubWing = %d",tb,ntb,kNSubWing));
Float_t h=tb;
 h/=ntb;
 h*=(kNSubWing-1);
return TMath::Nint(h);
}
