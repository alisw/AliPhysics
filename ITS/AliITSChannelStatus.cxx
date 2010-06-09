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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Implementation of the class for bad channel treatment in the tracker //
// Stores 1 status bit for each SPD pixel and SDD anode:                //
//  0 = bad channel                                                     //
//  1 = good channel                                                    //
// Dead and noisy channels are read from AliITSCalibration objects      //
// Origin: F.Prino, Torino, prino@to.infn.it                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliITSChannelStatus.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliCDBEntry.h"
#include "TMath.h"
#include "AliLog.h"

ClassImp(AliITSChannelStatus)


//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus():
TObject(),
fSPDChannelStatus(0),
fSDDChannelStatus(0),
fSSDChannelStatus(0)
{
  // default constructor 
  UInt_t nSPDchan=kSPDModules*kSPDNpxPerModule*kSPDNpzPerModule;
  fSPDChannelStatus=new TBits(nSPDchan);
  UInt_t nSDDchan=kSDDModules*kSDDAnodesPerModule;
  fSDDChannelStatus=new TBits(nSDDchan);
  UInt_t nSSDchan=kSSDModules*kSSDStripsPerModule;
  fSSDChannelStatus=new TBits(nSSDchan);
  InitDefaults();
}
//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus(AliCDBManager *cdb):
TObject(),
fSPDChannelStatus(0),
fSDDChannelStatus(0),
fSSDChannelStatus(0)
{
  AliCDBEntry* spdEntryD = cdb->Get("ITS/Calib/SPDDead");
  if (!spdEntryD) AliFatal("Cannot get CDB entry for SPDDead");
  TObjArray* deadArrSPD = (TObjArray*)spdEntryD->GetObject();
  if (!deadArrSPD) AliFatal("No object found in SPDDead file");

  AliCDBEntry* spdEntryN = cdb->Get("ITS/Calib/SPDNoisy");
  if (!spdEntryN) AliFatal("Cannot get CDB entry for SPDNoisy");
  TObjArray* noisArrSPD = (TObjArray*)spdEntryN->GetObject();
  if (!noisArrSPD) AliFatal("No object found in SPDNoisy file");

  AliCDBEntry* sddEntry = cdb->Get("ITS/Calib/CalibSDD");
  if (!sddEntry) AliFatal("Cannot get CDB entry for CalibSDD");
  TObjArray* calArrSDD = (TObjArray*)sddEntry->GetObject();
  if (!calArrSDD) AliFatal("No object found in CalibSDD file");

  AliCDBEntry* ssdEntry = cdb->Get("ITS/Calib/CalibSSD");
  if (!ssdEntry) AliFatal("Cannot get CDB entry for CalibSSD");
  TObjArray* calArrSSD = (TObjArray*)ssdEntry->GetObject();
  if (!calArrSSD) AliFatal("No object found in CalibSSD file");

  UInt_t nSPDchan=kSPDModules*kSPDNpxPerModule*kSPDNpzPerModule;
  fSPDChannelStatus=new TBits(nSPDchan);
  UInt_t nSDDchan=kSDDModules*kSDDAnodesPerModule;
  fSDDChannelStatus=new TBits(nSDDchan);
  UInt_t nSSDchan=kSSDModules*kSSDStripsPerModule;
  fSSDChannelStatus=new TBits(nSSDchan);
  InitFromOCDB(deadArrSPD,noisArrSPD,calArrSDD,calArrSSD);
}
//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus(const AliITSDetTypeRec *dtrec):
TObject(),
fSPDChannelStatus(0),
fSDDChannelStatus(0),
fSSDChannelStatus(0)
{
  UInt_t nSPDchan=kSPDModules*kSPDNpxPerModule*kSPDNpzPerModule;
  fSPDChannelStatus=new TBits(nSPDchan);
  
  UInt_t nSDDchan=kSDDModules*kSDDAnodesPerModule;
  fSDDChannelStatus=new TBits(nSDDchan);

  UInt_t nSSDchan=kSSDModules*kSSDStripsPerModule;
  fSSDChannelStatus=new TBits(nSSDchan);
  
  // SPD modules
  for(Int_t imod=0; imod<kSPDModules; imod++){
    for(Int_t ix=0; ix<kSPDNpxPerModule; ix++){
      for(Int_t iz=0; iz<kSPDNpzPerModule; iz++){
	Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
	fSPDChannelStatus->SetBitNumber(index,kTRUE);
      }
    }
    Int_t ix,iz;

    // Mask SPD dead pixels
    AliITSCalibrationSPD* deadspd=(AliITSCalibrationSPD*)dtrec->GetSPDDeadModel(imod);
    for(Int_t ipix=0; ipix<deadspd->GetNrBad();ipix++){
      deadspd->GetBadPixel(ipix,ix,iz);
      Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
      fSPDChannelStatus->SetBitNumber(index,kFALSE);      
    }
    // Mask SPD noisy pixels
    AliITSCalibrationSPD* noisspd=(AliITSCalibrationSPD*)dtrec->GetCalibrationModel(imod);
    for(Int_t ipix=0; ipix<noisspd->GetNrBad();ipix++){
      noisspd->GetBadPixel(ipix,ix,iz);
      Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
      fSPDChannelStatus->SetBitNumber(index,kFALSE);
    }
  }

  // SDD modules
  for(Int_t imod=0; imod<kSDDModules; imod++){
    AliITSCalibrationSDD* calsdd=(AliITSCalibrationSDD*)dtrec->GetCalibrationModel(imod+kSPDModules);
    for(Int_t ian=0; ian<kSDDAnodesPerModule; ian++){
      Bool_t cstatus=kTRUE;
      if(calsdd->IsBadChannel(ian)) cstatus=kFALSE;
      Int_t index=imod*kSDDAnodesPerModule+ian;
      fSDDChannelStatus->SetBitNumber(index,cstatus);
    }
  }

  // SSD modules
  for (Int_t imod = 0; imod < kSSDModules; imod++) {
    AliITSCalibrationSSD* calssd=(AliITSCalibrationSSD*)dtrec->GetCalibrationModel(kSSDFirstModule+imod);
    for(Int_t ip=0; ip<kSSDStripsPerModule; ip++) {
      Int_t index=imod*kSSDStripsPerModule+ip;
      Bool_t cstatus = kTRUE;
      if (ip < 768 && calssd->IsPChannelBad(ip))
	cstatus = kFALSE;
      if (ip >= 768 && calssd->IsNChannelBad(ip-768))
	cstatus = kFALSE;

      fSSDChannelStatus->SetBitNumber(index,cstatus);
    }
  }
}
//______________________________________________________________________
void  AliITSChannelStatus::InitDefaults(){
  // fill bitmaps setting all channels as good
  for(Int_t imod=0; imod<kSPDModules; imod++){
    for(Int_t ix=0; ix<kSPDNpxPerModule; ix++){
      for(Int_t iz=0; iz<kSPDNpzPerModule; iz++){
	Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
	fSPDChannelStatus->SetBitNumber(index,kTRUE);
      }
    }
  }
  for(Int_t imod=0; imod<kSDDModules; imod++){
    for(Int_t ian=0; ian<kSDDAnodesPerModule; ian++){
      Int_t index=imod*kSDDAnodesPerModule+ian;
      fSDDChannelStatus->SetBitNumber(index,kTRUE);
    }
  }
  for(Int_t imod=0; imod<kSSDModules; imod++){
    for(Int_t is=0; is<kSSDStripsPerModule; is++){
      Int_t index=imod*kSSDStripsPerModule+is;
      fSSDChannelStatus->SetBitNumber(index,kTRUE);
    }
  }
}
//______________________________________________________________________
void AliITSChannelStatus::InitFromOCDB(TObjArray* deadArrSPD, TObjArray* noisArrSPD, TObjArray* calArrSDD, TObjArray *calArrSSD){
// fills bitmaps from arrays of AliITSCalibrationSXD objects

  // SPD modules
  for(Int_t imod=0; imod<kSPDModules; imod++){
    for(Int_t ix=0; ix<kSPDNpxPerModule; ix++){
      for(Int_t iz=0; iz<kSPDNpzPerModule; iz++){
	Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
	fSPDChannelStatus->SetBitNumber(index,kTRUE);
      }
    }
    Int_t ix,iz;

    // Mask SPD dead pixels
    AliITSCalibrationSPD* deadspd=(AliITSCalibrationSPD*)deadArrSPD->At(imod);
    for(Int_t ipix=0; ipix<deadspd->GetNrBad();ipix++){
      deadspd->GetBadPixel(ipix,ix,iz);
      Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
      fSPDChannelStatus->SetBitNumber(index,kFALSE);      
    }

    // Mask SPD noisy pixels
    AliITSCalibrationSPD* noisspd=(AliITSCalibrationSPD*)noisArrSPD->At(imod);
    for(Int_t ipix=0; ipix<noisspd->GetNrBad();ipix++){
      noisspd->GetBadPixel(ipix,ix,iz);
      Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
      fSPDChannelStatus->SetBitNumber(index,kFALSE);
    }
  }

  // SDD modules
  for(Int_t imod=0; imod<kSDDModules; imod++){
    AliITSCalibrationSDD* calsdd=(AliITSCalibrationSDD*)calArrSDD->At(imod);
    for(Int_t ian=0; ian<kSDDAnodesPerModule; ian++){
      Bool_t cstatus=kTRUE;
      if(calsdd->IsBadChannel(ian)) cstatus=kFALSE;
      Int_t index=imod*kSDDAnodesPerModule+ian;
      fSDDChannelStatus->SetBitNumber(index,cstatus);
    }
  }

  // SSD modules
  for (Int_t imod = 0; imod < kSSDModules; imod++) {
    AliITSCalibrationSSD* calssd=(AliITSCalibrationSSD*)calArrSSD->At(imod);
    for(Int_t ip=0; ip<kSSDStripsPerModule; ip++) {
      Int_t index=imod*kSSDStripsPerModule+ip;
      Bool_t cstatus = kTRUE;
      if (ip < 768 && calssd->IsPChannelBad(ip))
	cstatus = kFALSE;
      if (ip >= 768 && calssd->IsNChannelBad(ip-768))
	cstatus = kFALSE;

      fSSDChannelStatus->SetBitNumber(index,cstatus);
    }
  }
}
//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus(const AliITSChannelStatus& cstatus):
TObject(),
fSPDChannelStatus(cstatus.fSPDChannelStatus),
fSDDChannelStatus(cstatus.fSDDChannelStatus),
fSSDChannelStatus(cstatus.fSSDChannelStatus)
{
  // copy constructor 
}
//______________________________________________________________________
AliITSChannelStatus& AliITSChannelStatus::operator=(const AliITSChannelStatus& cstatus)
{
  // assignment operator
  this->~AliITSChannelStatus();
  new(this) AliITSChannelStatus(cstatus);
  return *this;
}

//______________________________________________________________________
AliITSChannelStatus::~AliITSChannelStatus(){
  // destructor
  if(fSPDChannelStatus) delete fSPDChannelStatus;
  if(fSDDChannelStatus) delete fSDDChannelStatus;
  if(fSSDChannelStatus) delete fSSDChannelStatus;
}

//______________________________________________________________________
Bool_t AliITSChannelStatus::CheckBounds(Int_t imod, Int_t iz, Int_t ix) const {
  // check for out of bounds
  if(imod<0 || imod>=kSPDModules+kSDDModules+kSSDModules){
    AliError(Form("Module number out of range 0-%d",kSPDModules+kSDDModules));
    return kFALSE;
  }
  if(imod<kSPDModules){
    if(ix<0 || ix>=kSPDNpxPerModule || iz<0 || iz>=kSPDNpzPerModule){
      AliError("SPD: Pixel number out of range");
      return kFALSE;
    }
  }else if (imod<kSSDFirstModule){
    if(iz<0 || iz>=kSDDAnodesPerModule){
      AliError("SDD: anode number out of range");
      return kFALSE;
    }
  }
  else {
    if(iz<0 || iz>=kSSDStripsPerModule){
      AliError("SSD: strip number out of range");
      return kFALSE;
    }
  }
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSChannelStatus::GetChannelStatus(Int_t imod, Int_t iz, Int_t ix) const {
  // return status of inquired channel
  if(CheckBounds(imod,iz,ix)==kFALSE) return kFALSE;
  if(imod<kSPDModules){
    Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
    return fSPDChannelStatus->TestBitNumber(index);
  } else if (imod < kSSDFirstModule) {
    imod-=kSPDModules;
    Int_t index=imod*kSDDAnodesPerModule+iz;
    return fSDDChannelStatus->TestBitNumber(index);    
  }
  else { // SSD: iz is strip number 0-767 P-side, 768 - 1535 N-side
    imod-=kSSDFirstModule;
    Int_t index=imod*kSSDStripsPerModule+iz;
    return fSSDChannelStatus->TestBitNumber(index);    
  }
}
//______________________________________________________________________
void AliITSChannelStatus::SetChannelStatus(Bool_t cstatus, Int_t imod, Int_t iz, Int_t ix){
  // set status for given channel
  if(CheckBounds(imod,iz,ix)==kFALSE) return;
  if(imod<kSPDModules){
    Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
    fSPDChannelStatus->SetBitNumber(index,cstatus);
  }else if (imod < kSSDFirstModule) {
    imod-=kSPDModules;
    Int_t index=imod*kSDDAnodesPerModule+iz;
    fSDDChannelStatus->SetBitNumber(index,cstatus);
  }
  else { // SSD: iz is strip number 0-767 P-side, 768 - 1535 N-side
    imod-=kSSDFirstModule;
    Int_t index=imod*kSSDStripsPerModule+iz;
    return fSSDChannelStatus->SetBitNumber(index,cstatus);    
  }
}
//______________________________________________________________________
Bool_t AliITSChannelStatus::GetSDDLimits(Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& izmin, Int_t& izmax, Int_t& izmin2, Int_t& izmax2) const {
  // Returns min. and max. anode numbers from local coordindate
  AliITSsegmentationSDD *seg=new AliITSsegmentationSDD();
  Float_t dummySpeed=6.5; // to avoid warnings in SDD segmentation
  Float_t tolerance=0.9999;
  seg->SetDriftSpeed(dummySpeed);
  Float_t zHalfSize=0.5*seg->Dz()/10000.;
  zHalfSize*=tolerance;
  if(zlocmin<-zHalfSize) zlocmin=-zHalfSize;
  if(zlocmax>zHalfSize) zlocmax=zHalfSize;
  if(zlocmax<-zHalfSize || zlocmin>zHalfSize){
    AliWarning("Search region completely outside module");
    return kFALSE;
  }
  Float_t xHalfSize=seg->Dx()/10000.;
  xHalfSize*=tolerance;
  if(xlocmin<-xHalfSize) xlocmin=-xHalfSize;
  if(xlocmax>xHalfSize) xlocmax=xHalfSize;
  if(xlocmax<-xHalfSize || xlocmin>xHalfSize){
    AliWarning("Search region completely outside module");
    return kFALSE;
  }
  
  Int_t iSid1=seg->GetSideFromLocalX(xlocmin);
  Int_t iSid2=seg->GetSideFromLocalX(xlocmax);
  Int_t iz1,iz2,ixdummy;
  seg->LocalToDet(xlocmin,zlocmin,ixdummy,iz1);
  seg->LocalToDet(xlocmin,zlocmax,ixdummy,iz2);
  izmin=TMath::Min(iz1,iz2);
  izmax=TMath::Max(iz1,iz2);    
  if(iSid1==iSid2){
    izmax2=izmin2=-1;
  }else{
    seg->LocalToDet(xlocmax,zlocmin,ixdummy,iz1);
    seg->LocalToDet(xlocmax,zlocmax,ixdummy,iz2);
    izmin2=TMath::Min(iz1,iz2);
    izmax2=TMath::Max(iz1,iz2);    
  }
  delete seg;
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSChannelStatus::GetSPDLimits(Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& izmin, Int_t& izmax, Int_t& ixmin, Int_t& ixmax) const {
  // Returns min. and max. pixel numbers from local coordindate
  Float_t tolerance=0.9999;
  AliITSsegmentationSPD *seg=new AliITSsegmentationSPD();
  Float_t zHalfSize=0.5*seg->Dz()/10000.;
  zHalfSize*=tolerance;
  if(zlocmin<-zHalfSize) zlocmin=-zHalfSize;
  if(zlocmax>zHalfSize) zlocmax=zHalfSize;
  if(zlocmax<-zHalfSize || zlocmin>zHalfSize){
    AliWarning("Search region completely outside module");
    return kFALSE;
  }
  Float_t xHalfSize=0.5*seg->Dx()/10000.;
  xHalfSize*=tolerance;
  if(xlocmin<-xHalfSize) xlocmin=-xHalfSize;
  if(xlocmax>xHalfSize) xlocmax=xHalfSize;
  if(xlocmax<-xHalfSize || xlocmin>xHalfSize){
    AliWarning("Search region completely outside module");
    return kFALSE;
  }

  Int_t iz1,ix1,iz2,ix2;
  seg->LocalToDet(xlocmin,zlocmin,ix1,iz1);
  seg->LocalToDet(xlocmax,zlocmax,ix2,iz2);
  izmin=TMath::Min(iz1,iz2);
  izmax=TMath::Max(iz1,iz2);
  ixmin=TMath::Min(ix1,ix2);
  ixmax=TMath::Max(ix1,ix2);
  delete seg;
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSChannelStatus::GetSSDLimits(Int_t layer, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& iPmin, Int_t& iPmax, Int_t& iNmin, Int_t& iNmax) const {
  // Returns min, max strip for SSD, given a search window
  static AliITSsegmentationSSD seg;

  AliDebug(2,Form("xmin %f zmin %f xmax %f zmax %f\n",xlocmin,zlocmin,xlocmax,zlocmax));

  Int_t p,n;
  seg.SetLayer(layer);
  seg.GetPadIxz(xlocmin,zlocmin,p,n);
  iPmin = iPmax = p;
  iNmin = iNmax = n;
  AliDebug(5,Form("lay %d xmin, zmin p %d n %d\n",layer,p,n));
  seg.GetPadIxz(xlocmax,zlocmin,p,n);
  iPmin = TMath::Min(iPmin,p);
  iPmax = TMath::Max(iPmax,p);
  iNmin = TMath::Min(iNmin,n);
  iNmax = TMath::Max(iNmax,n);
  AliDebug(5,Form("lay %d xmax, zmin p %d n %d\n",layer,p,n));
  seg.GetPadIxz(xlocmax,zlocmax,p,n);
  iPmin = TMath::Min(iPmin,p);
  iPmax = TMath::Max(iPmax,p);
  iNmin = TMath::Min(iNmin,n);
  iNmax = TMath::Max(iNmax,n);
  AliDebug(5,Form("lay %d xmax, zmax p %d n %d\n",layer,p,n));
  seg.GetPadIxz(xlocmin,zlocmax,p,n);
  iPmin = TMath::Min(iPmin,p);
  iPmax = TMath::Max(iPmax,p);
  iNmin = TMath::Min(iNmin,n);
  iNmax = TMath::Max(iNmax,n);
  AliDebug(5,Form("lay %d xmin, zmax p %d n %d\n",layer,p,n));

  if (iPmin < 0)
    iPmin = 0;
  if (iNmin < 0)
    iNmin = 0;
  if (iPmax >= kSSDStripsPerSide)
    iPmax = kSSDStripsPerSide-1;
  if (iNmax >= kSSDStripsPerSide)
    iNmax = kSSDStripsPerSide-1;
  AliDebug(2,Form("lay %d p %d %d n %d %d\n",layer,iPmin,iPmax,iNmin,iNmax));
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSChannelStatus::AnyBadInRoad(Int_t imod, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax) const{
  // Checks if there is at least one bad channel in the search road
  // For SSD: return kTRUE if there is at least one bad strip on both sides
  //          if only bad strips on one side 1D clusters can be used
  AliDebug(2,Form("checking for module %d",imod));
  if(imod<kSPDModules){
    Int_t izmin,izmax,ixmin,ixmax;
    Bool_t retcode=GetSPDLimits(zlocmin,zlocmax,xlocmin,xlocmax,izmin,izmax,ixmin,ixmax);
    if(!retcode) return kFALSE;
    for(Int_t iz=izmin; iz<=izmax;iz++){
      for(Int_t ix=ixmin; ix<=ixmax;ix++){
	if(GetChannelStatus(imod,iz,ix)==kFALSE) return kTRUE;
      }
    }
  }else if (imod < kSSDFirstModule) {
    Int_t izmin,izmax,izmin2,izmax2;
    Bool_t retcode=GetSDDLimits(zlocmin,zlocmax,xlocmin,xlocmax,izmin,izmax,izmin2,izmax2);
    if(!retcode) return kFALSE;
    for(Int_t iz=izmin; iz<=izmax;iz++){
      if(GetChannelStatus(imod,iz,0)==kFALSE) { return kTRUE;}
    }
    if(izmin2!=-1 && izmax2!=-1){
      for(Int_t iz=izmin2; iz<=izmax2;iz++){
	if(GetChannelStatus(imod,iz,0)==kFALSE) { return kTRUE;}
      }
    }
  }
  else {
    Int_t layer = 5;
    if (imod > kSSDMaxModLay5) 
      layer = 6;
    Int_t iPmin,iPmax,iNmin,iNmax;
    Bool_t retcode=GetSSDLimits(layer,zlocmin,zlocmax,xlocmin,xlocmax,iPmin,iPmax,iNmin,iNmax);
    if(!retcode) return kFALSE;
    Int_t nPbad = 0;
    for(Int_t iP=iPmin; iP<=iPmax;iP++){
      if(GetChannelStatus(imod,iP,0)==kFALSE) nPbad++;
    }
    if (nPbad == 0)
      return kFALSE;
    for(Int_t iN=iNmin; iN<=iNmax;iN++){
      if(GetChannelStatus(imod,iN+768,0)==kFALSE) { return kTRUE; }
    }
  }
  return kFALSE;
}
//______________________________________________________________________
Float_t AliITSChannelStatus::FractionOfBadInRoad(Int_t imod, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax) const{
  // Calculate the fraction of bad channels in the road  
  // Note: SSD returns fraction of dead strips on 'best' side. This 
  //       is not proportional to dead surface area. 
  Float_t totChan=0.;
  Float_t badChan=0.;
  if(imod<kSPDModules){
    Int_t izmin,izmax,ixmin,ixmax;
    Bool_t retcode=GetSPDLimits(zlocmin,zlocmax,xlocmin,xlocmax,izmin,izmax,ixmin,ixmax);
    if(!retcode) return 0.;
    for(Int_t iz=izmin; iz<=izmax;iz++){
      for(Int_t ix=ixmin; ix<=ixmax;ix++){
	totChan+=1;
	if(GetChannelStatus(imod,iz,ix)==kFALSE) badChan+=1.;
      }
    }
  }else if (imod < kSSDFirstModule) {
    Int_t izmin,izmax,izmin2,izmax2;
    Bool_t retcode=GetSDDLimits(zlocmin,zlocmax,xlocmin,xlocmax,izmin,izmax,izmin2,izmax2);
    if(!retcode) return 0.;
    for(Int_t iz=izmin; iz<=izmax;iz++){
      totChan+=1;
      if(GetChannelStatus(imod,iz,0)==kFALSE) badChan+=1.;
    }
    if(izmin2!=-1 && izmax2!=-1){
      for(Int_t iz=izmin2; iz<=izmax2;iz++){
	totChan+=1;
	if(GetChannelStatus(imod,iz,0)==kFALSE) badChan+=1.;
      }
    }
  }
  else {
    Int_t layer = 5;
    if (imod > kSSDMaxModLay5) 
      layer = 6;
    Int_t iPmin,iPmax,iNmin,iNmax;
    Bool_t retcode=GetSSDLimits(layer,zlocmin,zlocmax,xlocmin,xlocmax,iPmin,iPmax,iNmin,iNmax);
    if(!retcode) return kFALSE;
    Int_t nPbad = 0;
    for(Int_t iP=iPmin; iP<=iPmax;iP++){
      if(GetChannelStatus(imod,iP,0)==kFALSE) nPbad++;
    }
    Float_t fracP = (Float_t) nPbad / (iPmax-iPmin+1);
    if (nPbad == 0)
      return 0;
    Int_t nNbad = 0;
    for(Int_t iN=iNmin; iN<=iNmax;iN++){
      if(GetChannelStatus(imod,iN+768,0)==kFALSE) nNbad++;
    }
    Float_t fracN = (Float_t) nNbad / (iPmax-iPmin+1);
    return TMath::Min(fracP,fracN);
  }
  if(totChan==0.) return 0.;
  else return badChan/totChan;
}

