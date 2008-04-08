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

/* $Id: $ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Implementation of the class for bad channel treatment in the tracker //
// Stores 1 status bit for each SPD pixel and SDD anode:                //
//  0 = bad channel                                                     //
//  1 = good channel                                                    //
// Dead and noisy channels are read from AliITSCalibration objects //
// Origin: F.Prino, Torino, prino@to.infn.it                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliITSChannelStatus.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSDD.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

ClassImp(AliITSChannelStatus)


//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus():
TObject(),
fSPDChannelStatus(0),
fSDDChannelStatus(0)
{
  // default constructor 
  UInt_t nSPDchan=kSPDModules*kSPDNpxPerModule*kSPDNpzPerModule;
  fSPDChannelStatus=new TBits(nSPDchan);
  UInt_t nSDDchan=kSDDModules*kSDDAnodesPerModule;
  fSDDChannelStatus=new TBits(nSDDchan);
  InitDefaults();
}
//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus(AliCDBManager *cdb):
TObject(),
fSPDChannelStatus(0),
fSDDChannelStatus(0)
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

  UInt_t nSPDchan=kSPDModules*kSPDNpxPerModule*kSPDNpzPerModule;
  fSPDChannelStatus=new TBits(nSPDchan);
  UInt_t nSDDchan=kSDDModules*kSDDAnodesPerModule;
  fSDDChannelStatus=new TBits(nSDDchan);
  InitFromOCDB(deadArrSPD,noisArrSPD,calArrSDD);
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
}
//______________________________________________________________________
void AliITSChannelStatus::InitFromOCDB(TObjArray* deadArrSPD, TObjArray* noisArrSPD, TObjArray* calArrSDD){
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
}
//______________________________________________________________________
AliITSChannelStatus::AliITSChannelStatus(const AliITSChannelStatus& cstatus):
TObject(),
fSPDChannelStatus(cstatus.fSPDChannelStatus),
fSDDChannelStatus(cstatus.fSDDChannelStatus)
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
}

//______________________________________________________________________
Bool_t AliITSChannelStatus::CheckBounds(Int_t imod, Int_t iz, Int_t ix) const {
  // check for out of bounds
  if(imod<0 || imod>=kSPDModules+kSDDModules){
    AliError(Form("Module number out of range 0-%d",kSPDModules+kSDDModules));
    return kFALSE;
  }
  if(imod<kSPDModules){
    if(ix<0 || ix>=kSPDNpxPerModule || iz<0 || iz>=kSPDNpzPerModule){
      AliError("SPD: Pixel number out of range");
      return kFALSE;
    }
  }else{
    if(iz<0 || iz>=kSDDAnodesPerModule){
      AliError("SDD: anode number out of range");
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
  }else{
    imod-=kSPDModules;
    Int_t index=imod*kSDDAnodesPerModule+iz;
    return fSDDChannelStatus->TestBitNumber(index);    
  }
}
//______________________________________________________________________
void AliITSChannelStatus::SetChannelStatus(Bool_t cstatus, Int_t imod, Int_t iz, Int_t ix){
  // set status for given channel
  if(CheckBounds(imod,iz,ix)==kFALSE) return;
  if(imod<kSPDModules){
    Int_t index=imod*kSPDNpxPerModule*kSPDNpzPerModule+ix*kSPDNpzPerModule+iz;
    fSPDChannelStatus->SetBitNumber(index,cstatus);
  }else{
    imod-=kSPDModules;
    Int_t index=imod*kSDDAnodesPerModule+iz;
    fSDDChannelStatus->SetBitNumber(index,cstatus);
  }
}
