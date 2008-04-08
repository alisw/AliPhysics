#ifndef ALIITSCHANNELSTATUS_H
#define ALIITSCHANNELSTATUS_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:  $ */

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Class  for bad channel treatment in the tracker                 //
// Stores 1 status bit for each SPD pixel and SDD anode:           //
//  0 = bad channel                                                //
//  1 = good channel                                               //
// Dead and noisy channels are read from AliITSCalibration objects //
// Origin: F.Prino, Torino, prino@to.infn.it                       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TBits.h>
#include "AliCDBManager.h"

class AliITSChannelStatus : public TObject {

 public:
  AliITSChannelStatus();
  AliITSChannelStatus(AliCDBManager *cdb);
  AliITSChannelStatus(const AliITSChannelStatus& cstatus);
  AliITSChannelStatus& operator=(const AliITSChannelStatus& cstatus);
  virtual ~AliITSChannelStatus();

  void SetChannelStatus(Bool_t cstatus, Int_t imod, Int_t iz, Int_t ix=0);

  Bool_t GetChannelStatus(Int_t imod, Int_t iz, Int_t ix=0) const;
  Int_t GetNSPDChannels()const {return fSPDChannelStatus->GetNbits();}
  Int_t GetNSDDChannels()const {return fSDDChannelStatus->GetNbits();}

 protected:
  void InitDefaults();
  void InitFromOCDB(TObjArray* deadArrSPD, TObjArray* noisArrSPD, TObjArray* calArrSDD);
  Bool_t CheckBounds(Int_t imod, Int_t iz, Int_t ix=0) const;
  enum {kSPDModules=240};
  enum {kSPDNpzPerModule=160};
  enum {kSPDNpxPerModule=256};
  enum {kSDDModules=260};
  enum {kSDDAnodesPerModule=512};

  TBits *fSPDChannelStatus;  // bit map with status of SPD channels
  TBits *fSDDChannelStatus;  // bit map with status of SDD channels

  ClassDef(AliITSChannelStatus,1);
};
#endif
