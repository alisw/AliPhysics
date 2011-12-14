#ifndef ALIITSCHANNELSTATUS_H
#define ALIITSCHANNELSTATUS_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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

class TString;
class AliCDBManager;
class AliITSDetTypeRec;

class AliITSChannelStatus : public TObject {

 public:
  AliITSChannelStatus();
  AliITSChannelStatus(TString config);
  AliITSChannelStatus(AliCDBManager *cdb);
  AliITSChannelStatus(const AliITSDetTypeRec *dtrec);
  AliITSChannelStatus(const AliITSChannelStatus& cstatus);
  AliITSChannelStatus& operator=(const AliITSChannelStatus& cstatus);
  virtual ~AliITSChannelStatus();

  void SetChannelStatus(Bool_t cstatus, Int_t imod, Int_t iz, Int_t ix=0);

  Bool_t GetChannelStatus(Int_t imod, Int_t iz, Int_t ix=0) const;

  Bool_t AnyBadInRoad(Int_t imod, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax) const;
  Float_t FractionOfBadInRoad(Int_t imod, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax) const;

  Int_t GetNSPDChannels()const {return fSPDChannelStatus->GetNbits();}
  Int_t GetNSDDChannels()const {return fSDDChannelStatus->GetNbits();}
  Int_t GetNSSDChannels()const {return fSSDChannelStatus->GetNbits();}
  
 protected:
  void CreateArrays();
  void InitDefaults();
  void InitFromOCDB(const TObjArray* deadArrSPD, const TObjArray* /* deadSparseArrSPD */, const TObjArray* noisArrSPD, const TObjArray* calArrSDD, const TObjArray* noisArrSSD);
  Bool_t CheckBounds(Int_t imod, Int_t iz, Int_t ix=0) const;
  Bool_t GetSPDLimits(Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& izmin, Int_t& izmax, Int_t& ixmin, Int_t& ixmax)  const;
  Bool_t GetSDDLimits(Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& izmin, Int_t& izmax, Int_t& izmin2, Int_t& izmax2) const;
  Bool_t GetSSDLimits(Int_t layer, Float_t zlocmin, Float_t zlocmax, Float_t xlocmin, Float_t xlocmax, Int_t& iPmin, Int_t& iPmax, Int_t& iNmin, Int_t& iNmax) const;
  enum {kSPDModules=240};
  enum {kSPDNpzPerModule=160};
  enum {kSPDNpxPerModule=256};
  enum {kSDDModules=260};
  enum {kSDDAnodesPerModule=512};
  enum {kSSDFirstModule=500};
  enum {kSSDMaxModLay5=1247};
  enum {kSSDModules=1698};
  enum {kSSDStripsPerModule=1536};
  enum {kSSDStripsPerSide=768};

  TBits *fSPDChannelStatus;  // bit map with status of SPD channels
  TBits *fSDDChannelStatus;  // bit map with status of SDD channels
  TBits *fSSDChannelStatus;  // bit map with status of SSD channels

  ClassDef(AliITSChannelStatus,2);
};
#endif

