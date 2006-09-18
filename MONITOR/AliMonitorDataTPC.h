#ifndef ALIMONITORDATATPC_H
#define ALIMONITORDATATPC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>


class AliMonitorDataTPC : public TObject {
public:
  AliMonitorDataTPC();
  AliMonitorDataTPC(Int_t size);
  virtual ~AliMonitorDataTPC();
  void     SetSize(Int_t size);
  void     SetNTracks(Int_t nTracks);
  void     SetData(Int_t i, Float_t pt, Float_t eta, Float_t phi);

private:
  AliMonitorDataTPC(const AliMonitorDataTPC& data);
  AliMonitorDataTPC& operator = (const AliMonitorDataTPC& data);

  Int_t    fNTracks;   // number of TPC tracks
  Float_t* fPt;        //[fNTracks]
  Float_t* fEta;       //[fNTracks]
  Float_t* fPhi;       //[fNTracks]

  Int_t    fSize;      //! size of the arrays

  ClassDef(AliMonitorDataTPC, 1)   // data structure for the TPC monitor tree branch
};
 

#endif
