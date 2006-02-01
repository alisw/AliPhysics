#ifndef MUONCDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// By Laurent Aphecetche

class TList;

TList* padList(Bool_t reset=kFALSE);

class Triplet : public TObject
{
public:
  Triplet(Int_t detElemId=0, Int_t manuId=0, Int_t manuChannel=0)
  : TObject(),fDetElemId(detElemId),fManuId(manuId),fManuChannel(manuChannel)
{}
  virtual ~Triplet() {}
  
  Int_t DetElemId() const { return fDetElemId; }
  Int_t ManuId() const { return fManuId; }
  Int_t ManuChannel() const { return fManuChannel; }
  
private:
    Int_t fDetElemId;
  Int_t fManuId;
  Int_t fManuChannel;
  ClassDef(Triplet,1)
};


#endif
