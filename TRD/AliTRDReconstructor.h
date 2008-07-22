#ifndef ALITRDRECONSTRUCTOR_H
#define ALITRDRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for TRD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"

class AliRawReader;

#include "AliDetectorRecoParam.h"
#include "AliTRDrecoParam.h"

class AliTRDReconstructor: public AliReconstructor 
{

 public:

  AliTRDReconstructor():AliReconstructor()                              { };
  AliTRDReconstructor(const AliTRDReconstructor &r):AliReconstructor(r) { };
  virtual ~AliTRDReconstructor();
  AliTRDReconstructor& operator = (const AliTRDReconstructor& /*r*/) 
                                                                 { return *this;            };

  virtual Bool_t           HasDigitConversion() const            { return kFALSE;           };
  virtual void             ConvertDigits(AliRawReader *rawReader, TTree *digitsTree) const;

  virtual void             Reconstruct(AliRawReader *rawReader, TTree *clusterTree) const;
  virtual void             Reconstruct(TTree *digitsTree, TTree *clusterTree) const;
  virtual AliTracker      *CreateTracker() const;

  virtual void             FillESD(AliRawReader */*rawReader*/, TTree *clusterTree, AliESDEvent *esd) const
                                                                 { FillESD((TTree * )NULL
                                                                 , clusterTree
                                                                 , esd);                    }
  virtual void             FillESD(TTree *digitsTree, TTree *clusterTree, AliESDEvent *esd) const;

  static const AliTRDrecoParam* GetRecoParam() { return dynamic_cast<const AliTRDrecoParam*>(AliReconstructor::GetRecoParam(2)); }

  ClassDef(AliTRDReconstructor,0)         //  Class for the TRD reconstruction

};

#endif
