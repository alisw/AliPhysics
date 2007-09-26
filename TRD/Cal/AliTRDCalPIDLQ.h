#ifndef ALITRDCALPIDLQ_H
#define ALITRDCALPIDLQ_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// PID distributions for the LQ method                                    //
//                                                                        //
// Author:                                                                //
// Alex Bercuci <A.Bercuci@gsi.de>                                        //
//                                                                        //          
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDCALPID_H
#include "AliTRDCalPID.h"
#endif

class AliTRDCalPIDLQ : public AliTRDCalPID
{

 public:

  enum {
    kNLength = 4
  };
	
  AliTRDCalPIDLQ();
  AliTRDCalPIDLQ(const Text_t *name, const Text_t *title);
  virtual        ~AliTRDCalPIDLQ();

  Bool_t          LoadReferences(Char_t* refFile);
  TObject*        GetModel(Int_t ip, Int_t iType, Int_t iPlane) const;
  static Double_t GetLength(Int_t il) { return (il<0 || il>=kNLength) ? -1. : fTrackSegLength[il]; }
         Double_t GetProbability(Int_t spec, Float_t mom, Float_t *dedx
                               , Float_t length, Int_t plane) const;

 private:

  AliTRDCalPIDLQ(const AliTRDCalPIDLQ& pd);
  AliTRDCalPIDLQ&   operator=(const AliTRDCalPIDLQ &c);

  void     Init();
  Int_t    GetModelID(Int_t mom, Int_t , Int_t) const;

 protected:

  static Float_t  fTrackSegLength[kNLength];  // Track segment lengths

  ClassDef(AliTRDCalPIDLQ, 1)                 // LQ PID reference manager

};
#endif

