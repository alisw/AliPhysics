#ifndef ALITRDCALPIDNN_H
#define ALITRDCALPIDNN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// PID distributions for the NN method                                    //
//                                                                        //
// Author:                                                                //
// Alex Wilk <wilka@uni-muenster.de>                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalPID.h"

class AliTRDCalPIDNN : public AliTRDCalPID
{
 public:

  enum{
    kMLPscale  = 16000            // scaling of the MLP input to be smaller than 1
  };

  AliTRDCalPIDNN();
  AliTRDCalPIDNN(const Text_t *name, const Text_t *title);
  virtual  ~AliTRDCalPIDNN();
  Bool_t    LoadReferences(Char_t *refFile);
  TObject  *GetModel(Int_t ip, Int_t iType, Int_t iPlane) const;
  static Int_t GetModelID(Int_t mom, Int_t ii, Int_t plane);
  Double_t  GetProbability(Int_t spec, Float_t mom
                         , const Float_t * const dedx
                         , Float_t length, Int_t plane) const;

 private:

  AliTRDCalPIDNN(const AliTRDCalPIDNN &pd);
  AliTRDCalPIDNN &operator=(const AliTRDCalPIDNN &c);
           
  void     Init();

  ClassDef(AliTRDCalPIDNN, 1) // NN PID reference manager

};
#endif
