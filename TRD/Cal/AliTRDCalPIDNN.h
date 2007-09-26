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

#ifndef ALITRDCALPID_H
#include "AliTRDCalPID.h"
#endif

class AliTRDCalPIDNN : public AliTRDCalPID
{

 public:

  AliTRDCalPIDNN();
  AliTRDCalPIDNN(const Text_t *name, const Text_t *title);
  virtual  ~AliTRDCalPIDNN();
  Bool_t    LoadReferences(Char_t *refFile);
  TObject  *GetModel(Int_t ip, Int_t iType, Int_t iPlane) const;
  Double_t  GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length, Int_t plane) const;

 private:

  AliTRDCalPIDNN(const AliTRDCalPIDNN &pd);
  AliTRDCalPIDNN &operator=(const AliTRDCalPIDNN &c);
           
  void     Init();
  Int_t    GetModelID(Int_t mom, Int_t , Int_t) const;

  ClassDef(AliTRDCalPIDNN, 1) // NN PID reference manager

};
#endif
