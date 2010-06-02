#ifndef ALITRDTENDERSUPPLY_H
#define ALITRDTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  TRD tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////



#include <AliTenderSupply.h>

class AliTRDTenderSupply: public AliTenderSupply {
  
public:
  AliTRDTenderSupply();
  AliTRDTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliTRDTenderSupply();

  void SetPIDmethod(Int_t pidMethod) { fPIDmethod = pidMethod; }
  
  virtual void              Init();
  virtual void              ProcessEvent();
  
  
private:
  enum{
    kNNpid = 0,
    k1DLQpid = 1,
    k2DLQpid = 2
  };
  AliESDpid          *fESDpid;       //! ESD PID object

  Int_t fPIDmethod;                  // PID method
  
  AliTRDTenderSupply(const AliTRDTenderSupply&c);
  AliTRDTenderSupply& operator= (const AliTRDTenderSupply&c);
  
  ClassDef(AliTRDTenderSupply, 1);  // TRD tender task
};


#endif

