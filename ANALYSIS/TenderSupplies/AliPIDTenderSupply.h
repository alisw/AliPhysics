#ifndef ALIPIDTENDERSUPPLY_H
#define ALIPIDTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  PID tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////



#include <AliTenderSupply.h>

class AliPIDTenderSupply: public AliTenderSupply {
  
public:
  AliPIDTenderSupply();
  AliPIDTenderSupply(const char *name, const AliTender *tender=NULL);
  
  virtual ~AliPIDTenderSupply(){;}
  
  virtual void              Init(){;}
  virtual void              ProcessEvent();

  void SetCachePID(Bool_t cachePID) { fCachePID=cachePID; }
private:
  Bool_t fCachePID;                    // Cache PID values in transient object
  
  AliPIDTenderSupply(const AliPIDTenderSupply&c);
  AliPIDTenderSupply& operator= (const AliPIDTenderSupply&c);
  
  ClassDef(AliPIDTenderSupply, 2);  // PID tender task
};


#endif

