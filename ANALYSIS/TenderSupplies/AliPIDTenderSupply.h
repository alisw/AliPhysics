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
  
private:
  
  AliPIDTenderSupply(const AliPIDTenderSupply&c);
  AliPIDTenderSupply& operator= (const AliPIDTenderSupply&c);
  
  ClassDef(AliPIDTenderSupply, 1);  // PID tender task
};


#endif

