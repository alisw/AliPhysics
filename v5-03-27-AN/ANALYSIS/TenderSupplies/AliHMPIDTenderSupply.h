#ifndef ALIHMPIDTENDERSUPPLY_H
#define ALIHMPIDTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  HMPID tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>
#include <AliLog.h>
#include <AliESDpid.h>

class AliESDpid;

class AliHMPIDTenderSupply: public AliTenderSupply {

public:
  AliHMPIDTenderSupply();
  AliHMPIDTenderSupply(const char *name, const AliTender *tender=NULL);

  virtual ~AliHMPIDTenderSupply(){;}

  virtual void              Init();
  virtual void              ProcessEvent();


private:

  AliHMPIDTenderSupply(const AliHMPIDTenderSupply&c);
  AliHMPIDTenderSupply& operator= (const AliHMPIDTenderSupply&c);

  ClassDef(AliHMPIDTenderSupply, 1);
};


#endif 
