#ifndef ALIHMPIDRECOPARAM_H
#define ALIHMPIDRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class to set HMPID reconstruction parameters (normal, HTA, UserCut ...    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
//Email: Levente.Molnar@ba.infn.it
//

#include "TNamed.h"

class AliHMPIDRecoParam : public TNamed
{
 public: 
  
  AliHMPIDRecoParam();                                                                  //ctor
  AliHMPIDRecoParam(const AliHMPIDRecoParam &p);                                            //copy ctor 
  AliHMPIDRecoParam& operator=(const AliHMPIDRecoParam &p);                                 // ass. op.
  virtual ~AliHMPIDRecoParam();                                                         //dtor
  
  Bool_t   GetRecoMode(            )       const  { return fRecoMode;          }
  Bool_t   GetUserCutMode(         )       const  { return fUserCutMode;       }
  Int_t   GetUserCut(Int_t iCh)            const  { return fUserCut[iCh];      }
  void    SetRecoMode(Bool_t recoMode)            { fRecoMode=recoMode;        }
  void    SetUserCut(Int_t iChamb,Int_t ucCh)     { fUserCut[iChamb]=ucCh;     }              //set user cut for a given chamber
  void    SetUserCutMode(Bool_t userCutMode)      { fUserCutMode=userCutMode;  }
  
  static   AliHMPIDRecoParam *GetUserModeParam();        // make reco parameters
  
 
  protected:
   
  Bool_t  fRecoMode;                                    //kTRUE = normal tracking reco, kFALSE = HTA
  Int_t   fUserCut[7];                                 //user cut for the 7 chambers
  Bool_t  fUserCutMode;                                 //kTRUE =  get user cut from OCDB, kFALSE = get user cut from AliHMPIDRecoParam
  

      ClassDef(AliHMPIDRecoParam, 1) 
};
#endif

