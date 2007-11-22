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
  
  Bool_t   GetRecoMode(            )       const  { return fRecoMode;       }
  inline  Int_t  GetUserCut(Int_t iCh    );
  Bool_t   GetUserCutMode(         )       const  { return fUserCutMode;    }

  void    SetRecoMode(Bool_t recoMode)           { fRecoMode=recoMode; }
  inline  void    SetUserCut(Int_t ucCh0,Int_t ucCh1,Int_t ucCh2,Int_t ucCh3,Int_t ucCh4,Int_t ucCh5,Int_t ucCh6);////set user cut
  void    SetUserCutMode(Bool_t userCutMode)      { fUserCutMode=userCutMode;   }
  
  static   AliHMPIDRecoParam *GetUserModeParam();        // make reco parameters
  
 
  protected:
   
  Bool_t  fRecoMode;                                    //kTRUE = normal tracking reco, kFALSE = HTA
  Int_t   fUserCut[7];                                 //user cut for the 7 chambers
  Bool_t  fUserCutMode;                                 //kTRUE =  get user cut from OCDB, kFALSE = get user cut from AliHMPIDRecoParam
  

      ClassDef(AliHMPIDRecoParam, 1) 
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void AliHMPIDRecoParam::SetUserCut(Int_t ucCh0,Int_t ucCh1,Int_t ucCh2,Int_t ucCh3,Int_t ucCh4,Int_t ucCh5,Int_t ucCh6)
{ 
  //
  // Set user cuts, ... make it nicer later ...
  //
  
  fUserCut[0]=ucCh0;                 
  fUserCut[1]=ucCh1;   
  fUserCut[2]=ucCh2;   
  fUserCut[3]=ucCh3;   
  fUserCut[4]=ucCh4;   
  fUserCut[5]=ucCh5;   
  fUserCut[6]=ucCh6;     
} //SetUserCut()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
Int_t  AliHMPIDRecoParam::GetUserCut(Int_t iCh)
{
  //
  //Return the UserCut for a given chamber
  //
  if( iCh < 0 || iCh > 6 ) return 3; // return a basic value actually does nothing , ADD AliError??? 
  else return fUserCut[iCh];   // return the actual user cut
  
}
#endif

