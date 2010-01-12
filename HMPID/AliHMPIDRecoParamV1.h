#ifndef ALIHMPIDRECOPARAMV1_H
#define ALIHMPIDRECOPARAMV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class to set HMPID reconstruction parameters (normal, HTA, UserCut ...    //
//                                                                           //
// This class is introduced as a copy of teh AliHMPIDRecoParam to solve      //
// the problem occured in v4-17-Rev-18                                       //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDetectorRecoParam.h"

class AliHMPIDRecoParamV1 : public AliDetectorRecoParam
{
 public: 
  
  AliHMPIDRecoParamV1();                                                                                           //ctor
  AliHMPIDRecoParamV1(const AliHMPIDRecoParamV1 &p);                                                               //copy ctor 
  AliHMPIDRecoParamV1& operator=(const AliHMPIDRecoParamV1 &p);                                                    // ass. op.
  virtual ~AliHMPIDRecoParamV1();                                                                                  //dtor

  
  static AliHMPIDRecoParamV1 *GetLowFluxParam();                                                                   // reco params for low flux env.
  static AliHMPIDRecoParamV1 *GetHighFluxParam();                                                                  // reco params for high flux env. 
  static AliHMPIDRecoParamV1 *GetCosmicParam();                                                                    // reco params for cosmic  
    
  Bool_t   GetHmpRecoMode(            )        const            { return fHmpRecoMode;          }                  //kTRUE = normal tracking reco, kFALSE = HTA     
  void     SetHmpRecoMode(Bool_t recoMode)                       { fHmpRecoMode=recoMode;        }                 //kTRUE = normal tracking reco, kFALSE = HTA   
  Int_t    GetHmpUserCut(Int_t iCh)            const             { return fHmpUserCut[iCh];      }                 //user cut for the 7 chambers
  void     SetHmpUserCut(Int_t iChamb,Int_t ucCh)     { fHmpUserCut[iChamb]=ucCh; Printf("fUserCut[%d]=%d",iChamb,ucCh);    }       //set user cut (DAQ Sigma) for a given chamber
  Bool_t   IsFixedDistCut()                    const             { return fHmpFixedDistCut;      }                  //if kTRUE the track matching  distance is a fix number, if kFALSE the distance depends on momentum
  void     SetIsFixedDistCut(Bool_t isFix)                       { fHmpFixedDistCut=isFix;       }                  //Change from fix distance cut to parameterized
  Double_t GetHmpTrackMatchingDist()           const             { return fHmpTrackMatchingDist; }                  //Get distance between the MIP cluster
  void     SetHmpTrackMatchingDist(Double_t dist)                { fHmpTrackMatchingDist=dist;   }                  //Set distance between the MIP cluster
  Double_t GetHmpTrackMatchingDistParam(Int_t par) const         {return fHmpTrackMatchingDistParas[par];}          //Prevision to get  momentum dependen track matching parameters
  void     SetHmpTrackMatchingDistParam(Int_t par, Double_t val) {fHmpTrackMatchingDistParas[par]=val;}             //Prevision to set  momentum dependen track matching parameters

  virtual void PrintParameters() const;
  
 
  protected:
   
  Bool_t   fHmpRecoMode;                                    //kTRUE = normal tracking reco, kFALSE = HTA
  Int_t    fHmpUserCut[7];                                  //user cut for the 7 chambers
  Bool_t   fHmpFixedDistCut;                                //if kTRUE the track matching  distance is a fix number, if kFALSE the distance depends on momentum
  Double_t fHmpTrackMatchingDist;                           //distance between the MIP cluster
  Double_t fHmpTrackMatchingDistParas[5];                   //Prevision for momentum dependen track matching
  

  ClassDef(AliHMPIDRecoParamV1, 1) 
};
#endif

