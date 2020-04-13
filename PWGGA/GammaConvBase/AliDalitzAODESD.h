#ifndef ALIDALITZAODESD_H
#define ALIDALITZAODESD_H
/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Edgar Dominguez Rosas, Pedro Gonzalez Zamora                   *
* Version 1                                                              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// This class is to support AOD and ESD data
 
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDInputHandler.h"

class AliESDtrack;
class AliESDVertex;
class AliAODtrack;
class AliAODVertex;

class AliVEvent;
class AliVTrack;
class AliESDEvent;
class AliESDtrack;
class AliAODEvent;
class AliAODtrack;


//class AliDalitzAODESD{
class AliDalitzAODESD : public TObject{
public:
        AliDalitzAODESD();
      // virtual ~AliDalitzAODESD();
    
    AliDalitzAODESD(AliESDtrack* lESDtrack):
        TObject(),
	fIsESD(kTRUE),
	fESDtrack(lESDtrack),
	fAODtrack(NULL)
        {
      //  fIsESD=kTRUE;
      //  fESDtrack=lESDtrack;
    };
            AliDalitzAODESD(AliAODTrack* lAODtrack):
        TObject(),
	fIsESD(kFALSE),
	fESDtrack(NULL),
	fAODtrack(lAODtrack)
    {
        //fIsESD=kFALSE;
        //fAODtrack=lAODtrack;
    };
     virtual ~AliDalitzAODESD() {}
    //ESD
     AliDalitzAODESD(const AliDalitzAODESD & );//Copy constructor
     AliDalitzAODESD & operator = (const AliDalitzAODESD & );//Overwrite
     void ComputeImpactParameter();
     void ComputeImpactParameter(const AliVVertex* vx,Double_t bmag);
     Double_t GetPtG();
     Double_t GetPxG();
     Double_t GetPyG();
     Double_t GetPzG();
     Double_t GetPhiG();
     Double_t GetConstrainedParamPhiG();
     Int_t GetLabelG();
     Bool_t GetConstrainedPxPyPzG(Double_t* p) const;
     void GetParamG(const AliVVertex* vx,Double_t bmag, Double_t vector[4]) const;
     const AliExternalTrackParam* GetParamG(const AliVVertex* vx,Double_t bmag);
     void GetImpactParametersG(Float_t* p,Float_t* cov) const;
     Double_t GetSignG();
     Double_t GetEtaG();
     Double_t GetTPCNclsFG();
     Double_t GetNclsG();
     Double_t GetTPCCrossedRowsG();
     Int_t    GetITSclsG();
     Bool_t   HasPointOnITSLayerG(Int_t i);
     Double_t GetDCAxy();
     Double_t GetDCAz();
     Bool_t TestFilterBitG(UInt_t bit) const;
     Bool_t   GetIsESD(){
         return fIsESD;
     }
     AliVTrack* GetDalitzVTrack(){
        if (fIsESD==kTRUE) return dynamic_cast<AliVTrack*>(fESDtrack);
        else return dynamic_cast<AliVTrack*>(fAODtrack);
     };
        AliESDtrack* GetDalitzESDTrack(){
        return fESDtrack;
     };
        AliAODTrack* GetDalitzAODTrack(){
        return fAODtrack;
     };
     
private:
    
    Bool_t fIsESD;
    AliESDtrack* fESDtrack;
    AliAODTrack* fAODtrack;
    Double_t b[2];
    Double_t bCov[3];
    ClassDef( AliDalitzAODESD, 1 );
    };


#endif // ALIDALITZAODESD_H

// // 
