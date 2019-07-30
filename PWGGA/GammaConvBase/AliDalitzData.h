#ifndef ALIDALITZDATA_H
#define ALIDALITZDATA_H
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
#include <memory>
#include <utility>
#include <iostream>
//using  std::unique_ptr;
using namespace std;
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliDalitzAODESD.h"

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

class AliDalitzData{
public:

    AliDalitzData();
    virtual ~AliDalitzData();
    AliDalitzData(const AliDalitzData & );//Copy constructor
    AliDalitzData & operator = (const AliDalitzData & );//Overwrite
    
    AliDalitzData(AliESDEvent* lESDEvent);
    AliDalitzData(AliAODEvent* lAODEvent);
        //fIsESD(kFALSE),
       // fAODEvent(lAODEvent)
    //ESD
       void SetInputEvent(AliESDEvent* lESDEvent){
            fIsESD=kTRUE;
            fESDEvent=lESDEvent;
        };
       void SetInputEvent(AliAODEvent* lAODEvent){
        fIsESD=kFALSE;
        fAODEvent=lAODEvent;
       };
       AliVEvent* GetInputEvent(){
        if(fIsESD) return fESDEvent;
        else return fAODEvent;
       }
    AliDalitzAODESD* GetTrack(Int_t i);
   // Int_t GetTrack(Int_t i);
    Int_t GetNumberOfTrackletsG();
    Int_t GetNumberOfTracks();
    Int_t GetNumberOfITSClustersG(Int_t i);
   const AliVVertex* GetPrimaryVertex();
     Int_t GetNumberOfV0s();
     Bool_t GetIsESD(){
      return fIsESD;
     }

     
     
private:
    
    //Bool_t fIsESD;
    AliESDEvent* fESDEvent;
    AliAODEvent* fAODEvent;
    Bool_t fIsESD;
ClassDef( AliDalitzData, 1 );
};

#endif // ALIDALITZDATA_H

