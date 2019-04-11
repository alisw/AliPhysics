#ifndef ALIDALITZEVENTMC_H
#define ALIDALITZEVENTMC_H
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
//#include "AliAODEvent.h"
#include "AliDalitzAODESDMC.h"
#include <memory>
#include <utility>
//using  std::unique_ptr;
using namespace std;

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


class AliDalitzEventMC{
public:

    
    AliDalitzEventMC();
    virtual ~AliDalitzEventMC();
    AliDalitzEventMC(const AliDalitzEventMC & );//Copy constructor
    AliDalitzEventMC & operator = (const AliDalitzEventMC & );//Overwrite
    
    AliDalitzEventMC(AliMCEvent* lESDMCEvent);
    AliDalitzEventMC(AliAODEvent* lAODMCEvent);
    
    //ESD
       void SetInputEvent(AliMCEvent* lESDMCEvent){
        fIsESDMC=kTRUE;
        fESDEventMC=lESDMCEvent;
    };
    //AOD
       void SetInputEvent(AliAODEvent* lAODMCEvent){
        fIsESDMC=kFALSE;
        fAODEvent=lAODMCEvent;
       };
       AliDalitzAODESDMC* Particle(Int_t i);
       
       
        //    AliDalitzAODESD* GetTrack(Int_t i);
        //  Int_t GetNumberOfTracks();
        //  const AliVVertex* GetPrimaryVertex();
        // Int_t GetNumberOfV0s();
     
private:
    
    AliMCEvent* fESDEventMC;
    AliAODEvent* fAODEvent;
    Bool_t fIsESDMC;
ClassDef( AliDalitzEventMC, 1 );
};

#endif // ALIDALITZEVENTMC_H
