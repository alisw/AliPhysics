#ifndef ALIDALITZAODESDMC_H
#define ALIDALITZAODESDMC_H
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

// This class is to support AOD and ESD MonteCarlo data
#include "AliVTrack.h"
#include "AliDalitzAODESDMC.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"

class AliMCEvent;
class AliMCParticle;
class AliAODMCParticle;


class AliDalitzAODESDMC: public TObject{
public:
    
    AliDalitzAODESDMC(AliMCParticle* lalimcparticle):
	TObject(),
	fIsESDMC(kTRUE),
    alimcparticle(lalimcparticle),
	aliaodmcparticle(NULL)
    {
//        fIsESDMC=kTRUE;
  //      alimcparticle=lalimcparticle;
    };
    AliDalitzAODESDMC(AliAODMCParticle* laliaodmcparticle):
	TObject(),
        fIsESDMC(kFALSE),
        alimcparticle(NULL),
        aliaodmcparticle(laliaodmcparticle)
    {
    //    fIsESDMC=kFALSE;
      //  aliaodmcparticle=laliaodmcparticle;
    };
    virtual ~AliDalitzAODESDMC();
    AliDalitzAODESDMC(const AliDalitzAODESDMC & );//Copy constructor
    AliDalitzAODESDMC & operator = (const AliDalitzAODESDMC & );//Overwrite
    Int_t GetMotherG();
    Int_t GetFirstDaughterG();
    Int_t GetLastDaughterG();
    Int_t GetNDaughtersG();
    Int_t GetPdgCodeG();
    Double_t GetRatioVxyG();
    Double_t PtG();
    Double_t PxG();
    Double_t PyG();
    Double_t PzG();
    Double_t EtaG();
    Int_t GetUniqueIDG();
    Double_t VertexOnZ();
    Double_t EnergyG();
     
     Bool_t   GetIsESDMC(){
         return fIsESDMC;
     }
        AliMCParticle* GetDalitzESDparticle(){
        return alimcparticle;
     };
        AliAODMCParticle* GetDalitzAODparticle(){
        return aliaodmcparticle;
     };
     
private:
    
    Bool_t fIsESDMC;
    AliMCParticle* alimcparticle;
    AliAODMCParticle* aliaodmcparticle;
    ClassDef( AliDalitzAODESDMC, 2 )
    };


#endif // ALIDALITZAODESDMC_H

// // 
