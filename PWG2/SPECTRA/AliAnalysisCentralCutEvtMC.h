#ifndef ALIANALYSISCENTRALCUTEVTMC_H
#define ALIANALYSISCENTRALCUTEVTMC_H


/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

// ---------------------------------------------------
//  MC event level cuts for azimuthal isotropic
//  expansion in highly central collisions analysis
//  author: Cristian Andrei
//          acristian@niham.nipne.ro
// ----------------------------------------------------



#include "AliAnalysisCuts.h"


class TObject;
class AliMCEvent;
class AliStack;


class AliAnalysisCentralCutEvtMC: public AliAnalysisCuts {
 public:
	AliAnalysisCentralCutEvtMC(const char *name="AliAnalysisCentralCutEvtMC", const char *title="MC_cuts");
    virtual ~AliAnalysisCentralCutEvtMC();

    Bool_t  IsSelected(TObject* obj);
    Bool_t  IsSelected(TList* /*list*/) {return kTRUE;}

    void SetMultiplicityRange(Int_t r1=0, Int_t r2=1000000){fReqMult = kTRUE; fMultMin=r1;  fMultMax=r2;}
    void SetDirectivityRange(Float_t r1=-1e10, Float_t r2=1e10) {fReqDir = kTRUE; fDirMin=r1; fDirMax=r2;}
    void SetDirUnitRange(Float_t r1=-1e10, Float_t r2=1e10) {fReqDirUnit = kTRUE; fDirUMin=r1; fDirUMax=r2;}


 private:
	AliAnalysisCentralCutEvtMC(const AliAnalysisCentralCutEvtMC& ref);
	AliAnalysisCentralCutEvtMC& operator=(const AliAnalysisCentralCutEvtMC& ref);

    Bool_t fReqMult, fReqDir, fReqDirUnit; //set whether to compute multiplicity, directivity or dir unity
    Double_t fMultMin, fMultMax; //stores the multiplicity cut interval
    Double_t fDirMin, fDirMax; //stores the directivity cut interval
    Double_t fDirUMin, fDirUMax; //stores the directivity unity cut interval


    Int_t CalcMult(AliMCEvent* const mcEv);
    Double_t CalcDir(AliMCEvent* const mcEv);
    Double_t CalcDirUnit(AliMCEvent* const mcEv);


    ClassDef(AliAnalysisCentralCutEvtMC, 1);
};

#endif
