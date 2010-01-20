#ifndef ALIANALYSISCENTRALCUTEVTESD_H
#define ALIANALYSISCENTRALCUTEVTESD_H

/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

// ---------------------------------------------------
//  ESD event level cuts for azimuthal isotropic
//  expansion in highly central collisions analysis
//  author: Cristian Andrei
//           acristian@niham.nipne.ro
//  --------------------------------------------------


#include "AliAnalysisCuts.h"


class TObject;
class TObjArray;
class TList;

class AliESDEvent;


class AliAnalysisCentralCutEvtESD: public AliAnalysisCuts {
public:
    AliAnalysisCentralCutEvtESD(const char *name="AliAnalysisCentralCutEvtESD", const char *title="ESD_Event_cuts");
    virtual ~AliAnalysisCentralCutEvtESD();

    Bool_t  IsSelected(TObject* obj);
    Bool_t  IsSelected(TList* /*list*/) {return kTRUE;}

    void SetMultiplicityRange(Int_t r1 = 0, Int_t r2 = 1000000) {fReqMult = kTRUE; fMultMin = r1; fMultMax = r2;}
    void SetDirectivityRange(Double_t r1 = 0.0, Double_t r2 = 1e10) {fReqDir = kTRUE; fDirMin = r1; fDirMax = r2;}
    void SetSPDMultiplicityRange(Int_t r1 = 0, Int_t r2 = 1000000) {fReqSPDMult = kTRUE; fSPDMultMin = r1; fSPDMultMax = r2;}
    void SetSPDDirectivityRange(Double_t r1 = 0.0, Double_t r2 = 1e10) {fReqSPDDir = kTRUE; fSPDDirMin = r1; fSPDDirMax = r2;}

private:
	AliAnalysisCentralCutEvtESD(const AliAnalysisCentralCutEvtESD& ref);
	AliAnalysisCentralCutEvtESD& operator=(const AliAnalysisCentralCutEvtESD& ref);

    Bool_t fReqMult, fReqDir, fReqSPDMult, fReqSPDDir; //set whether to compute multiplicity, directivity or dir unity
    Double_t fMultMin, fMultMax; //stores the multiplicity cut interval
    Double_t fDirMin, fDirMax;  //stores the directivity cut interval
    Double_t fSPDMultMin, fSPDMultMax;  //stores the SPD multiplicity cut interval
    Double_t fSPDDirMin, fSPDDirMax;  //stores the directivity unity cut interval

    TObjArray *fCutsList[10];  //list containing the cuts

    void InitCuts(); //used to initialize the cuts

    Bool_t CheckIntCuts(Int_t no, TObject *obj) const; //used to check if a track/particle is selected

    Int_t CalcMult(AliESDEvent* const esdEv);
    Double_t CalcDir(AliESDEvent* const esdEv);
    Int_t CalcSPDMult(AliESDEvent* const esdEv);
    Double_t CalcSPDDir(AliESDEvent* const esdEv);

    ClassDef(AliAnalysisCentralCutEvtESD, 1);
};

#endif
