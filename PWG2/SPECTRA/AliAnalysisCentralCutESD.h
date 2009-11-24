#ifndef ALIANALYSISCENTRALCUTESD_H
#define ALIANALYSISCENTRALCUTESD_H
/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

//  ***************************************************
//  * ESD particle level cuts for azimuthal isotropic *
//  * expansion in highly central collisions analysis *
//  * author: Cristian Andrei                         *
//  *         acristian@niham.nipne.ro                *
//  ***************************************************

#include <TPDGCode.h>

#include "AliAnalysisCuts.h"

class TObject;
class TList;
class TF1;
class TString;

class AliESDtrack;


class AliAnalysisCentralCutESD: public AliAnalysisCuts {
public:
    AliAnalysisCentralCutESD(const char *name="AliAnalysisCentralCutESD", const char *title="ESD_cuts");
    virtual ~AliAnalysisCentralCutESD();

    Bool_t  IsSelected(TObject* obj);
    Bool_t  IsSelected(TList* /*list*/) {return kTRUE;}

    void SetPartType(PDG_t type) {fReqPID = kTRUE; fPartType = type;}
    void SetPIDtype(TString type) {fPIDtype = type;}
    void SetPriorFunctions(Bool_t pfunc) {fPriorsFunc = pfunc;}
    
    void SetReqIsCharged() {fReqCharge = kTRUE;} 

private:
  AliAnalysisCentralCutESD(const AliAnalysisCentralCutESD &ref);
  AliAnalysisCentralCutESD& operator=(const AliAnalysisCentralCutESD &ref);

    Bool_t fReqPID;  //kTRUE -> run the PID

    Bool_t fReqCharge; //kTRUE -> only charged particles selected

    PDG_t fPartType; //can be kPi*, kK* or kProton

    TString fPIDtype; //PID method -> can be Custom or Bayesian
    Bool_t fPriorsFunc; // if kTRUE -> use function priors 

    Double_t fPartPriors[10]; //prior probabilities
    TF1  *fElectronFunction; //momentum dependence of the prior probs
    TF1  *fMuonFunction; //momentum dependence of the prior probs
    TF1  *fPionFunction; //momentum dependence of the prior probs
    TF1  *fKaonFunction; //momentum dependence of the prior probs
    TF1  *fProtonFunction; //momentum dependence of the prior probs

    Bool_t IsA(AliESDtrack *track, PDG_t type);

    Bool_t IsCharged  (AliESDtrack* const track) const;

    Double_t GetPriors(Int_t i, Double_t P);

    ClassDef(AliAnalysisCentralCutESD, 1);
};

#endif
