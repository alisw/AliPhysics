#ifndef ALIEMCALHADRONCORRECTIONV1_H
#define ALIEMCALHADRONCORRECTIONV1_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//                  
//*-- Author: Mark Horner (LBL/UCT)
//
#include "AliEMCALHadronCorrection.h"


#define HCPARAMETERS    6
#define HCPARAMETERSETS 2

class AliEMCALGeometry;

class AliEMCALHadronCorrectionv1: public AliEMCALHadronCorrection {
 public:
    static  AliEMCALHadronCorrectionv1* Instance();
    virtual Double_t GetEnergy(const Double_t pmom,const Double_t eta,const Int_t gid); 
    Double_t GetEnergy(const Double_t pmom, const Double_t eta) 
	{return GetEnergy(pmom,eta,7);}
    
    void SetGeometry(TString name); 
    void SetGeometry(AliEMCALGeometry *geometry); 
    virtual ~AliEMCALHadronCorrectionv1() {}
 protected:
    AliEMCALHadronCorrectionv1(const char *name="HadronCorrectionv1", const char *title="Hadron Correction");
    
//    AliEMCALHadronCorrectionv1(const char *name="HadronCorrectionv1", const char *title="Hadron Correction",AliEMCALGeometry *geometry = NULL);
 private:
    void SetParameters(TString name = "") {;}
    
    static AliEMCALHadronCorrectionv1* fHadrCorr;
    Double_t fPar[6];
    
    ClassDef(AliEMCALHadronCorrectionv1,1) // Hadron correction for EMC (version for MDC)
};

#endif // ALIEMCALHADRONCORRECTIONV1_H
