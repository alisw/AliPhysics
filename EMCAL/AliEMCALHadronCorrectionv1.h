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
    virtual Double_t GetEnergy(Double_t pmom, Double_t eta, Int_t gid); 
    Double_t GetEnergy(Double_t pmom, Double_t eta) 
	{return GetEnergy(pmom,eta,7);}
    
    void SetGeometry(TString name, Double_t fs = 1.); 
    virtual ~AliEMCALHadronCorrectionv1() {}
 protected:
    AliEMCALHadronCorrectionv1(const char *name="HadronCorrectionv1", const char *title="Hadron Correction");
    
//    AliEMCALHadronCorrectionv1(const char *name="HadronCorrectionv1", const char *title="Hadron Correction",AliEMCALGeometry *geometry = NULL);
    void SetGeometry(AliEMCALGeometry *geometry);
    
 private:
    void SetParameters(TString name = "") {Warning("SetParameter","Dummy method with argument %s",name.Data());}
    
    static AliEMCALHadronCorrectionv1* fgHadrCorr;  // Pointer to global instance (singleton)
    static Double_t fgParLookup[HCPARAMETERS][HCPARAMETERSETS]; // Global array with parameters for hadronic response
    Double_t fPar[6];
    Float_t  fSamplingFraction;  // Sampling fraction
    
    
    ClassDef(AliEMCALHadronCorrectionv1,2) // Hadron correction for EMC (version for MDC)
};

#endif // ALIEMCALHADRONCORRECTIONV1_H
