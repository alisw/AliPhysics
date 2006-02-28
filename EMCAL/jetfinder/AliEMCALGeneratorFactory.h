#ifndef ALIEMCALGENERATORFACTORY_H
#define ALIEMCALGENERATORFACTORY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Class for generator factory which used in production for EMCAL.
//  Based on Config.C file 
//*-- Author: Aleksei Pavlinov (WSU)
#include <TObject.h>
#include <TMath.h>
#include "AliPDG.h"
// 23-aug-04 for kGamma and kPi0
#include <TPDGCode.h>

class AliGenerator;
class AliGenHijing;
class AliGenPythia;
class AliGenBox;
class TString;

enum PprRunFact_t 
{
    kTest50,
    kParam_8000,   kParam_4000,  kParam_2000, 
    kHijing_cent1, kHijing_cent2, 
    kHijing_per1,  kHijing_per2, kHijing_per3, kHijing_per4,  kHijing_per5,
    kHijing_jj25,  kHijing_jj50, kHijing_jj75, kHijing_jj100, kHijing_jj125, 
    kHijing_gj25,  kHijing_gj50, kHijing_gj75, kHijing_gj100, kHijing_gj125, 
    kJetPlusBg,    kGammaPlusBg, 
    kParam_8000_Ecal, kParam_4000_Ecal, 
    kJets_50,       kJets_75,      kJets_100,      kJets_200,
    kJets_100RadOn, 
    kGammaJets_50,  kGammaJets_75, kGammaJets_100, kGammaJets_200,
    kGammaJets_250, kGammaJets_300,
    kGammaGun, kGammaBox,
    kGammaBoxOne, kPi0BoxOne 
};

enum PprRadFact_t
{ // Concern only HIJING 
    kGluonRadiation, kNoGluonRadiation
};

class AliEMCALGeneratorFactory : public TObject{

 public:
  explicit AliEMCALGeneratorFactory
  (PprRunFact_t run=kJets_50, PprRadFact_t rad = kGluonRadiation);
  explicit AliEMCALGeneratorFactory(PprRunFact_t run, Float_t p);
  AliGenHijing* HijingStandard();
  AliGenPythia* PythiaJets(Float_t energy);
  AliGenPythia* PythiaGamma(Float_t energy);
  AliGenBox*    OneParticleWithFixedEnergy(Int_t type=kGamma, Float_t p=1.0);

  AliGenerator* GetGenerator() {return fGenerator;}
  AliGenerator* GetBgGenerator() {return fBgGenerator;}
  AliGenerator* GetSignalGenerator() {return fSignalGenerator;}
  PprRunFact_t  GetRunType()   {return fRunType;}
  PprRadFact_t  GetRadiation() {return fRadiation;}
 
  TString*      GetComment() {return fComment;}
  static Float_t EtaToTheta(Float_t arg) {return (180./TMath::Pi())*2.*atan(exp(-arg));}

protected:
  AliGenerator* fGenerator;        //! 
  AliGenerator* fBgGenerator;      //! 
  AliGenerator* fSignalGenerator;  //! 
  PprRunFact_t  fRunType;
  PprRadFact_t  fRadiation;
  Float_t       fMomentum;
  TString      *fComment;          //!

  ClassDef(AliEMCALGeneratorFactory,1) // Generator Factory for EMCAL production

};
#endif // ALIEMCALGENERATORFACTORY_H
