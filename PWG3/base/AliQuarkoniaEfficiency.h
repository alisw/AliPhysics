#ifndef ALIQUARKONIAEFFICIENCY_H
#define ALIQUARKONIAEFFICIENCY_H

/* $Id$ */ 

//===================================================================
//  Class AliQuarkoniaEfficiency                               
//
//  This class will provide the quarkonia reconstruction efficiency 
//  in ALICE without acceptance consideration
//  for different resonances :
//    kJpsi
//    kPsiP
//    kUpsilon
//    kUpsilonP
//    kUpsilonPP
//  for some vector mesons :
//    kPhi
//    kOmega
//  different decay channels:
//    kDimuon
//    kDielectron
//  different trigger configurations:
//    kSinglePlusLpt, kSinglePlusHpt, kSinglePlusApt,
//    kSingleMinusLpt, kSingleMinusHpt, kSingleMinusApt,
//    kSingleUndefLpt, kSingleUndefHpt, kSingleUndefApt,
//    kPairUnlikeLpt, kPairUnlikeHpt, kPairUnlikeApt,
//    kPairLikeLpt, kPairLikeHpt, kPairLikeApt
//  different parameterizations:
//    kFlat
//    kCDFscaled
//    kCDFscaledPP
//
//
//  Reconstruction efficiency has been evaluated by means of a flat
//  y and pt distribution of quarkonia in -4 < y < -2.5, 
//  0 < pt < 20 GeV/c. Weights have been used to evaluate the
//  reconstruction efficiency in different parameterizations.
//                                                              
//  Subatech 2006
//===================================================================

#include "TNamed.h"

class TH2F;
class TString;

class AliQuarkoniaEfficiency : public TNamed{

 public:

  enum quarkonia{kJpsi, kPsiP, kUpsilon, kUpsilonP, kUpsilonPP, kOmega, kPhi};
  enum decay{kDimuon, kDielectron};
  enum parameterization{kFlat, kCDFscaled, kCDFscaledPP};
  enum trigger{kSinglePlusLpt, kSinglePlusHpt, kSinglePlusApt,
	       kSingleMinusLpt, kSingleMinusHpt, kSingleMinusApt,
	       kSingleUndefLpt, kSingleUndefHpt, kSingleUndefApt,
	       kPairUnlikeLpt, kPairUnlikeHpt, kPairUnlikeApt,
	       kPairLikeLpt, kPairLikeHpt, kPairLikeApt };


  AliQuarkoniaEfficiency(Int_t quarkoniaResonance=kJpsi, Int_t decayChannel=kDimuon, Int_t simParameterization=kCDFscaledPP);
  virtual ~AliQuarkoniaEfficiency();    
  void   Init(); 

  TH2F*  GetEfficiencyHisto() const;
  void   GetEfficiency(Float_t rap, Float_t pT, Double_t & eff, Double_t & error); 

  void  SetEfficiencyFileName(char * efficiencyFileName) { fEfficiencyFileName = efficiencyFileName; }
  void  SetQuarkoniaResonance(Int_t quarkoniaResonance = kJpsi) { fQuarkoniaResonance= quarkoniaResonance;}
  void  SetDecayChannel(Int_t decayChannel = kDimuon) { fDecayChannel = decayChannel;}
  void  SetSimulatedParameterization(Int_t simParameterization = kCDFscaledPP) { 
    fParameterization = simParameterization;
  }
  void  SetTrigger(Bool_t trig = kFALSE, Int_t triggerType = kPairUnlikeApt){ 
    fTrigger = trig; fTriggerType = triggerType;
  }
  
 private:

  AliQuarkoniaEfficiency(const AliQuarkoniaEfficiency& rhs);
  AliQuarkoniaEfficiency& operator=(const AliQuarkoniaEfficiency& rhs);

  TString        fEfficiencyFileName;      // Name of the efficiency root file
  Int_t          fQuarkoniaResonance;      // Resonance Acceptance
  Int_t          fDecayChannel;            // Studied decay channel
  Int_t          fParameterization;        // Quarkonia simulated parameterization  
  Int_t          fTriggerType;             // Trigger type to be considered  
  Bool_t         fTrigger;                 // Boolean to decide if consider or not trigger
  TH2F *         fEfficiency;              // Efficiency histogram

  ClassDef(AliQuarkoniaEfficiency,1)
    };

#endif
