#ifndef ALIANAPARTICLEISOLATION_H
#define ALIANAPARTICLEISOLATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaParticleIsolation.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________

// Class for the analysis of particle isolation
// Input is selected particles put in AOD branch (AliAODPWG4ParticleCorrelation)
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
class AliAODPWG4Particle;
class AliAODPWG4ParticleCorrelation ;


class AliAnaParticleIsolation : public AliAnaPartCorrBaseClass {

 public:   
  AliAnaParticleIsolation() ; // default ctor
  virtual ~AliAnaParticleIsolation() ; //virtual dtor

 private:
  AliAnaParticleIsolation(const AliAnaParticleIsolation & g) ; // cpy ctor
  AliAnaParticleIsolation & operator = (const AliAnaParticleIsolation & g) ;//cpy assignment

 public:

  Bool_t CheckInvMass(const Int_t icalo,const AliAODPWG4Particle * ph) const ;
  
  TList *  GetCreateOutputObjects();
  
  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 
  
  void MakeSeveralICAnalysis(AliAODPWG4ParticleCorrelation * ph); 
  
  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
  
  Int_t   GetNCones()              const {return fNCones ; }
  Int_t   GetNPtThresFrac()        const {return fNPtThresFrac ; }
  Float_t GetConeSizes(Int_t i)    const {return fConeSizes[i] ; }
  Float_t GetPtThresholds(Int_t i) const {return fPtThresholds[i] ; }
  Float_t GetPtFractions(Int_t i)  const {return fPtFractions[i] ; }
  
  void InitParameters();
  
  void SetNCones(Int_t ncs)                  {fNCones = ncs ; }
  void SetNPtThresFrac(Int_t npt)            {fNPtThresFrac = npt; }
  void SetConeSizes(Int_t i, Float_t r)      {fConeSizes[i] = r ; }
  void SetPtThresholds(Int_t i, Float_t pt)  {fPtThresholds[i] = pt; }
  void SetPtFractions(Int_t i, Float_t pt)   {fPtFractions[i] = pt; } 
  
  Bool_t IsReIsolationOn() const {return fReMakeIC ; }
  void SwitchOnReIsolation()  { fReMakeIC = kTRUE;}
  void SwitchOffReIsolation() { fReMakeIC = kFALSE;}
  
  Bool_t IsSeveralIsolationOn() const {return fMakeSeveralIC ; }
  void SwitchOnSeveralIsolation()  { fMakeSeveralIC = kTRUE;}
  void SwitchOffSeveralIsolation() { fMakeSeveralIC = kFALSE;}
  
  Bool_t IsInvariantMassOn() const {return fMakeInvMass ; }
  void SwitchOnInvariantMass()  { fMakeInvMass = kTRUE;}
  void SwitchOffInvariantMass() { fMakeInvMass = kFALSE;}
  
  //Histogrammes setters and getters
  virtual void SetHistoPtSumRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtSumBins = n ;
    fHistoPtSumMax = max ;
    fHistoPtSumMin = min ;
  }
  
  Int_t   GetHistoNPtSumBins() const { return fHistoNPtSumBins ; }
  Float_t GetHistoPtSumMin()   const { return fHistoPtSumMin ; }
  Float_t GetHistoPtSumMax()   const { return fHistoPtSumMax ; }
  
  virtual void SetHistoPtInConeRangeAndNBins(Float_t min, Float_t max, Int_t n) {
    fHistoNPtInConeBins = n ;
    fHistoPtInConeMax = max ;
    fHistoPtInConeMin = min ;
  }
  
  Int_t   GetHistoNPtInConeBins() const { return fHistoNPtInConeBins ; }
  Float_t GetHistoPtInConeMin()   const { return fHistoPtInConeMin ; }
  Float_t GetHistoPtInConeMax()   const { return fHistoPtInConeMax ; }
  
  
 private:
  
  TString  fCalorimeter ;   // Calorimeter where neutral particles in cone for isolation are;
  Bool_t   fReMakeIC ;      // Do isolation analysis
  Bool_t   fMakeSeveralIC ; // Do analysis for different IC
  Bool_t   fMakeInvMass;    // Select candidate if no pair from decay
  
  //Histograms  
  
  TH1F * fhPtIso ;     //! Number of isolated particles
  TH2F * fhPhiIso ;    //! Phi of isolated particles
  TH2F * fhEtaIso ;    //! eta of isolated particles
  TH2F * fhConeSumPt ; //! Sum Pt in the cone
  TH2F * fhPtInCone ;  //! Particle Pt in the cone
  
  //Prompt photon analysis data members for multiple cones and pt thresholds 
  Int_t       fNCones ;          //! Number of cone sizes to test
  Int_t       fNPtThresFrac ;    //! Number of ptThres and ptFrac to test
  
  Float_t   fConeSizes[5] ;    //! Array with cones to test
  Float_t   fPtThresholds[5] ; //! Array with pt thresholds to test
  Float_t   fPtFractions[5] ;  //! Array with pt thresholds to test
  
  TH1F* fhPtThresIsolated[5][5] ; //! Isolated particle with pt threshold 
  TH1F* fhPtFracIsolated[5][5] ;  //! Isolated particle with pt threshold 
  TH2F* fhPtSumIsolated[5] ;      //! Isolated particle with threshold on cone pt sum
  
  //MC
  TH1F * fhPtIsoPrompt;   //! Number of isolated prompt gamma 
  TH2F * fhPhiIsoPrompt;  //! Phi of isolated prompt gamma
  TH2F * fhEtaIsoPrompt;  //! eta of isolated prompt gamma
  TH1F * fhPtThresIsolatedPrompt[5][5];   //! Isolated prompt gamma with pt threshold 
  TH1F * fhPtFracIsolatedPrompt[5][5];    //! Isolated prompt gamma with pt frac
  TH2F * fhPtSumIsolatedPrompt[5];        //! Isolated prompt gamma with threshold on cone pt sume
  TH1F * fhPtIsoFragmentation;   //! Number of isolated fragmentation gamma 
  TH2F * fhPhiIsoFragmentation;  //! Phi of isolated fragmentation gamma
  TH2F * fhEtaIsoFragmentation;  //! eta of isolated fragmentation gamma
  TH1F * fhPtThresIsolatedFragmentation[5][5];  //! Isolated fragmentation gamma with pt threshold 
  TH1F * fhPtFracIsolatedFragmentation[5][5];   //! Isolated fragmentation gamma with pt frac
  TH2F * fhPtSumIsolatedFragmentation[5];       //! Isolated fragmentation gamma with threshold on cone pt sume
  TH1F * fhPtIsoPi0Decay;   //!Number of isolated Pi0Decay gamma 
  TH2F * fhPhiIsoPi0Decay;  //! Phi of isolated Pi0Decay gamma
  TH2F * fhEtaIsoPi0Decay;  //! eta of isolated Pi0Decay gamma
  TH1F * fhPtThresIsolatedPi0Decay[5][5];  //! Isolated Pi0Decay gamma with pt threshold 
  TH1F * fhPtFracIsolatedPi0Decay[5][5];   //! Isolated Pi0Decay gamma with pt frac
  TH2F * fhPtSumIsolatedPi0Decay[5];       //! Isolated Pi0Decay gamma with threshold on cone pt sume
  TH1F * fhPtIsoOtherDecay;   //! Number of isolated OtherDecay gamma 
  TH2F * fhPhiIsoOtherDecay;  //! Phi of isolated OtherDecay gamma
  TH2F * fhEtaIsoOtherDecay;  //! eta of isolated OtherDecay gamma
  TH1F * fhPtThresIsolatedOtherDecay[5][5];  //! Isolated OtherDecay gamma with pt threshold 
  TH1F * fhPtFracIsolatedOtherDecay[5][5];   //! Isolated OtherDecay gamma with pt frac
  TH2F * fhPtSumIsolatedOtherDecay[5];       //! Isolated OtherDecay gamma with threshold on cone pt sume	
  TH1F * fhPtIsoConversion;   //! Number of isolated Conversion gamma 
  TH2F * fhPhiIsoConversion;  //! Phi of isolated Conversion gamma
  TH2F * fhEtaIsoConversion;  //! eta of isolated Conversion gamma
  TH1F * fhPtThresIsolatedConversion[5][5];  //! Isolated Conversion gamma with pt threshold 
  TH1F * fhPtFracIsolatedConversion[5][5];   //! Isolated Conversion gamma with pt frac
  TH2F * fhPtSumIsolatedConversion[5];       //! Isolated Conversion gamma with threshold on cone pt sume
  TH1F * fhPtIsoUnknown;   //! Number of isolated Unknown 
  TH2F * fhPhiIsoUnknown;  //! Phi of isolated Unknown
  TH2F * fhEtaIsoUnknown;  //! eta of isolated Unknown
  TH1F * fhPtThresIsolatedUnknown[5][5];  //! Isolated Unknown gamma with pt threshold 
  TH1F * fhPtFracIsolatedUnknown[5][5];   //! Isolated Unknown gamma with pt frac
  TH2F * fhPtSumIsolatedUnknown[5];       //! Isolated Unknown gamma with threshold on cone pt sume
  
  //Histograms settings
  Int_t   fHistoNPtSumBins; // Number of bins in PtSum histograms
  Float_t fHistoPtSumMax;   // PtSum maximum in histogram
  Float_t fHistoPtSumMin;	  // PtSum minimum in histogram
  Int_t   fHistoNPtInConeBins; // Number of bins in PtInCone histogram
  Float_t fHistoPtInConeMax;   // PtInCone maximum in histogram
  Float_t fHistoPtInConeMin;   // PtInCone maximum in histogram 
  
  ClassDef(AliAnaParticleIsolation,3)
    } ;


#endif //ALIANAPARTICLEISOLATION_H



