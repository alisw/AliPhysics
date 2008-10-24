#ifndef AliAnaGammaDirect_H
#define AliAnaGammaDirect_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaGammaDirect.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________

// Class for the analysis of prompt gamma, isolation cut. 
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h>  
#include <TH2F.h>
#include <TString.h>

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
class AliAODParticleCorrelations ;

class TList ;

class AliAnaGammaDirect : public AliAnaPartCorrBaseClass {

public: 

  AliAnaGammaDirect() ; // default ctor
  AliAnaGammaDirect(const AliAnaGammaDirect & g) ; // cpy ctor
  AliAnaGammaDirect & operator = (const AliAnaGammaDirect & g) ;//cpy assignment
  virtual ~AliAnaGammaDirect() ; //virtual dtor

  enum mcTypes {kPrompt, kFragmentation, kPi0Decay, kEtaDecay, kOtherDecay, kPi0, kEta, kElectron, kConversion, kUnknown};

  Bool_t CheckInvMass(const Int_t icalo,const TLorentzVector mom, Double_t *v, TClonesArray * pl);
  Int_t CheckOrigin(const Int_t label);
  
  TList *  GetCreateOutputObjects();

  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 

  void MakeSeveralICAnalysis(AliAODParticleCorrelation * ph, Double_t v[3]); 
  
  void Print(const Option_t * opt)const;
  
  TString GetDetector()   const {return fDetector ; }
  void SetDetector(TString det)    {fDetector = det ; }

  Int_t    GetNCones()                  const {return fNCones ; }
  Int_t    GetNPtThresFrac()                const {return fNPtThresFrac ; }
  Float_t GetConeSizes(Int_t i)      const {return fConeSizes[i] ; }
  Float_t GetPtThresholds(Int_t i)  const {return fPtThresholds[i] ; }
  Float_t GetPtFractions(Int_t i)  const {return fPtFractions[i] ; }

  void InitParameters();
 
  void SetNCones(Int_t ncs)              {fNCones = ncs ; }
  void SetNPtThresFrac(Int_t npt)        {fNPtThresFrac = npt; }
  void SetConeSizes(Int_t i, Float_t r)         {fConeSizes[i] = r ; }
  void SetPtThresholds(Int_t i, Float_t pt)   {fPtThresholds[i] = pt; }
  void SetPtFractions(Int_t i, Float_t pt)   {fPtFractions[i] = pt; } 

  Bool_t IsIsolationOn() {return fMakeIC ; }
  void SwitchOnIsolation() { fMakeIC = kTRUE;}
  void SwitchOffIsolation() { fMakeIC = kFALSE;}

  Bool_t IsReIsolationOn() {return fReMakeIC ; }
  void SwitchOnReIsolation() { fReMakeIC = kTRUE;}
  void SwitchOffReIsolation() { fReMakeIC = kFALSE;}

  Bool_t IsSeveralIsolationOn() {return fMakeSeveralIC ; }
  void SwitchOnSeveralIsolation() { fMakeSeveralIC = kTRUE;}
  void SwitchOffSeveralIsolation() { fMakeSeveralIC = kFALSE;}

  Bool_t IsInvariantMassOn() {return fMakeInvMass ; }
  void SwitchOnInvariantMass() { fMakeInvMass = kTRUE;}
  void SwitchOffInvariantMass() { fMakeInvMass = kFALSE;}

  Bool_t SelectCluster(AliAODCaloCluster * calo, Double_t vertex[3], TLorentzVector & mom);

  private:
 
  TString fDetector ; // Detector where the gamma is searched;
  Bool_t fMakeIC ; //Do isolation analysis
  Bool_t fReMakeIC ; //Do isolation analysis
  Bool_t fMakeSeveralIC ; //Do analysis for different IC
  Bool_t fMakeInvMass; //Select candidate if no pair from decay

  //Histograms  
  TH1F * fhPtGamma    ;  //!Number of identified (isolated) gamma 
  TH2F * fhPhiGamma  ; //! Phi of identified  (isolated) gamma
  TH2F * fhEtaGamma  ; //! eta of identified  (isolated) gamma
  TH2F * fhConeSumPt ; //! Sum Pt in the cone

  //Prompt photon analysis data members for multiple cones and pt thresholds 
  Int_t         fNCones   ; //!Number of cone sizes to test
  Int_t         fNPtThresFrac ; //!Number of ptThres and ptFrac to test
  Float_t     fConeSizes[5] ; //! Array with cones to test
  Float_t     fPtThresholds[5] ; //! Array with pt thresholds to test
  Float_t     fPtFractions[5] ; //! Array with pt thresholds to test

  TH1F* fhPtThresIsolated[5][5]; //! Isolated gamma with pt threshold 
  TH1F* fhPtFracIsolated[5][5]; //! Isolated gamma with pt threshold 
  TH2F* fhPtSumIsolated[5] ;  //!  Isolated gamma with threshold on cone pt sume

  //MC
  TH1F * fhPtPrompt; //!Number of identified (isolated) prompt gamma 
  TH2F * fhPhiPrompt;  //! Phi of identified  (isolated) prompt gamma
  TH2F * fhEtaPrompt;  //! eta of identified  (isolated) prompt gamma
  TH1F * fhPtThresIsolatedPrompt[5][5];  //! Isolated prompt gamma with pt threshold 
  TH1F * fhPtFracIsolatedPrompt[5][5];    //! Isolated prompt gamma with pt frac
  TH2F * fhPtSumIsolatedPrompt[5]; //!  Isolated prompt gamma with threshold on cone pt sume
  TH1F * fhPtFragmentation; //!Number of identified (isolated) fragmentation gamma 
  TH2F * fhPhiFragmentation;  //! Phi of identified  (isolated) fragmentation gamma
  TH2F * fhEtaFragmentation;  //! eta of identified  (isolated) fragmentation gamma
  TH1F * fhPtThresIsolatedFragmentation[5][5];  //! Isolated fragmentation gamma with pt threshold 
  TH1F * fhPtFracIsolatedFragmentation[5][5];    //! Isolated fragmentation gamma with pt frac
  TH2F * fhPtSumIsolatedFragmentation[5]; //!  Isolated fragmentation gamma with threshold on cone pt sume
  TH1F * fhPtPi0Decay; //!Number of identified (isolated) Pi0Decay gamma 
  TH2F * fhPhiPi0Decay;  //! Phi of identified  (isolated) Pi0Decay gamma
  TH2F * fhEtaPi0Decay;  //! eta of identified  (isolated) Pi0Decay gamma
  TH1F * fhPtThresIsolatedPi0Decay[5][5];  //! Isolated Pi0Decay gamma with pt threshold 
  TH1F * fhPtFracIsolatedPi0Decay[5][5];    //! Isolated Pi0Decay gamma with pt frac
  TH2F * fhPtSumIsolatedPi0Decay[5]; //!  Isolated Pi0Decay gamma with threshold on cone pt sume
  TH1F * fhPtOtherDecay; //!Number of identified (isolated) OtherDecay gamma 
  TH2F * fhPhiOtherDecay;  //! Phi of identified  (isolated) OtherDecay gamma
  TH2F * fhEtaOtherDecay;  //! eta of identified  (isolated) OtherDecay gamma
  TH1F * fhPtThresIsolatedOtherDecay[5][5];  //! Isolated OtherDecay gamma with pt threshold 
  TH1F * fhPtFracIsolatedOtherDecay[5][5];    //! Isolated OtherDecay gamma with pt frac
  TH2F * fhPtSumIsolatedOtherDecay[5]; //!  Isolated OtherDecay gamma with threshold on cone pt sume	
  TH1F * fhPtConversion; //!Number of identified (isolated) Conversion gamma 
  TH2F * fhPhiConversion;  //! Phi of identified  (isolated) Conversion gamma
  TH2F * fhEtaConversion;  //! eta of identified  (isolated) Conversion gamma
  TH1F * fhPtThresIsolatedConversion[5][5];  //! Isolated Conversion gamma with pt threshold 
  TH1F * fhPtFracIsolatedConversion[5][5];    //! Isolated Conversion gamma with pt frac
  TH2F * fhPtSumIsolatedConversion[5]; //!  Isolated Conversion gamma with threshold on cone pt sume
  TH1F * fhPtUnknown; //!Number of identified (isolated) Unknown gamma 
  TH2F * fhPhiUnknown;  //! Phi of identified  (isolated) Unknown gamma
  TH2F * fhEtaUnknown;  //! eta of identified  (isolated) Unknown gamma
  TH1F * fhPtThresIsolatedUnknown[5][5];  //! Isolated Unknown gamma with pt threshold 
  TH1F * fhPtFracIsolatedUnknown[5][5];    //! Isolated Unknown gamma with pt frac
  TH2F * fhPtSumIsolatedUnknown[5]; //!  Isolated Unknown gamma with threshold on cone pt sume
						
  ClassDef(AliAnaGammaDirect,1)
} ;
 

#endif //AliAnaGammaDirect_H



