#ifndef ALIANAPHOTON_H
#define ALIANAPHOTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaPhoton.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
//

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F ;
class TH1F;
class TString ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
//#include "AliStack.h"
//#include "TParticle.h"
class AliStack;
class TParticle;

class TList ;

class AliAnaPhoton : public AliAnaPartCorrBaseClass {

 public: 
  AliAnaPhoton() ; // default ctor
  virtual ~AliAnaPhoton() ; //virtual dtor
 private:
  AliAnaPhoton(const AliAnaPhoton & g) ; // cpy ctor
  AliAnaPhoton & operator = (const AliAnaPhoton & g) ;//cpy assignment

 public:
	
  TObjString * GetAnalysisCuts();
  TList      * GetCreateOutputObjects();

  void Init();

  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 
  
  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  void SetCalorimeter(TString det)    {fCalorimeter = det ; }

  Bool_t IsTrackMatchRejectionOn()   const {return fRejectTrackMatch ; }
  void SwitchOnTrackMatchRejection()  {fRejectTrackMatch = kTRUE ; }
  void SwitchOffTrackMatchRejection() {fRejectTrackMatch = kFALSE ; }  

  Bool_t IsCheckConversionOn()   const {return fCheckConversion ; }
  void SwitchOnConversionChecker()  {fCheckConversion = kTRUE ; }
  void SwitchOffConversionChecker() {fCheckConversion = kFALSE ; }  
	
  Bool_t AreConvertedPairsInAOD()   const {return fAddConvertedPairsToAOD ; }
  void SwitchOnAdditionConvertedPairsToAOD()  {fAddConvertedPairsToAOD = kTRUE ; }
  void SwitchOffAdditionConvertedPairsToAOD() {fAddConvertedPairsToAOD = kFALSE ; }  
	
  void InitParameters();
 
  void SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
    fMinDist = m1;
    fMinDist2 = m2;
    fMinDist3 = m3;
  }
	
  Float_t GetMassCut()    const {return fMassCut ; }
  void SetMassCut(Float_t m)    {fMassCut = m ; }
	
  void SetTimeCut(Double_t min, Double_t max) {fTimeCutMin = min; fTimeCutMax = max;}
  Double_t GetTimeCutMin() const {return fTimeCutMin;}
  Double_t GetTimeCutMax() const {return fTimeCutMax;}	
	
  void SetNCellCut(Int_t n) {fNCellsCut = n;}
  Double_t GetNCellCut() const {return fNCellsCut;}
	
	
  private:
 
  TString fCalorimeter ; // Calorimeter where the gamma is searched;
  Float_t fMinDist ;     // Minimal distance to bad channel to accept cluster
  Float_t fMinDist2;     // Cuts on Minimal distance to study acceptance evaluation
  Float_t fMinDist3;     // One more cut on distance used for acceptance-efficiency study
  Bool_t  fRejectTrackMatch ;      //If PID on, reject clusters which have an associated TPC track
  Bool_t  fCheckConversion;        // Combine pairs of clusters with mass close to 0
  Bool_t  fAddConvertedPairsToAOD; // Put Converted pairs in AOD
  Float_t fMassCut;                // Mass cut for the conversion pairs selection
  Double_t fTimeCutMin  ;    // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;    // Remove clusters/cells with time larger than this value, in ns
  Int_t fNCellsCut ;     // Accept for the analysis clusters with more than fNCellsCut cells
	
  //Histograms  
  TH1F * fhPtPhoton   ; //! Number of identified photon vs transerse momentum 
  TH2F * fhPhiPhoton  ; //! Azimuthal angle of identified  photon vs transerse momentum 
  TH2F * fhEtaPhoton  ; //! Pseudorapidity of identified  photon vs transerse momentum 
  		
  //MC
  TH1F * fhDeltaE  ; //! MC-Reco E distribution      
  TH1F * fhDeltaPt ; //! MC-Reco pT distribution
  TH1F * fhRatioE  ; //! Reco/MC E distribution      
  TH1F * fhRatioPt ; //! Reco/MC pT distribution
  TH2F * fh2E  ; //! E distribution, Reco vs MC
  TH2F * fh2Pt ; //! pT distribution, Reco vs MC
  
  TH1F * fhPtMCPhoton;   //! Number of identified gamma 
  TH2F * fhPhiMCPhoton;  //! Phi of identified gamma
  TH2F * fhEtaMCPhoton;  //! eta of identified gamma	
	
  TH1F * fhPtPrompt;   //! Number of identified prompt gamma 
  TH2F * fhPhiPrompt;  //! Phi of identified  prompt gamma
  TH2F * fhEtaPrompt;  //! eta of identified  prompt gamma

  TH1F * fhPtFragmentation;   //! Number of identified fragmentation gamma 
  TH2F * fhPhiFragmentation;  //! Phi of identified  fragmentation gamma
  TH2F * fhEtaFragmentation;  //! eta of identified  fragmentation gamma

  TH1F * fhPtISR;   //! Number of identified initial state radiation gamma 
  TH2F * fhPhiISR;  //! Phi of identified initial state radiation gamma
  TH2F * fhEtaISR;  //! eta of identified initial state radiation gamma

  TH1F * fhPtPi0Decay;   //! Number of identified Pi0Decay gamma 
  TH2F * fhPhiPi0Decay;  //! Phi of identified  Pi0Decay gamma
  TH2F * fhEtaPi0Decay;  //! eta of identified  Pi0Decay gamma

  TH1F * fhPtOtherDecay;   //! Number of identified OtherDecay gamma 
  TH2F * fhPhiOtherDecay;  //! Phi of identified  OtherDecay gamma
  TH2F * fhEtaOtherDecay;  //! eta of identified  OtherDecay gamma

  TH1F * fhPtConversion;   //! Number of identified Conversion gamma 
  TH2F * fhPhiConversion;  //! Phi of identified  Conversion gamma
  TH2F * fhEtaConversion;  //! eta of identified  Conversion gamma

  TH1F * fhPtUnknown;   //! Number of identified Unknown gamma 
  TH2F * fhPhiUnknown;  //! Phi of identified  Unknown gamma
  TH2F * fhEtaUnknown;  //! eta of identified  Unknown gamma

   ClassDef(AliAnaPhoton,8)

} ;
 

#endif//ALIANAPHOTON_H



