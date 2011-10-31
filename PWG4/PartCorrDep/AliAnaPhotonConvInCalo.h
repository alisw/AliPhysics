#ifndef ALIANAPHOTONINCALO_H
#define ALIANAPHOTONINCALO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Conversions pairs analysis
// Check if cluster comes from a conversion in the material in front of the calorimeter
// Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
// Input are selected clusters with AliAnaPhoton
//
//-- Author: Gustavo Conesa (LPSC-IN2P3-CNRS)

// --- ROOT system ---
class TH2F;
class TH1F;
class TString ;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

class AliAnaPhotonConvInCalo : public AliAnaPartCorrBaseClass {

 public: 
  AliAnaPhotonConvInCalo() ;              // default ctor
  virtual ~AliAnaPhotonConvInCalo() { ; } // virtual dtor
 private:
  AliAnaPhotonConvInCalo(const AliAnaPhotonConvInCalo & g) ;               // cpy ctor
  AliAnaPhotonConvInCalo & operator = (const AliAnaPhotonConvInCalo & g) ; // cpy assignment

 public:
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         InitParameters();

  void         MakeAnalysisFillAOD()  ;

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const ;
    
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------
    
  Float_t      GetMassCut()                     const { return fMassCut            ; }
  void         SetMassCut(Float_t m)                  { fMassCut    = m            ; }
	
  Bool_t       AreConvertedPairsInAOD()         const { return fAddConvertedPairsToAOD   ; }
  void         SwitchOnAdditionConvertedPairsToAOD()  { fAddConvertedPairsToAOD = kTRUE  ; }
  void         SwitchOffAdditionConvertedPairsToAOD() { fAddConvertedPairsToAOD = kFALSE ; }  
	
  Bool_t       AreConvertedPairsRemoved()       const { return fRemoveConvertedPair      ; }
  void         SwitchOnConvertedPairsRemoval()        { fRemoveConvertedPair  = kTRUE    ; }
  void         SwitchOffConvertedPairsRemoval()       { fRemoveConvertedPair  = kFALSE   ; }    
  
  void         SetConvAsymCut(Float_t c)              { fConvAsymCut = c           ; }
  Float_t      GetConvAsymCut()                 const { return fConvAsymCut        ; }
  
  void         SetConvDEtaCut(Float_t c)              { fConvDEtaCut = c           ; }
  Float_t      GetConvDEtaCut()                 const { return fConvDEtaCut        ; }
  
  void         SetConvDPhiCut(Float_t min, Float_t max)  { fConvDPhiMinCut = min   ;  
                                                           fConvDPhiMaxCut = max   ; }
  Float_t      GetConvDPhiMinCut()              const { return fConvDPhiMinCut     ; }
  Float_t      GetConvDPhiMaxCut()              const { return fConvDPhiMaxCut     ; }
  
  private:
 
  Bool_t   fRemoveConvertedPair;          // Remove conversion pairs
  Bool_t   fAddConvertedPairsToAOD;       // Put Converted pairs in AOD
  Float_t  fMassCut;                      // Mass cut for the conversion pairs selection  
  Float_t  fConvAsymCut;                  // Select conversion pairs when asymmetry is smaller than cut
  Float_t  fConvDEtaCut;                  // Select conversion pairs when deta of pair smaller than cut
  Float_t  fConvDPhiMinCut;               // Select conversion pairs when dphi of pair lager than cut
  Float_t  fConvDPhiMaxCut;               // Select conversion pairs when dphi of pair smaller than cut

  // Histograms
  TH1F * fhPtPhotonConv   ;               //! Number of identified photon vs transerse momentum 
  TH2F * fhEtaPhiPhotonConv  ;            //! Pseudorapidity vs Phi of identified  photon for transerse momentum > 0.5, for converted
  TH2F * fhEtaPhi05PhotonConv  ;          //! Pseudorapidity vs Phi of identified  photon for transerse momentum < 0.5, for converted
  TH2F * fhConvDeltaEta;                  //! Small mass photons, correlation in eta
  TH2F * fhConvDeltaPhi;                  //! Small mass photons, correlation in phi
  TH2F * fhConvDeltaEtaPhi;               //! Small mass photons, correlation in phi and eta
  TH2F * fhConvAsym;                      //! Small mass photons, correlation in energy asymmetry
  TH2F * fhConvPt;                        //! Small mass photons, pT of pair
  
  //Vertex distance
  TH2F * fhConvDistEta;                   //! Approx distance to vertex vs cluster Eta 
  TH2F * fhConvDistEn;                    //! Approx distance to vertex vs Energy
  TH2F * fhConvDistMass;                  //! Approx distance to vertex vs Mass
  TH2F * fhConvDistEtaCutEta;             //! Approx distance to vertex vs cluster Eta, dEta < 0.05 
  TH2F * fhConvDistEnCutEta;              //! Approx distance to vertex vs Energy, dEta < 0.05
  TH2F * fhConvDistMassCutEta;            //! Approx distance to vertex vs Mass, dEta < 0.05
  TH2F * fhConvDistEtaCutMass;            //! Approx distance to vertex vs cluster Eta, dEta < 0.05, m < 10 MeV 
  TH2F * fhConvDistEnCutMass;             //! Approx distance to vertex vs Energy, dEta < 0.05, m < 10 MeV
  TH2F * fhConvDistEtaCutAsy;             //! Approx distance to vertex vs cluster Eta, dEta < 0.05, m < 10 MeV, A < 0.1
  TH2F * fhConvDistEnCutAsy;              //! Approx distance to vertex vs energy, dEta < 0.05, m < 10 MeV, A < 0.1

  //Conversion pairs analysis histograms
  TH1F * fhPtConversionTagged;            //! Number of identified gamma from Conversion , tagged as conversion 
  TH1F * fhPtAntiNeutronTagged;           //! Number of identified gamma from AntiNeutrons gamma, tagged as conversion 
  TH1F * fhPtAntiProtonTagged;            //! Number of identified gamma from AntiProtons gamma, tagged as conversion 
  TH1F * fhPtUnknownTagged;               //! Number of identified gamma from unknown, tagged as conversion 
  
  TH2F * fhConvDeltaEtaMCConversion;      //! Small mass cluster pairs, correlation in eta, origin of both clusters is conversion
  TH2F * fhConvDeltaPhiMCConversion;      //! Small mass cluster pairs, correlation in phi, origin of both clusters is conversion
  TH2F * fhConvDeltaEtaPhiMCConversion;   //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is conversion
  TH2F * fhConvAsymMCConversion;          //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is conversion
  TH2F * fhConvPtMCConversion;            //! Small mass cluster pairs, pt of pair, origin of both clusters is conversion
  TH2F * fhConvDispersionMCConversion;    //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2
  TH2F * fhConvM02MCConversion;           //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2 

  TH2F * fhConvDeltaEtaMCAntiNeutron;     //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti neutron
  TH2F * fhConvDeltaPhiMCAntiNeutron;     //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti neutron
  TH2F * fhConvDeltaEtaPhiMCAntiNeutron;  //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti neutron
  TH2F * fhConvAsymMCAntiNeutron;         //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti neutron
  TH2F * fhConvPtMCAntiNeutron;           //! Small mass cluster pairs, pt of pair, origin of both clusters is anti neutron
  TH2F * fhConvDispersionMCAntiNeutron;   //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti neutron
  TH2F * fhConvM02MCAntiNeutron;          //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti neutron

  TH2F * fhConvDeltaEtaMCAntiProton;      //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti proton
  TH2F * fhConvDeltaPhiMCAntiProton;      //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti proton
  TH2F * fhConvDeltaEtaPhiMCAntiProton;   //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti proton
  TH2F * fhConvAsymMCAntiProton;          //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti proton
  TH2F * fhConvPtMCAntiProton;            //! Small mass cluster pairs, pt of pairs, origin of both clusters is anti proton
  TH2F * fhConvDispersionMCAntiProton;    //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti proton
  TH2F * fhConvM02MCAntiProton;           //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti proton

  TH2F * fhConvDeltaEtaMCString;          //! Small mass cluster pairs, correlation in eta, origin of both clusters is string
  TH2F * fhConvDeltaPhiMCString;          //! Small mass cluster pairs, correlation in phi, origin of both clusters is string
  TH2F * fhConvDeltaEtaPhiMCString;       //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is string
  TH2F * fhConvAsymMCString;              //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is string
  TH2F * fhConvPtMCString;                //! Small mass cluster pairs, pt of pairs, origin of both clusters is string
  TH2F * fhConvDispersionMCString;        //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvM02MCString;               //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvDistMCConversion;          //! Calculated conversion distance vs real distance to vertex       
  TH2F * fhConvDistMCConversionCuts;      //! Calculated conversion distance vs real distance to vertex        

   ClassDef(AliAnaPhotonConvInCalo,1)

} ;
 
#endif//ALIANAPHOTONINCALO_H



