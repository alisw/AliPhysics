#ifndef ALIANAPHOTONINCALO_H
#define ALIANAPHOTONINCALO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPhotonConvInCalo
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Conversions pairs clusters analysis.
///
/// Conversions pairs clusters analysis.
/// Check if cluster comes from a conversion in the material in front of the calorimeter
/// Do invariant mass of all pairs, if mass is close to 0, then it is conversion.
/// Input are selected clusters with AliAnaPhoton
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaPhotonInCalo).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
class TH2F;
class TH1F;
class TString ;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPhotonConvInCalo : public AliAnaCaloTrackCorrBaseClass {

 public: 
               AliAnaPhotonConvInCalo() ;

  /// Virtual destructor.
  virtual     ~AliAnaPhotonConvInCalo() { ; }
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         InitParameters();

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const ;
    
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------
    
  Float_t      GetMassCut()                     const { return fMassCut            ; }
  void         SetMassCut(Float_t m)                  { fMassCut    = m            ; }

  Float_t      GetMassCutTight()                const { return fMassCutTight       ; }
  void         SetMassCutTight(Float_t m)             { fMassCutTight    = m       ; }
  
  Bool_t       AreConvertedPairsInAOD()         const { return fAddConvertedPairsToAOD   ; }
  void         SwitchOnAdditionConvertedPairsToAOD()  { fAddConvertedPairsToAOD = kTRUE  ; }
  void         SwitchOffAdditionConvertedPairsToAOD() { fAddConvertedPairsToAOD = kFALSE ; }  
  
  Bool_t       AreConvertedPairsRemoved()       const { return fRemoveConvertedPair      ; }
  void         SwitchOnConvertedPairsRemoval()        { fRemoveConvertedPair  = kTRUE    ; }
  void         SwitchOffConvertedPairsRemoval()       { fRemoveConvertedPair  = kFALSE   ; }    
  
  void         SwitchOnClusterConvDistHistograms()    { fFillClusterConvDistHisto= kTRUE ; }
  void         SwitchOffClusterConvDistHistograms()   { fFillClusterConvDistHisto= kFALSE; }  
  
  void         SetConvAsymCut(Float_t c)              { fConvAsymCut = c           ; }
  Float_t      GetConvAsymCut()                 const { return fConvAsymCut        ; }
  
  void         SetConvDEtaCut(Float_t c)              { fConvDEtaCut = c           ; }
  Float_t      GetConvDEtaCut()                 const { return fConvDEtaCut        ; }
  
  void         SetConvDPhiCut(Float_t min, Float_t max)  { fConvDPhiMinCut = min   ;  
                                                           fConvDPhiMaxCut = max   ; }
  Float_t      GetConvDPhiMinCut()              const { return fConvDPhiMinCut     ; }
  Float_t      GetConvDPhiMaxCut()              const { return fConvDPhiMaxCut     ; }
  
  private:
 
  Bool_t   fRemoveConvertedPair;          ///<  Remove conversion pairs
  Bool_t   fAddConvertedPairsToAOD;       ///<  Put Converted pairs in AOD
  Bool_t   fFillClusterConvDistHisto;     ///<  Fill histograms with calculated conversion distance with data clusters
  
  Float_t  fMassCut;                      ///<  Mass cut for the conversion pairs selection
  Float_t  fMassCutTight;                 ///<  Mass cut for the conversion pairs selection, tighter
  Float_t  fConvAsymCut;                  ///<  Select conversion pairs when asymmetry is smaller than cut
  Float_t  fConvDEtaCut;                  ///<  Select conversion pairs when deta of pair smaller than cut
  Float_t  fConvDPhiMinCut;               ///<  Select conversion pairs when dphi of pair lager than cut
  Float_t  fConvDPhiMaxCut;               ///<  Select conversion pairs when dphi of pair smaller than cut

  TLorentzVector fMomentum ;              //!<! Cluster momentum
  TVector3       fProdVertex;             //!<! Production vertex
  
  // Histograms
  TH1F * fhPtPhotonConv   ;               //!<! Number of identified photon vs transverse momentum 
  TH2F * fhEtaPhiPhotonConvPaired[6]  ;   //!<! Pseudorapidity vs Phi of identified  photon conv leg  for 6 transverse momentum bins
  TH2F * fhEtaPhiPhotonConv[6]  ;         //!<! Pseudorapidity vs Phi of identified  photon conv pair for 6 transverse momentum bins
  
  TH2F * fhConvDeltaEta;                  //!<! Small mass photons, correlation in eta
  TH2F * fhConvDeltaPhi;                  //!<! Small mass photons, correlation in phi
  TH2F * fhConvDeltaEtaPhi;               //!<! Small mass photons, correlation in phi and eta
  TH2F * fhConvAsym;                      //!<! Small mass photons, correlation in energy asymmetry
  TH2F * fhConvPt;                        //!<! Small mass photons, pT of pair
//  TH2F * fhConvPtRcut[6];                 //!<! Small mass photons, pT of pair, for difference conversion distances
  
  // Vertex distance
  TH2F * fhConvDistEta;                   //!<! Approx distance to vertex vs rapidity 
  TH2F * fhConvDistPhi;                   //!<! Approx distance to vertex vs azimuth
  TH2F * fhConvDistEn;                    //!<! Approx distance to vertex vs Energy
  TH2F * fhConvDistMass;                  //!<! Approx distance to vertex vs Mass
  
  TH2F * fhConvDistEtaCutEta;             //!<! Approx distance to vertex vs rapidity, dEta < fConvDEtaCut 
  TH2F * fhConvDistPhiCutEta;             //!<! Approx distance to vertex vs azimuth
  TH2F * fhConvDistEnCutEta;              //!<! Approx distance to vertex vs Energy,   dEta < fConvDEtaCut
  TH2F * fhConvDistMassCutEta;            //!<! Approx distance to vertex vs Mass,     dEta < fConvDEtaCut
  
  TH2F * fhConvDistEtaCutMass;            //!<! Approx distance to vertex vs rapidity, M < 20 MeV/c^2 
  TH2F * fhConvDistPhiCutMass;            //!<! Approx distance to vertex vs azimuth,  M < 20 MeV/c^2
  TH2F * fhConvDistEnCutMass;             //!<! Approx distance to vertex vs Energy,   M < 20 MeV/c^2
  
  TH2F * fhConvDistEtaCutAsy;             //!<! Approx distance to vertex vs rapidity, A < fConvAsymCut
  TH2F * fhConvDistPhiCutAsy;             //!<! Approx distance to vertex vs azimuth,  A < fConvAsymCut
  TH2F * fhConvDistEnCutAsy;              //!<! Approx distance to vertex vs energy,   A < fConvAsymCut
  TH2F * fhConvDistMassCutAsy;            //!<! Approx distance to vertex vs mass,     A < fConvAsymCut

  TH2F * fhConvDistEtaCutAll;             //!<! Approx distance to vertex vs rapidity, dEta < fConvDEtaCut, M < 20 MeV/c^2, A < fConvAsymCut
  TH2F * fhConvDistPhiCutAll;             //!<! Approx distance to vertex vs azimuth,  dEta < fConvDEtaCut, M < 20 MeV/c^2, A < fConvAsymCut
  TH2F * fhConvDistEnCutAll;              //!<! Approx distance to vertex vs energy,   dEta < fConvDEtacut, M < 20 MeV/c^2, A < fConvAsymCut  
  
  // Conversion pairs analysis histograms
  TH1F * fhPtConversionTagged;            //!<! Number of identified gamma from Conversion , tagged as conversion 
  TH1F * fhPtAntiNeutronTagged;           //!<! Number of identified gamma from AntiNeutrons gamma, tagged as conversion 
  TH1F * fhPtAntiProtonTagged;            //!<! Number of identified gamma from AntiProtons gamma, tagged as conversion 
  TH1F * fhPtUnknownTagged;               //!<! Number of identified gamma from unknown, tagged as conversion 
  
  TH2F * fhConvDeltaEtaMCConversion;      //!<! Small mass cluster pairs, correlation in eta, origin of both clusters is conversion
  TH2F * fhConvDeltaPhiMCConversion;      //!<! Small mass cluster pairs, correlation in phi, origin of both clusters is conversion
  TH2F * fhConvDeltaEtaPhiMCConversion;   //!<! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is conversion
  TH2F * fhConvAsymMCConversion;          //!<! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is conversion
  TH2F * fhConvPtMCConversion;            //!<! Small mass cluster pairs, pt of pair, origin of both clusters is conversion
  TH2F * fhConvPtMCConversionRcut[6];     //!<! Small mass cluster pairs, pt of pair, origin of both clusters is conversion, for different production vertices
//TH2F * fhConvDispersionMCConversion;    //!<! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2
  TH2F * fhConvM02MCConversion;           //!<! Small mass cluster pairs, m02 of cluster 1 vs cluster 2 
  
  TH2F * fhConvDeltaEtaMCAntiNeutron;     //!<! Small mass cluster pairs, correlation in eta, origin of both clusters is anti neutron
  TH2F * fhConvDeltaPhiMCAntiNeutron;     //!<! Small mass cluster pairs, correlation in phi, origin of both clusters is anti neutron
  TH2F * fhConvDeltaEtaPhiMCAntiNeutron;  //!<! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti neutron
  TH2F * fhConvAsymMCAntiNeutron;         //!<! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti neutron
  TH2F * fhConvPtMCAntiNeutron;           //!<! Small mass cluster pairs, pt of pair, origin of both clusters is anti neutron
//TH2F * fhConvDispersionMCAntiNeutron;   //!<! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti neutron
  TH2F * fhConvM02MCAntiNeutron;          //!<! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti neutron

  TH2F * fhConvDeltaEtaMCAntiProton;      //!<! Small mass cluster pairs, correlation in eta, origin of both clusters is anti proton
  TH2F * fhConvDeltaPhiMCAntiProton;      //!<! Small mass cluster pairs, correlation in phi, origin of both clusters is anti proton
  TH2F * fhConvDeltaEtaPhiMCAntiProton;   //!<! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti proton
  TH2F * fhConvAsymMCAntiProton;          //!<! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti proton
  TH2F * fhConvPtMCAntiProton;            //!<! Small mass cluster pairs, pt of pairs, origin of both clusters is anti proton
//TH2F * fhConvDispersionMCAntiProton;    //!<! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti proton
  TH2F * fhConvM02MCAntiProton;           //!<! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti proton

  TH2F * fhConvDeltaEtaMCString;          //!<! Small mass cluster pairs, correlation in eta, origin of both clusters is string
  TH2F * fhConvDeltaPhiMCString;          //!<! Small mass cluster pairs, correlation in phi, origin of both clusters is string
  TH2F * fhConvDeltaEtaPhiMCString;       //!<! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is string
  TH2F * fhConvAsymMCString;              //!<! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is string
  TH2F * fhConvPtMCString;                //!<! Small mass cluster pairs, pt of pairs, origin of both clusters is string
//TH2F * fhConvDispersionMCString;        //!<! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvM02MCString;               //!<! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvDistMCConversion;          //!<! Calculated conversion distance vs real distance to vertex       
  TH2F * fhConvDistMCConversionCuts;      //!<! Calculated conversion distance vs real distance to vertex        

  /// Copy constructor not implemented.
  AliAnaPhotonConvInCalo(              const AliAnaPhotonConvInCalo & g) ;
    
  /// Assignment operator not implemented.
  AliAnaPhotonConvInCalo & operator = (const AliAnaPhotonConvInCalo & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPhotonConvInCalo,3) ;
  /// \endcond

} ;
 
#endif//ALIANAPHOTONINCALO_H



