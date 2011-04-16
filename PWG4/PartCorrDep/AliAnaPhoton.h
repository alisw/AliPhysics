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
// Produces input for other analysis classes like AliAnaPi0, 
// AliAnaParticleHadronCorrelation ... 
//

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F ;
class TH1F;
class TH3D;
class TString ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"
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
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------

  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();

  void         Init();

  void         InitParameters();

  void         MakeAnalysisFillAOD()  ;

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const;
  
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------

  TString GetCalorimeter()                  const {return fCalorimeter ; }
  void    SetCalorimeter(TString  & det)          {fCalorimeter = det  ; }
  
  // ** Cluster selection methods **
  
  void    SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
    fMinDist = m1;
    fMinDist2 = m2;
    fMinDist3 = m3;
  }

  void     SetTimeCut(Double_t min, Double_t max) {fTimeCutMin = min; fTimeCutMax = max ; }
  Double_t GetTimeCutMin()                  const {return fTimeCutMin ; }
  Double_t GetTimeCutMax()                  const {return fTimeCutMax ; }	
	
  void     SetNCellCut(Int_t n)                   {fNCellsCut = n    ; }
  Double_t GetNCellCut()                    const {return fNCellsCut ; }
  
  Bool_t   IsTrackMatchRejectionOn()        const {return fRejectTrackMatch   ; }
  void     SwitchOnTrackMatchRejection()          {fRejectTrackMatch = kTRUE  ; }
  void     SwitchOffTrackMatchRejection()         {fRejectTrackMatch = kFALSE ; }  

  // ** Conversion pair analysis **
  
  Float_t  GetMassCut()                     const { return fMassCut           ; }
  void     SetMassCut(Float_t m)                  { fMassCut    = m           ; }
  
  Bool_t   IsCheckConversionOn()            const { return fCheckConversion   ; }
  void     SwitchOnConversionChecker()            { fCheckConversion = kTRUE  ; }
  void     SwitchOffConversionChecker()           { fCheckConversion = kFALSE ; }  
	
  Bool_t   AreConvertedPairsInAOD()         const { return fAddConvertedPairsToAOD   ; }
  void     SwitchOnAdditionConvertedPairsToAOD()  { fAddConvertedPairsToAOD = kTRUE  ; fCheckConversion = kTRUE ; }
  void     SwitchOffAdditionConvertedPairsToAOD() { fAddConvertedPairsToAOD = kFALSE ; }  
	
  Bool_t   AreConvertedPairsRemoved()       const { return fRemoveConvertedPair      ; }
  void     SwitchOnConvertedPairsRemoval()        { fRemoveConvertedPair  = kTRUE    ; fCheckConversion = kTRUE ; }
  void     SwitchOffConvertedPairsRemoval()       { fRemoveConvertedPair  = kFALSE   ; }    
  
  void     SetConvAsymCut(Float_t c)              { fConvAsymCut = c    ; }
  Float_t  GetConvAsymCut()                 const { return fConvAsymCut ; }
  
  void     SetConvDEtaCut(Float_t c)              { fConvDEtaCut = c    ; }
  Float_t  GetConvDEtaCut()                 const { return fConvDEtaCut ; }
  
  void     SetConvDPhiCut(Float_t min, Float_t max)  { fConvDPhiMinCut = min ;  
                                                       fConvDPhiMaxCut = max ; }
  Float_t  GetConvDPhiMinCut()              const { return fConvDPhiMinCut ; }
  Float_t  GetConvDPhiMaxCut()              const { return fConvDPhiMaxCut ; }

  private:
 
  TString  fCalorimeter ;                // Calorimeter where the gamma is searched;
  Float_t  fMinDist ;                    // Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                    // Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                    // One more cut on distance used for acceptance-efficiency study
  Bool_t   fRejectTrackMatch ;           // If PID on, reject clusters which have an associated TPC track
  Double_t fTimeCutMin  ;                // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                // Remove clusters/cells with time larger than this value, in ns
  Int_t    fNCellsCut ;                  // Accept for the analysis clusters with more than fNCellsCut cells
  
  //Conversion pairs selection cuts
  Bool_t   fCheckConversion;             // Combine pairs of clusters with mass close to 0
  Bool_t   fRemoveConvertedPair;         // Combine pairs of clusters with mass close to 0
  Bool_t   fAddConvertedPairsToAOD;      // Put Converted pairs in AOD
  Float_t  fMassCut;                     // Mass cut for the conversion pairs selection  
  Float_t  fConvAsymCut;                 // Select conversion pairs when asymmetry is smaller than cut
	Float_t  fConvDEtaCut;                 // Select conversion pairs when deta of pair smaller than cut
  Float_t  fConvDPhiMinCut;              // Select conversion pairs when dphi of pair lager than cut
  Float_t  fConvDPhiMaxCut;              // Select conversion pairs when dphi of pair smaller than cut

  //Histograms 
  TH2F * fhNtraNclu;                     //! track multiplicity distribution vs cluster multiplicity
  TH2F * fhNCellsPt;                     //! number of cells in cluster vs pt 
  TH1F * fhEPhoton    ;                  //! Number of identified photon vs energy
  TH1F * fhPtPhoton   ;                  //! Number of identified photon vs transerse momentum 
  TH2F * fhPhiPhoton  ;                  //! Azimuthal angle of identified  photon vs transerse momentum 
  TH2F * fhEtaPhoton  ;                  //! Pseudorapidity of identified  photon vs transerse momentum 
  TH2F * fhEtaPhiPhoton  ;               //! Pseudorapidity vs Phi of identified  photon for transerse momentum > 0.5
  TH2F * fhEtaPhi05Photon  ;             //! Pseudorapidity vs Phi of identified  photon for transerse momentum < 0.5
  
  //Conversion pairs
  TH1F * fhPtPhotonConv   ;              //! Number of identified photon vs transerse momentum 
  TH2F * fhEtaPhiPhotonConv  ;           //! Pseudorapidity vs Phi of identified  photon for transerse momentum > 0.5, for converted
  TH2F * fhEtaPhi05PhotonConv  ;         //! Pseudorapidity vs Phi of identified  photon for transerse momentum < 0.5, for converted
  TH2F * fhConvDeltaEta;                 //! Small mass photons, correlation in eta
  TH2F * fhConvDeltaPhi;                 //! Small mass photons, correlation in phi
  TH2F * fhConvDeltaEtaPhi;              //! Small mass photons, correlation in phi and eta
  TH2F * fhConvAsym;                     //! Small mass photons, correlation in energy asymmetry
  TH2F * fhConvPt;                       //! Small mass photons, pT of pair

  //Fill MC dependent histograms
  TH1F * fhDeltaE  ;                     //! MC-Reco E distribution      
  TH1F * fhDeltaPt ;                     //! MC-Reco pT distribution
  TH1F * fhRatioE  ;                     //! Reco/MC E distribution      
  TH1F * fhRatioPt ;                     //! Reco/MC pT distribution
  TH2F * fh2E  ;                         //! E distribution, Reco vs MC
  TH2F * fh2Pt ;                         //! pT distribution, Reco vs MC
  
  //Origin of this cluster is ...
  TH1F * fhPtMCPhoton;                   //! Number of identified gamma 
  TH2F * fhPhiMCPhoton;                  //! Phi of identified gamma
  TH2F * fhEtaMCPhoton;                  //! eta of identified gamma	
	
  TH1F * fhPtPrompt;                     //! Number of identified prompt gamma 
  TH2F * fhPhiPrompt;                    //! Phi of identified  prompt gamma
  TH2F * fhEtaPrompt;                    //! eta of identified  prompt gamma

  TH1F * fhPtFragmentation;              //! Number of identified fragmentation gamma 
  TH2F * fhPhiFragmentation;             //! Phi of identified  fragmentation gamma
  TH2F * fhEtaFragmentation;             //! eta of identified  fragmentation gamma

  TH1F * fhPtISR;                        //! Number of identified initial state radiation gamma 
  TH2F * fhPhiISR;                       //! Phi of identified initial state radiation gamma
  TH2F * fhEtaISR;                       //! eta of identified initial state radiation gamma

  TH1F * fhPtPi0Decay;                   //! Number of identified Pi0Decay gamma 
  TH2F * fhPhiPi0Decay;                  //! Phi of identified  Pi0Decay gamma
  TH2F * fhEtaPi0Decay;                  //! eta of identified  Pi0Decay gamma

  TH1F * fhPtOtherDecay;                 //! Number of identified OtherDecay gamma 
  TH2F * fhPhiOtherDecay;                //! Phi of identified  OtherDecay gamma
  TH2F * fhEtaOtherDecay;                //! eta of identified  OtherDecay gamma

  TH1F * fhPtConversion;                 //! Number of identified Conversion gamma 
  TH2F * fhPhiConversion;                //! Phi of identified  Conversion gamma
  TH2F * fhEtaConversion;                //! eta of identified  Conversion gamma

  TH2F * fhEtaPhiConversion  ;           //! Pseudorapidity vs Phi of identified  photon for transerse momentum > 0.5
  TH2F * fhEtaPhi05Conversion  ;         //! Pseudorapidity vs Phi of identified  photon for transerse momentum < 0.5

  TH1F * fhPtAntiNeutron;                //! Origin are AntiNeutrons or AntiProtons
  TH2F * fhPhiAntiNeutron;               //! Origin are AntiNeutrons or AntiProtons
  TH2F * fhEtaAntiNeutron;               //! Origin are AntiNeutrons or AntiProtons
    
  TH1F * fhPtAntiProton;                 //! Origin are AntiNeutrons or AntiProtons
  TH2F * fhPhiAntiProton;                //! Origin are AntiNeutrons or AntiProtons
  TH2F * fhEtaAntiProton;                //! Origin are AntiNeutrons or AntiProtons
  
  TH1F * fhPtUnknown;                    //! Number of identified Unknown gamma 
  TH2F * fhPhiUnknown;                   //! Phi of identified  Unknown gamma
  TH2F * fhEtaUnknown;                   //! eta of identified  Unknown gamma
  
  //Conversion pairs analysis histograms
  TH1F * fhPtConversionTagged;           //! Number of identified gamma from Conversion , tagged as conversion 
  TH1F * fhPtAntiNeutronTagged;          //! Number of identified gamma from AntiNeutrons gamma, tagged as conversion 
  TH1F * fhPtAntiProtonTagged;           //! Number of identified gamma from AntiProtons gamma, tagged as conversion 
  TH1F * fhPtUnknownTagged;              //! Number of identified gamma from unknown, tagged as conversion 
  
  TH2F * fhConvDeltaEtaMCConversion;     //! Small mass cluster pairs, correlation in eta, origin of both clusters is conversion
  TH2F * fhConvDeltaPhiMCConversion;     //! Small mass cluster pairs, correlation in phi, origin of both clusters is conversion
  TH2F * fhConvDeltaEtaPhiMCConversion;  //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is conversion
  TH2F * fhConvAsymMCConversion;         //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is conversion
  TH2F * fhConvPtMCConversion;           //! Small mass cluster pairs, pt of pair, origin of both clusters is conversion
  TH2F * fhConvDispersionMCConversion;   //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2
  TH2F * fhConvM02MCConversion;          //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2

  TH2F * fhConvDeltaEtaMCAntiNeutron;    //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti neutron
  TH2F * fhConvDeltaPhiMCAntiNeutron;    //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti neutron
  TH2F * fhConvDeltaEtaPhiMCAntiNeutron; //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti neutron
  TH2F * fhConvAsymMCAntiNeutron;        //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti neutron
  TH2F * fhConvPtMCAntiNeutron;          //! Small mass cluster pairs, pt of pair, origin of both clusters is anti neutron
  TH2F * fhConvDispersionMCAntiNeutron;  //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti neutron
  TH2F * fhConvM02MCAntiNeutron;         //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti neutron

  TH2F * fhConvDeltaEtaMCAntiProton;     //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti proton
  TH2F * fhConvDeltaPhiMCAntiProton;     //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti proton
  TH2F * fhConvDeltaEtaPhiMCAntiProton;  //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti proton
  TH2F * fhConvAsymMCAntiProton;         //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti proton
  TH2F * fhConvPtMCAntiProton;           //! Small mass cluster pairs, pt of pairs, origin of both clusters is anti proton
  TH2F * fhConvDispersionMCAntiProton;   //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti proton
  TH2F * fhConvM02MCAntiProton;          //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti proton

  TH2F * fhConvDeltaEtaMCString;         //! Small mass cluster pairs, correlation in eta, origin of both clusters is string
  TH2F * fhConvDeltaPhiMCString;         //! Small mass cluster pairs, correlation in phi, origin of both clusters is string
  TH2F * fhConvDeltaEtaPhiMCString;      //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is string
  TH2F * fhConvAsymMCString;             //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is string
  TH2F * fhConvPtMCString;               //! Small mass cluster pairs, pt of pairs, origin of both clusters is string
  TH2F * fhConvDispersionMCString;       //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvM02MCString;              //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is string

   ClassDef(AliAnaPhoton,11)

} ;
 

#endif//ALIANAPHOTON_H



