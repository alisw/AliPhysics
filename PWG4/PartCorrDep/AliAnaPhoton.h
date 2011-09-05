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
    
  
  // Analysis methods
  
  Bool_t       ClusterSelected(AliVCluster* cl, TLorentzVector mom) ;
  
  void         FillAcceptanceHistograms();

  //           Fill Shower Shape histograms
  void         FillShowerShapeHistograms( AliVCluster* cluster, const Int_t mcTag) ;
  
  void         SwitchOnFillShowerShapeHistograms()    { fFillSSHistograms = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()   { fFillSSHistograms = kFALSE ; }  
  
  
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------
  
  TString      GetCalorimeter()                 const { return fCalorimeter        ; }
  void         SetCalorimeter(TString  & det)         { fCalorimeter = det         ; }
    
  // ** Cluster selection methods **
  
  void         SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3; }

  void         SetTimeCut(Double_t min, Double_t max) { fTimeCutMin = min; 
                                                        fTimeCutMax = max          ; }
  Double_t     GetTimeCutMin()                  const { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()                  const { return fTimeCutMax         ; }	
	
  void         SetNCellCut(Int_t n)                   { fNCellsCut = n             ; }
  Double_t     GetNCellCut()                    const { return fNCellsCut          ; }
  
  Bool_t       IsTrackMatchRejectionOn()        const { return fRejectTrackMatch   ; }
  void         SwitchOnTrackMatchRejection()          { fRejectTrackMatch = kTRUE  ; }
  void         SwitchOffTrackMatchRejection()         { fRejectTrackMatch = kFALSE ; }  
  
  // ** Conversion pair analysis **
  
  Float_t      GetMassCut()                     const { return fMassCut            ; }
  void         SetMassCut(Float_t m)                  { fMassCut    = m            ; }
  
  Bool_t       IsCheckConversionOn()            const { return fCheckConversion    ; }
  void         SwitchOnConversionChecker()            { fCheckConversion = kTRUE   ; }
  void         SwitchOffConversionChecker()           { fCheckConversion = kFALSE  ; }  
	
  Bool_t       AreConvertedPairsInAOD()         const { return fAddConvertedPairsToAOD   ; }
  void         SwitchOnAdditionConvertedPairsToAOD()  { fAddConvertedPairsToAOD = kTRUE  ; 
                                                        fCheckConversion        = kTRUE  ; }
  void         SwitchOffAdditionConvertedPairsToAOD() { fAddConvertedPairsToAOD = kFALSE ; }  
	
  Bool_t       AreConvertedPairsRemoved()       const { return fRemoveConvertedPair      ; }
  void         SwitchOnConvertedPairsRemoval()        { fRemoveConvertedPair  = kTRUE    ; 
                                                        fCheckConversion      = kTRUE    ; }
  void         SwitchOffConvertedPairsRemoval()       { fRemoveConvertedPair  = kFALSE   ; }    
  
  void         SetConvAsymCut(Float_t c)              { fConvAsymCut = c           ; }
  Float_t      GetConvAsymCut()                 const { return fConvAsymCut        ; }
  
  void         SetConvDEtaCut(Float_t c)              { fConvDEtaCut = c           ; }
  Float_t      GetConvDEtaCut()                 const { return fConvDEtaCut        ; }
  
  void         SetConvDPhiCut(Float_t min, Float_t max)  { fConvDPhiMinCut = min   ;  
                                                           fConvDPhiMaxCut = max   ; }
  Float_t      GetConvDPhiMinCut()              const { return fConvDPhiMinCut     ; }
  Float_t      GetConvDPhiMaxCut()              const { return fConvDPhiMaxCut     ; }
  
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
    if(n > 14) fNOriginHistograms = 14; }
  void         FillNPrimaryHistograms(Int_t n)        { fNPrimaryHistograms= n ;
    if(n > 7)  fNPrimaryHistograms = 7; }

  // For histograms in arrays, index in the array, corresponding to a particle
  enum mcTypes    { mcPhoton = 0,        mcPi0Decay = 1,       mcOtherDecay = 2,  
                    mcPi0 = 3,           mcEta = 4,            mcElectron = 5,       
                    mcConversion = 6,    mcOther = 7,          mcAntiNeutron = 8,    
                    mcAntiProton = 9,    mcPrompt = 10,        mcFragmentation = 11, 
                    mcISR = 12,          mcString = 13                               };  

  enum mcPTypes   { mcPPhoton = 0,       mcPPi0Decay = 1,       mcPOtherDecay = 2,  mcPOther = 3,
                    mcPPrompt = 4,       mcPFragmentation = 5,  mcPISR = 6           };  
  
  enum mcssTypes  { mcssPhoton = 0,      mcssOther = 1,       mcssPi0 = 2,         
                    mcssEta = 3,         mcssConversion = 4,  mcssElectron = 5       };  
  
  private:
 
  TString  fCalorimeter ;                // Calorimeter where the gamma is searched;
  Float_t  fMinDist ;                    // Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                    // Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                    // One more cut on distance used for acceptance-efficiency study
  Bool_t   fRejectTrackMatch ;           // If PID on, reject clusters which have an associated TPC track
  Double_t fTimeCutMin  ;                // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                // Remove clusters/cells with time larger than this value, in ns
  Int_t    fNCellsCut ;                  // Accept for the analysis clusters with more than fNCellsCut cells
  Bool_t   fFillSSHistograms ;           // Fill shower shape histograms
  Int_t    fNOriginHistograms;           // Fill only NOriginHistograms of the 14 defined types
  Int_t    fNPrimaryHistograms;          // Fill only NPrimaryHistograms of the 7 defined types

  //Conversion pairs selection cuts
  Bool_t   fCheckConversion;             // Combine pairs of clusters with mass close to 0
  Bool_t   fRemoveConvertedPair;         // Remove conversion pairs
  Bool_t   fAddConvertedPairsToAOD;      // Put Converted pairs in AOD
  Float_t  fMassCut;                     // Mass cut for the conversion pairs selection  
  Float_t  fConvAsymCut;                 // Select conversion pairs when asymmetry is smaller than cut
  Float_t  fConvDEtaCut;                 // Select conversion pairs when deta of pair smaller than cut
  Float_t  fConvDPhiMinCut;              // Select conversion pairs when dphi of pair lager than cut
  Float_t  fConvDPhiMaxCut;              // Select conversion pairs when dphi of pair smaller than cut

  //Histograms 
  TH2F * fhNCellsE;                      //! number of cells in cluster vs E 
  TH2F * fhMaxCellDiffClusterE;          //! Fraction of energy carried by cell with maximum energy
  
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

  //Shower shape
  
  TH2F * fhDispE;                         //! cluster dispersion vs E
  TH2F * fhLam0E;                         //! cluster lambda0 vs  E
  TH2F * fhLam1E;                         //! cluster lambda1 vs  E  
  TH2F * fhdDispE;                        //! cluster dispersion/Ncells vs E
  TH2F * fhdLam0E;                        //! cluster lambda0/Ncells vs  E
  TH2F * fhdLam1E;                        //! cluster lambda1/Ncells vs  E
  
  TH2F * fhDispETRD;                      //! cluster dispersion vs E, SM covered by TRD
  TH2F * fhLam0ETRD;                      //! cluster lambda0 vs  E, SM covered by TRD
  TH2F * fhLam1ETRD;                      //! cluster lambda1 vs  E, SM covered by TRD 
  TH2F * fhdDispETRD;                     //! cluster dispersion/Ncells vs E, SM covered by TRD
  TH2F * fhdLam0ETRD;                     //! cluster lambda0/Ncells vs  E, SM covered by TRD
  TH2F * fhdLam1ETRD;                     //! cluster lambda1/Ncells vs  E, SM covered by TRD    
  
  TH2F * fhNCellsLam0LowE;                //! number of cells in cluster vs lambda0
  TH2F * fhNCellsLam1LowE;                //! number of cells in cluster vs lambda1
  TH2F * fhNCellsDispLowE;                //! number of cells in cluster vs dispersion
  TH2F * fhNCellsLam0HighE;               //! number of cells in cluster vs lambda0, E>2
  TH2F * fhNCellsLam1HighE;               //! number of cells in cluster vs lambda1, E>2
  TH2F * fhNCellsDispHighE;               //! number of cells in cluster vs dispersion, E>2
  
  TH2F * fhNCellsdLam0LowE;               //! number of cells in cluster vs lambda0/ncells
  TH2F * fhNCellsdLam1LowE;               //! number of cells in cluster vs lambda1/ncells
  TH2F * fhNCellsdDispLowE;               //! number of cells in cluster vs dispersion/ncells
  TH2F * fhNCellsdLam0HighE;              //! number of cells in cluster vs lambda0/ncells, E>2
  TH2F * fhNCellsdLam1HighE;              //! number of cells in cluster vs lambda1/ncells, E>2
  TH2F * fhNCellsdDispHighE;              //! number of cells in cluster vs dispersion/ncells, E>2  
  
  TH2F * fhEtaLam0LowE;                   //! cluster eta vs lambda0, E<2
  TH2F * fhPhiLam0LowE;                   //! cluster phi vs lambda0, E<2
  TH2F * fhEtaLam0HighE;                  //! cluster eta vs lambda0, E>2
  TH2F * fhPhiLam0HighE;                  //! cluster phi vs lambda0, E>2
  TH2F * fhLam0DispLowE;                  //! cluster lambda0 vs dispersion, E<2
  TH2F * fhLam0DispHighE;                 //! cluster lambda0 vs dispersion, E>2
  TH2F * fhLam1Lam0LowE;                  //! cluster lambda1 vs lambda0, E<2
  TH2F * fhLam1Lam0HighE;                 //! cluster lambda1 vs lambda0, E>2
  TH2F * fhDispLam1LowE;                  //! cluster disp vs lambda1, E<2
  TH2F * fhDispLam1HighE;                 //! cluster disp vs lambda1, E>2
  
  TH2F * fhEtadLam0LowE;                  //! cluster eta vs lambda0/ncells, E<2
  TH2F * fhPhidLam0LowE;                  //! cluster phi vs lambda0/ncells, E<2
  TH2F * fhEtadLam0HighE;                 //! cluster eta vs lambda0/ncells, E>2
  TH2F * fhPhidLam0HighE;                 //! cluster phi vs lambda0/ncells, E>2
  TH2F * fhdLam0dDispLowE;                //! cluster lambda0/ncells vs dispersion/ncells, E<2
  TH2F * fhdLam0dDispHighE;               //! cluster lambda0/ncells vs dispersion/ncells, E>2
  TH2F * fhdLam1dLam0LowE;                //! cluster lambda1/ncells vs lambda0/ncells, E<2
  TH2F * fhdLam1dLam0HighE;               //! cluster lambda1/ncells vs lambda0/ncells, E>2
  TH2F * fhdDispdLam1LowE;                //! cluster disp/ncells vs lambda1/ncells, E<2
  TH2F * fhdDispdLam1HighE;               //! cluster disp/ncells vs lambda1/ncells, E>2
  
  
  //Fill MC dependent histograms
  TH1F * fhDeltaE  ;                          //! MC-Reco E distribution      
  TH1F * fhDeltaPt ;                          //! MC-Reco pT distribution
  TH1F * fhRatioE  ;                          //! Reco/MC E distribution      
  TH1F * fhRatioPt ;                          //! Reco/MC pT distribution
  TH2F * fh2E  ;                              //! E distribution, Reco vs MC
  TH2F * fh2Pt ;                              //! pT distribution, Reco vs MC
  
  //Origin of this cluster is ...
  TH1F * fhMCE[14];                           //! Number of identified photon vs cluster energy coming from MC particle
  TH1F * fhPtMC[14];                          //! Number of identified photon vs cluster pT     coming from MC particle
  TH2F * fhPhiMC[14];                         //! Phi of identified photon coming from MC particle
  TH2F * fhEtaMC[14];                         //! eta of identified photon coming from MC particle

  TH1F * fhEPrimMC[7];                        //! Number of generated photon vs energy
  TH1F * fhPtPrimMC[7];                       //! Number of generated photon vs pT   
  TH2F * fhPhiPrimMC[7];                      //! Phi of generted photon
  TH2F * fhYPrimMC[7];                        //! Rapidity of generated photon 
  
  TH1F * fhEPrimMCAcc[7];                     //! Number of generated photon vs energy, in calorimeter acceptance
  TH1F * fhPtPrimMCAcc[7];                    //! Number of generated photon vs pT, in calorimeter acceptance   
  TH2F * fhPhiPrimMCAcc[7];                   //! Phi of generted photon, in calorimeter acceptance
  TH2F * fhYPrimMCAcc[7];                     //! Rapidity of generated photon, in calorimeter acceptance   
  
  //Conversion pairs analysis histograms
  TH1F * fhPtConversionTagged;                //! Number of identified gamma from Conversion , tagged as conversion 
  TH1F * fhPtAntiNeutronTagged;               //! Number of identified gamma from AntiNeutrons gamma, tagged as conversion 
  TH1F * fhPtAntiProtonTagged;                //! Number of identified gamma from AntiProtons gamma, tagged as conversion 
  TH1F * fhPtUnknownTagged;                   //! Number of identified gamma from unknown, tagged as conversion 
  
  TH2F * fhEtaPhiConversion  ;                //! Pseudorapidity vs Phi for transerse momentum > 0.5, for MC converted
  TH2F * fhEtaPhi05Conversion  ;              //! Pseudorapidity vs Phi for transerse momentum < 0.5, for MC converted
  
  TH2F * fhConvDeltaEtaMCConversion;          //! Small mass cluster pairs, correlation in eta, origin of both clusters is conversion
  TH2F * fhConvDeltaPhiMCConversion;          //! Small mass cluster pairs, correlation in phi, origin of both clusters is conversion
  TH2F * fhConvDeltaEtaPhiMCConversion;       //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is conversion
  TH2F * fhConvAsymMCConversion;              //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is conversion
  TH2F * fhConvPtMCConversion;                //! Small mass cluster pairs, pt of pair, origin of both clusters is conversion
  TH2F * fhConvDispersionMCConversion;        //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2
  TH2F * fhConvM02MCConversion;               //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2 

  TH2F * fhConvDeltaEtaMCAntiNeutron;         //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti neutron
  TH2F * fhConvDeltaPhiMCAntiNeutron;         //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti neutron
  TH2F * fhConvDeltaEtaPhiMCAntiNeutron;      //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti neutron
  TH2F * fhConvAsymMCAntiNeutron;             //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti neutron
  TH2F * fhConvPtMCAntiNeutron;               //! Small mass cluster pairs, pt of pair, origin of both clusters is anti neutron
  TH2F * fhConvDispersionMCAntiNeutron;       //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti neutron
  TH2F * fhConvM02MCAntiNeutron;              //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti neutron

  TH2F * fhConvDeltaEtaMCAntiProton;          //! Small mass cluster pairs, correlation in eta, origin of both clusters is anti proton
  TH2F * fhConvDeltaPhiMCAntiProton;          //! Small mass cluster pairs, correlation in phi, origin of both clusters is anti proton
  TH2F * fhConvDeltaEtaPhiMCAntiProton;       //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is anti proton
  TH2F * fhConvAsymMCAntiProton;              //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is anti proton
  TH2F * fhConvPtMCAntiProton;                //! Small mass cluster pairs, pt of pairs, origin of both clusters is anti proton
  TH2F * fhConvDispersionMCAntiProton;        //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is anti proton
  TH2F * fhConvM02MCAntiProton;               //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is anti proton

  TH2F * fhConvDeltaEtaMCString;              //! Small mass cluster pairs, correlation in eta, origin of both clusters is string
  TH2F * fhConvDeltaPhiMCString;              //! Small mass cluster pairs, correlation in phi, origin of both clusters is string
  TH2F * fhConvDeltaEtaPhiMCString;           //! Small mass cluster pairs, correlation in eta-phi, origin of both clusters is string
  TH2F * fhConvAsymMCString;                  //! Small mass cluster pairs, correlation in energy asymmetry, origin of both clusters is string
  TH2F * fhConvPtMCString;                    //! Small mass cluster pairs, pt of pairs, origin of both clusters is string
  TH2F * fhConvDispersionMCString;            //! Small mass cluster pairs, dispersion of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvM02MCString;                   //! Small mass cluster pairs, m02 of cluster 1 vs cluster 2, origin of both clusters is string
  TH2F * fhConvDistMCConversion;              //! Calculated conversion distance vs real distance to vertex       
  TH2F * fhConvDistMCConversionCuts;          //! Calculated conversion distance vs real distance to vertex        

  // Shower Shape MC

  TH2F * fhMCELambda0[6] ;                    //! E vs Lambda0     from MC particle
  TH2F * fhMCEdLambda0[6];                    //! E vs dLambda0    from MC particle
  TH2F * fhMCELambda1[6] ;                    //! E vs Lambda1     from MC particle
  TH2F * fhMCEdLambda1[6];                    //! E vs dLambda1    from MC particle
  TH2F * fhMCEDispersion[6] ;                 //! E vs Dispersion  from MC particle
  TH2F * fhMCEdDispersion[6];                 //! E vs dDispersion from MC particle
  
  TH2F * fhMCPhotonELambda0NoOverlap ;        //! E vs Lambda0     from MC photons, no overlap
  TH2F * fhMCPhotonELambda0TwoOverlap ;       //! E vs Lambda0     from MC photons, 2 particles overlap
  TH2F * fhMCPhotonELambda0NOverlap ;         //! E vs Lambda0     from MC photons, N particles overlap
  TH2F * fhMCPhotonEdLambda0NoOverlap ;       //! E vs dLambda0    from MC photons, no overlap
  TH2F * fhMCPhotonEdLambda0TwoOverlap ;      //! E vs dLambda0    from MC photons, 2 particles overlap
  TH2F * fhMCPhotonEdLambda0NOverlap ;        //! E vs dLambda0    from MC photons, N particles overlap
  
  TH2F * fhMCLambda0vsClusterMaxCellDiffE0[6];  //! Lambda0 vs fraction of energy of max cell for E < 2 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE2[6];  //! Lambda0 vs fraction of energy of max cell for 2< E < 6 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE6[6];  //! Lambda0 vs fraction of energy of max cell for E > 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE0[6];   //! NCells  vs fraction of energy of max cell for E < 2
  TH2F * fhMCNCellsvsClusterMaxCellDiffE2[6];   //! NCells  vs fraction of energy of max cell for 2 < E < 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE6[6];   //! NCells  vs fraction of energy of max cell for E > 6
  TH2F * fhMCNCellsE[6];                        //! NCells per cluster vs energy
  TH2F * fhMCMaxCellDiffClusterE[6];            //! Fraction of energy carried by cell with maximum energy

  //Embedding
  TH2F * fhEmbeddedSignalFractionEnergy ;     //! Fraction of photon energy of embedded signal vs cluster energy
  
  TH2F * fhEmbedPhotonELambda0FullSignal ;    //!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPhotonEdLambda0FullSignal ;   //! dLambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPhotonELambda0MostlySignal ;  //!  Lambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPhotonEdLambda0MostlySignal ; //! dLambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPhotonELambda0MostlyBkg ;     //!  Lambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPhotonEdLambda0MostlyBkg ;    //! dLambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPhotonELambda0FullBkg ;       //!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  TH2F * fhEmbedPhotonEdLambda0FullBkg ;      //! dLambda0 vs E for embedded photons with less than 10% of the cluster energy
  
  TH2F * fhEmbedPi0ELambda0FullSignal ;       //!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPi0EdLambda0FullSignal ;      //! dLambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPi0ELambda0MostlySignal ;     //!  Lambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPi0EdLambda0MostlySignal ;    //! dLambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPi0ELambda0MostlyBkg ;        //!  Lambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPi0EdLambda0MostlyBkg ;       //! dLambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPi0ELambda0FullBkg ;          //!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  TH2F * fhEmbedPi0EdLambda0FullBkg ;         //! dLambda0 vs E for embedded photons with less than 10% of the cluster energy  
  
   ClassDef(AliAnaPhoton,16)

} ;
 

#endif//ALIANAPHOTON_H



