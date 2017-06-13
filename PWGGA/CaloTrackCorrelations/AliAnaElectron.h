#ifndef ALIANAELECTRON_H
#define ALIANAELECTRON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaElectron
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Selection of electron clusters in calorimeter.
///
/// Class for the electron identification,
/// Clusters from calorimeters are identified as electrons
/// and kept in the AOD, possibility to select hadronic clusters.
/// Few histograms produced.
/// Copy of AliAnaPhoton just add electron id.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaElectron).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
class TH2F ;
class TH1F;
class TH3D;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class TList ;

class AliAnaElectron : public AliAnaCaloTrackCorrBaseClass {

 public: 
  
  AliAnaElectron() ;
  
  /// Virtual destructor.
  virtual ~AliAnaElectron() { ; }
	
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
  
  Bool_t       ClusterSelected(AliVCluster* cl, Int_t nMaxima) ;
  
  void         FillShowerShapeHistograms( AliVCluster* cluster, Int_t mcTag , Int_t pidTag) ;
  
  void         SwitchOnFillShowerShapeHistograms()    { fFillSSHistograms = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()   { fFillSSHistograms = kFALSE ; }  
    
  void         WeightHistograms(AliVCluster *clus);
  
  void         SwitchOnFillWeightHistograms()         { fFillWeightHistograms = kTRUE  ; }
  void         SwitchOffFillWeightHistograms()        { fFillWeightHistograms = kFALSE ; }  
  
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------
  
  // ** Cluster selection methods **
  
  void         SetdEdxCut(Double_t min, Double_t max) { fdEdxMin    = min ; 
                                                        fdEdxMax    = max          ; }
  
  void         SetEOverP(Double_t min, Double_t max)  { fEOverPMin  = min ; 
                                                        fEOverPMax  = max          ; }

  
  void         SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3; }

  void         SetTimeCut(Double_t min, Double_t max) { fTimeCutMin = min; 
                                                        fTimeCutMax = max          ; }
  Double_t     GetTimeCutMin()                  const { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()                  const { return fTimeCutMax         ; }	
	
  void         SetNCellCut(Int_t n)                   { fNCellsCut = n             ; }
  Double_t     GetNCellCut()                    const { return fNCellsCut          ; }
  
  void         SetNLMCut(Int_t min, Int_t max)        { fNLMCutMin = min;
                                                        fNLMCutMax = max                ; }
  Int_t        GetNLMCutMin()                   const { return fNLMCutMin               ; }
  Int_t        GetNLMCutMax()                   const { return fNLMCutMax               ; }
  
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
                                                        if(n > 10) fNOriginHistograms = 10; }

  
  void         FillAODWithElectrons()                 { fAODParticle = AliCaloPID::kElectron      ; }
  void         FillAODWithHadrons()                   { fAODParticle = AliCaloPID::kChargedHadron ; }
  void         FillAODWithAny()                       { fAODParticle = 0 ; }

  void         SwitchOnOnlySimpleSSHistoFill()        { fFillOnlySimpleSSHisto = kTRUE  ; }
  void         SwitchOffOnlySimpleHistoFill()         { fFillOnlySimpleSSHisto = kFALSE ; }
  
  /// For histograms in arrays, index in the array, corresponding to the originating particle of the cluster
  enum mcTypes    { kmcPhoton = 0,        kmcPi0Decay = 1,       kmcOtherDecay = 2,  
                    kmcPi0 = 3,           kmcEta = 4,            kmcElectron = 5,       
                    kmcConversion = 6,    kmcOther = 7,          kmcAntiNeutron = 8,    
                    kmcAntiProton = 9                                                 };    
  
  /// For histograms in arrays, dependent on shower shape, index in the array,
  /// corresponding to the originating particle of the cluster
  enum mcssTypes  { kmcssPhoton = 0,      kmcssOther = 1,       kmcssPi0 = 2,         
                    kmcssEta = 3,         kmcssConversion = 4,  kmcssElectron = 5       };  
  
  private:
 
  Float_t  fMinDist ;                           ///<  Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                           ///<  Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                           ///<  One more cut on distance used for acceptance-efficiency study
  Double_t fTimeCutMin  ;                       ///<  Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                       ///<  Remove clusters/cells with time larger than this value, in ns
  Int_t    fNCellsCut ;                         ///<  Accept for the analysis clusters with more than fNCellsCut cells
  Int_t    fNLMCutMin  ;                        ///<  Remove clusters/cells with number of local maxima smaller than this value
  Int_t    fNLMCutMax  ;                        ///<  Remove clusters/cells with number of local maxima larger than this value
  Bool_t   fFillSSHistograms ;                  ///<  Fill shower shape histograms
  Bool_t   fFillOnlySimpleSSHisto;              ///<  Fill selected cluster histograms, selected SS histograms
  Bool_t   fFillWeightHistograms ;              ///<  Fill weigth histograms
  Int_t    fNOriginHistograms;                  ///<  Fill only NOriginHistograms of the 14 defined types

  Float_t  fdEdxMin;                            ///<  Max dEdx for electrons
  Float_t  fdEdxMax;                            ///<  Min dEdx for electrons
  Float_t  fEOverPMin;                          ///<  Max E/p for electrons, after dEdx cut
  Float_t  fEOverPMax;                          ///<  Min E/p for electrons, after dEdx cut

  Int_t    fAODParticle;                        ///<  Select the type of particle to put in AODs for other analysis
  
  TLorentzVector fMomentum;                     //!<! cluster momentum
  TLorentzVector fMomentumMC;                   //!<! mc particle momentum
  TVector3       fProdVertex;                   //!<! mc particle production vertex

  //Histograms
  TH2F * fhdEdxvsE;                             //!<! Matched track dEdx vs cluster E
  TH2F * fhdEdxvsP;                             //!<! Matched track dEdx vs track P
  TH2F * fhEOverPvsE;                           //!<! Matched track E cluster over P track vs cluster E, after dEdx cut
  TH2F * fhEOverPvsP;                           //!<! Matched track E cluster over P track vs track P, after dEdx cut

  TH2F * fhdEdxvsECutM02;                       //!<! Matched track dEdx vs cluster E, mild M02 cut
  TH2F * fhdEdxvsPCutM02;                       //!<! Matched track dEdx vs track P, mild M02 cut
  TH2F * fhEOverPvsECutM02;                     //!<! Matched track E cluster over P track vs cluster E, after dEdx cut, mild M02 cut
  TH2F * fhEOverPvsPCutM02;                     //!<! Matched track E cluster over P track vs track P, after dEdx cut, mild M02 cut
  
  TH2F * fhdEdxvsECutEOverP;                    //!<! Matched track dEdx vs cluster E , cut on EOverP
  TH2F * fhdEdxvsPCutEOverP;                    //!<! Matched track dEdx vs track P, cut on EOverP
  TH2F * fhEOverPvsECutM02CutdEdx;              //!<! Matched track E cluster over P track vs cluster E, after dEdx cut and mild M02 cut
  TH2F * fhEOverPvsPCutM02CutdEdx;              //!<! Matched track E cluster over P track vs track P, after dEdx cut and mild M02 cut

  TH2F * fhMCdEdxvsE[10];                       //!<! Matched track dEdx vs cluster E, coming from MC particle
  TH2F * fhMCdEdxvsP[10];                       //!<! Matched track dEdx vs track P, coming from MC particle
  TH2F * fhMCEOverPvsE[10];                     //!<! Matched track E cluster over P track vs cluster E, after dEdx cut, coming from MC particle
  TH2F * fhMCEOverPvsP[10];                     //!<! Matched track E cluster over P track vs track P, after dEdx cut, coming from MC particle
  
  TH2F * fhNCellsE[2];                          //!<! Number of cells in cluster vs E
  TH2F * fhNLME[2];                             //!<! Number of local maxima in cluster vs E
  TH2F * fhMaxCellDiffClusterE[2];              //!<! Fraction of energy carried by cell with maximum energy
  TH2F * fhTimeE[2];                            //!<! E vs Time of selected cluster

  TH1F * fhE[2]    ;                            //!<! Number of identified electron vs energy
  TH1F * fhPt[2]   ;                            //!<! Number of identified electron vs transerse momentum
  TH2F * fhPhi[2]  ;                            //!<! Azimuthal angle of identified  electron vs transerse momentum
  TH2F * fhEta[2]  ;                            //!<! Pseudorapidity of identified  electron vs transerse momentum
  TH2F * fhEtaPhi[2]  ;                         //!<! Pseudorapidity vs Phi of identified  electron for transerse momentum > 0.5
  TH2F * fhEtaPhi05[2]  ;                       //!<! Pseudorapidity vs Phi of identified  electron for transerse momentum < 0.5
  
  //Shower shape
  
  TH2F * fhDispE[2];                            //!<! cluster dispersion vs E
  TH2F * fhLam0E[2];                            //!<! cluster lambda0 vs  E
  TH2F * fhLam1E[2];                            //!<! cluster lambda1 vs  E

  TH2F * fhDispETRD[2];                         //!<! cluster dispersion vs E, SM covered by TRD
  TH2F * fhLam0ETRD[2];                         //!<! cluster lambda0 vs  E, SM covered by TRD
  TH2F * fhLam1ETRD[2];                         //!<! cluster lambda1 vs  E, SM covered by TRD
  
  TH2F * fhNCellsLam0LowE[2];                   //!<! cluster N cells vs lambda0, E<2
  TH2F * fhNCellsLam0HighE[2];                  //!<! cluster N Cells vs lambda0, E>2

  TH2F * fhEtaLam0LowE[2];                      //!<! cluster eta vs lambda0, E<2
  TH2F * fhPhiLam0LowE[2];                      //!<! cluster phi vs lambda0, E<2
  TH2F * fhEtaLam0HighE[2];                     //!<! cluster eta vs lambda0, E>2
  TH2F * fhPhiLam0HighE[2];                     //!<! cluster phi vs lambda0, E>2
    
  TH2F * fhDispEtaE[2] ;                        //!<! shower dispersion in eta direction
  TH2F * fhDispPhiE[2] ;                        //!<! shower dispersion in phi direction
  TH2F * fhSumEtaE[2] ;                         //!<! shower dispersion in eta direction
  TH2F * fhSumPhiE[2] ;                         //!<! shower dispersion in phi direction
  TH2F * fhSumEtaPhiE[2] ;                      //!<! shower dispersion in eta and phi direction
  TH2F * fhDispEtaPhiDiffE[2] ;                 //!<! shower dispersion eta - phi
  TH2F * fhSphericityE[2] ;                     //!<! shower sphericity in eta vs phi
  TH2F * fhDispEtaDispPhiEBin[2][5] ;           //!<! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]

 //  Weight studies
  
  TH2F * fhECellClusterRatio;                   //!<! E cell / e cluster vs e cluster for selected electrons
  TH2F * fhECellClusterLogRatio;                //!<! log (E cell / E cluster)  vs E cluster for selected electrons
  TH2F * fhEMaxCellClusterRatio;                //!<! E max cell / E cluster vs E cluster for selected electrons
  TH2F * fhEMaxCellClusterLogRatio;             //!<! log (e max cell / e cluster) vs e cluster for selected electrons
  TH2F * fhLambda0ForW0[14];                    //!<! L0 for 7 defined w0= 3, 3.5 ... 6 for selected electrons
//TH2F * fhLambda1ForW0[14];                    //!<! L1 for 7 defined w0= 3, 3.5 ... 6 for selected electrons
  
  // Fill MC dependent histograms, Origin of this cluster is ...

  TH2F * fhMCDeltaE[2][10]  ;                   //!<! MC-Reco E distribution coming from MC particle
  TH2F * fhMC2E[2][10]  ;                       //!<! E distribution, Reco vs MC coming from MC particle
  
  TH1F * fhMCE[2][10];                          //!<! Number of identified electron vs cluster energy coming from MC particle
  TH1F * fhMCPt[2][10];                         //!<! Number of identified electron vs cluster energy coming from MC particle
  TH2F * fhMCPhi[2][10];                        //!<! Phi of identified electron coming from MC particle
  TH2F * fhMCEta[2][10];                        //!<! eta of identified electron coming from MC particle
  
  // Shower Shape MC

  TH2F * fhMCELambda0[2][6] ;                   //!<! E vs Lambda0 from MC particle
  
  TH2F * fhMCEDispEta[2][6] ;                   //!<! shower dispersion in eta direction from MC particle
  TH2F * fhMCEDispPhi[2][6] ;                   //!<! shower dispersion in phi direction from MC particle
  TH2F * fhMCESumEtaPhi[2][6] ;                 //!<! shower dispersion in eta vs phi direction from MC particle
  TH2F * fhMCEDispEtaPhiDiff[2][6] ;            //!<! shower dispersion in eta -phi direction from MC particle
  TH2F * fhMCESphericity[2][6] ;                //!<! shower sphericity, eta vs phi from MC particle

  TH2F * fhMCElectronELambda0NoOverlap ;        //!<! E vs Lambda0 from MC electrons, no overlap
  TH2F * fhMCElectronELambda0TwoOverlap ;       //!<! E vs Lambda0 from MC electrons, 2 particles overlap
  TH2F * fhMCElectronELambda0NOverlap ;         //!<! E vs Lambda0 from MC electrons, N particles overlap
  
  //Embedding
  TH2F * fhEmbeddedSignalFractionEnergy ;       //!<! Fraction of electron energy of embedded signal vs cluster energy
  
  TH2F * fhEmbedElectronELambda0FullSignal ;    //!<! Lambda0 vs E for embedded electrons with more than 90% of the cluster energy
  TH2F * fhEmbedElectronELambda0MostlySignal ;  //!<! Lambda0 vs E for embedded electrons with 90%<fraction<50%
  TH2F * fhEmbedElectronELambda0MostlyBkg ;     //!<! Lambda0 vs E for embedded electrons with 50%<fraction<10%
  TH2F * fhEmbedElectronELambda0FullBkg ;       //!<! Lambda0 vs E for embedded electrons with less than 10% of the cluster energy
  
  /// Copy constructor not implemented.
  AliAnaElectron(              const AliAnaElectron & el) ;
    
  /// Assignment operator not implemented.
  AliAnaElectron & operator = (const AliAnaElectron & el) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaElectron,5) ;
  /// \endcond

} ;
 

#endif//ALIANAELECTRON_H



