#ifndef ALIANAELECTRON_H
#define ALIANAELECTRON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaElectron.h 27413 2008-07-18 13:28:12Z gconesab $ */

//_________________________________________________________________________
//
// Class for the electron identification, 
// Clusters from calorimeters are identified as electrons
// and kept in the AOD. Few histograms produced.
// Copy of AliAnaPhoton just add electron id.
//

//-- Author: Gustavo Conesa (LPSC-IN2P3-CNRS)

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

class AliAnaElectron : public AliAnaPartCorrBaseClass {

 public: 
  
  AliAnaElectron() ;                                       // default ctor
  
  virtual ~AliAnaElectron() { ; }                          // virtual dtor
  
 private:
  
  AliAnaElectron(const AliAnaElectron & g) ;               // cpy ctor
  
  AliAnaElectron & operator = (const AliAnaElectron & g) ; // cpy assignment

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
  
  void         FillShowerShapeHistograms( AliVCluster* cluster, const Int_t mcTag , const Int_t pidTag) ;
  
  void         SwitchOnFillShowerShapeHistograms()    { fFillSSHistograms = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()   { fFillSSHistograms = kFALSE ; }  
  
  void         RecalibrateCellAmplitude(Float_t  & amp,  const Int_t absId);
  
  void         WeightHistograms(AliVCluster *clus);
  
  void         SwitchOnFillWeightHistograms()         { fFillWeightHistograms = kTRUE  ; }
  void         SwitchOffFillWeightHistograms()        { fFillWeightHistograms = kFALSE ; }  
  
  //---------------------------------------
  // Analysis parameters setters getters
  //---------------------------------------
  
  TString      GetCalorimeter()                 const { return fCalorimeter        ; }
  void         SetCalorimeter(TString  & det)         { fCalorimeter = det         ; }
    
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
    
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
    if(n > 10) fNOriginHistograms = 10; }

  // For histograms in arrays, index in the array, corresponding to a particle
  enum mcTypes    { mcPhoton = 0,        mcPi0Decay = 1,       mcOtherDecay = 2,  
                    mcPi0 = 3,           mcEta = 4,            mcElectron = 5,       
                    mcConversion = 6,    mcOther = 7,          mcAntiNeutron = 8,    
                    mcAntiProton = 9                                                 };    
  
  enum mcssTypes  { mcssPhoton = 0,      mcssOther = 1,       mcssPi0 = 2,         
                    mcssEta = 3,         mcssConversion = 4,  mcssElectron = 5       };  
  
  private:
 
  TString  fCalorimeter ;                      // Calorimeter where the gamma is searched;
  Float_t  fMinDist ;                          // Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                          // Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                          // One more cut on distance used for acceptance-efficiency study
  Double_t fTimeCutMin  ;                      // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                      // Remove clusters/cells with time larger than this value, in ns
  Int_t    fNCellsCut ;                        // Accept for the analysis clusters with more than fNCellsCut cells
  Bool_t   fFillSSHistograms ;                 // Fill shower shape histograms
  Bool_t   fFillWeightHistograms ;             // Fill weigth histograms
  Int_t    fNOriginHistograms;                 // Fill only NOriginHistograms of the 14 defined types

  Float_t  fdEdxMin;                           // Max dEdx for electrons
  Float_t  fdEdxMax;                           // Min dEdx for electrons
  Float_t  fEOverPMin;                         // Max E/p for electrons, after dEdx cut
  Float_t  fEOverPMax;                         // Min E/p for electrons, after dEdx cut

  //Histograms 
  TH2F * fhdEdxvsE;                            //! matched track dEdx vs cluster E 
  TH2F * fhdEdxvsP;                            //! matched track dEdx vs track P
  TH2F * fhEOverPvsE;                          //! matched track E cluster over P track vs cluster E, after dEdx cut 
  TH2F * fhEOverPvsP;                          //! matched track E cluster over P track vs track P, after dEdx cut 

  TH2F * fhNCellsE[2];                         //! number of cells in cluster vs E 
  TH2F * fhMaxCellDiffClusterE[2];             //! Fraction of energy carried by cell with maximum energy
  TH2F * fhTimeE[2];                           //! E vs Time of selected cluster 

  TH1F * fhE[2]    ;                           //! Number of identified electron vs energy
  TH1F * fhPt[2]   ;                           //! Number of identified electron vs transerse momentum 
  TH2F * fhPhi[2]  ;                           //! Azimuthal angle of identified  electron vs transerse momentum 
  TH2F * fhEta[2]  ;                           //! Pseudorapidity of identified  electron vs transerse momentum 
  TH2F * fhEtaPhi[2]  ;                        //! Pseudorapidity vs Phi of identified  electron for transerse momentum > 0.5
  TH2F * fhEtaPhi05[2]  ;                      //! Pseudorapidity vs Phi of identified  electron for transerse momentum < 0.5
  
  //Shower shape
  
  TH2F * fhDispE[2];                           //! cluster dispersion vs E
  TH2F * fhLam0E[2];                           //! cluster lambda0 vs  E
  TH2F * fhLam1E[2];                           //! cluster lambda1 vs  E  

  TH2F * fhDispETRD[2];                        //! cluster dispersion vs E, SM covered by TRD
  TH2F * fhLam0ETRD[2];                        //! cluster lambda0 vs  E, SM covered by TRD
  TH2F * fhLam1ETRD[2];                        //! cluster lambda1 vs  E, SM covered by TRD 
  
  TH2F * fhNCellsLam0LowE[2];                  //! cluster N cells vs lambda0, E<2
  TH2F * fhNCellsLam0HighE[2];                 //! cluster N Cells vs lambda0, E>2

  TH2F * fhEtaLam0LowE[2];                     //! cluster eta vs lambda0, E<2
  TH2F * fhPhiLam0LowE[2];                     //! cluster phi vs lambda0, E<2
  TH2F * fhEtaLam0HighE[2];                    //! cluster eta vs lambda0, E>2
  TH2F * fhPhiLam0HighE[2];                    //! cluster phi vs lambda0, E>2
    
  // Weight studies
  
  TH2F * fhECellClusterRatio;                  //! e cell / e cluster vs e cluster for selected electrons
  TH2F * fhECellClusterLogRatio;               //! log (e cell / e cluster)  vs e cluster for selected electrons
  TH2F * fhEMaxCellClusterRatio;               //! e max cell / e cluster vs e cluster for selected electrons
  TH2F * fhEMaxCellClusterLogRatio;            //! log (e max cell / e cluster) vs e cluster for selected electrons
  TH2F * fhLambda0ForW0[14];                    //! L0 for 7 defined w0= 3, 3.5 ... 6 for selected electrons
  //TH2F * fhLambda1ForW0[14];                    //! L1 for 7 defined w0= 3, 3.5 ... 6 for selected electrons
  
  //Fill MC dependent histograms, Origin of this cluster is ...

  TH2F * fhMCDeltaE[2][10]  ;                  //! MC-Reco E distribution coming from MC particle     
  TH2F * fhMC2E[2][10]  ;                      //! E distribution, Reco vs MC coming from MC particle
  
  TH1F * fhMCE[2][10];                          //! Number of identified electron vs cluster energy coming from MC particle
  TH1F * fhMCPt[2][10];                         //! Number of identified electron vs cluster energy coming from MC particle
  TH2F * fhMCPhi[2][10];                        //! Phi of identified electron coming from MC particle
  TH2F * fhMCEta[2][10];                        //! eta of identified electron coming from MC particle
  
  // Shower Shape MC

  TH2F * fhMCELambda0[2][6] ;                   //! E vs Lambda0     from MC particle
  
  TH2F * fhMCElectronELambda0NoOverlap ;        //! E vs Lambda0     from MC electrons, no overlap
  TH2F * fhMCElectronELambda0TwoOverlap ;       //! E vs Lambda0     from MC electrons, 2 particles overlap
  TH2F * fhMCElectronELambda0NOverlap ;         //! E vs Lambda0     from MC electrons, N particles overlap
  
  //Embedding
  TH2F * fhEmbeddedSignalFractionEnergy ;       //! Fraction of electron energy of embedded signal vs cluster energy
  
  TH2F * fhEmbedElectronELambda0FullSignal ;    //!  Lambda0 vs E for embedded electrons with more than 90% of the cluster energy
  TH2F * fhEmbedElectronELambda0MostlySignal ;  //!  Lambda0 vs E for embedded electrons with 90%<fraction<50% 
  TH2F * fhEmbedElectronELambda0MostlyBkg ;     //!  Lambda0 vs E for embedded electrons with 50%<fraction<10% 
  TH2F * fhEmbedElectronELambda0FullBkg ;       //!  Lambda0 vs E for embedded electrons with less than 10% of the cluster energy
  
   ClassDef(AliAnaElectron,2)

} ;
 

#endif//ALIANAELECTRON_H



