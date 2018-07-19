#ifndef ALIANALYSISTASKANTIHE4_CXX
#define ALIANALYSISTASKANTIHE4_CXX

// anti-alpha analysis
// Authors: Alexander Kalweit and Nicole Martin

class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDVertex;

#include "AliAnalysisTaskSE.h"
#include "THn.h"
#include "TH3F.h"
#include "TGraph.h"
#include "AliStack.h"

class AliAnalysisTaskAntiHe4 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAntiHe4();
  AliAnalysisTaskAntiHe4(const char *name);
  virtual ~AliAnalysisTaskAntiHe4() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  Int_t  Initialize();
  Int_t  SetupEvent();
  void   ResetEvent();

  
 private:
  AliInputEventHandler *fEventHandler;                              //  for ESDs or AODs
  AliESDEvent          *fESD;                                       //  ESD object
  TH1F                 *fHistCentralityClass10;                     //! histo to look at the centrality distribution    
  TH1F                 *fHistCentralityPercentile;                  //! histo to look at the centrality distribution 
  TH1F                 *fHistTriggerStat;                           //! Trigger statistics
  TH1F                 *fHistTriggerStatAfterEventSelection;        //! Trigger statistics
  TH3F                 *fHistDEDx;                                  //! final histograms for anti-alpha analysis
  TH3F                 *fHistTOF3D;                                 //! final histograms for anti-alpha analysis
  TH1F                 *fHistAlpha;                                 //! Alpha plot TOF mass
  TH1F                 *fHistAlphaSignal;                           //! Alpha plot TOF mass only signal candidates
  TGraph               *fGraphAlphaSignal;                          //! TGraph with alphas / dE/dx vs. mom
  Int_t                 fNCounter;                                  //! counts alphas

  TH2F                 *fHistDeDx;                                  //! histo for a dE/dx     
  TH3F                 *fHistDeDxRegion;                            //! histo for a dE/dx per Region 
  TH2F                 *fHistDeDxSharp;                             //! histo for a dE/dx with sharp cuts

  TH2F                 *fHistTOF2D;                                 //! histo for a TOF
  TH2F                 *fHistTOFnuclei;                             //! histo for a TOF nuclei  

  Int_t                 fNTriggers;                                 //! N Triggers used
  Double_t              fBBParametersLightParticles[5];             //! Bethe Bloch paramters for light paritcles
  Double_t              fBBParametersNuclei[5];                     //! Bethe Bloch paramters for nuclei
  Bool_t                fMCtrue;                                    //! flag if real data or MC is processed
  Bool_t                fTriggerFired[6];                           //! TriggerFired 0: INT7 | 1: MB | 2: CE | 3: SC | 4: EJE | 5: EGA
  //
  AliESDtrackCuts      *fESDtrackCuts;                              //  basic cut variables
  AliESDtrackCuts      *fESDtrackCutsSharp;                         //  sharp cut variables for final results
  AliESDpid            *fESDpid;                                    //  basic TPC object for n-sigma cuts
  THnF                 *fAntiAlpha;                                 //! histogram for particle ratios as a function of momentum: (0.) dca, (1.) sign, (2.) particle Type, (3.) p_tot

  TH1F                 *fHistHelium4PtGen;                          //! for MC
  TH1F                 *fHistHelium4PtGenPrim;                      //! for MC
  TH1F                 *fHistHelium4PtGenSec;                       //! for MC
  TH1F                 *fHistHelium4PtGenEta;                       //! for MC 
  TH1F                 *fHistHelium4PtGenPrimEta;                   //! for MC
  TH1F                 *fHistAntiHelium4PtGen;                      //! for MC
  TH1F                 *fHistAntiHelium4PtGenPrim;                  //! for MC                        
  TH1F                 *fHistAntiHelium4PtGenSec;                   //! for MC
  TH1F                 *fHistAntiHelium4PtGenEta;                   //! for MC
  TH1F                 *fHistHelium4PtAso;                          //! for MC
  TH1F                 *fHistHelium4PtAsoPrim;                      //! for MC
  TH1F                 *fHistHelium4PtAsoSec;                       //! for MC
  TH1F                 *fHistAntiHelium4PtAso;                      //! for MC


  void                  BinLogAxis(const THn *h, Int_t axisNumber); //  define function for log axis for search for Anti-Alpha candidates
  void                  BinLogAxis(const TH3 *h, Int_t axisNumber); 
  void                  BinLogAxis(const TH1 *h);
  Bool_t                IsTriggered();
  void                  SetBBParameters(Int_t runnumber);
  void                  MCGenerated(AliStack* stack); 

  //
  // output containers
  //
  TTree *fTree;  
  TObjArray * fOutputContainer; // ! output data container
  //
  // tree variables
  //
  Char_t fName[1000];  
  Int_t  fEvnt;
  Char_t  fFileName[1000];
  Int_t   fEventNumber[1000];
  // 
  Int_t  fItrk;
  //
  Double_t fEta[1000];
  Int_t    fKinkIndex[1000];
  Float_t    fCentrality[1000];

  //
  UShort_t   fTPCNsignal[1000];
  UShort_t   fTPCnCluster[1000];
  Double_t   fChi2PerClusterTPC[1000];
  Bool_t  fTPCRefit[1000];
  Double_t fTPCsignal0[1000];
  Double_t fTPCsignal1[1000];
  Double_t fTPCsignal2[1000];
  Double_t fTPCsignal3[1000];
  Int_t   fTPCSharedClusters[1000];
  UShort_t   fTPCNclsIter1[1000];
  //
  Double_t fITSsignal[1000];
  Int_t   fITSnCluster[1000];
  Double_t   fChi2PerClusterITS[1000];
  Bool_t  fITSRefit[1000];
  //
  Bool_t  fTOFRefit[1000];
  Bool_t  fTOFtime[1000];
  Bool_t  fTOFout[1000];
  Double_t fTOFsignalDz[1000];
  Double_t fTOFsignalDx[1000];
  //
  Bool_t fTRDin[1000];
  //
  Float_t fDCAZ[1000];
  Float_t fDCAXY[1000];
  //
  Double_t fTrkPtot[1000];
  Double_t fTPCPtot[1000];
  Double_t fTrackPt[1000];
  Double_t fDeDx[1000];
  Double_t fSign[1000];
  Float_t  fMass[1000];
  Float_t  fTime[1000];
  Float_t  fLength[1000];
  Double_t fSigmaQP[1000];
  //
  Bool_t  fAssociated[1000];
 
  //
  //
  AliAnalysisTaskAntiHe4(const AliAnalysisTaskAntiHe4&); // not implemented
  AliAnalysisTaskAntiHe4& operator=(const AliAnalysisTaskAntiHe4&); // not implemented
  
  ClassDef(AliAnalysisTaskAntiHe4, 1); // example of analysis
};

#endif
