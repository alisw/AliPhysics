#ifndef AliHighPtTreeAnalysis_h
#define AliHighPtTreeAnalysis_h




// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxfileName = 1;
const Int_t kMaxvtxESD = 1;
const Int_t kMaxesdTrack = 1;
const Int_t kMaxextTPCInnerC = 1;
const Int_t kMaxextInnerParamC = 1;
const Int_t kMaxextInnerParam = 1;
const Int_t kMaxextInnerParamRef = 1;

class TTree;

#include "TTreeStream.h"
#include "TObject.h"



class AliHighPtTreeAnalysis : public TObject {
public :
  AliHighPtTreeAnalysis();
  virtual ~AliHighPtTreeAnalysis();
  void InitAnalysis( TString file );
  TTreeSRedirector * GetStreamer() { return fStreamer;}
  void InitStreamer( const char * streamerName="AliHighPtTreeAnalysisStreamer.root") { fStreamer=new TTreeSRedirector(streamerName,"recreate");}
  void CloseStreamer(  ) { if (fStreamer) delete fStreamer;}


  virtual Int_t         GetEntry(Long64_t entry);
  virtual Long64_t      LoadTree(Long64_t entry);
  virtual Long64_t      LoadV0Tree(Long64_t entry);
  virtual void          Init(TTree *tree);
  virtual void          InitV0tree(TTree *tree);
  virtual void          Loop();
  virtual Bool_t        BaseCut(); // returns true if track passes baseCut
   virtual void         BookHistos();
  virtual void          FillHistos();
  virtual void          Terminate();
  //
  virtual TGraphErrors* Calc2DProfileContent(TH3D *h1, const char *projAxisName);
  //
  virtual void          MakeDCArPullFitsMI();


  virtual void          MakeDCArPullFits();
  virtual void          MakeDCArResFits();
  virtual void          MakePhiFits();
  virtual void          Make1pTresCovFits();
  virtual void          MakeTPCITSMatchingEff();
  virtual void          MakeK0trends();
  virtual Double_t      qoverptCorr(Double_t trEta, Double_t trPhi, Int_t type);
  virtual Bool_t        ConnectGenericHistos( const char *genericHistoFile );
  virtual void          RunPeriod();
  virtual void          SetApplyCorrections( const char *correctionFile );
  virtual void          MakePowerFit(Int_t entries);  // make power fit
  virtual Bool_t        GetK0TrendFitFunction(TF1 *fLinearFitK0sShift, TF1 *fLinearFitK0sSigma, Int_t Type, Int_t Charge);
  virtual void          MakedcaRTrends();
  virtual void          MakeDeltaPhiTrends();
  virtual void          MakeEfficiencyTrends();
  
  virtual void          Plot1D(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName);
   virtual void          Plot2D(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName);
  virtual void          Plot2DK0s(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName);
  virtual void          Plot1PtRes(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX,  const char *xaxisName, const char *plotName);
  virtual void          PlotEff(TH3D *hTPCITS, TH3D *hTPC, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName ,const char *plotName);
  virtual void          SetHistoProperties(TH1D *hist, Int_t marker, Int_t color, Double_t yMin, Double_t yMax);
  virtual void          MakeAllPlots();
  
  virtual void          SetMakePlots(Bool_t makeAllPlots);
  virtual void          SetHighPtTree(TTree *HighPtTree) { Init(HighPtTree); };
  virtual void          SetV0Tree(TTree *V0Tree) { InitV0tree(V0Tree); fMakeV0s = kTRUE; };
  virtual void		 SetRunV0s(Bool_t bV0s )  { fMakeV0s = bV0s; };
  virtual void          SetPeriodName( const char *ch ) {  fPeriodName = new TString( ch ); };
  virtual void		 SetMakeFitPerfomancePlots( Bool_t bFPP ) { fMakeFitPerfomancePlots = bFPP; }
  
public:  
  TTree          *fChain;                     //!pointer to the analyzed TTree or TChain
  TTree          *fV0Chain;                   //!pointer to the V0 tree    
  TTree          *OutTree;                    // output tree
  TTreeSRedirector *fStreamer;                //! fStreamer
  //
  // Declaration of leaf types for the highPt tree
  //
  TObjString            *fChunkName;          // highPt - chunk name  -  never used -  TO REMOVE    MI!!!
  Double_t               runNumber;           // highPt - run number as double - should be carefully removed MI!!! 
  Int_t                  fRunNumberInt;       // highPt - run number as Ingeter
  Int_t                  fMult;
  TObjString            *triggerClass;
  Double_t               Bz;
  Int_t                  BzInt;
  AliESDVertex          *vtxESD;
  AliESDtrack           *esdTrack;
  TParticle             *particle;
  AliExternalTrackParam *extTPCInnerC;
  AliExternalTrackParam *extInnerParamC;
  AliExternalTrackParam *extInnerParam;
  AliExternalTrackParam *extInnerParamRef;
  Double_t               chi2TPCInnerC;
  Double_t               chi2InnerC;
  //
  TString        *fPeriodName;                // LHC period name
  Bool_t          fMakePeriod;                // switch to run per period ExperQA - who is setting ????      
  Bool_t          fHasMC;                     // swith has MC information - automatically got from the tree
  Double_t        fPtCut;                     // pt cut lowe rvalue of pt ranges for some histograms - name is misleading as it is not general
  Int_t           fBfield;                    // cache value for B field -- TO FIX - use abasolute field, currently only sign e.g -0.5 T, 0.5 T
  //
  //   
  Bool_t          fApplyCorrections;          // switch to apply q/pt correction
  Bool_t          fMakePlots;                 // switch to run  MakePlot and produce figure 
  Bool_t	  fMakeV0s;                   // switch to make plots/trending for the V0s  (currently harrdwired kTRUE)
  Bool_t	  fMakeFitPerfomancePlots;    // swith to make performance plot
  // set apply q/pt correction - here we stroe parameters of corrections - factors com from external files set by ::SetApplyCorrections( const char *correctionFile )
  Double_t       *fCorrectionAside;           // delta(q/pt) per sector on A side  - combined tracks
  Double_t       *fCorrectionCside;           // delta(q/pt) per sector on C side  - combined tracks
  Double_t       *fCorrectionAsideTPCInner;   // delta(q/pt) per sector on A side  - TPC only
  Double_t       *fCorrectionCsideTPCInner;   // delta(q/pt) per sector on C side  - TPC only
  Double_t       *fCorrectionAsideTPCInnerC;  // delta(q/pt) per sector on A side  - TPC constrained 
  Double_t       *fCorrectionCsideTPCInnerC;  // delta(q/pt) per sector on C side  - TPC constrained 
  
  Int_t fNtracks_TPCLowPt;                    // 
  Int_t fNtracks_TPCHighPt;
  Int_t fNtracks_TPCITSLowPt;
  Int_t fNtracks_TPCITSHighPt;
   
  // Declaration of V0 data members
  AliESDv0              *v0;                  
  AliESDtrack           *v0track0;
  AliESDtrack           *v0track1;
  // Histos
  TH3D *hPulldcaRTPConly_vs_eta_1pT;          
  TH3D *hPulldcaRcomb_vs_eta_1pT;
  TH3D *hResdcaRTPConly_vs_eta_1pT;           
  TH3D *hResdcaRcomb_vs_eta_1pT;
  TH3D *hphiPull_vs_eta_1pT;                  
  TH3D *hphiRes_vs_eta_1pT; 
  TH3D *hPulldcaR_vs_eta_pT_Aside;            
  TH3D *hPulldcaR_vs_eta_pT_Cside;
  TH3D *hPulldcaRTPCInner_vs_eta_pT_Aside;    
  TH3D *hPulldcaRTPCInner_vs_eta_pT_Cside;
  TH3D *hResdcaR_vs_eta_pT_Aside;             
  TH3D *hResdcaR_vs_eta_pT_Cside;
  TH3D *hResdcaRTPCInner_vs_eta_pT_Aside;     
  TH3D *hResdcaRTPCInner_vs_eta_pT_Cside;
  TH3D *hphiPull_vs_eta_pT_Aside;             
  TH3D *hphiPull_vs_eta_pT_Cside;
  TH3D *hphiRes_vs_eta_pT_Aside;              
  TH3D *hphiRes_vs_eta_pT_Cside;
  TH3D *hPulldcaR_vs_phi_pT_Aside;            
  TH3D *hPulldcaR_vs_phi_pT_Cside;
  TH3D *hPulldcaRTPCInner_vs_phi_pT_Aside;    
  TH3D *hPulldcaRTPCInner_vs_phi_pT_Cside;
  TH3D *hResdcaR_vs_phi_pT_Aside;             
  TH3D *hResdcaR_vs_phi_pT_Cside;
  TH3D *hResdcaRTPCInner_vs_phi_pT_Aside;     
  TH3D *hResdcaRTPCInner_vs_phi_pT_Cside;
  TH3D *hphiPull_vs_phi_pT_Aside;             
  TH3D *hphiPull_vs_phi_pT_Cside;
  TH3D *hphiRes_vs_phi_pT_Aside;              
  TH3D *hphiRes_vs_phi_pT_Cside;
  TH3D *heta_phi_pT;
  TH3D *hphi_vs_eta_pT_cutTPC;                
  TH3D *hphi_vs_eta_pT_cutTPCITS;
  // histogram for 1/pt shift calculation
  TH3F *h1pt_vs_eta_phi;
  TH3D *h1ptRes_vs_phi_pT_Aside;        
  TH3D *h1ptRes_vs_phi_pT_Cside;        // 1/pT resolution from cov. matrix
  TH3D *h1ptRes_vs_mult_pT_Aside;       
  TH3D *h1ptRes_vs_mult_pT_Cside;           // 1/pT resolution from cov. matrix vs mult.
  TH3D *h1ptSigma_vs_phi_pT_Aside;      
  TH3D *h1ptSigma_vs_phi_pT_Cside;      // sigma 1/pT from cov. matrix
  TH3D *h1ptSigma_vs_mult_pT_Aside;     
  TH3D *h1ptSigma_vs_mult_pT_Cside;      // sigma 1/pT from cov. matrix vs mult.
  //
  // histogram for 1/pt shift calculation for TPCInnerC
  //
  TH3F *h1ptTPCInnerC_vs_eta_phi;
  TH3D *h1ptResTPCInnerC_vs_phi_pT_Aside;    
  TH3D *h1ptResTPCInnerC_vs_phi_pT_Cside;  // 1/pT resolution from cov. matrix TPCInnerC
  TH3D *h1ptResTPCInnerC_vs_mult_pT_Aside;   
  TH3D *h1ptResTPCInnerC_vs_mult_pT_Cside;  // 1/pT resolution from cov. matrix vs mult. TPCInnerC
  TH3D *h1ptSigmaTPCInnerC_vs_phi_pT_Aside;  
  TH3D *h1ptSigmaTPCInnerC_vs_phi_pT_Cside; // 1/pT sigma from cov. matrix TPCInnerC
  TH3D *h1ptSigmaTPCInnerC_vs_mult_pT_Aside; 
  TH3D *h1ptSigmaTPCInnerC_vs_mult_pT_Cside; // 1/pT sigma from cov. matrix vs mult. TPCInnerC
  //
  // histogram for 1/pt shift calculation for TPCInner
  TH3F *h1ptTPCInner_vs_eta_phi; 
  TH3D *h1ptResTPCInner_vs_phi_pT_Aside;    
  TH3D *h1ptResTPCInner_vs_phi_pT_Cside;     // 1/pT resolution from cov. matrix TPCInner
  TH3D *h1ptResTPCInner_vs_mult_pT_Aside;  
  TH3D *h1ptResTPCInner_vs_mult_pT_Cside;   // 1/pT resolution from cov. matrix vs mult. TPCInner
  TH3D *h1ptSigmaTPCInner_vs_phi_pT_Aside;  
  TH3D *h1ptSigmaTPCInner_vs_phi_pT_Cside;   // 1/pT sigma from cov. matrix TPCInner
  TH3D *h1ptSigmaTPCInner_vs_mult_pT_Aside; 
  TH3D *h1ptSigmaTPCInner_vs_mult_pT_Cside;   // 1/pT sigma from cov. matrix vs mult. TPCInner
  //
  // Histogramm for V0s
  //
  TH3D *hK0sPull_vs_alpha_1pT_pos;
  TH3D *hK0sRes_vs_alpha_1pT_pos;
  TH3D *hK0sPull_vs_alpha_1pT_neg;
  TH3D *hK0sRes_vs_alpha_1pT_neg;
  // MC info
  TH3D *hptPull_vs_eta_pT;
  TH3D *hptRes_vs_eta_pT;
  TH3D *hptPullTPCInnerC_vs_eta_pT;
  TH3D *hptResTPCInnerC_vs_eta_pT;
  TH3D *hptPullTPCInner_vs_eta_pT;
  TH3D *hptResTPCInner_vs_eta_pT;

  TH3D *hptPull_vs_phi_pT;
  TH3D *hptRes_vs_phi_pT;
  TH3D *hptPullTPCInnerC_vs_phi_pT;
  TH3D *hptResTPCInnerC_vs_phi_pT;
  TH3D *hptPullTPCInner_vs_phi_pT;
  TH3D *hptResTPCInner_vs_phi_pT;
  
  ClassDef(AliHighPtTreeAnalysis, 2); // example of analysis
};

#endif
