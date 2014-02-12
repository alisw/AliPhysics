//////////////////////////////////////////////////////////

#ifndef AliHighPtTreeAnalysis_h
#define AliHighPtTreeAnalysis_h

#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TPaveStats.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TObjString.h>
#include <TObject.h>
#include "AliESDVertex.h"
#include <TNamed.h>
#include "AliVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliExternalTrackParam.h"
#include <TBits.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TTreeStream.h>
#include "TString.h"
#include "TSystem.h"
#include "TMinuit.h"
#include "TCanvas.h"


// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxfileName = 1;
const Int_t kMaxvtxESD = 1;
const Int_t kMaxesdTrack = 1;
const Int_t kMaxextTPCInnerC = 1;
const Int_t kMaxextInnerParamC = 1;
const Int_t kMaxextInnerParam = 1;
const Int_t kMaxextInnerParamRef = 1;

class AliHighPtTreeAnalysis {
public :

  
   AliHighPtTreeAnalysis( );
   AliHighPtTreeAnalysis( TString file );
   virtual ~AliHighPtTreeAnalysis();
   virtual Int_t         GetEntry(Long64_t entry);
   virtual Long64_t      LoadTree(Long64_t entry);
   virtual Long64_t      LoadV0Tree(Long64_t entry);
   virtual void          Init(TTree *tree);
   virtual void          InitV0tree(TTree *tree);
   virtual void          Loop();
   virtual Bool_t        BaseCut(); // returns true if track passes baseCut
   virtual void          BookHistos();
   virtual void          FillHistos();
   virtual void          Terminate();
   virtual TGraphErrors* Calc2DProfileContent(TH3D *h1, const char *projAxisName);
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
   virtual void          SetV0Tree(TTree *V0Tree) { InitV0tree(V0Tree); fV0s = kTRUE; };
   virtual void		 SetRunV0s(Bool_t bV0s )  { fV0s = bV0s; };
   virtual void          SetPeriodName( const char *ch ) {  fPeriodName = new TString( ch ); };
   virtual void		 SetMakeFitPerfomancePlots( Bool_t bFPP ) { fMakeFitPerfomancePlots = bFPP; }
   
private:
    
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *fV0Chain;
   
   TTree          *OutTree;
    
   Bool_t          fApplyCorrections;
   Bool_t          fMakePlots;
   Bool_t	   fV0s;
   Bool_t          fZipIn;
   Bool_t	   fMakeFitPerfomancePlots;
   
   Double_t       *fCorrectionAside;
   Double_t       *fCorrectionCside;
   Double_t       *fCorrectionAsideTPCInner;
   Double_t       *fCorrectionCsideTPCInner;
   Double_t       *fCorrectionAsideTPCInnerC;
   Double_t       *fCorrectionCsideTPCInnerC;
   
   Int_t fNtracks_TPCLowPt;   
   Int_t fNtracks_TPCHighPt;
   Int_t fNtracks_TPCITSLowPt;
   Int_t fNtracks_TPCITSHighPt;
    
   TString        *fPeriodName;
   Bool_t          fPeriod;
   Bool_t          hasMC;
   Double_t        pTcut;
   Int_t           fBfield;
    // Declaration of leaf types
  TObjString            *fileName;
  Double_t               runNumber;
  Int_t                  runNumberInt;
  Int_t                  mult;
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

  // Declaration of V0 data members
  AliESDv0              *v0;
  AliESDtrack           *v0track0;
  AliESDtrack           *v0track1;
  
  
   // Histos
 
  TH3D *hPulldcaRTPConly_vs_eta_1pT;          TH3D *hPulldcaRcomb_vs_eta_1pT;
  TH3D *hResdcaRTPConly_vs_eta_1pT;           TH3D *hResdcaRcomb_vs_eta_1pT;
  
  TH3D *hphiPull_vs_eta_1pT;                  
  TH3D *hphiRes_vs_eta_1pT; 
 
  TH3D *hPulldcaR_vs_eta_pT_Aside;            TH3D *hPulldcaR_vs_eta_pT_Cside;
  TH3D *hPulldcaRTPCInner_vs_eta_pT_Aside;    TH3D *hPulldcaRTPCInner_vs_eta_pT_Cside;
  TH3D *hResdcaR_vs_eta_pT_Aside;             TH3D *hResdcaR_vs_eta_pT_Cside;
  TH3D *hResdcaRTPCInner_vs_eta_pT_Aside;     TH3D *hResdcaRTPCInner_vs_eta_pT_Cside;
  TH3D *hphiPull_vs_eta_pT_Aside;             TH3D *hphiPull_vs_eta_pT_Cside;
  TH3D *hphiRes_vs_eta_pT_Aside;              TH3D *hphiRes_vs_eta_pT_Cside;
  TH3D *hPulldcaR_vs_phi_pT_Aside;            TH3D *hPulldcaR_vs_phi_pT_Cside;
  TH3D *hPulldcaRTPCInner_vs_phi_pT_Aside;    TH3D *hPulldcaRTPCInner_vs_phi_pT_Cside;
  TH3D *hResdcaR_vs_phi_pT_Aside;             TH3D *hResdcaR_vs_phi_pT_Cside;
  TH3D *hResdcaRTPCInner_vs_phi_pT_Aside;     TH3D *hResdcaRTPCInner_vs_phi_pT_Cside;
  TH3D *hphiPull_vs_phi_pT_Aside;             TH3D *hphiPull_vs_phi_pT_Cside;
  TH3D *hphiRes_vs_phi_pT_Aside;              TH3D *hphiRes_vs_phi_pT_Cside;
  TH3D *heta_phi_pT;
  TH3D *hphi_vs_eta_pT_cutTPC;                TH3D *hphi_vs_eta_pT_cutTPCITS;


// histogram for 1/pt shift calculation
  TH3F *h1pt_vs_eta_phi;
  TH3D *h1ptRes_vs_phi_pT_Aside;        TH3D *h1ptRes_vs_phi_pT_Cside;        // 1/pT resolution from cov. matrix
  TH3D *h1ptRes_vs_mult_pT_Aside;       TH3D *h1ptRes_vs_mult_pT_Cside;           // 1/pT resolution from cov. matrix vs mult.
  TH3D *h1ptSigma_vs_phi_pT_Aside;      TH3D *h1ptSigma_vs_phi_pT_Cside;      // sigma 1/pT from cov. matrix
  TH3D *h1ptSigma_vs_mult_pT_Aside;     TH3D *h1ptSigma_vs_mult_pT_Cside;      // sigma 1/pT from cov. matrix vs mult.

// histogram for 1/pt shift calculation for TPCInnerC
  TH3F *h1ptTPCInnerC_vs_eta_phi;
  TH3D *h1ptResTPCInnerC_vs_phi_pT_Aside;    TH3D *h1ptResTPCInnerC_vs_phi_pT_Cside;  // 1/pT resolution from cov. matrix TPCInnerC
  TH3D *h1ptResTPCInnerC_vs_mult_pT_Aside;   TH3D *h1ptResTPCInnerC_vs_mult_pT_Cside;  // 1/pT resolution from cov. matrix vs mult. TPCInnerC
  TH3D *h1ptSigmaTPCInnerC_vs_phi_pT_Aside;  TH3D *h1ptSigmaTPCInnerC_vs_phi_pT_Cside; // 1/pT sigma from cov. matrix TPCInnerC
  TH3D *h1ptSigmaTPCInnerC_vs_mult_pT_Aside; TH3D *h1ptSigmaTPCInnerC_vs_mult_pT_Cside; // 1/pT sigma from cov. matrix vs mult. TPCInnerC

// histogram for 1/pt shift calculation for TPCInner
  TH3F *h1ptTPCInner_vs_eta_phi; 
  TH3D *h1ptResTPCInner_vs_phi_pT_Aside;    TH3D *h1ptResTPCInner_vs_phi_pT_Cside;     // 1/pT resolution from cov. matrix TPCInner
  TH3D *h1ptResTPCInner_vs_mult_pT_Aside;   TH3D *h1ptResTPCInner_vs_mult_pT_Cside;   // 1/pT resolution from cov. matrix vs mult. TPCInner
  TH3D *h1ptSigmaTPCInner_vs_phi_pT_Aside;  TH3D *h1ptSigmaTPCInner_vs_phi_pT_Cside;   // 1/pT sigma from cov. matrix TPCInner
  TH3D *h1ptSigmaTPCInner_vs_mult_pT_Aside; TH3D *h1ptSigmaTPCInner_vs_mult_pT_Cside;   // 1/pT sigma from cov. matrix vs mult. TPCInner
// Histogramm for V0s
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

  

};

#endif

#ifdef AliHighPtTreeAnalysis_cxx

AliHighPtTreeAnalysis::AliHighPtTreeAnalysis( ) :
  fChain(0),fV0Chain(0),fApplyCorrections(kFALSE),fMakePlots(kTRUE),fV0s(kFALSE),fZipIn(kFALSE),fMakeFitPerfomancePlots(kFALSE),fCorrectionAside(0),fCorrectionCside(0),fCorrectionAsideTPCInner(0),fCorrectionCsideTPCInner(0),
  fCorrectionAsideTPCInnerC(0),fCorrectionCsideTPCInnerC(0),fNtracks_TPCLowPt(0),fNtracks_TPCHighPt(0),fNtracks_TPCITSLowPt(0),fNtracks_TPCITSHighPt(0),fPeriodName(0),
  fPeriod(kFALSE),hasMC(0),pTcut(3.),fBfield(0),fileName(0),runNumber(0),runNumberInt(0),mult(0),triggerClass(0),Bz(0),BzInt(0),vtxESD(0),
  esdTrack(0),particle(0),extTPCInnerC(0),extInnerParamC(0),extInnerParam(0),extInnerParamRef(0),chi2TPCInnerC(0),chi2InnerC(0),v0(0),v0track0(0),v0track1(0),hPulldcaRTPConly_vs_eta_1pT(0),
  hPulldcaRcomb_vs_eta_1pT(0),hResdcaRTPConly_vs_eta_1pT(0), hResdcaRcomb_vs_eta_1pT(0),hphiPull_vs_eta_1pT(0),hphiRes_vs_eta_1pT(0),
  hPulldcaR_vs_eta_pT_Aside(0),hPulldcaR_vs_eta_pT_Cside(0),hPulldcaRTPCInner_vs_eta_pT_Aside(0),hPulldcaRTPCInner_vs_eta_pT_Cside(0),hResdcaR_vs_eta_pT_Aside(0),
  hResdcaR_vs_eta_pT_Cside(0),hResdcaRTPCInner_vs_eta_pT_Aside(0),hResdcaRTPCInner_vs_eta_pT_Cside(0),hphiPull_vs_eta_pT_Aside(0),hphiPull_vs_eta_pT_Cside(0),hphiRes_vs_eta_pT_Aside(0),
  hphiRes_vs_eta_pT_Cside(0),hPulldcaR_vs_phi_pT_Aside(0),hPulldcaR_vs_phi_pT_Cside(0),hPulldcaRTPCInner_vs_phi_pT_Aside(0),hPulldcaRTPCInner_vs_phi_pT_Cside(0),hResdcaR_vs_phi_pT_Aside(0),
  hResdcaR_vs_phi_pT_Cside(0),hResdcaRTPCInner_vs_phi_pT_Aside(0),hResdcaRTPCInner_vs_phi_pT_Cside(0),hphiPull_vs_phi_pT_Aside(0),hphiPull_vs_phi_pT_Cside(0),
  hphiRes_vs_phi_pT_Aside(0),hphiRes_vs_phi_pT_Cside(0),heta_phi_pT(0),hphi_vs_eta_pT_cutTPC(0),hphi_vs_eta_pT_cutTPCITS(0),h1pt_vs_eta_phi(0),h1ptRes_vs_phi_pT_Aside(0),h1ptRes_vs_phi_pT_Cside(0),
  h1ptRes_vs_mult_pT_Aside(0),h1ptRes_vs_mult_pT_Cside(0),h1ptSigma_vs_phi_pT_Aside(0),h1ptSigma_vs_phi_pT_Cside(0),h1ptSigma_vs_mult_pT_Aside(0),h1ptSigma_vs_mult_pT_Cside(0),h1ptTPCInnerC_vs_eta_phi(0),
  h1ptResTPCInnerC_vs_phi_pT_Aside(0),h1ptResTPCInnerC_vs_phi_pT_Cside(0),h1ptResTPCInnerC_vs_mult_pT_Aside(0),h1ptResTPCInnerC_vs_mult_pT_Cside(0),h1ptSigmaTPCInnerC_vs_phi_pT_Aside(0),
  h1ptSigmaTPCInnerC_vs_phi_pT_Cside(0),h1ptSigmaTPCInnerC_vs_mult_pT_Aside(0),h1ptSigmaTPCInnerC_vs_mult_pT_Cside(0),h1ptTPCInner_vs_eta_phi(0),h1ptResTPCInner_vs_phi_pT_Aside(0),
  h1ptResTPCInner_vs_phi_pT_Cside(0),h1ptResTPCInner_vs_mult_pT_Aside(0),h1ptResTPCInner_vs_mult_pT_Cside(0),h1ptSigmaTPCInner_vs_phi_pT_Aside(0),h1ptSigmaTPCInner_vs_phi_pT_Cside(0),
  h1ptSigmaTPCInner_vs_mult_pT_Aside(0),h1ptSigmaTPCInner_vs_mult_pT_Cside(0),hK0sPull_vs_alpha_1pT_pos(0),hK0sRes_vs_alpha_1pT_pos(0),hK0sPull_vs_alpha_1pT_neg(0),hK0sRes_vs_alpha_1pT_neg(0)
{
}

AliHighPtTreeAnalysis::AliHighPtTreeAnalysis( TString file ) :
  fChain(0),
  fV0Chain(0),
  fApplyCorrections(kFALSE),
  fMakePlots(kTRUE),
  fV0s(kFALSE),
  fZipIn(kFALSE),
  fMakeFitPerfomancePlots(kFALSE),
  fCorrectionAside(0),
  fCorrectionCside(0),
  fCorrectionAsideTPCInner(0),
  fCorrectionCsideTPCInner(0),
  fCorrectionAsideTPCInnerC(0),
  fCorrectionCsideTPCInnerC(0),
  fNtracks_TPCLowPt(0),fNtracks_TPCHighPt(0),fNtracks_TPCITSLowPt(0),fNtracks_TPCITSHighPt(0),
  fPeriodName(0),
  fPeriod(kFALSE),
  hasMC(0),
  pTcut(3.),
  fBfield(0),
  fileName(0),
  runNumber(0),
  runNumberInt(0),
  mult(0),
  triggerClass(0),
  Bz(0),
  BzInt(0),
  vtxESD(0),
  esdTrack(0),
  particle(0),
  extTPCInnerC(0),
  extInnerParamC(0),
  extInnerParam(0),
  extInnerParamRef(0),
  chi2TPCInnerC(0),
  chi2InnerC(0),
  v0(0),
  v0track0(0),
  v0track1(0),
  hPulldcaRTPConly_vs_eta_1pT(0), hPulldcaRcomb_vs_eta_1pT(0),
  hResdcaRTPConly_vs_eta_1pT(0), hResdcaRcomb_vs_eta_1pT(0),
  hPulldcaR_vs_eta_pT_Aside(0), hPulldcaR_vs_eta_pT_Cside(0),
  hPulldcaRTPCInner_vs_eta_pT_Aside(0), hPulldcaRTPCInner_vs_eta_pT_Cside(0),
  hResdcaR_vs_eta_pT_Aside(0), hResdcaR_vs_eta_pT_Cside(0),
  hResdcaRTPCInner_vs_eta_pT_Aside(0), hResdcaRTPCInner_vs_eta_pT_Cside(0),
  hphiPull_vs_eta_pT_Aside(0), hphiPull_vs_eta_pT_Cside(0),
  hphiRes_vs_eta_pT_Aside(0), hphiRes_vs_eta_pT_Cside(0),
  hPulldcaR_vs_phi_pT_Aside(0), hPulldcaR_vs_phi_pT_Cside(0),
  hPulldcaRTPCInner_vs_phi_pT_Aside(0), hPulldcaRTPCInner_vs_phi_pT_Cside(0),
  hResdcaR_vs_phi_pT_Aside(0), hResdcaR_vs_phi_pT_Cside(0),
  hResdcaRTPCInner_vs_phi_pT_Aside(0), hResdcaRTPCInner_vs_phi_pT_Cside(0),
  hphiPull_vs_phi_pT_Aside(0), hphiPull_vs_phi_pT_Cside(0),
  hphiRes_vs_phi_pT_Aside(0), hphiRes_vs_phi_pT_Cside(0),
  heta_phi_pT(0),
  hphi_vs_eta_pT_cutTPC(0),
  hphi_vs_eta_pT_cutTPCITS(0),
  h1pt_vs_eta_phi(0),
  h1ptRes_vs_phi_pT_Aside(0),
  h1ptRes_vs_phi_pT_Cside (0),
  h1ptRes_vs_mult_pT_Aside(0),
  h1ptRes_vs_mult_pT_Cside(0),
  h1ptSigma_vs_phi_pT_Aside(0),
  h1ptSigma_vs_phi_pT_Cside(0),
  h1ptSigma_vs_mult_pT_Aside(0),
  h1ptSigma_vs_mult_pT_Cside(0),
  h1ptTPCInnerC_vs_eta_phi(0),
  h1ptResTPCInnerC_vs_phi_pT_Aside(0),
  h1ptResTPCInnerC_vs_phi_pT_Cside(0),
  h1ptResTPCInnerC_vs_mult_pT_Aside(0),
  h1ptResTPCInnerC_vs_mult_pT_Cside(0),
  h1ptSigmaTPCInnerC_vs_phi_pT_Aside(0),
  h1ptSigmaTPCInnerC_vs_phi_pT_Cside(0),
  h1ptSigmaTPCInnerC_vs_mult_pT_Aside(0),
  h1ptSigmaTPCInnerC_vs_mult_pT_Cside(0),
  h1ptTPCInner_vs_eta_phi(0),
  h1ptResTPCInner_vs_phi_pT_Aside(0),
  h1ptResTPCInner_vs_phi_pT_Cside(0),
  h1ptResTPCInner_vs_mult_pT_Aside(0),
  h1ptResTPCInner_vs_mult_pT_Cside(0),
  h1ptSigmaTPCInner_vs_phi_pT_Aside(0),
  h1ptSigmaTPCInner_vs_phi_pT_Cside(0),
  h1ptSigmaTPCInner_vs_mult_pT_Aside(0),
  h1ptSigmaTPCInner_vs_mult_pT_Cside(0),
  hK0sPull_vs_alpha_1pT_pos(0),
  hK0sRes_vs_alpha_1pT_pos(0),
  hK0sPull_vs_alpha_1pT_neg(0),
  hK0sRes_vs_alpha_1pT_neg(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
      TString V0file(0);
      
      if( file == 0 ){
        std::cout << "You have to specify the Input file" << std::endl;
        exit(1);
      }
      if( file.Contains("root_archive.zip") ){
        file += "#FilterEvents_Trees.root";
        fZipIn = kTRUE;
      }
      else{
      V0file = file.Copy();
      V0file.Replace(V0file.Length() - 11, 11, "V0s.root",8);
      }
      

      TTree *tree = NULL;
      TTree *V0tree = NULL;
      if(fZipIn){
      TFile *f = TFile::Open(file);
      f->GetObject("highPt",tree);
      f->GetObject("V0s",V0tree);  
      if(V0tree->GetEntries() > 1) fV0s = kTRUE;
     // if(!V0tree || V0tree->GetEntries() < 1) V0log
      }
      if(!fZipIn){
      TFile *f   = TFile::Open(file);
      if(f) f->GetObject("highPt",tree);
      TFile *fV0 = TFile::Open(V0file);
      if(fV0) fV0->GetObject("V0s",V0tree);
      if(V0tree) fV0s = kTRUE;
      }

 //   fMakePlots = kTRUE;
    if(tree)   Init(tree);
    if(V0tree) InitV0tree(V0tree);
    
}


AliHighPtTreeAnalysis::~AliHighPtTreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AliHighPtTreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t AliHighPtTreeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;

   return centry;
}

Long64_t AliHighPtTreeAnalysis::LoadV0Tree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fV0Chain) return -5;
   Long64_t centry = fV0Chain->LoadTree(entry);
   if (centry < 0) return centry;

   return centry;
}

void AliHighPtTreeAnalysis::InitV0tree(TTree *tree)
{
   // Set branch addresses and branch pointers

  if(!tree) return;
  fV0Chain = tree;
 
  fV0Chain->SetBranchAddress("v0.", &v0);
  fV0Chain->SetBranchAddress("track0.", &v0track0);
  fV0Chain->SetBranchAddress("track1.", &v0track1);
  
}

void AliHighPtTreeAnalysis::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  TString str(fChain->GetBranch("runNumber")->GetTitle());
            if(str[str.Length()-1]=='I')
              fChain->SetBranchAddress("runNumber", &runNumberInt);
            if(str[str.Length()-1]=='D')
              fChain->SetBranchAddress("runNumber", &runNumber);

          str = fChain->GetBranch("Bz")->GetTitle();
            if(str[str.Length()-1]=='I')
              fChain->SetBranchAddress("Bz", &BzInt);
            if(str[str.Length()-1]=='D')
              fChain->SetBranchAddress("Bz", &Bz);
  

   fChain->SetBranchAddress("esdTrack.", &esdTrack);
   fChain->SetBranchAddress("vtxESD.", &vtxESD);
//   fChain->SetBranchAddress("runNumber", &runNumber);
   fChain->SetBranchAddress("extTPCInnerC.", &extTPCInnerC);
   fChain->SetBranchAddress("chi2TPCInnerC", &chi2TPCInnerC);
   fChain->SetBranchAddress("mult", &mult);
//   fChain->SetBranchAddress("Bz", &Bz);

   if( fChain->GetBranchStatus("particle.") ){
     hasMC = kTRUE;
     fChain->SetBranchAddress("particle.", &particle);
   }
   else
     hasMC = kFALSE;

}

Bool_t AliHighPtTreeAnalysis::ConnectGenericHistos( const char *genericHistoFile )
{
    TFile *f = TFile::Open(genericHistoFile);
    if(!f) return kFALSE;
    if(strstr(genericHistoFile,"Bpos")) fBfield = 1;
      else if(strstr(genericHistoFile,"Bneg")) fBfield = -1;
             else fBfield = 0;
    
    
    heta_phi_pT                        = (TH3D*) f->Get("heta_phi_pT");
    hphi_vs_eta_pT_cutTPC              = (TH3D*) f->Get("hphi_vs_eta_pT_cutTPC");
    hphi_vs_eta_pT_cutTPCITS           = (TH3D*) f->Get("hphi_vs_eta_pT_cutTPCITS");
    
    h1pt_vs_eta_phi                    = (TH3F*) f->Get("h1pt_vs_eta_phi");
    h1ptTPCInner_vs_eta_phi            = (TH3F*) f->Get("h1ptTPCInner_vs_eta_phi");
    h1ptTPCInnerC_vs_eta_phi           = (TH3F*) f->Get("h1ptTPCInnerC_vs_eta_phi");
    
    hPulldcaR_vs_eta_pT_Aside          = (TH3D*) f->Get("hPulldcaR_vs_eta_pT_Aside");           hPulldcaR_vs_eta_pT_Cside          = (TH3D*) f->Get("hPulldcaR_vs_eta_pT_Cside");
    hPulldcaR_vs_phi_pT_Aside          = (TH3D*) f->Get("hPulldcaR_vs_phi_pT_Aside");           hPulldcaR_vs_phi_pT_Cside          = (TH3D*) f->Get("hPulldcaR_vs_phi_pT_Cside");
    hPulldcaRTPCInner_vs_eta_pT_Aside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_eta_pT_Aside");   hPulldcaRTPCInner_vs_eta_pT_Cside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_eta_pT_Cside");
    hPulldcaRTPCInner_vs_phi_pT_Aside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_phi_pT_Aside");   hPulldcaRTPCInner_vs_phi_pT_Cside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_phi_pT_Cside");
    hResdcaR_vs_eta_pT_Aside           = (TH3D*) f->Get("hResdcaR_vs_eta_pT_Aside");            hResdcaR_vs_eta_pT_Cside           = (TH3D*) f->Get("hResdcaR_vs_eta_pT_Cside");
    hResdcaR_vs_phi_pT_Aside           = (TH3D*) f->Get("hResdcaR_vs_phi_pT_Aside");            hResdcaR_vs_phi_pT_Cside           = (TH3D*) f->Get("hResdcaR_vs_phi_pT_Cside");
    hResdcaRTPCInner_vs_eta_pT_Aside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_eta_pT_Aside");    hResdcaRTPCInner_vs_eta_pT_Cside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_eta_pT_Cside");
    hResdcaRTPCInner_vs_phi_pT_Aside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_phi_pT_Aside");    hResdcaRTPCInner_vs_phi_pT_Cside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_phi_pT_Cside");
    
    hphiPull_vs_phi_pT_Aside           = (TH3D*) f->Get("hphiPull_vs_phi_pT_Aside");            hphiPull_vs_phi_pT_Cside           = (TH3D*) f->Get("hphiPull_vs_phi_pT_Cside");
    hphiRes_vs_phi_pT_Aside            = (TH3D*) f->Get("hphiRes_vs_phi_pT_Aside");             hphiRes_vs_phi_pT_Cside            = (TH3D*) f->Get("hphiRes_vs_phi_pT_Cside");
     
    h1ptRes_vs_phi_pT_Aside            = (TH3D*) f->Get("h1ptRes_vs_phi_pT_Aside");             h1ptRes_vs_phi_pT_Cside            = (TH3D*) f->Get("h1ptRes_vs_phi_pT_Cside");
    h1ptRes_vs_mult_pT_Aside           = (TH3D*) f->Get("h1ptRes_vs_mult_pT_Aside");            h1ptRes_vs_mult_pT_Cside           = (TH3D*) f->Get("h1ptRes_vs_mult_pT_Cside");
    h1ptSigma_vs_phi_pT_Aside          = (TH3D*) f->Get("h1ptSigma_vs_phi_pT_Aside");           h1ptSigma_vs_phi_pT_Cside          = (TH3D*) f->Get("h1ptSigma_vs_phi_pT_Cside");
    h1ptSigma_vs_mult_pT_Aside         = (TH3D*) f->Get("h1ptSigma_vs_mult_pT_Aside");          h1ptSigma_vs_mult_pT_Cside         = (TH3D*) f->Get("h1ptSigma_vs_mult_pT_Cside");
    h1ptResTPCInnerC_vs_phi_pT_Aside   = (TH3D*) f->Get("h1ptResTPCInnerC_vs_phi_pT_Aside");    h1ptResTPCInnerC_vs_phi_pT_Cside   = (TH3D*) f->Get("h1ptResTPCInnerC_vs_phi_pT_Cside");
    h1ptResTPCInnerC_vs_mult_pT_Aside  = (TH3D*) f->Get("h1ptResTPCInnerC_vs_mult_pT_Aside");   h1ptResTPCInnerC_vs_mult_pT_Cside  = (TH3D*) f->Get("h1ptResTPCInnerC_vs_mult_pT_Cside");
    h1ptSigmaTPCInnerC_vs_phi_pT_Aside = (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_phi_pT_Aside");  h1ptSigmaTPCInnerC_vs_phi_pT_Cside = (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_phi_pT_Cside");
    h1ptSigmaTPCInnerC_vs_mult_pT_Aside= (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_mult_pT_Aside"); h1ptSigmaTPCInnerC_vs_mult_pT_Cside= (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_mult_pT_Cside");
    h1ptResTPCInner_vs_phi_pT_Aside    = (TH3D*) f->Get("h1ptResTPCInner_vs_phi_pT_Aside");     h1ptResTPCInner_vs_phi_pT_Cside    = (TH3D*) f->Get("h1ptResTPCInner_vs_phi_pT_Cside");
    h1ptResTPCInner_vs_mult_pT_Aside   = (TH3D*) f->Get("h1ptResTPCInner_vs_mult_pT_Aside");    h1ptResTPCInner_vs_mult_pT_Cside   = (TH3D*) f->Get("h1ptResTPCInner_vs_mult_pT_Cside");
    h1ptSigmaTPCInner_vs_phi_pT_Aside  = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_phi_pT_Aside");   h1ptSigmaTPCInner_vs_phi_pT_Cside  = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_phi_pT_Cside");
    h1ptSigmaTPCInner_vs_mult_pT_Aside = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_mult_pT_Aside");  h1ptSigmaTPCInner_vs_mult_pT_Cside = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_mult_pT_Cside");
    
    hK0sPull_vs_alpha_1pT_pos = (TH3D*) f->Get("hK0sPull_vs_alpha_1pT_pos");
    hK0sRes_vs_alpha_1pT_pos  = (TH3D*) f->Get("hK0sRes_vs_alpha_1pT_pos");
    hK0sPull_vs_alpha_1pT_neg = (TH3D*) f->Get("hK0sPull_vs_alpha_1pT_neg");
    hK0sRes_vs_alpha_1pT_neg  = (TH3D*) f->Get("hK0sRes_vs_alpha_1pT_neg");
    
    hPulldcaRTPConly_vs_eta_1pT  = (TH3D*) f->Get("hPulldcaRTPConly_vs_eta_1pT");
    hPulldcaRcomb_vs_eta_1pT     = (TH3D*) f->Get("hPulldcaRcomb_vs_eta_1pT");
    hResdcaRTPConly_vs_eta_1pT   = (TH3D*) f->Get("hResdcaRTPConly_vs_eta_1pT");
    hResdcaRcomb_vs_eta_1pT      = (TH3D*) f->Get("hResdcaRcomb_vs_eta_1pT");

    hphiPull_vs_eta_1pT          = (TH3D*) f->Get("hphiPull_vs_eta_1pT");
    hphiRes_vs_eta_1pT           = (TH3D*) f->Get("hphiRes_vs_eta_1pT");
    
    if(hK0sPull_vs_alpha_1pT_pos && hK0sRes_vs_alpha_1pT_pos && hK0sPull_vs_alpha_1pT_neg && hK0sRes_vs_alpha_1pT_neg) fV0s = kTRUE;
    
    return kTRUE;
}

void AliHighPtTreeAnalysis::SetMakePlots(Bool_t makeAllPlots){ fMakePlots = makeAllPlots; }

void AliHighPtTreeAnalysis::SetApplyCorrections( const char *correctionFile )
{
    fApplyCorrections = kTRUE;
    
    TFile *fCorr = TFile::Open( correctionFile );
    TTree *CorrectionTree;
    fCorr->GetObject("PeriodTree",CorrectionTree);
  
  
    fCorrectionAside          = new Double_t[18];
    fCorrectionCside          = new Double_t[18];
    fCorrectionAsideTPCInner  = new Double_t[18];
    fCorrectionCsideTPCInner  = new Double_t[18];
    fCorrectionAsideTPCInnerC = new Double_t[18];
    fCorrectionCsideTPCInnerC = new Double_t[18];

    TH1D *hCorrtmp = new TH1D("hname","htitle",100,-1.,1.);
    for(Int_t i = 0; i != 18; ++i){
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAside.fElements[%d]",i));
      fCorrectionAside[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();
      
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCside.fElements[%d]",i));
      fCorrectionCside[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();
      
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAsideTPCInner.fElements[%d]",i));
      fCorrectionAsideTPCInner[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();
      
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCsideTPCInner.fElements[%d]",i));
      fCorrectionCsideTPCInner[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();
      
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAsideTPCInnerC.fElements[%d]",i));
      fCorrectionAsideTPCInnerC[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();
      
      CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCsideTPCInnerC.fElements[%d]",i));
      fCorrectionCsideTPCInnerC[i] = hCorrtmp->GetMean();
      hCorrtmp->Reset();     
    }
    

    
}
#endif // #ifdef AliHighPtTreeAnalysis_cxx
