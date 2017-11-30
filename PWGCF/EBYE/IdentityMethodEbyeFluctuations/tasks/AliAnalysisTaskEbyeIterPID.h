#ifndef ALIANALYSISEBYERATIOS_H
#define ALIANALYSISEBYERATIOS_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class THn;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliPIDResponse;
class AliHeader;
class AliESDpid;
class fPIDCombined;


#include "AliAnalysisTaskSE.h"
#include "AliPIDCombined.h"
#include "AliTPCdEdxInfo.h"
#include "THnSparse.h"
#include "THn.h"
#include "TCutG.h"
#include "TTreeStream.h"
#include "AliESDv0Cuts.h"

// class AliAnalysisTaskPIDetaTreeElectrons : public AliAnalysisTaskPIDV0base {
class AliAnalysisTaskEbyeIterPID : public AliAnalysisTaskSE {
 public:
   
// ---------------------------------------------------------------------------------
//                           Constructor and Destructor
// ---------------------------------------------------------------------------------
   
  AliAnalysisTaskEbyeIterPID(const char *name);
  AliAnalysisTaskEbyeIterPID();
  virtual ~AliAnalysisTaskEbyeIterPID();
  
  enum momentType {kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8,kLa=9,kLaLa=10,kCh=11,kChCh=12};
  enum momentTypeUnlike {
    kPiPosPiNeg=0,
    kPiPosKaNeg=1,
    kPiPosPrNeg=2,
    kKaPosPiNeg=3,
    kKaPosKaNeg=4,
    kKaPosPrNeg=5,
    kPrPosPiNeg=6,
    kPrPosKaNeg=7,
    kPrPosPrNeg=8,
    kLaPosLaNeg=9,
    kChPosChNeg=10,
  };
  enum trackCutBit {
    kNCrossedRowsTPC60=0,
    kNCrossedRowsTPC80=1,
    kNCrossedRowsTPC100=2,
    kMaxChi2PerClusterTPC3=3,
    kMaxChi2PerClusterTPC4=4,
    kMaxChi2PerClusterTPC5=5,
    kMaxDCAToVertexXYPtDepSmall=6,
    kMaxDCAToVertexXYPtDep=7,
    kMaxDCAToVertexXYPtDepLarge=8,
    kRequireITSRefit=9,
    kClusterRequirementITS=10,
    kVertexZSmall=11,
    kVertexZ=12,
    kVertexZLarge=13,
  };

// ---------------------------------------------------------------------------------
//                                    Methods
// ---------------------------------------------------------------------------------
  
  virtual void   UserCreateOutputObjects();            // create output objects
  virtual void   UserExec(Option_t *option);           // run over event-by-event and fill output objects
  virtual void   Terminate(Option_t *);                // run only once and terminate 
  
// ---------------------------------------------------------------------------------
//                                    Settings
// ---------------------------------------------------------------------------------

  void   SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void   SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
  void   Initialize();
  
   // Some boolian settings
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)          {fIncludeITS          = ifITSCuts;}
  void   SetFillArmPodTree(const Bool_t ifArmpodTree = kTRUE)       {fFillArmPodTree      = ifArmpodTree;}
  void   SetTightCuts(const Bool_t ifTightCuts = kFALSE)            {fTightCuts           = ifTightCuts;}
  void   SetDeDxCheck(const Bool_t ifDeDxCheck = kFALSE)            {fdEdxCheck           = ifDeDxCheck;}
  void   SetEffMatrix(const Bool_t ifEffMatrix = kFALSE)            {fEffMatrix           = ifEffMatrix;}
  void   SetCleanSamplesOnly(const Bool_t ifSamplesOnly = kFALSE)   {fCleanSamplesOnly    = ifSamplesOnly;}
  void   SetFillBayesianProb(const Bool_t ifBayesProb = kFALSE)     {fFillBayes           = ifBayesProb;}
  void   SetFillAllCutVariables(const Bool_t ifAllCuts = kFALSE)    {fFillCuts            = ifAllCuts;}
  void   SetFillDeDxTree(const Bool_t ifDeDxTree = kFALSE)          {fFillDeDxTree        = ifDeDxTree;}
  void   SetRunFastSimulation(const Bool_t ifFastSimul = kFALSE)    {fRunFastSimulation   = ifFastSimul;}
  void   SetRunFastHighMomentCal(const Bool_t ifFastHighMom = kFALSE)  {fRunFastHighMomentCal   = ifFastHighMom;}
  void   SetFillDnchDeta(const Bool_t ifDnchDetaCal = kFALSE)       {fFillDnchDeta        = ifDnchDetaCal;}
  void   SetIncludeTOF(const Bool_t ifIncludeTOF = kFALSE)          {fIncludeTOF          = ifIncludeTOF;}
  void   SetUseThnSparse(const Bool_t ifUseThnSparse = kFALSE)    {fUseThnSparse        = ifUseThnSparse;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)            {fUseCouts            = ifUseCouts;}  
  void   SetWeakAndMaterial(const Bool_t ifWeakAndMaterial = kFALSE){fWeakAndMaterial     = ifWeakAndMaterial;}
  void   SetFillTIdenTrees(const Bool_t ifTIdentity = kFALSE)       {fTIdentity           = ifTIdentity;}

  // Setters for the systematic uncertainty checks
  void   SetSystCentEstimator(const Int_t systCentEstimator = 0)  {fSystCentEstimatetor = systCentEstimator;}
  void   SetSystDCAxy(const Int_t systDCAxy = 0)                  {fSystDCAxy           = systDCAxy;}
  void   SetSystNCrossedRows(const Int_t systNCrossedRows = 0)    {fSystCrossedRows     = systNCrossedRows;}
  void   SetSystTPCChi2(const Int_t systTPCChi2 = 0)              {fSystChi2            = systTPCChi2;}
  void   SetSystVz(const Int_t systVz = 0)                        {fSystVz              = systVz;}

  // Setters for the eta momentum dEdx and centrality bins 
  void   SetSampleDeDxUpperEdge(const Float_t dEdxCleanUp = 140.) {fdEdxCleanUp         = dEdxCleanUp;}
  void   SetDeDxBinWidth(const Float_t dEdxBinWidth = 2.5)        {fdEdxBinWidth        = dEdxBinWidth;}
  void   SetDeDxLowerEdge(const Float_t dEdxLowerEdge = 20.)      {fdEdxDown            = dEdxLowerEdge;}
  void   SetDeDxUpperEdge(const Float_t dEdxUpperEdge = 1020.)    {fdEdxUp              = dEdxUpperEdge;}
  void   SetNEtabins(const Int_t nEtaBins = 16)                   {fnEtaBins            = nEtaBins;}
  void   SetPercentageOfEvents(const Int_t nPercentageOfEvents = 0) {fPercentageOfEvents = nPercentageOfEvents;}

  
  void   SetEtaLowerEdge(const Float_t etaLowerEdge = -0.8)       {fEtaDown             = etaLowerEdge;}
  void   SetEtaUpperEdge(const Float_t etaUpperEdge = 0.8)        {fEtaUp               = etaUpperEdge;}
  void   SetNMomBins(const Int_t nMombins = 150)                  {fnMomBins            = nMombins;}
  void   SetMomLowerEdge(const Float_t momLowerEdge = 0.2)        {fMomDown             = momLowerEdge;}
  void   SetMomUpperEdge(const Float_t momUpperEdge = 3.2)        {fMomUp               = momUpperEdge;}
  
  // Set the binning of centrality
  void   SetCentralityBinning(const Int_t tmpCentbins, Float_t tmpfxCentBins[])
  {
    // Create the histograms to be used in the binning of eta, cent and momentum
    std::cout << " Info::marsland: !!!!!! Centrality binning is being set !!!!!!! " << std::endl;
    fhEta  =  new TH1F("fhEta" ,"Eta Bins"       ,fnEtaBins        ,fEtaDown, fEtaUp );
    fhPtot =  new TH1F("fhPtot","Momentum Bins"  ,fnMomBins        ,fMomDown, fMomUp ); 
    fhCent =  new TH1F("fhCent","Centrality Bins",tmpCentbins-1    ,tmpfxCentBins );
    // ==========================================
    // prepare real data centrality bins
    fnCentbinsData = tmpCentbins; 
    fnCentBinsMC   = tmpCentbins-1;
    fxCentBins = new Float_t[fnCentbinsData]; 
    for (Int_t i=0; i<fnCentbinsData; i++) fxCentBins[i] =  tmpfxCentBins[i]; 
    fcentDownArr = new Float_t[fnCentBinsMC]; 
    fcentUpArr   = new Float_t[fnCentBinsMC]; 
    for (Int_t i=0; i<fnCentbinsData-1; i++) fcentDownArr[i] =  tmpfxCentBins[i];
    for (Int_t i=1; i<fnCentbinsData; i++)   fcentUpArr[i-1] =  tmpfxCentBins[i];
  }
  
  void SetMCEtaScanArray(const Int_t tmpEtaBinsMC, Float_t tmpetaDownArr[], Float_t tmpetaUpArr[])
  {    
    // set MC eta values to scan 
    std::cout << " Info::marsland: !!!!!! SetMCEtaScanArray is being set !!!!!!! " << std::endl;
    fnEtaWinBinsMC = tmpEtaBinsMC;
    fetaDownArr = new Float_t[fnEtaWinBinsMC]; 
    fetaUpArr   = new Float_t[fnEtaWinBinsMC]; 
    for (Int_t i=0; i<fnEtaWinBinsMC; i++) {
      fetaDownArr[i] =  tmpetaDownArr[i]; 
      fetaUpArr[i]   =  tmpetaUpArr[i]; 
    } 
  }
  
  void SetMCResonanceArray(const Int_t tmpNRes, TString tmpResArr[])
  {    
    // set MC eta values to scan 
    std::cout << " Info::marsland: !!!!!! SetMCResonanceArray is being set !!!!!!! " << std::endl;
    fnResBins = tmpNRes;
    fResonances = new TString[fnResBins]; 
    for (Int_t i=0; i<fnResBins; i++) fResonances[i] = tmpResArr[i]; 
     
  }
  
  void SetMCMomScanArray(const Int_t tmpMomBinsMC, Float_t tmppDownArr[], Float_t tmppUpArr[])
  {       
    // set MC momentum values to scan
    std::cout << " Info::marsland: !!!!!! SetMCMomScanArray is being set !!!!!!! " << std::endl;
    fnMomBinsMC = tmpMomBinsMC;
    fpDownArr = new Float_t[fnMomBinsMC]; 
    fpUpArr   = new Float_t[fnMomBinsMC]; 
    for (Int_t i=0; i<fnMomBinsMC; i++) {
      fpDownArr[i] =  tmppDownArr[i]; 
      fpUpArr[i]   =  tmppUpArr[i]; 
    }
  }
  
  void SetLookUpTableFirstMoments(TTree *lookUpTree, Int_t partType, Float_t pArr[],Float_t centArr[],Float_t etaArr[],const Int_t tmpMomBinsMC, const Int_t tmpCentbins, const Int_t tmpEtaBinsMC)
  {    
    // set MC eta values to scan 
    std::cout << " Info::marsland: !!!!!! SetLookUpTableFirstMoments is being set !!!!!!!   " << std::endl;
    //
    // fill arrays from lookup table
    TH1D *h=NULL, *h1=NULL;
    for (Int_t imom=0; imom<tmpMomBinsMC; imom++){
        for (Int_t icent=0; icent<tmpCentbins; icent++){
            for (Int_t ieta=0; ieta<tmpEtaBinsMC; ieta++){
                //
                // with resonances
                lookUpTree->Draw(Form("momentPos.fElements[%d]-momentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
                h= (TH1D*)lookUpTree->GetHistogram()->Clone(); h-> SetName("Res");
                if (partType==0)  fPiFirstMoments[0][imom][icent][ieta] = h->GetMean();
                if (partType==1)  fKaFirstMoments[0][imom][icent][ieta] = h->GetMean();
                if (partType==2)  fPrFirstMoments[0][imom][icent][ieta] = h->GetMean();
                if (partType==9)  fLaFirstMoments[0][imom][icent][ieta] = h->GetMean();
                if (partType==11) fChFirstMoments[0][imom][icent][ieta] = h->GetMean();
                delete h;
                //
                // without resonances
                lookUpTree->Draw(Form("noResmomentPos.fElements[%d]-noResmomentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
                h1= (TH1D*)lookUpTree->GetHistogram()->Clone(); h1-> SetName("noRes");
                if (partType==0)  fPiFirstMoments[1][imom][icent][ieta] = h1->GetMean();
                if (partType==1)  fKaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
                if (partType==2)  fPrFirstMoments[1][imom][icent][ieta] = h1->GetMean();
                if (partType==9)  fLaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
                if (partType==11) fChFirstMoments[1][imom][icent][ieta] = h1->GetMean();
                delete h1;
                
            }
        }
    }

  }
  
 private:
   
   AliAnalysisTaskEbyeIterPID(const AliAnalysisTaskEbyeIterPID&); 
   AliAnalysisTaskEbyeIterPID& operator=(const AliAnalysisTaskEbyeIterPID&); 

   // ---------------------------------------------------------------------------------
   //                                   Functions
   // ---------------------------------------------------------------------------------

   void  FillTPCdEdxReal();                   // Main function to fill all info + TIden 
   void  FillTPCdEdxCheck();                  // Quick check for the TPC dEdx 
   void  FillTPCdEdxMC();                     // Fill all info + TIdenMC from MC to do MC closure test
   void  FastGen();                           // Run over galice.root for Fastgen
   void  CalculateFastGenHigherMoments();     // Run over galice.root for Fastgen and calculate higher moments
   void  WeakAndMaterial();                   // Look full acceptance, weak decay and material 
   void  FillDnchDeta();                      // Fill dnch/deta values for each cent and eta bin  
   void  FillTPCdEdxMCEffMatrix();            // Prepare efficiency matrix
   void  FillCleanElectrons();                // Fill Clean Electrons 
   void  FillCleanPions();                    // Fill Clean Pions
   void  SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1); 
   Bool_t  ApplyDCAcutIfNoITSPixel(AliESDtrack *track);
   Int_t CountEmptyEvents(Int_t counterBin);  // Just count if there is empty events
   void  BinLogAxis(TH1 *h);
  
  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------
  
  
  AliPIDResponse   * fPIDResponse;            //! PID response object
  AliESDEvent      * fESD;                    //! ESD object
  TList            * fListHist;               //! list for histograms
  AliESDtrackCuts  * fESDtrackCuts;           // basic cut variables
  AliESDv0Cuts     * fESDtrackCutsV0;         // basic cut variables for V0
  AliESDtrackCuts  * fESDtrackCutsCleanSamp;  // basic cut variables for clean pion and electron form V0s
  AliPIDCombined   * fPIDCombined;            //! combined PID object
  AliTPCdEdxInfo   * fTPCdEdxInfo;            // detailed dEdx info
  AliStack         * fMCStack;                  // stack object to get Mc info
  
  TTree            * fTree;                   // data Tree for real Data
  TTree            * fIdenTree;               // data tree for TIdentity
  TTree            * fIdenTreeMC;             // data tree for TIdentity
  TTree            * fArmPodTree;             // Tree for clean pion and proton selection
  TTreeSRedirector * fTreeSRedirector;        //! temp tree to dump output
  TTree            * fTreeMCrec;              // tree for reconstructed moments
  TTree            * fTreeMCgen;              // tree for reconstructed moments
  TTree            * fTreeDnchDeta;           // tree for dnch/deta calculation
  TTree            * fTreeMC;                 // tree for mc samples
  TTree            * fTreedEdxCheck;          // tree to check dEdx performance for a small data sample 
  TTree            * fTreeBayes;              // tree to save bayesian probabilities
  TTree            * fTreeCuts;               // tree to save all variables for control plots
  TTree            * fTreeMCFullAcc;          // tree with full acceptance filled with MC
  TTree            * fTreeResonance;          // tree with full acceptance filled with MC
  TTree            * fTreeMCgenMoms;          // tree with higher moment calculations

  TH1F             * fhEta;                   // helper histogram for TIdentity tree
  TH1F             * fhCent;                  // helper histogram for TIdentity tree
  TH1F             * fhPtot;                  // helper histogram for TIdentity tree
  THnSparseF       * fhndEdx;                 // histogram which hold all dEdx info
  THnSparseF       * fhnExpected;             // histogram which hold all PIDresponse info
  THnSparseF       * fhnCleanEl;              // histogram which hold Clean Electrons
  THnSparseF       * fhnCleanKa;              // histogram which hold Clean Kaons
  THnSparseF       * fhnCleanDe;              // histogram which hold Clean Deuterons
 
  UInt_t            fTrackCutBits;           // integer which hold all cut variations as bits
  Int_t             myBin[3];                // binning array to be used for TIdentity module
  Int_t             myBinMC[3];              // binning array to be used for MC TIdentity module
  Double_t          fEtaDown;
  Double_t          fEtaUp;
  Int_t             fnEtaBins;
  Int_t             fPercentageOfEvents;     // when only a fPercentageOfEvents is enough
      
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fTIdentity;              // flag if tidentity trees are to be filled
  Bool_t            fWeakAndMaterial;        // flag for the Weak and Material analysis
  Bool_t            fEffMatrix;              // flag for efficiency matrix filling
  Bool_t            fdEdxCheck;              // flag to check only the dEdx performance
  Bool_t            fCleanSamplesOnly;       // flag for only clean sample production
  Bool_t            fTightCuts;              // addtional cuts from jens and Marian
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillBayes;              // switch whether to use bayesian PID or not
  Bool_t            fFillCuts;               // switch whether to fill all cut variables
  Bool_t            fFillDeDxTree;           // switch whether to fill dEdx tree
  Bool_t            fFillArmPodTree;         // switch whether to fill clean sample tree
  Bool_t            fRunFastSimulation;      // when running over galice.root do not fill other objects
  Bool_t            fRunFastHighMomentCal;   // when running over galice.root do not fill other objects
  
  Bool_t            fFillDnchDeta;           // switch on calculation of the dncdeta for fastgens
  Bool_t            fIncludeTOF;             // Include TOF information to investigate the efficiency loss effects on observable
  Bool_t            fUseThnSparse;           // in case thnsparse is filled
  Bool_t            fUseCouts;               // for debugging

  Int_t             fnMomBins;               // number of mombins --> for 20MeV slice 150 and 10MeV 300
  Float_t           fMomDown;                // bottom limit for the momentum range (default 0.2)
  Float_t           fMomUp;                  // uppper limit for the momentum range (default 3.2)
  Float_t           fdEdxBinWidth;           // bin width for the dEdx histograms (default 2.5)
  Float_t           fdEdxUp;                 // bottom limit for dEdx histogram (default 20) 
  Float_t           fdEdxDown;               // upper limit for dEdx histogram (default 1020)  
  Float_t           fdEdxCleanUp;            // upper limit for dEdx histogram of clean kaons and electrons (default 140)
  
  Float_t           fArmPodTPCSignal;
  Float_t           fArmPodptot;
  Float_t           fArmPodEta;
  Float_t           fArmPodCentrality;
  Float_t           fQt;
  Float_t           fAlfa;
  Float_t           fPiNSigmasTOF;           // TOF N sigma for Pion
  Float_t           fPrNSigmasTOF;           // TOF N sigma for Proton 
  
  Float_t           fdEdxEl;                 // Expected Electron dEdx
  Float_t           fdEdxKa;                 // Expected Kaon dEdx
  Float_t           fdEdxPi;                 // Expected Pion dEdx
  Float_t           fdEdxPr;                 // Expected Proton dEdx
  Float_t           fdEdxDe;                 // Expected Proton dEdx

  Float_t           fSigmaEl;                // Expected Electron sigma
  Float_t           fSigmaKa;                // Expected Kaon sigma
  Float_t           fSigmaPi;                // Expected Pion sigma
  Float_t           fSigmaPr;                // Expected Proton sigma
  Float_t           fSigmaDe;                // Expected Proton sigma

  Float_t           fNSigmasElTPC;           // TOF N sigma for Electron 
  Float_t           fNSigmasPiTPC;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTPC;           // TOF N sigma for Kaon 
  Float_t           fNSigmasPrTPC;           // TOF N sigma for Proton 
  Float_t           fNSigmasDeTPC;           // TOF N sigma for Proton 

  Float_t           fTPCSignalMC;
  Float_t           fptotMC;
  Float_t           fptotMCtruth;
  Float_t           fpTMC;
  Float_t           fEtaMC;
  Float_t           fCentralityMC;
  Int_t             fSignMC;
  
  Float_t           fPxMC;                     // x component of momentum
  Float_t           fPyMC;                     // y component of momentum
  Float_t           fPzMC;                     // z component of momentum
  
  Float_t           fElMC;
  Float_t           fPiMC;
  Float_t           fKaMC;
  Float_t           fPrMC;
  Float_t           fDeMC;
  Float_t           fMuMC;
  Float_t           fLaMC;
  
  Float_t           fptotMCgen;
  Float_t           fpTMCgen;
  Float_t           fEtaMCgen;
  Float_t           fCentralityMCgen;
  Int_t             fSignMCgen;
  Double_t          fMCImpactParameter;
  
  Float_t           fElMCgen;
  Float_t           fPiMCgen;
  Float_t           fKaMCgen;
  Float_t           fPrMCgen;
  Float_t           fDeMCgen;
  Float_t           fMuMCgen;
  Float_t           fLaMCgen;
   
  Float_t           fPx;                     // x component of momentum
  Float_t           fPy;                     // y component of momentum
  Float_t           fPz;                     // z component of momentum
  Float_t           fptot;                   // TPC momentum
  Float_t           fpT;                     // Transverse momentum
  Float_t           fY;                      // rapidity

  Int_t              fMultiplicity;           // Multiplicity in case of PbPb
  Int_t              fMultiplicityMC;
  Float_t            fCentrality;             // centrality information
  Double_t           fvZ;                     // Vertex position
  Int_t              fEventGID;               // global Event Id
  Int_t              fEventGIDMC;             // global MC event id
  Int_t              fEventCountInFile;       // event count per job
  Int_t              fEvent;                  // Event counter for Christian
  Int_t              fEventMC;                // Event id for MC data
  Int_t              fEventMCgen;             // Event id for MC generated
  
  Float_t            fTPCSignal;              // Measured dE/dx
  Double_t           myDeDx;                  // dEdx for TIdentity module    
  Int_t              signNew;                 // Sign Info for TIdentity module  
  Double_t           myDeDxMC;                  // dEdx for TIdentity module    
  Int_t              signNewMC;                 // Sign Info for TIdentity module  
  Float_t            fEta;                    // pseudo rapidity
  Double_t           fNContributors;          // Ntracks 
  Double_t           fTheta;                  // theta
  Double_t           fPhi;                    // azimut angle                 
  Int_t              fSign;                   // sign of the particle
  Int_t              fTPCShared;              // number of shared clusters
  Int_t              fNcl;                    // number of points used for dEdx

  Int_t              fnResBins;
  Int_t              fnEtaWinBinsMC;
  Int_t              fnMomBinsMC;
  Int_t              fnCentBinsMC;
  Int_t              fnResModeMC;
  Int_t              fnCentbinsData;
  Float_t            fMissingCl;
  
  // Additional cuts from marian
  Bool_t             fIsITSpixel01;           // if track has hits in innermost 2 pixels of ITS
  Int_t              fnITSclusters;           // number of ITS clusters
  Float_t            fPrimRestriction;        // prim vertex cut recommended by marian
  Float_t            fTPCvZ;                  // TPC vertex

  //   CleanSample cuts
  Bool_t             fCleanPionsFromK0;
  Bool_t             fCleanPion0FromK0;
  Bool_t             fCleanPion1FromK0;
  Bool_t             fCleanPion0FromLambda;
  Bool_t             fCleanPion1FromLambda;
  Bool_t             fCleanProton0FromLambda;
  Bool_t             fCleanProton1FromLambda;
  Bool_t             fHasTrack0FirstITSlayer;
  Bool_t             fHasTrack1FirstITSlayer;
  Bool_t             fHasV0FirstITSlayer;
  
  TCutG              *fPionCutG;
  TCutG              *fAntiProtonCutG;
  TCutG              *fProtonCutG;


  
  //  Variables for systematic uncertainty checks
  //  B field configurations -->  use default settings and analyse the following set of runs
  //  ***********************************************
  //  Field (++)  --> run interval is [137161, 138275]
  //  Field (--)  --> run interval is [138364, 139510]
  //  ***********************************************
  Int_t              fSystCentEstimatetor;   // 0 --> "V0M"   ||| -1 -->  "TRK" ||| +1 --> "CL1" 
  Int_t              fSystCrossedRows;       // 0 -->  80     ||| -1 -->   60   ||| +1 -->  100 
  Int_t              fSystDCAxy;             // 0 --> default ||| -1 --> -sigma ||| +1 --> +sigma 
  Int_t              fSystChi2;              // 0 -->  4      ||| -1 -->    3   ||| +1 -->   5
  Int_t              fSystVz;                // 0 -->  10     ||| -1 -->    8   ||| +1 -->   12
  Float_t            fPiFirstMoments[2][4][20][20];    //[fnResModeMC][fnMomBinsMC][fnCentBinsMC][fnEtaWinBinsMC]
  Float_t            fKaFirstMoments[2][4][20][20];    //[fnResModeMC][fnMomBinsMC][fnCentBinsMC][fnEtaWinBinsMC]
  Float_t            fPrFirstMoments[2][4][20][20];    //[fnResModeMC][fnMomBinsMC][fnCentBinsMC][fnEtaWinBinsMC]
  Float_t            fLaFirstMoments[2][4][20][20];    //[fnResModeMC][fnMomBinsMC][fnCentBinsMC][fnEtaWinBinsMC]
  Float_t            fChFirstMoments[2][4][20][20];    //[fnResModeMC][fnMomBinsMC][fnCentBinsMC][fnEtaWinBinsMC]
  Float_t            *fetaDownArr;           //[fnEtaWinBinsMC]
  Float_t            *fetaUpArr;             //[fnEtaWinBinsMC]
  Float_t            *fcentDownArr;          //[fnCentBinsMC]
  Float_t            *fcentUpArr;            //[fnCentBinsMC]
  Float_t            *fpDownArr;             //[fnMomBinsMC]
  Float_t            *fpUpArr;               //[fnMomBinsMC]
  Float_t            *fxCentBins;            //[fnCentbinsData]
  TString            *fResonances;           //[fnResBins]

  //
  // control and QA histograms
  //
  THnF             * fHistPosEffMatrixRec;       // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixRec;       // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixGen;       // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixGen;       // histogram efficiency matrix --> generated pions
 
  TH1F             * fHistEmptyEvent;         // control histogram for empty event
  TH1F             * fHistCentrality;         // control histogram for centrality
  TH1F             * fHistVertex;             // control histogram for vertexZ
  THnF             * fHistdEdxTPC;            // 5D hist of dEdx from all TPC
  TH2F             * fHistArmPod;             // control histogram for Armanteros Podolanski plot
   
  ClassDef(AliAnalysisTaskEbyeIterPID, 2);
  
};

#endif
