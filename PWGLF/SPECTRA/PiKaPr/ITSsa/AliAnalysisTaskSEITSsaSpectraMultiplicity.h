#ifndef ALIANALYSISTASKSEITSSASPECTRAMULTIPLICITY_H
#define ALIANALYSISTASKSEITSSASPECTRAMULTIPLICITY_H

///////////////////////////////////////////////////////////////////////////
/// AliAnalysisTaskSE for the extraction of the various histograms to
/// study the pt spectra of identified hadrons vs multiplicity:
/// - multiplicity estimated with Reference Multiplicity and "V0M" percentiles
/// - log(dEdx)-log(dEdxBB) distributions for pions, kaons and protons in pt bins
/// - Pt distributions of pions, kaons and protons with nSigma PID
/// Authors: 
/// E. Biolcati, biolcati@to.infn.it
/// L. Milano,   milano@to.infn.it
/// F. Prino,    prino@to.infn.it
/// N. Jacazio,  jacazio@to.infn.it
///////////////////////////////////////////////////////////////////////////

/* $Id$ */

class TString;
class TTree;
class TH1F;
class TF1;
class TH2F;
class TRandom3;
class AliESDEvent;
class TNtuple;
class AliITSPIDResponse;
class THnSparse;//per la vettorizzazione
class AliESDtrackCuts;
class AliPPVsMultUtils;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"
#include <iostream>

using std::cout;
using std::endl;

#ifdef LOG_NO_INFO
#define errormsg(msg) do { } while (false)
#define warningmsg(msg) do { } while (false)
#define infomsg(msg) do { } while (false)
#else
#define errormsg(msg) AliErrorF("%s%s%s",redTxt,msg,normalTxt)
#define warningmsg(msg) AliWarningF("%s%s%s",magentaTxt,msg,normalTxt)
#define infomsg(msg) AliInfoF("%s%s%s",blueTxt,msg,normalTxt)
#endif

const char redTxt[] =     { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
const char greenTxt[] =   { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
const char yellowTxt[] =  { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
const char blueTxt[] =    { 0x1b, '[', '1', ';', '3', '4', 'm', 0 };
const char magentaTxt[] = { 0x1b, '[', '1', ';', '3', '5', 'm', 0 };
const char cyanTxt[] =    { 0x1b, '[', '1', ';', '3', '6', 'm', 0 };
const char normalTxt[] =  { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };

class AliAnalysisTaskSEITSsaSpectraMultiplicity : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEITSsaSpectraMultiplicity();
  virtual ~AliAnalysisTaskSEITSsaSpectraMultiplicity();
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);
  
  void SetMaxChi2Clu(Double_t chi=2.5){
    fMaxChi2Clu=chi;
  }
  void SetRapidityRange(Double_t dy=0.5){
    fMaxY=dy;
  }
  void SetMinNSigma(Double_t ns=1.5){
    fMinNSigma=ns;
  }
  void SetMindEdx(Double_t mind=0.){
    fMindEdx=mind;
  }
  void SetMinSPDPoints(Int_t np=1){
    fMinSPDPts=np;
  }
  void SetMinNdEdxSamples(Int_t np=3){
    fMinNdEdxSamples=np;
  }
  void SetDCACuts(Double_t nsxy=7., Double_t nsz=7.){
    fNSigmaDCAxy=nsxy;
    fNSigmaDCAz=nsz; 
  }
  void SetMultBin(Double_t LowBin=-1,Double_t UpBin=-1){
    fLowMult=LowBin;
    fUpMult=UpBin;
  }
  void SetCentralityCut(Float_t low, Float_t up){
    if((up>low)&&(!(low<0.0))&&(!(up>100.0))){
      fLowCentrality=low; fUpCentrality=up;
    }
  }
  void SetHImode(){fHImode=kTRUE;}
  void SetUseTrackMultiplicityEstimator(){fMultEstimator=0;}
  void SetUseTrackletsMultiplicityEstimator(){fMultEstimator=1;}
  void SetUseClustersSPD1MultiplicityEstimator(){fMultEstimator=2;}
  void SetPileupContributors(Int_t cont = 3){fPileupContributors = cont;}
  void SetPileupDistance(Double_t distance = 0.8 ){fPileupDistance = distance;}
  
  void SetEtaMax(Double_t maxeta){
    fEtaRange=maxeta;
  }
  
  void SetReadMC(Bool_t flag = kTRUE) {fMC = flag;}
  void SetFillTree(Bool_t fill=kTRUE) {fFillTree=fill;}
  void SetLowEnergypp(Bool_t opt=kTRUE) {fLowEnergypp=opt;}
  void SetV0Estimator(Bool_t flag = kTRUE) {fuseV0 = flag;}
  void SetRejectPileUp(Bool_t flag = kTRUE) {fPileupRej = flag;}
  
  void SetSmearMC(Double_t smearp, Double_t smeardedx){
    fSmearMC=kTRUE;
    fSmearP=smearp;
    fSmeardEdx=smeardedx;
  }
  
  Double_t CookdEdx(Double_t *s) const; 
  Double_t Eta2y(Double_t pt, Double_t m, Double_t eta) const;
  Bool_t DCAcut(Double_t impactXY, Double_t impactZ, Double_t pt) const;
  Bool_t DCAcutXY(Double_t impactXY, Double_t pt) const;
  Bool_t DCAcutZ(Double_t impactZ, Double_t pt) const;
  
  void PrintStatus(){
    Printf("Printing Status");
    cout<<"Using V0 Estimator "<<fuseV0<<endl;
    cout<<"Using MC "<<fMC<<endl;
    cout<<"Using fLowMult "<<fLowMult<<endl;
    cout<<"Using fUpMult "<<fUpMult<<endl;
  }
  
private:
  AliAnalysisTaskSEITSsaSpectraMultiplicity(const AliAnalysisTaskSEITSsaSpectraMultiplicity &source); 
  AliAnalysisTaskSEITSsaSpectraMultiplicity& operator=(const AliAnalysisTaskSEITSsaSpectraMultiplicity &source);
  
  AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
  AliPPVsMultUtils *fPPVsMultUtils; //V0 multiplicity estimator
  AliAnalysisUtils *fUtils;         // analysis utils (for MV pileup selection)
  
  enum {kNbins=22};//Number of pt bins
  enum {kIsNewEvent,kWithSPDVtx,kWithVtxCont,kVtxInRange,kVtxRes,kGoodZVtx,kPassPileUp,kMCVtxInRange,kbitEv=1*8};//Event selection bitmask
  enum {kIsPositive,kIsPion,kIsKaon,kIsProton,kIsPhysicalPrimary,kIsFromStrangeness,kIsFromMaterial};//MC information
  enum {kHasNoSelection=1,kIsITSsa,kIsITSrefit,kIsNotNeutralParticle,kHasOneSPD,kHasTwoSPD,kPassSPD,kHas3PIDcls,kHas4PIDcls,kPassPIDcls,kPassChi2Ncls,kIsInEta,kPassdEdx,kPassDCAzcut,kPassDCAxycut};//!track selection criteria
  AliESDEvent *fESD; //ESD object
  
  //////////////////////////////
  /////////////Lists////////////
  //////////////////////////////
  TList *fOutput; //! tlist with output
  TList *fListCuts; // list of functions storing DCA cut 
  
  //////////////////////////////
  /////////////Histos///////////
  //////////////////////////////
  TH1F *fHistNEvents; //! histo with number of events
  TH1F *fHistMCSampSel; //! histo for normalization to the number of events that passed sample selection
  TH1F *fHistMult; //! histo with multiplicity of the events or histo with V0 multiplicity of the events depending on fuseV0 flag
  TH1F *fHistMultAftEvSel; //! histo for normalization to the number of events that passed event selection and are in the multiplicity bin
  THnSparse *fHistEventMultiplicity; //! histo for normalization to the number of events that passed event selection and are in the multiplicity bin
    
  TH1F *fHistBefPileUp; //! histo with distribution of the difference between the primary vertex and the secondaries before pileup selection
  TH1F *fHistTaggedPileUp; //! histo with distribution of the difference between the primary vertex and the secondaries after pileup selection
  
  THnSparse *fHistSparse;//thsparse with multiplicity, pt, hipotesis particle asym, hipotesis particle sym, DCA , MC Hipotesis Prim;SecSt;SetMat on pt, MCpt
  THnSparse *fHistSparseBefEvSel;//thsparse with information from MC before event selection  multiplicity, MCpt, hipotesis particle (with primary Y and eta cut information) and event selection
  
  THnSparse *fHistPtCorr;//! Correlation between the pt reconstructed and the pt from MC
  
  TH1F *fHistCharge[4]; //! histo with charge distribution to check the calibration 
  
  TH1F *fHistCen; //! histo with multiplicity of the events
  TH1F *fHistNTracks; //! histo with number of tracks
  TH1F *fHistNTracksPos; //! histo with number of tracks
  TH1F *fHistNTracksNeg; //! histo with number of tracks
  
  TH2F *fHistDEDX; //! histo with dedx versus momentum
  TH2F *fHistDEDXdouble; //! histo with dedx versus signed momentum  
  THnSparse *fHistDEDXMulti; //! histo with dedx versus momentum versus multiplicity
  THnSparse *fHistDEDXdoubleMulti; //! histo with dedx versus signed momentum versus multiplicity
  TH2F *fHistPosNSigmaSep[3]; //! histo nsigma separation vs momentum
  TH2F *fHistNegNSigmaSep[3]; //! histo nsigma separation vs momentum
  
  TH1F *fHistBeforeEvSel; //! histo with pt distribution before my event selection
  TH1F *fHistAfterEvSel; //! histo with pt distribution after my event selection
  
  //Utility Histos
  TH1F *fHistITSchi2ncls;  //! histo with the chi2 / nclusters for each track
    
  //////////////////////////////
  /////////////MultiHistos//////
  //////////////////////////////
  
  TH1F *fHistPrimMCpos[3]; //! histo with spectra of primaries from the MC truth (positive)
  TH1F *fHistPrimMCneg[3]; //! histo with spectra of primaries from the MC truth (negative)  
  TH1F *fHistPrimMCposEtaY[3]; //! histo with spectra of primaries from the MC truth (positive)
  TH1F *fHistPrimMCnegEtaY[3]; //! histo with spectra of primaries from the MC truth (negative)
  TH1F *fHistSecStrMCpos[3]; //! histo with spectra of strange decays from the MC truth (positive)
  TH1F *fHistSecStrMCneg[3]; //! histo with spectra of strange decays from the MC truth (negative)
  TH1F *fHistSecMatMCpos[3]; //! histo with spectra of sec. from material from the MC truth (positive)
  TH1F *fHistSecMatMCneg[3]; //! histo with spectra of sec. from material from the MC truth (negative)
  //Before Event Selection
  TH1F *fHistPrimMCposBefEvSel[3]; //! histo with spectra of primaries from the MC truth (positive)
  TH1F *fHistPrimMCnegBefEvSel[3]; //! histo with spectra of primaries from the MC truth (negative)
  TH1F *fHistPrimMCposBefEvSelEtaY[3]; //! histo with spectra of primaries from the MC truth with both cut on y and eta(positive)
  TH1F *fHistPrimMCnegBefEvSelEtaY[3]; //! histo with spectra of primaries from the MC truth with both cut on y and eta(negative)
  TH1F *fHistPrimMCposBefEvSelEta[3]; //! histo with spectra of primaries from the MC truth with cut on eta only(positive)
  TH1F *fHistPrimMCnegBefEvSelEta[3]; //! histo with spectra of primaries from the MC truth with cut on eta only(negative)
  TH1F *fHistSecStrMCposBefEvSel[3]; //! histo with spectra of strange decays from the MC truth (positive)
  TH1F *fHistSecStrMCnegBefEvSel[3]; //! histo with spectra of strange decays from the MC truth (negative)
  TH1F *fHistSecMatMCposBefEvSel[3]; //! histo with spectra of sec. from material from the MC truth (positive)
  TH1F *fHistSecMatMCnegBefEvSel[3]; //! histo with spectra of sec. from material from the MC truth (negative)
  //Reconstructed with MC Truth
  TH1F *fHistPrimMCposReco[3]; //! histo with spectra of primaries from the MC truth (positive)
  TH1F *fHistPrimMCnegReco[3]; //! histo with spectra of primaries from the MC truth (negative)
  TH1F *fHistPrimMCposRecoEtaY[3]; //! histo with spectra of primaries from the MC truth with both cut on y and eta(positive)
  TH1F *fHistPrimMCnegRecoEtaY[3]; //! histo with spectra of primaries from the MC truth with both cut on y and eta(negative)
  TH1F *fHistSecStrMCposReco[3]; //! histo with spectra of strange decays from the MC truth (positive)
  TH1F *fHistSecStrMCnegReco[3]; //! histo with spectra of strange decays from the MC truth (negative)
  TH1F *fHistSecMatMCposReco[3]; //! histo with spectra of sec. from material from the MC truth (positive)
  TH1F *fHistSecMatMCnegReco[3]; //! histo with spectra of sec. from material from the MC truth (negative)
  
  //dEdx distributions 
  TH1F *fHistPosPi[kNbins]; //! histo with dedx distibution in the pions hypotesis (positive)
  TH1F *fHistPosK[kNbins]; //! histo with dedx distibution in the kaons hypotesis (positive)
  TH1F *fHistPosP[kNbins]; //! histo with dedx distibution in the protons hypotesis (positive)
  TH1F *fHistNegPi[kNbins]; //! histo with dedx distibution in the pions hypotesis (negative)
  TH1F *fHistNegK[kNbins]; //! histo with dedx distibution in the kaons hypotesis (negative)
  TH1F *fHistNegP[kNbins]; //! histo with dedx distibution in the protons hypotesis (negative)
  
  //DCA distributions
  TH1F *fHistDCAPosPi[kNbins]; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH1F *fHistDCAPosK[kNbins]; //! histo with DCA distibution in the kaons hypotesis (positive)
  TH1F *fHistDCAPosP[kNbins]; //! histo with DCA distibution in the protons hypotesis (positive)
  TH1F *fHistDCANegPi[kNbins]; //! histo with DCA distibution in the pions hypotesis (negative)
  TH1F *fHistDCANegK[kNbins]; //! histo with DCA distibution in the kaons hypotesis (negative)
  TH1F *fHistDCANegP[kNbins]; //! histo with DCA distibution in the protons hypotesis (negative)
  
  //DCA Templates
  TH1F *fHistMCPrimDCAPosPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCPrimDCAPosK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCPrimDCAPosP[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCPrimDCANegPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCPrimDCANegK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCPrimDCANegP[kNbins]; //! histo with DCA distibution, MC truth
  
  TH1F *fHistMCSecStDCAPosPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecStDCAPosK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecStDCAPosP[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecStDCANegPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecStDCANegK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecStDCANegP[kNbins]; //! histo with DCA distibution, MC truth
  
  TH1F *fHistMCSecMatDCAPosPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecMatDCAPosK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecMatDCAPosP[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecMatDCANegPi[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecMatDCANegK[kNbins]; //! histo with DCA distibution, MC truth 
  TH1F *fHistMCSecMatDCANegP[kNbins]; //! histo with DCA distibution, MC truth
  
  //dEdx distributions for MC
  TH1F *fHistMCPosOtherHypPion[kNbins]; //! histo with dedx using the MC truth 
  TH1F *fHistMCPosOtherHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosOtherHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosElHypPion[kNbins]; //! histo with dedx using the MC truth 
  TH1F *fHistMCPosElHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosElHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosPiHypPion[kNbins]; //! histo with dedx using the MC truth 
  TH1F *fHistMCPosPiHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosPiHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosKHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosKHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosKHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosPHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosPHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCPosPHypProton[kNbins]; //! histo with dedx using the MC truth
  
  TH1F *fHistMCNegOtherHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegOtherHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegOtherHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegElHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegElHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegElHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPiHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPiHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPiHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegKHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegKHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegKHypProton[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPHypPion[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPHypKaon[kNbins]; //! histo with dedx using the MC truth
  TH1F *fHistMCNegPHypProton[kNbins]; //! histo with dedx using the MC truth
  
  //Raw Spectra MEAN
  TH1F *fHistPosNSigmaMean[3];       //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaMCMean[3];       //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaPrimMean[3];   //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaPrimMCMean[3]; //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaMean[3];       //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaMCMean[3];       //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaPrimMean[3];   //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaPrimMCMean[3]; //! NSigma histos for 6 species
  
  //Raw Spectra
  TH1F *fHistPosNSigma[3];       //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaMC[3];       //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaPrim[3];   //! NSigma histos for 6 species
  TH1F *fHistPosNSigmaPrimMC[3]; //! NSigma histos for 6 species
  TH1F *fHistNegNSigma[3];       //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaMC[3];       //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaPrim[3];   //! NSigma histos for 6 species
  TH1F *fHistNegNSigmaPrimMC[3]; //! NSigma histos for 6 species
  
  //parameterizations functions
  TF1* fDCAxyCutFunc;  // function with DCAz cut vs. pt
  TF1* fDCAzCutFunc;   // function with DCAxy cut vs. pt
  
  //////////////////////////////
  /////////////Parameters///////
  //////////////////////////////
  Double_t fPtBinLimits[kNbins+1]; // limits of Pt Bins
  AliITSPIDResponse* fITSPIDResponse; //! class with BB parameterizations
  Int_t    fMinSPDPts;       // minimum number of SPD Points
  Int_t    fMinNdEdxSamples; // minimum number of SDD+SSD points
  Double_t fMindEdx;         // minimum dE/dx value in a layer (to cut noise)
  Double_t fMinNSigma;       // minimum number of sigmas
  Double_t fMaxY;            // maximum rapidity
  Double_t fMaxChi2Clu;      // maximum cluster
  Double_t fNSigmaDCAxy;     // DCA cut in bend. plane
  Double_t fNSigmaDCAz;      // DCA cut along z
  Double_t fEtaRange;        // limits in pseudorap
  Double_t fvZRange;         // limits in vertex Z position
  Double_t fLowMult;         // Lower Multiplicity bin
  Double_t fUpMult;          // Upper Multiplicity bin
  Float_t fLowCentrality;    //low Centrality cut
  Float_t fUpCentrality;     //up  Centrality cut
  Int_t fMultEstimator;      // multiplicty estimator 
  Bool_t fHImode;            //use spd2 as mulestimator 
  Int_t fPileupContributors; // Number of contributors requested for the event to be tagged as pileup
  Double_t fPileupDistance;  // Distance of the identified pileup vertex for the event to be tagged as pileup
  Int_t fYear;               // Year (2009, 2010)
  Bool_t   fMC;              //flag to switch on the MC analysis for the efficiency estimation
  Bool_t   fSmearMC;         // flag to apply extra smearing on MC 
  Double_t fSmearP;          // extra relative smearing on simulated momentum
  Double_t fSmeardEdx;       // extra relative smearing on simulated dE/dx
  TRandom3* fRandGener;      // generator for smearing
  Bool_t   fFillTree;        // fill ntuple  
  Bool_t   fLowEnergypp;     // flag for counrint ev. in p-p 2.76
  Bool_t   fuseV0;           // flag to use the V0 multiplicity estimator instead of the reference multiplicity
  Bool_t   fPileupRej;       // flag to reject pileup using the IsPileUpFromSPD
  Bool_t   fUseThnSparse;    // flag to use THnSparse
  TTree *fTreeNSigma;        //! output TRee
  TTree *fTreeMC;            //! output MC TRee
  
  //Values in fTReeNSigma
  Int_t   fMult;
  Float_t fV0Mult;
  Float_t fdedx;
  Float_t fP;
  Char_t  fTrkSign;
  UChar_t fClumap;
  Float_t fEta;
  UChar_t fevSelmask;
  Float_t fImpactXY;
  Float_t fImpactZ;
  //MC
  UChar_t fParticleMap;
  Float_t fptMC;
  Float_t fChi2;
  
  ClassDef(AliAnalysisTaskSEITSsaSpectraMultiplicity, 12);
};

#endif
