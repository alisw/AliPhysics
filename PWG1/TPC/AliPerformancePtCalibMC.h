#ifndef ALIPERFORMANCEPTCALIBMC_H
#define ALIPERFORMANCEPTCALIBMC_H
//----------------------------------------------------------------------------------------------------
// Class to study systematic shifts in pt and charge/pt respectively. Furthermore a comparison between
// either ESD or TPC and MC track momenta is included.
// Track cuts and a user defined shift in 1/pt can be switched on and off by user.
//
// Analysis with class AliPerfAnalyzeInvPt via AliPerformancePtCalibMC::Analyse().:
// Projection of 1/pt vs theta and vs phi resp. Histograms will be fitted with either
// polynomial or gaussian fit function to extract minimum position of 1/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set* functions of AliPerformancePtCalibMC when running AliPerformancePtCalibMC::Analyse().
//
// Author: S. Schuchmann 11/13/2009 
//----------------------------------------------------------------------------------------------------

class TString;
class TNamed;
class TCanvas;
class TH1F;
class TH2F;
class TList;

class AliESDVertex;
class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliTrackReference;
class AliESDEvent; 
class AliESDfriend; 
class AliESDfriendTrack; 
class AliMCParticle;
class AliMCInfoCuts;
class AliRecInfoCuts;
class AliESDtrackCuts;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformancePtCalibMC : public AliPerformanceObject {
public:
   AliPerformancePtCalibMC();
   AliPerformancePtCalibMC(const char *name, const char *title);
   virtual ~AliPerformancePtCalibMC() ;

   // Init data members
   virtual void  Init();

   // Execute analysis
   virtual void  Exec(AliMCEvent* const mcEvent, AliESDEvent *const esdEvent, AliESDfriend *const esdFriend, const Bool_t bUseMC, const Bool_t bUseESDfriend);

   // Merge output objects (needed by PROOF) 
   virtual Long64_t Merge(TCollection* const list);

   // Analyse output histograms
   virtual void Analyse();

   // Get analysis folder
   virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

   // Options for track cuts
   
   void SetReadTPCTracks(const Bool_t readTPC)        {fOptTPC   = readTPC;}//read only ESD tracks
   void SetEtaRange(const Double_t eta)               {fEtaAcceptance =  eta ;}//sets eta window
   
   void SetAliESDtrackCuts( AliESDtrackCuts* esdTrackCuts) { fESDTrackCuts = esdTrackCuts;fESDcuts=kTRUE;}//neu
    
   //user defined shift in charge/pt
   void SetPtShift(const Double_t shiftVal); //sets user defined shift in charge/pt
   
   // setters for analysis with AliPerformancePtCalibMC::Analyse()
   void SetProjBinsPhi(const Double_t *pBins,const Int_t sizep,const Double_t minTheta, const Double_t maxTheta);// set phi bins for projection and theta range selection (rad)
   void SetProjBinsTheta(const Double_t *tBins, const Int_t sizet,const Double_t minPhi, const Double_t maxPhi);// set theta bins for projection and phi range selection (rad)
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );// set fit options
   void SetDoRebin(const Int_t rebin){if(rebin) {fDoRebin = kTRUE; fRebin = rebin;}}
   void SetAnaMCOff() {fAnaMC = kFALSE;} // switch analysis of MC true tracks off
   const TList *GetHistoList() {return fList;} // get list of histograms
   
 
   // Create folder for analysed histograms
   TFolder *CreateFolder(TString folder = "folderPtCalib",TString title = "Analysed PtCalib histograms");

   // Export objects to folder
   TFolder *ExportToFolder(TObjArray * array=0);

   // Selection cuts
   void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}   
   void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}
    
   AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
   AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

protected:
   // variables for fitting in Analyse() function
   Double_t fThetaBins[100];// array of theta bins for projection of charge/pt vs theta
   Double_t fPhiBins[100]; // array of phi bins for projection of charge/pt vs theta

   Int_t fNThetaBins;// sets number of theta bins
   Int_t fNPhiBins ;// sets number of phi bins
   Double_t fMaxPhi;// max phi for 2D projection on theta and charge/pt axis
   Double_t fMinPhi;// min phi for 2D projection on theta and charge/pt axis
   Double_t fMaxTheta;// max theta for 2D projection on phi and charge/pt axis
   Double_t fMinTheta;// min theta for 2D projection on phi and charge/pt axis
   Double_t fRange;// sets fit range
   Double_t fExclRange ;// sets range of rejection of points around 0
   Bool_t fFitGaus ;// flag for usage of gaussian fit function
   Bool_t fDoRebin;// flag for rebin 1D histos before fitting
   Int_t fRebin;// number of bins for rebin
   Bool_t  fAnaMC;// flag for analysis of MC tracks
   
    
private:
   // option for user defined shift in charge/pt
   Bool_t fShift;//flag for shift in charge/pt
   Double_t fDeltaInvP;// shift value of charge/pt
   
   //options for cuts
   Bool_t fOptTPC;// flag for reading of TPC tracks in Exec
   Bool_t fESDcuts;//flag for usage of esd track cuts
   
   //ESD track cut values
   Double_t fEtaAcceptance;//sets value of eta window
   AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
   AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks

   
   TList       *fList;// list of histograms
   
   // histograms and THnSparse
   THnSparseF  *fHistInvPtPtThetaPhi;// is filled with charge/pt, pt, theta, phi for ESD or TPC

   TH1F        *fHistPtShift0;//if shift in charge/pt is set by user, this histogram shows pt wihtout shift
   TH1F        *fHistPrimaryVertexPosX;// primary vertex position x          
   TH1F        *fHistPrimaryVertexPosY;// primary vertex position y        
   TH1F        *fHistPrimaryVertexPosZ;// primary vertex position z        
   TH1F        *fHistTrackMultiplicity; // track multiplicity         
   TH1F        *fHistTrackMultiplicityCuts;//track multiplicity after all cuts are applied

   TH2F        *fHistTPCMomentaPosP;//TPC p vs global esd track p for positive tracks
   TH2F        *fHistTPCMomentaNegP;//TPC p vs global esd track p for negative tracks
   TH2F        *fHistTPCMomentaPosPt;//TPC pt vs global esd track p positive tracks
   TH2F        *fHistTPCMomentaNegPt;//TPC pt vs global esd track p for negative tracks

   THnSparseF *fHistInvPtPtThetaPhiMC;// is filled with charge/pt, pt, theta, phi for MC true

   TH2F        *fHistInvPtMCESD;// charge/pt of ESD vs MC
   TH2F        *fHistInvPtMCTPC;// charge/pt of TPC vs MC
   TH2F        *fHistPtMCESD;//pt of ESD vs MC
   TH2F        *fHistPtMCTPC;//pt of TPC vs MC
   TH2F        *fHistMomresMCESD;   //(pt ESD - pt MC)/ptMC vs pt MC
   TH2F        *fHistMomresMCTPC;   //(pt TPC - pt MC)/ptMC vs pt MC
   TH2F        *fHistTPCMomentaPosInvPtMC;//TPC-MC of 1/pt vs global ESD-MC of charge/pt of positive tracks
   TH2F        *fHistTPCMomentaNegInvPtMC;//TPC-MC of 1/pt vs global ESD-MC of charge/pt of negative tracks
   TH2F        *fHistTPCMomentaPosPtMC;//TPC-MC of pt vs global ESD-MC of pt of positive tracks
   TH2F        *fHistTPCMomentaNegPtMC;//TPC-MC of pt vs global ESD-MC of pt of negative tracks
  
   TH1F        *fHistESDMomentaPosInvPtMC;//ESD-MC of charge/pt of positive tracks
   TH1F        *fHistESDMomentaNegInvPtMC;//ESD-MC of charge/pt of negative tracks
   TH1F        *fHistESDMomentaPosPtMC;//ESD-MC of pt of positive tracks
   TH1F        *fHistESDMomentaNegPtMC;//ESD-MC of pt of negative tracks
   
   TH1F        *fHistUserPtShift;// shows the shift value if set by user
   
   AliESDtrackCuts* fESDTrackCuts;// esd track cuts
  
   // analysis folder 
   TFolder *fAnalysisFolder; // folder for analysed histograms

   AliPerformancePtCalibMC(const AliPerformancePtCalibMC&);            // not implemented 
   AliPerformancePtCalibMC& operator=(const AliPerformancePtCalibMC&); // not implemented 

   ClassDef(AliPerformancePtCalibMC, 1); 
};

#endif
