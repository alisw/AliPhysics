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
class TH3F;
class TList;

class AliESDVertex;
class AliESDtrack;
class AliMCEvent;
class AliStack;
class AliTrackReference;
class AliESDEvent; 
class AliESDfriend; 
class AliESDfriendTrack; 
class AliMCEvent;
class AliMCParticle;
class AliMCInfoCuts;
class AliRecInfoCuts;
class AliESDtrackCuts;

#include "AliPerformanceObject.h"

class AliPerformancePtCalibMC : public AliPerformanceObject {
public:
   AliPerformancePtCalibMC();
   AliPerformancePtCalibMC(const char *name, const char *title);//, Int_t analysisMode, Bool_t hptGenerator);
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
   Bool_t AddTPCcuts(const AliESDtrack *esdTrack);// applies TPC cuts
   Bool_t AddITScuts(const AliESDtrack *esdTrack);// applies ITS cuts
   Bool_t AddDCAcuts(const AliESDtrack *esdTrack);// applies DCA cuts
 
   void SetReadTPCTracks(const Bool_t readTPC)        {fOptTPC   = readTPC;}//read only ESD tracks
   void SetTPCRefit(const Bool_t refitTPC)            {fRefitTPC = refitTPC;} //switch TPC refit flag on/off
   void SetITSRefit(const Bool_t refitITS)            {fRefitITS = refitITS;} //switch ITS refit flag on/off
   void SetESDCuts(const Bool_t esdCuts)              {fESDcuts  = esdCuts;} //switch ESD track cuts on/off
   void SetDCACuts(const Bool_t dcaCut)               {fDCAcut   = dcaCut;} //switch DCA cut off
   void SetAcceptKinkDaughters(const Bool_t kink)     {fAcceptKinkDaughters = kink;} //switch accept kink daughters on/off
   void SetRequireSigmaToVertex(const Bool_t sigmaTo) {fRequireSigmaToVertex = sigmaTo;}//switch require SigmaToVertex on/off
   void SetfDCAToVertex2D(const Bool_t dcaTo)         {fDCAToVertex2D = dcaTo;}//switch DCA to vertex2D on/off
   void SetEtaRange(const Double_t eta)               {fEtaAcceptance =  eta ;}//sets eta window
   void SetESDcutValues(const Double_t * esdCutValues);// set ESD track cut values as array of size 6 according to:
   //    fMinPt                = esdCutValues[0]; 
   //    fMaxPt                = esdCutValues[1];
   //    fMinNClustersTPC      = esdCutValues[2];
   //    fMaxChi2PerClusterTPC = esdCutValues[3];
   //    fMaxDCAtoVertexXY     = esdCutValues[4];
   //    fMaxDCAtoVertexZ      = esdCutValues[5];

   //user defined shift in charge/pt
   void SetPtShift(const Double_t shiftVal); //shift in 1/pt
   

   // for analysis with AliPerformancePtCalibMC::Analyse()
   void SetProjBinsPhi(const Double_t *pBins,Int_t sizep);// set phi bins and nr of phi bins for projection
   void SetProjBinsTheta(const Double_t *tBins, Int_t sizet);//  set theta bins and nr of phi bins for projection
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );// set fit options
   void SetAnaMCOff() {fAnaMC = kFALSE;} // switch analysis of MC tracks off
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
   
   Double_t fThetaBins[100];// array of theta bins for projection of 1/pt vs theta
   Double_t fPhiBins[100]; // array of phi bins for projection of 1/pt vs theta
   
   Int_t fNThetaBins;// sets number of theta bins
   Int_t fNPhiBins ;// sets number of phi bins
   Double_t fRange;// sets fit range
   Double_t fExclRange ;// sets range of rejection of points around 0
   Bool_t fFitGaus ;// flag for usage of gaussian fit function
   Bool_t  fAnaMC;// flag for analysis of MC tracks
   
    
private:
   Bool_t fShift;//flag for shift in charge/pt
   Double_t fDeltaInvP;// shift value of charge/pt
   
   //options for cuts
   Bool_t fOptTPC;// flag for reading of TPC tracks in Exec
   Bool_t fESDcuts;//flag for usage of esd track cuts
   Bool_t fRefitTPC;//flag for TPC refit
   Bool_t fRefitITS;// flag for ITS refit
   Bool_t fDCAcut;//flag for usage of DCA cut
   
   Double_t fEtaAcceptance;//sets value of eta window
   Double_t fMinPt;//sets minimum pt for esd track cuts
   Double_t fMaxPt;//sets maximum pt for esd track cuts
   Double_t fMinNClustersTPC;// set minimum number of clusters in TPC for esd track cuts
   Double_t fMaxChi2PerClusterTPC;//set maximum of chi2 per cluster in TPC for esd track cuts
   Double_t fMaxDCAtoVertexXY;//set maximum of dca to vertex in xy direction for esd track cuts
   Double_t fMaxDCAtoVertexZ;//set maximum of dca to vertex in z for esd track cuts

   Bool_t fAcceptKinkDaughters;// flag for acception of kink daughters
   Bool_t fRequireSigmaToVertex;// flag for requirering sigma to vertex
   Bool_t fDCAToVertex2D;//flag for dca to vertex in 2d cut
   
   AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
   AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks

   
   TList       *fList;// list of histograms
   TH2F        *fHistInvPtTheta;//theta vs charge/pt
   TH2F        *fHistInvPtPhi;//phi vs charge/pt
   TH2F        *fHistPtTheta;//theta vs pt
   TH2F        *fHistPtPhi;//phi vs pt

   TH1F        *fHistPtShift0;//if shift in 1/pt is set by user, this histogram shows pt wihtout shift
   TH1F        *fHistPrimaryVertexPosX;// primary vertex position x          
   TH1F        *fHistPrimaryVertexPosY;// primary vertex position y        
   TH1F        *fHistPrimaryVertexPosZ;// primary vertex position z        
   TH1F        *fHistTrackMultiplicity; // track multiplicity         
   TH1F        *fHistTrackMultiplicityCuts;//track multiplicity after all cuts are applied

   TH2F        *fHistTPCMomentaPosP;//TPC p vs global esd track p for positive tracks
   TH2F        *fHistTPCMomentaNegP;//TPC p vs global esd track p for negative tracks
   TH2F        *fHistTPCMomentaPosPt;//TPC pt vs global esd track p positive tracks
   TH2F        *fHistTPCMomentaNegPt;//TPC pt vs global esd track p for negative tracks
   
   TH2F        *fHistInvPtThetaMC;//theta vs charge/pt for MC tracks
   TH2F        *fHistInvPtPhiMC;//phi vs charge/pt for MC tracks
   TH2F        *fHistPtThetaMC;//theta vs pt for MC tracks
   TH2F        *fHistPtPhiMC; //phi vs pt for MC tracks
   TH2F        *fHistInvPtMCESD;// 1/pt of ESD vs MC
   TH2F        *fHistInvPtMCTPC;// 1/pt of TPC vs MC
   TH2F        *fHistPtMCESD;//pt of ESD vs MC
   TH2F        *fHistPtMCTPC;//pt of TPC vs MC
   TH2F        *fHistMomresMCESD;   //(pt ESD - pt MC)/ptMC vs pt MC
   TH2F        *fHistMomresMCTPC;   //(pt TPC - pt MC)/ptMC vs pt MC
   TH2F        *fHistTPCMomentaPosInvPtMC;//TPC-MC of 1/pt vs global ESD-MC of 1/pt of positive tracks
   TH2F        *fHistTPCMomentaNegInvPtMC;//TPC-MC of 1/pt vs global ESD-MC of 1/pt of negative tracks
   TH2F        *fHistTPCMomentaPosPtMC;//TPC-MC of pt vs global ESD-MC of pt of positive tracks
   TH2F        *fHistTPCMomentaNegPtMC;//TPC-MC of pt vs global ESD-MC of pt of negative tracks
  
   TH1F        *fHistESDMomentaPosInvPtMC;//ESD-MC of 1/pt of positive tracks
   TH1F        *fHistESDMomentaNegInvPtMC;//ESD-MC of 1/pt of negative tracks
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
