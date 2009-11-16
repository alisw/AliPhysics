
#ifndef ALIPERFORMANCEPTCALIB_H
#define ALIPERFORMANCEPTCALIB_H
//----------------------------------------------------------------------------------------------------
// Class to study systematic shifts in pt and charge/pt respectively.Furthermore a comparison between
// ESD and TPC track momenta is included.
// Track cuts and a user defined shift in 1/pt can be switched on and off by user.
//
// Analysis with class AliPerfAnalyzeInvPt via AliPerformancePtCalib::Analyse(). :
// Projection of charge/pt vs theta and vs phi resp. Histograms will be fitted with either
// polynomial or gaussian fit function to extract minimum position of 1/pt.
// Fit options and theta, phi bins can be set by user.
// Attention: use the Set* functions of AliPerformancePtCalib when running AliPerformancePtCalib::Analyse().
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
class AliAnalysisTask;

#include "AliPerformanceObject.h"

class AliPerformancePtCalib : public AliPerformanceObject {
public:
   AliPerformancePtCalib(); 
   AliPerformancePtCalib(Char_t* name, Char_t* title);//, Int_t analysisMode, Bool_t hptGenerator);
   virtual ~AliPerformancePtCalib();

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
   
   const Bool_t AddTPCcuts(const AliESDtrack *esdTrack);// applies TPC cuts
   const Bool_t AddITScuts(const AliESDtrack *esdTrack);// applies ITS cuts
   const Bool_t AddDCAcuts(const AliESDtrack *esdTrack);// applies DCA cuts
 
   void SetReadTPCTracks(const Bool_t readTPC)        {fOptTPC   = readTPC;}//read only ESD tracks
   void SetTPCRefit(const Bool_t refitTPC)            {fRefitTPC = refitTPC;} //switch TPC refit flag on/off
   void SetITSRefit(const Bool_t refitITS)            {fRefitITS = refitITS;} //switch ITS refit flag on/off
   void SetESDCuts(const Bool_t esdCuts)              {fESDcuts  = esdCuts;} //switch ESD track cuts on/off
   void SetDCACuts(const Bool_t dcaCut)               {fDCAcut   = dcaCut;} //switch DCA cut off
   void SetAcceptKinkDaughters(const Bool_t kink)     {fAcceptKinkDaughters = kink;} //switch accept kink daughters on/off
   void SetRequireSigmaToVertex(const Bool_t sigmaTo) {fRequireSigmaToVertex = sigmaTo;}//switch require SigmaToVertex on/off
   void SetfDCAToVertex2D(const Bool_t dcaTo)         {fDCAToVertex2D = dcaTo;}//switch DCA to vertex2D on/off
   void SetEtaRange(const Double_t eta)               {fEtaAcceptance =  eta ;}//sets eta window
   void SetESDcutValues(const Double_t * esdCutValues);// set ESD cut values as an array with size 6
   //  fMinPt                = esdCutValues[0]; 
   //  fMaxPt                = esdCutValues[1];
   //  fMinNClustersTPC      = esdCutValues[2];
   //  fMaxChi2PerClusterTPC = esdCutValues[3];
   //  fMaxDCAtoVertexXY     = esdCutValues[4];
   //  fMaxDCAtoVertexZ      = esdCutValues[5];

   //user defined shift in charge/pt
   void SetPtShift(const Double_t shiftVal); // set user defined shift in charge/pt
   
   // for analysis with AliPerformancePtCalibMC::Analyse()  
   void SetProjBinsPhi(const Double_t *pBins,const Int_t sizep);// set phi bins for projection
   void SetProjBinsTheta(const Double_t *tBins, const Int_t sizet);// set theta bins for projection
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );//set fit options
   const TList *GetHistoList() {return fList;}// get list of histograms for analysis
   
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

   Bool_t fAcceptKinkDaughters; // flag for acception of kink daughters
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
   

   TH1F        *fHistUserPtShift;// shows the shift value if set by user
   
   AliESDtrackCuts* fESDTrackCuts;// esd track cuts
  
   // analysis folder 
   TFolder *fAnalysisFolder; // folder for analysed histograms

   AliPerformancePtCalib(const AliPerformancePtCalib&);            // not implemented 
   AliPerformancePtCalib& operator=(const AliPerformancePtCalib&); // not implemented 

   ClassDef(AliPerformancePtCalib, 1); 
};

#endif
