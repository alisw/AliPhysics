
#ifndef ALIPERFORMANCEPTCALIB_H
#define ALIPERFORMANCEPTCALIB_H
//----------------------------------------------------------------------------------------------------
// Class to study systematic shifts in pt and charge/pt respectively.Furthermore a comparison between
// ESD and TPC track momenta is included.
// Track cuts and a user defined shift in charge/pt can be switched on and off by user.
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
class TList;

class AliESDVertex;
class AliESDtrack;
class AliMCEvent;
class AliTrackReference;
class AliESDEvent; 
class AliESDfriend; 
class AliESDfriendTrack; 
class AliMCInfoCuts;
class AliRecInfoCuts;
class AliESDtrackCuts;
class AliESDpid;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformancePtCalib : public AliPerformanceObject {
public:
  AliPerformancePtCalib(const Char_t * name="AliPerformancePtCalib",const Char_t* title ="AliPerformancePtCalib");
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
 
   void SetReadTPCTracks(const Bool_t readTPC)        {fOptTPC   = readTPC;}//read only ESD tracks
   void SetEtaRange(const Double_t eta)               {fEtaAcceptance =  eta ;}//sets eta window
  
   void SetAliESDtrackCuts( AliESDtrackCuts* esdTrackCuts) { fESDTrackCuts = esdTrackCuts;fESDcuts=kTRUE;}//esd track cuts

   void SetAnalysePions(const Bool_t anaPions) {fPions = anaPions;}
   void SetPtShift(const Double_t shiftVal); // set user defined shift in charge/pt

   // setters for analysis with AliPerformancePtCalib::Analyse()  
   void SetProjBinsPhi(const Double_t *pBins,const Int_t sizep,const Double_t minTheta, const Double_t maxTheta);// set phi bins for projection and theta range selection (rad)
   void SetProjBinsTheta(const Double_t *tBins, const Int_t sizet,const Double_t minPhi, const Double_t maxPhi);// set theta bins for projection and phi range selection (rad)
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );//set fit options
   void SetDoRebin(const Int_t rebin){if(rebin) {fDoRebin = kTRUE; fRebin = rebin;}}
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
   // variables for fitting in Analyse() function
   Double_t fThetaBins[100];// array of theta bins for projection of 1/pt vs theta
   Double_t fPhiBins[100]; // array of phi bins for projection of 1/pt vs theta

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
    
private:
   // option for user defined shift in charge/pt
   Bool_t fShift;//flag for shift in charge/pt
   Double_t fDeltaInvP;// shift value of charge/pt
   
   //options for cuts
   Bool_t fOptTPC;// flag for reading of TPC tracks in Exec
   Bool_t fESDcuts;//flag for usage of esd track cuts
   Bool_t fPions;// flag for analzsing pions instead of all charged particles

   //ESD track cut values
   Double_t fEtaAcceptance;//sets value of eta window

   AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
   AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks
  

    
   TList       *fList;// list of histograms
   
   //histograms and THnSparse
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
   

   TH1F        *fHistUserPtShift;// shows the shift value if set by user
   TH2F        *fHistdedxPions;// dEdx vs cahrge*pt
   
   AliESDtrackCuts* fESDTrackCuts;// esd track cuts
   //pid
   AliESDpid *fESDpid;
   // analysis folder 
   TFolder *fAnalysisFolder; // folder for analysed histograms

   AliPerformancePtCalib(const AliPerformancePtCalib&);            // not implemented 
   AliPerformancePtCalib& operator=(const AliPerformancePtCalib&); // not implemented 

   ClassDef(AliPerformancePtCalib, 1); 
};

#endif
