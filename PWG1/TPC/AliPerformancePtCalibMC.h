#ifndef ALIPERFORMANCEPTCALIBMC_H
#define ALIPERFORMANCEPTCALIBMC_H

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
class AliMCEvent;
class AliMCParticle;
class AliMCInfoCuts;
class AliRecInfoCuts;
class AliESDtrackCuts;

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

   Bool_t AddTPCcuts(AliESDtrack *ESDTrack);
   Bool_t AddITScuts(AliESDtrack *ESDTrack);
   Bool_t AddDCAcuts(AliESDtrack *ESDTrack);

   void SetReadTPCTracksOff(){fOptTPC   = kFALSE;}// function  can be used in the macro to set the variables
   void SetTPCRefitOn()      {fRefitTPC = kTRUE;} // function  can be used in the macro to set the variables
   void SetITSRefitOn()      {fRefitITS = kTRUE;} // function  can be used in the macro to set the variables
   void SetESDcutsOff()      {fESDcuts  = kFALSE;}// function  can be used in the macro to set the variables
   void SetDCAcutsOff()      {fDCAcut   = kFALSE;}// function  can be used in the macro to set the variables

   //void SetPtShift(Double_t shiftVal ){if(!(shiftVal==0)) {fShift=kTRUE; fDeltaInvP = shiftVal;} };
   void SetPtShift(Double_t shiftVal); //shift in 1/pt

   void SetESDcutValues(Double_t * esdCutValues){// functions can be used in the macro to set the variables
      fMinPt                = esdCutValues[0]; 
      fMaxPt                = esdCutValues[1];
      fMinNClustersTPC      = esdCutValues[2];
      fMaxChi2PerClusterTPC = esdCutValues[3];
      fMaxDCAtoVertexXY     = esdCutValues[4];
      fMaxDCAtoVertexZ      = esdCutValues[5];
   }
   
   void SetProjBinsPhi(const Double_t *pBins,Int_t sizep);
   void SetProjBinsTheta(const Double_t *tBins, Int_t sizet);
   void SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR );
   void SetAnaMCOff() {fAnaMC = kFALSE;}
   TList *GetHistoList() {return fList;}
   
 
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
   
  Double_t fThetaBins[100];
  Double_t fPhiBins[100];
   
   Int_t fNThetaBins;
   Int_t fNPhiBins ;
   Double_t fRange;
   Double_t fExclRange ;
   Bool_t fFitGaus ;
   Bool_t  fAnaMC;
   
    
 private:
  Bool_t fShift;
  Double_t fDeltaInvP;
  //options for cuts
  Bool_t fOptTPC;
  Bool_t fESDcuts;
  Bool_t fRefitTPC;
  Bool_t fRefitITS;
  Bool_t fDCAcut;
  
  Double_t fMinPt;
  Double_t fMaxPt;
  Double_t fMinNClustersTPC;
  Double_t fMaxChi2PerClusterTPC;
  Double_t fMaxDCAtoVertexXY;
  Double_t fMaxDCAtoVertexZ;

  AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks

   
  TList       *fList;
  TH2F        *fHistInvPtTheta;
  TH2F        *fHistInvPtPhi;
  TH2F        *fHistPtTheta;
  TH2F        *fHistPtPhi;

  TH1F        *fHistPtShift0;
  TH1F        *fHistPrimaryVertexPosX;          
  TH1F        *fHistPrimaryVertexPosY;          
  TH1F        *fHistPrimaryVertexPosZ;          
  TH1F        *fHistTrackMultiplicity;          
  TH1F        *fHistTrackMultiplicityCuts;

  TH2F        *fHistTPCMomentaPosP;
  TH2F        *fHistTPCMomentaNegP;
  TH2F        *fHistTPCMomentaPosPt;
  TH2F        *fHistTPCMomentaNegPt;
   
  TH2F        *fHistInvPtThetaMC;
  TH2F        *fHistInvPtPhiMC;
  TH2F        *fHistPtThetaMC;
  TH2F        *fHistPtPhiMC;   
  TH2F        *fHistInvPtMCESD;
  TH2F        *fHistInvPtMCTPC;
  TH2F        *fHistPtMCESD;
  TH2F        *fHistPtMCTPC;
  TH2F        *fHistMomresMCESD;   
  TH2F        *fHistMomresMCTPC;   
  TH2F        *fHistTPCMomentaPosInvPtMC;
  TH2F        *fHistTPCMomentaNegInvPtMC;
  TH2F        *fHistTPCMomentaPosPtMC;
  TH2F        *fHistTPCMomentaNegPtMC;
  
  TH1F        *fHistESDMomentaPosInvPtMC;
  TH1F        *fHistESDMomentaNegInvPtMC;
  TH1F        *fHistESDMomentaPosPtMC;
   TH1F        *fHistESDMomentaNegPtMC;
   TH1F        *fHistUserPtShift;
   
  AliESDtrackCuts* fESDTrackCuts;
  
  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformancePtCalibMC(const AliPerformancePtCalibMC&);            // not implemented 
  AliPerformancePtCalibMC& operator=(const AliPerformancePtCalibMC&); // not implemented 

  ClassDef(AliPerformancePtCalibMC, 1); 
};

#endif
