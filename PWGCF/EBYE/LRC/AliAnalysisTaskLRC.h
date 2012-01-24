#ifndef ALIANALYSISTASKLRC_H
#define ALIANALYSISTASKLRC_H

// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCProcess objects that are processing LRC analysis
// for a given Eta window 

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

/*  See cxx source for full Copyright notice */

#include <AliAnalysisTaskSE.h>
class AliLRCProcess;
class AliLRCBase;
class AliESDtrackCuts;
class TH1D;
class TH2D;

class AliAnalysisTaskLRC : public AliAnalysisTaskSE {

public:
   //Constructors 
  AliAnalysisTaskLRC(const char *name = "AliAnalysisTaskLRC",Bool_t runKine=kFALSE);
  virtual ~AliAnalysisTaskLRC() {}
  
  //AliAnalysisTaskSE overloading
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //----------------------------------  
  
  void AddLRCProcess(AliLRCBase *newProc); //Adds new AliLRCProcess to analysis task
  
  // Setters
  
  void SetMaxPtLimit(Double_t MaxPtLimit);   //Sets  Max Pt filter
  void SetMinPtLimit(Double_t MinPtLimit);   //Sets  Min Pt filter 
  void SetCheckForkVtx(Bool_t CheckForkVtx){fCheckForkVtx=CheckForkVtx;} // Accept only events with veretex
  void SetCheckForVtxPosition(Bool_t CheckForVtxPosition ){fCheckForVtxPosition=CheckForVtxPosition;} //Accept only events with veretex in slected range
  void SetTrackCuts(AliESDtrackCuts* const cuts)  { fEsdTrackCuts = cuts; }
  void SetShowEventStats(Bool_t ShowEventStats)  {fShowEventStats= ShowEventStats;}
  void SetShowPerTrackStats(Bool_t ShowPerTrackStats) {fShowPerTrackStats=ShowPerTrackStats;}
  void SetVtxDiamond(Double_t Vx, Double_t Vy, Double_t Vz) {fVxMax = Vx;fVyMax =Vy;fVzMax = Vz;}
  void SetNchCuts(Int_t minNch, Int_t maxNch){fMinAceptedTracksCut=minNch; fMaxAceptedTracksCut=maxNch;}
  
// Getters
  TList* GetListOfProcessors() { return &fLRCproc;} // Returns list of included 
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  AliLRCBase * Proc(Int_t index);// Get Processor i 
  
  
  
protected:
// Track cuts

  AliESDtrackCuts *fEsdTrackCuts;               // esd track cuts



// Aceptance cuts

  Double_t fMaxPtLimit;  //Max Pt filter
  Double_t fMinPtLimit;  // Min Pt filter 

// Nch cuts
  Int_t fMinAceptedTracksCut;   //Minimum number of accepted tracks in event
  Int_t fMaxAceptedTracksCut;   //Maximum number of accepted tracks in event
  
// Vtx cuts
  Bool_t fCheckForkVtx;		// Check for vertex
  Bool_t fCheckForVtxPosition;  // Check if vertex position in range
  Double_t fVxMax;	// X vrtx max
  Double_t fVyMax;	// Y vrtx max
  Double_t fVzMax;	// Z vrtx max


  TList fLRCproc;       //  AliLRCProcess objects list
  TList* fOutList;      //! Task Output data container 
  
  Bool_t fRunKine;      // ESD/AOD  - KINE switch
  Bool_t fShowEventStats; //  Allows per event debug output (trigger Nch, cuts etc)
  Bool_t fShowPerTrackStats; // Allows per track debug output
      

// QA histos 

  TH1D *fHistEventCutStats;  //! Event cut statistics
  TH1D *fHistTrackCutStats;  //! Track cut statistics

  TH1D *fHistVx;  //!Vx hist
  TH1D *fHistVy;  //!Vy hist
  TH1D *fHistVz;  //!Vz hist

  TH2D *fHistEtaVsZvCoverage; //! Statistics on tracks Zv and Eta for all tracks
  TH2D *fHistEtaVsZvCoverageAcepted; //!  Statistics on tracks Zv and Eta for acepted tracks

  TH1D *fHistAceptedMult;   //! Number of acepted tracks histo



  AliAnalysisTaskLRC(const AliAnalysisTaskLRC&); // not implemented
  AliAnalysisTaskLRC& operator=(const AliAnalysisTaskLRC&); // not implemented
  
  ClassDef(AliAnalysisTaskLRC, 1); 
};

#endif
