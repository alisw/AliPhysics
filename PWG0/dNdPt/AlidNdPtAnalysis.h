#ifndef ALIDEDPTANALYSIS_H
#define ALIDEDPTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtAnalysis class. 
// 
// a. functionality:
// - fills analysis control histograms
// - fills generic correction matrices 
// - generates correction matrices 
//
// b. data members:
// - generic correction matrices
// - control histograms
//
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------

class TProfile;
class TFolder;
class TObjArray;
class TString;

class AliESDtrackCuts;
class AliVertexerTracks;
class AliESD;
class AliESDfriend;
class AliESDfriendTrack;

#include "THnSparse.h"
#include "AlidNdPt.h"
#include "AlidNdPtHelper.h"

class AlidNdPtAnalysis : public AlidNdPt {
public :
  AlidNdPtAnalysis(); 
  AlidNdPtAnalysis(Char_t* name, Char_t* title);
  ~AlidNdPtAnalysis();

  // Init data members
  virtual void Init();

  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms 
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * array=0);

  // Get analysis folder
  TFolder* GetAnalysisFolder() {return fAnalysisFolder;}

  // Fill control histograms
  void SetHistogramsOn(Bool_t histOn=kTRUE) {fHistogramsOn = histOn;}
  Bool_t IsHistogramsOn() {return fHistogramsOn;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj);
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj);
  void FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec);

  // Getters
  THnSparseF *GetEventMultCorrelationMatrix() {return fEventMultCorrelationMatrix;}
  THnSparseF *GetTrackPtCorrelationMatrix()   {return fTrackPtCorrelationMatrix;}

  //
  THnSparseF *GetGenEventMatrix() {return fGenEventMatrix;}
  THnSparseF *GetGenSDEventMatrix() {return fGenSDEventMatrix;}
  THnSparseF *GetGenDDEventMatrix() {return fGenDDEventMatrix;}
  THnSparseF *GetGenNDEventMatrix() {return fGenNDEventMatrix;}
  THnSparseF *GetGenNSDEventMatrix() {return fGenNSDEventMatrix;}

  THnSparseF *GetTriggerEventMatrix() {return fTriggerEventMatrix;}
  THnSparseF *GetTriggerSDEventMatrix() {return fTriggerSDEventMatrix;}
  THnSparseF *GetTriggerDDEventMatrix() {return fTriggerDDEventMatrix;}
  THnSparseF *GetTriggerNDEventMatrix() {return fTriggerNDEventMatrix;}
  THnSparseF *GetTriggerNSDEventMatrix() {return fTriggerNSDEventMatrix;}

  THnSparseF *GetRecEventMatrix() {return fRecEventMatrix;}
  THnSparseF *GetRecSDEventMatrix() {return fRecSDEventMatrix;}
  THnSparseF *GetRecDDEventMatrix() {return fRecDDEventMatrix;}
  THnSparseF *GetRecNDEventMatrix() {return fRecNDEventMatrix;}
  THnSparseF *GetRecNSDEventMatrix() {return fRecNSDEventMatrix;}

  // 
  THnSparseF *GetGenTrackEventMatrix() {return fGenTrackEventMatrix;}
  THnSparseF *GetGenTrackSDEventMatrix() {return fGenTrackSDEventMatrix;}
  THnSparseF *GetGenTrackDDEventMatrix() {return fGenTrackDDEventMatrix;}
  THnSparseF *GetGenTrackNDEventMatrix() {return fGenTrackNDEventMatrix;}
  THnSparseF *GetGenTrackNSDEventMatrix() {return fGenTrackNSDEventMatrix;}

  THnSparseF *GetTriggerTrackEventMatrix() {return fTriggerTrackEventMatrix;}
  THnSparseF *GetTriggerTrackSDEventMatrix() {return fTriggerTrackSDEventMatrix;}
  THnSparseF *GetTriggerTrackDDEventMatrix() {return fTriggerTrackDDEventMatrix;}
  THnSparseF *GetTriggerTrackNDEventMatrix() {return fTriggerTrackNDEventMatrix;}
  THnSparseF *GetTriggerTrackNSDEventMatrix() {return fTriggerTrackNSDEventMatrix;}

  THnSparseF *GetRecTrackEventMatrix() {return fRecTrackEventMatrix;}
  THnSparseF *GetRecTrackSDEventMatrix() {return fRecTrackSDEventMatrix;}
  THnSparseF *GetRecTrackDDEventMatrix() {return fRecTrackDDEventMatrix;}
  THnSparseF *GetRecTrackNDEventMatrix() {return fRecTrackNDEventMatrix;}
  THnSparseF *GetRecTrackNSDEventMatrix() {return fRecTrackNSDEventMatrix;}

  //
  THnSparseF *GetGenPrimTrackMatrix() {return fGenPrimTrackMatrix;}
  THnSparseF *GetRecPrimTrackMatrix() {return fRecPrimTrackMatrix;}

  THnSparseF *GetRecTrackMatrix() {return fRecTrackMatrix;}
  THnSparseF *GetRecSecTrackMatrix() {return fRecSecTrackMatrix;}
  THnSparseF *GetRecMultTrackMatrix() {return fRecMultTrackMatrix;}

  //
  // control histograms
  //
  THnSparseF *GetMCEventHist1() {return fMCEventHist1;}
  THnSparseF *GetRecEventHist1() {return fRecEventHist1;}
  THnSparseF *GetRecEventHist2() {return fRecEventHist2;}
  THnSparseF *GetRecMCEventHist1() {return fRecMCEventHist1;}
  THnSparseF *GetRecMCEventHist2() {return fRecMCEventHist2;}
  THnSparseF *GetRecMCEventHist3() {return fRecMCEventHist3;}

  THnSparseF *GetMCTrackHist1(Int_t i) {return fMCTrackHist1[i];}
  THnSparseF *GetMCPrimTrackHist1(Int_t i) {return fMCPrimTrackHist1[i];}
  THnSparseF *GetMCSecTrackHist1(Int_t i) {return fMCSecTrackHist1[i];}

  THnSparseF *GetRecTrackHist1(Int_t i) {return fRecTrackHist1[i];}
  THnSparseF *GetRecTrackMultHist1(Int_t i) {return fRecTrackMultHist1[i];}

  THnSparseF *GetRecMCTrackHist1() {return fRecMCTrackHist1;}
  THnSparseF *GetMCMultRecTrackHist1() {return fMCMultRecTrackHist1;}

private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  // switch on/off filling of control histograms
  Bool_t fHistogramsOn; 

  // 
  // correlation matrices (histograms)
  //

  // event rec. track vs true track multiplicity correlation matrix 
  THnSparseF *fEventMultCorrelationMatrix; //-> mult:mult_true_tracks

  // rec. track pt vs true track pt correlation matrix for given eta
  THnSparseF *fTrackPtCorrelationMatrix; //-> Pt:mcPt:mcEta

  //
  // event level correction 
  //

  // all genertated
  THnSparseF *fGenEventMatrix; //-> mcZv:mult (inelastic)
  THnSparseF *fGenSDEventMatrix; //-> mcZv:mult (single diffractive)
  THnSparseF *fGenDDEventMatrix; //-> mcZv:mult (single diffractive)
  THnSparseF *fGenNDEventMatrix; //-> mcZv:mult (non diffractive)
  THnSparseF *fGenNSDEventMatrix; //-> mcZv:mult (non single diffractive)

  // trigger bias corrections (fTriggerEventMatrix / fGenEventMatrix)
  THnSparseF *fTriggerEventMatrix; //-> mcZv:mult
  THnSparseF *fTriggerSDEventMatrix; //-> mcZv:mult
  THnSparseF *fTriggerDDEventMatrix; //-> mcZv:mult
  THnSparseF *fTriggerNDEventMatrix; //-> mcZv:mult
  THnSparseF *fTriggerNSDEventMatrix; //-> mcZv:mult

  // event vertex rec. eff correction (fRecEventMatrix / fTriggerEventMatrix)
  THnSparseF *fRecEventMatrix; //-> mcZv:mult 
  THnSparseF *fRecSDEventMatrix; //-> mcZv:mult
  THnSparseF *fRecDDEventMatrix; //-> mcZv:mult
  THnSparseF *fRecNDEventMatrix; //-> mcZv:mult
  THnSparseF *fRecNSDEventMatrix; //-> mcZv:mult

  //
  // track-event level correction 
  //

  THnSparseF *fGenTrackEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fGenTrackSDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fGenTrackDDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fGenTrackNDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fGenTrackNSDEventMatrix; //-> mcZv:mcPt:mcEta

  // trigger bias corrections (fTriggerTrackEventMatrix / fGenTrackEventMatrix)
  THnSparseF *fTriggerTrackEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fTriggerTrackSDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fTriggerTrackDDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fTriggerTrackNDEventMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fTriggerTrackNSDEventMatrix; //-> mcZv:mcPt:mcEta

  // event vertex rec. corrections (fRecTrackEventMatrix / fTriggerTrackEventMatrix)
  THnSparseF *fRecTrackEventMatrix; //-> mcZv:Pt:mcEta
  THnSparseF *fRecTrackSDEventMatrix; //-> mcZv:Pt:mcEta
  THnSparseF *fRecTrackDDEventMatrix; //-> mcZv:Pt:mcEta
  THnSparseF *fRecTrackNDEventMatrix; //-> mcZv:Pt:mcEta
  THnSparseF *fRecTrackNSDEventMatrix; //-> mcZv:Pt:mcEta

  //
  // track level correction 
  //

  // track rec. efficiency correction (fRecPrimTrackMatrix / fGenPrimTrackMatrix)
  THnSparseF *fGenPrimTrackMatrix; //-> mcZv:mcPt:mcEta
  THnSparseF *fRecPrimTrackMatrix; //-> mcZv:mcPt:mcEta

  // secondary track contamination correction (fRecSecTrackMatrix / fRecTrackMatrix)
  THnSparseF *fRecTrackMatrix;    //-> mcZv:mcPt:mcEta
  THnSparseF *fRecSecTrackMatrix; //-> mcZv:mcPt:mcEta

  // multiple rec. track corrections (fRecMultTrackMatrix / fRecTrackMatrix)
  THnSparseF *fRecMultTrackMatrix; //-> mcZv:Pt:mcEta

  //
  // ESD and MC control analysis histograms
  //

  // THnSparse event histograms
  THnSparseF *fMCEventHist1;  //-> mcXv:mcYv:mcZv
  THnSparseF *fRecEventHist1; //-> Xv:Yv:Zv
  THnSparseF *fRecEventHist2; //-> Zv:multMB
  THnSparseF *fRecMCEventHist1; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv
  THnSparseF *fRecMCEventHist2; //-> Xv-mcXv:Zv-mcZv:Mult
  THnSparseF *fRecMCEventHist3; //-> Mult:EventType (ND, DD, SD)

  // THnSparse track histograms
  // [0] - after charged track selection, [1] - after acceptance cuts, [2] - after esd track cuts

  THnSparseF *fMCTrackHist1[AlidNdPtHelper::kCutSteps];     //-> mcPt:mcEta:mcPhi
  THnSparseF *fMCPrimTrackHist1[AlidNdPtHelper::kCutSteps]; //-> mcPt:mcEta:pid:mech:mother
  THnSparseF *fMCSecTrackHist1[AlidNdPtHelper::kCutSteps];  //-> mcPt:mcEta:pid:mech:mother

  THnSparseF *fRecTrackHist1[AlidNdPtHelper::kCutSteps];     //-> Pt:Eta:Phi
  THnSparseF *fRecTrackMultHist1[AlidNdPtHelper::kCutSteps]; //-> Pt:Mult

  THnSparseF *fRecMCTrackHist1; //-> mcPt:mcEta:(Pt-mcPt)/mcPt:(Eta-mcEta)

  //multple reconstructed tracks
  THnSparseF *fMCMultRecTrackHist1; //-> mcPt:mcEta:pid

  AlidNdPtAnalysis(const AlidNdPtAnalysis&); // not implemented
  AlidNdPtAnalysis& operator=(const AlidNdPtAnalysis&); // not implemented

  ClassDef(AlidNdPtAnalysis,1);
};

#endif
