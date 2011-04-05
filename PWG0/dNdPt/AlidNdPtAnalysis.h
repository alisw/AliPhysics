#ifndef ALIDNDPTANALYSIS_H
#define ALIDNDPTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtAnalysis class used for dNdPt analysis. 
// 
// Author: J.Otwinowski 04/11/2008 
// last change: 2011-04-04 by M.Knichel
//------------------------------------------------------------------------------

class iostream;

class TFile;
class TCint;
class TProfile;
class TFolder;
class TObjArray;
class TString;
class THnSparse;

class AliESDtrackCuts;
class AliVertexerTracks;
class AliESD;
class AliESDfriend;
class AliESDfriendTrack;
class AlidNdPtHelper;

#include "AlidNdPt.h"

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
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms 
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * const array=0);

  // Get analysis folder
  TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Fill control histograms
  void SetHistogramsOn(const Bool_t histOn=kTRUE) {fHistogramsOn = histOn;}
  Bool_t IsHistogramsOn() const {return fHistogramsOn;}

  // Define 0-multiplicity bin for LHC
  // background calculation
  static Bool_t IsBinZeroSPDvtx(const AliESDEvent* esdEvent);
  static Bool_t IsBinZeroTrackSPDvtx(const AliESDEvent* esdEvent);
    
  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");

  // Set binning for Histograms (if not set default binning is used)
  void SetBinsMult(Int_t nbins, Double_t* edges) { fMultNbins = nbins; fBinsMult = edges; }
  void SetBinsPt(Int_t nbins, Double_t* edges) { fPtNbins = nbins; fBinsPt = edges; }
  void SetBinsPtCorr(Int_t nbins, Double_t* edges) { fPtCorrNbins = nbins; fBinsPtCorr = edges; }
  void SetBinsEta(Int_t nbins, Double_t* edges) { fEtaNbins = nbins; fBinsEta = edges; }
  void SetBinsZv(Int_t nbins, Double_t* edges) { fZvNbins = nbins; fBinsZv = edges; }  


  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Double_t zv, AlidNdPtHelper::TrackObject trackObj);
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj);
  void FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec);

  // Getters
  THnSparseF *GetEventMultCorrelationMatrix() const {return fEventMultCorrelationMatrix;}
  THnSparseF *GetTrackPtCorrelationMatrix()   const {return fTrackPtCorrelationMatrix;}

  //
  THnSparseF *GetGenEventMatrix() const {return fGenEventMatrix;}
  THnSparseF *GetGenSDEventMatrix() const {return fGenSDEventMatrix;}
  THnSparseF *GetGenDDEventMatrix() const {return fGenDDEventMatrix;}
  THnSparseF *GetGenNDEventMatrix() const {return fGenNDEventMatrix;}
  THnSparseF *GetGenNSDEventMatrix() const {return fGenNSDEventMatrix;}

  THnSparseF *GetTriggerEventMatrix() const {return fTriggerEventMatrix;}
  THnSparseF *GetTriggerSDEventMatrix() const {return fTriggerSDEventMatrix;}
  THnSparseF *GetTriggerDDEventMatrix() const {return fTriggerDDEventMatrix;}
  THnSparseF *GetTriggerNDEventMatrix() const {return fTriggerNDEventMatrix;}
  THnSparseF *GetTriggerNSDEventMatrix() const {return fTriggerNSDEventMatrix;}

  THnSparseF *GetRecEventMatrix() const {return fRecEventMatrix;}
  THnSparseF *GetRecSDEventMatrix() const {return fRecSDEventMatrix;}
  THnSparseF *GetRecDDEventMatrix() const {return fRecDDEventMatrix;}
  THnSparseF *GetRecNDEventMatrix() const {return fRecNDEventMatrix;}
  THnSparseF *GetRecNSDEventMatrix() const {return fRecNSDEventMatrix;}

  //
  THnSparseF *GetRecCandleEventMatrix() const {return fRecCandleEventMatrix;}

  // 
  THnSparseF *GetGenTrackEventMatrix() const {return fGenTrackEventMatrix;}
  THnSparseF *GetGenTrackSDEventMatrix() const {return fGenTrackSDEventMatrix;}
  THnSparseF *GetGenTrackDDEventMatrix() const {return fGenTrackDDEventMatrix;}
  THnSparseF *GetGenTrackNDEventMatrix() const {return fGenTrackNDEventMatrix;}
  THnSparseF *GetGenTrackNSDEventMatrix() const {return fGenTrackNSDEventMatrix;}

  THnSparseF *GetTriggerTrackEventMatrix() const {return fTriggerTrackEventMatrix;}
  THnSparseF *GetTriggerTrackSDEventMatrix() const {return fTriggerTrackSDEventMatrix;}
  THnSparseF *GetTriggerTrackDDEventMatrix() const {return fTriggerTrackDDEventMatrix;}
  THnSparseF *GetTriggerTrackNDEventMatrix() const {return fTriggerTrackNDEventMatrix;}
  THnSparseF *GetTriggerTrackNSDEventMatrix() const {return fTriggerTrackNSDEventMatrix;}

  THnSparseF *GetRecTrackEventMatrix() const {return fRecTrackEventMatrix;}
  THnSparseF *GetRecTrackSDEventMatrix() const {return fRecTrackSDEventMatrix;}
  THnSparseF *GetRecTrackDDEventMatrix() const {return fRecTrackDDEventMatrix;}
  THnSparseF *GetRecTrackNDEventMatrix() const {return fRecTrackNDEventMatrix;}
  THnSparseF *GetRecTrackNSDEventMatrix() const {return fRecTrackNSDEventMatrix;}

  //
  THnSparseF *GetGenTrackMatrix() const {return fGenTrackMatrix;}
  THnSparseF *GetGenPrimTrackMatrix() const {return fGenPrimTrackMatrix;}
  THnSparseF *GetRecPrimTrackMatrix() const {return fRecPrimTrackMatrix;}

  THnSparseF *GetRecTrackMatrix() const {return fRecTrackMatrix;}
  THnSparseF *GetRecSecTrackMatrix() const {return fRecSecTrackMatrix;}
  THnSparseF *GetRecMultTrackMatrix() const {return fRecMultTrackMatrix;}

  //
  // control histograms
  //
  THnSparseF *GetMCEventHist1() const {return fMCEventHist1;}
  THnSparseF *GetRecEventHist1() const {return fRecEventHist1;}
  THnSparseF *GetRecEventHist2() const {return fRecEventHist2;}
  THnSparseF *GetRecMCEventHist1() const {return fRecMCEventHist1;}
  THnSparseF *GetRecMCEventHist2() const {return fRecMCEventHist2;}
  THnSparseF *GetRecMCEventHist3() const {return fRecMCEventHist3;}

  THnSparseF *GetMCTrackHist1(Int_t i) const {return fMCTrackHist1[i];}
  THnSparseF *GetMCPrimTrackHist1(Int_t i) const {return fMCPrimTrackHist1[i];}
  THnSparseF *GetMCPrimTrackHist2(Int_t i) const {return fMCPrimTrackHist2[i];}
  THnSparseF *GetMCSecTrackHist1(Int_t i) const {return fMCSecTrackHist1[i];}

  THnSparseF *GetRecTrackHist1(Int_t i) const {return fRecTrackHist1[i];}
  THnSparseF *GetRecTrackMultHist1(Int_t i) const {return fRecTrackMultHist1[i];}

  THnSparseF *GetRecMCTrackHist1() const {return fRecMCTrackHist1;}
  THnSparseF *GetMCMultRecTrackHist1() const {return fMCMultRecTrackHist1;}

  THnSparseF *GetRecTrackHist2() const {return fRecTrackHist2;}


  //
  // Generic histograms to be corrected
  //
  THnSparseF *GetRecEventHist() const {return fRecEventHist;} 
  THnSparseF *GetRecTrackHist() const {return fRecTrackHist;} 

private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  Bool_t fHistogramsOn; // switch on/off filling of control histograms 

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
  THnSparseF *fGenEventMatrix; //-> mcZv:multMB (inelastic)
  THnSparseF *fGenSDEventMatrix; //-> mcZv:multMB (single diffractive)
  THnSparseF *fGenDDEventMatrix; //-> mcZv:multMB (single diffractive)
  THnSparseF *fGenNDEventMatrix; //-> mcZv:multMB (non diffractive)
  THnSparseF *fGenNSDEventMatrix; //-> mcZv:multMB (non single diffractive)

  // trigger bias corrections (fTriggerEventMatrix / fGenEventMatrix)
  THnSparseF *fTriggerEventMatrix; //-> mcZv:multMB
  THnSparseF *fTriggerSDEventMatrix; //-> mcZv:multMB
  THnSparseF *fTriggerDDEventMatrix; //-> mcZv:multMB
  THnSparseF *fTriggerNDEventMatrix; //-> mcZv:multMB
  THnSparseF *fTriggerNSDEventMatrix; //-> mcZv:multMB

  // event vertex rec. eff correction (fRecEventMatrix / fTriggerEventMatrix)
  THnSparseF *fRecEventMatrix; //-> mcZv:multMB 
  THnSparseF *fRecSDEventMatrix; //-> mcZv:multMB
  THnSparseF *fRecDDEventMatrix; //-> mcZv:multMB
  THnSparseF *fRecNDEventMatrix; //-> mcZv:multMB
  THnSparseF *fRecNSDEventMatrix; //-> mcZv:multMB


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
  THnSparseF *fGenTrackMatrix; //-> mcZv:mcPt:mcEta
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
  THnSparseF *fRecEventHist2; //-> Zv:multMB:mult
  THnSparseF *fRecMCEventHist1; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv
  THnSparseF *fRecMCEventHist2; //-> Xv-mcXv:Zv-mcZv:mult
  THnSparseF *fRecMCEventHist3; //-> mult:EventType (ND, DD, SD)

  // THnSparse track histograms
  // [0] - after charged track selection, [1] - after acceptance cuts, [2] - after esd track cuts

  THnSparseF *fMCTrackHist1[AlidNdPtHelper::kCutSteps];     //-> mcPt:mcEta:mcPhi
  THnSparseF *fMCPrimTrackHist1[AlidNdPtHelper::kCutSteps]; //-> mcPt:mcEta:pid:mech:mother
  THnSparseF *fMCPrimTrackHist2[AlidNdPtHelper::kCutSteps]; //-> pdg:mech:mother
  THnSparseF *fMCSecTrackHist1[AlidNdPtHelper::kCutSteps];  //-> mcPt:mcEta:pid:mech:mother

  THnSparseF *fRecTrackHist1[AlidNdPtHelper::kCutSteps];     //-> Pt:Eta:Phi
  THnSparseF *fRecTrackMultHist1[AlidNdPtHelper::kCutSteps]; //-> Pt:mult

  THnSparseF *fRecMCTrackHist1; //-> mcPt:mcEta:(Pt-mcPt)/mcPt:(Eta-mcEta)

  //multple reconstructed tracks
  THnSparseF *fMCMultRecTrackHist1; //-> mcPt:mcEta:pid

  // track control histograms
  THnSparseF *fRecTrackHist2;  //-> nclust:chi2:Pt:Eta:Phi

  //
  // Generic histograms to be corrected
  //
  THnSparseF *fRecEventHist; //-> Zv:multMB
  THnSparseF *fRecTrackHist; //-> Zv:pT:eta

  //
  // candle events track corrections
  //
  THnSparseF *fRecCandleEventMatrix; //-> Zv:multMB

  //binning for THNsparse
  Int_t	fMultNbins;
  Int_t fPtNbins;
  Int_t fPtCorrNbins;
  Int_t fEtaNbins;
  Int_t fZvNbins;
  Double_t* fBinsMult;
  Double_t* fBinsPt;
  Double_t* fBinsPtCorr;
  Double_t* fBinsEta;
  Double_t* fBinsZv;
  
  Bool_t fIsInit;

  AlidNdPtAnalysis(const AlidNdPtAnalysis&); // not implemented
  AlidNdPtAnalysis& operator=(const AlidNdPtAnalysis&); // not implemented

  ClassDef(AlidNdPtAnalysis,5);
};

#endif
