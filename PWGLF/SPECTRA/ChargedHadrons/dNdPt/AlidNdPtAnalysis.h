#ifndef ALIDNDPTANALYSIS_H
#define ALIDNDPTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtAnalysis class used for dNdPt analysis. 
// 
// Author: J.Otwinowski 04/11/2008 
// last change: 2016-02-11 by E. Perez Lezama
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

class AliTRDTriggerAnalysis;
class AliESDVertex;
class AliTRDSensorArray;

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
  void SetBinsMult(Int_t nbins, Double_t* edges) { if (CanChangeBins()) {  fMultNbins = nbins; fBinsMult = CloneArray(fMultNedges = nbins+1,edges); } }
  void SetBinsPt(Int_t nbins, Double_t* edges) { if (CanChangeBins()) {  fPtNbins = nbins; fBinsPt = CloneArray(fPtNedges = nbins+1,edges); } }
  void SetBinsPtCorr(Int_t nbins, Double_t* edges) { if (CanChangeBins()) {  fPtCorrNbins = nbins; fBinsPtCorr = CloneArray(fPtCorrNedges = nbins+1,edges); } }
  void SetBinsEta(Int_t nbins, Double_t* edges) { if (CanChangeBins()) {  fEtaNbins = nbins; fBinsEta = CloneArray(fEtaNedges = nbins+1,edges); } }
  void SetBinsZv(Int_t nbins, Double_t* edges) { if (CanChangeBins()) {  fZvNbins = nbins; fBinsZv = CloneArray(fZvNedges = nbins+1,edges); } }
  

  void SetTRDTriggerHQU(const Bool_t require = kFALSE) { fTRDTriggerRequiredHQU = require;}
  void SetTRDTriggerHJT(const Bool_t require = kFALSE) { fTRDTriggerRequiredHJT = require;}
  void SetTRDTriggerHSE(const Bool_t require = kFALSE) { fTRDTriggerRequiredHSE = require;}
  

  
  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Double_t zv, AlidNdPtHelper::TrackObject trackObj, Int_t multMB);
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj);
  void FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec);

  // PYTHIA hard pT bins
  void                        SetIsPythia(const Bool_t i)    { fIsPythia          = i                              ; }
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fTrials, Int_t &pthard); 

  
  // Getters
  THnSparseF *GetEventMultCorrelationMatrix() const {return fEventMultCorrelationMatrix;}
  THnSparseF *GetTrackPtCorrelationMatrix()   const {return fTrackPtCorrelationMatrix;}

  //
  THnSparseF *GetGenEventMatrix() const {return fGenEventMatrix;}
  THnSparseF *GetGenSDEventMatrix() const {return fGenSDEventMatrix;}
  THnSparseF *GetGenDDEventMatrix() const {return fGenDDEventMatrix;}
  THnSparseF *GetGenNDEventMatrix() const {return fGenNDEventMatrix;}
  THnSparseF *GetGenNSDEventMatrix() const {return fGenNSDEventMatrix;}
  THnSparseF *GetGenINEL0EventMatrix() const {return fGenINEL0EventMatrix;}

  THnSparseF *GetTriggerEventMatrix() const {return fTriggerEventMatrix;}
  THnSparseF *GetTriggerSDEventMatrix() const {return fTriggerSDEventMatrix;}
  THnSparseF *GetTriggerDDEventMatrix() const {return fTriggerDDEventMatrix;}
  THnSparseF *GetTriggerNDEventMatrix() const {return fTriggerNDEventMatrix;}
  THnSparseF *GetTriggerNSDEventMatrix() const {return fTriggerNSDEventMatrix;}
  THnSparseF *GetTriggerINEL0EventMatrix() const {return fTriggerINEL0EventMatrix;}

  THnSparseF *GetRecEventMatrix() const {return fRecEventMatrix;}
  THnSparseF *GetRecSDEventMatrix() const {return fRecSDEventMatrix;}
  THnSparseF *GetRecDDEventMatrix() const {return fRecDDEventMatrix;}
  THnSparseF *GetRecNDEventMatrix() const {return fRecNDEventMatrix;}
  THnSparseF *GetRecNSDEventMatrix() const {return fRecNSDEventMatrix;}
  THnSparseF *GetRecINEL0EventMatrix() const {return fRecINEL0EventMatrix;}

  //
  THnSparseF *GetRecCandleEventMatrix() const {return fRecCandleEventMatrix;}

  // 
  THnSparseF *GetGenTrackEventMatrix() const {return fGenTrackEventMatrix;}
  THnSparseF *GetGenTrackSDEventMatrix() const {return fGenTrackSDEventMatrix;}
  THnSparseF *GetGenTrackDDEventMatrix() const {return fGenTrackDDEventMatrix;}
  THnSparseF *GetGenTrackNDEventMatrix() const {return fGenTrackNDEventMatrix;}
  THnSparseF *GetGenTrackNSDEventMatrix() const {return fGenTrackNSDEventMatrix;}
  THnSparseF *GetGenTrackINEL0EventMatrix() const {return fGenTrackINEL0EventMatrix;}

  THnSparseF *GetTriggerTrackEventMatrix() const {return fTriggerTrackEventMatrix;}
  THnSparseF *GetTriggerTrackSDEventMatrix() const {return fTriggerTrackSDEventMatrix;}
  THnSparseF *GetTriggerTrackDDEventMatrix() const {return fTriggerTrackDDEventMatrix;}
  THnSparseF *GetTriggerTrackNDEventMatrix() const {return fTriggerTrackNDEventMatrix;}
  THnSparseF *GetTriggerTrackNSDEventMatrix() const {return fTriggerTrackNSDEventMatrix;}
  THnSparseF *GetTriggerTrackINEL0EventMatrix() const {return fTriggerTrackINEL0EventMatrix;}

  THnSparseF *GetRecTrackEventMatrix() const {return fRecTrackEventMatrix;}
  THnSparseF *GetRecTrackSDEventMatrix() const {return fRecTrackSDEventMatrix;}
  THnSparseF *GetRecTrackDDEventMatrix() const {return fRecTrackDDEventMatrix;}
  THnSparseF *GetRecTrackNDEventMatrix() const {return fRecTrackNDEventMatrix;}
  THnSparseF *GetRecTrackNSDEventMatrix() const {return fRecTrackNSDEventMatrix;}
  THnSparseF *GetRecTrackINEL0EventMatrix() const {return fRecTrackINEL0EventMatrix;}


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

  THnSparseF *GetRecPrimTrackHist() const {return fRecPrimTrackHist;}
  THnSparseF *GetGenPrimTrackHist() const {return fGenPrimTrackHist;}

  
  THnSparseF *GetRecMCTrackHist1() const {return fRecMCTrackHist1;}
  THnSparseF *GetMCMultRecTrackHist1() const {return fMCMultRecTrackHist1;}

  THnSparseF *GetRecTrackHist2() const {return fRecTrackHist2;}
  THnSparseF *GetEventCount() const {return fEventCount;}  
  TH1D       *GetPileUpCount() const {return fPileUpCount;}

  //
  // Generic histograms to be corrected
  //
  THnSparseF *GetRecEventHist() const {return fRecEventHist;} 
  THnSparseF *GetRecTrackHist() const {return fRecTrackHist;} 
  THnSparseF* GetMCPrimTrackHist() const {return fMCPrimTrackHist;}   

//   Bool_t IsRequiredCompleteDAQ() const  { return f2015IsIncompleteDAQ;}
 

  TH1D       *GetHistoTrials() const {return fHistoTrials;}
  TH1D       *GetHEvents() const {return fHEvents;}
  TProfile   *GetProfileXsection() const {return fProfXsection;}
  
  
private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  Bool_t fHistogramsOn; // switch on/off filling of control histograms 

  // 
  // correlation matrices (histograms)
  //

  // event rec. track vs true track multiplicity correlation matrix 
  THnSparseF *fEventMultCorrelationMatrix; //-> multRecMult:multTrueMC:multMB

  // rec. track pt vs true track pt correlation matrix for given eta
  THnSparseF *fTrackPtCorrelationMatrix; //-> Pt:mcPt:mcEta

  //
  // event level correction 
  //

  // all genertated
  THnSparseF *fGenEventMatrix; //-> mcZv:multTrueMC (inelastic)
  THnSparseF *fGenSDEventMatrix; //-> mcZv:multTrueMC (single diffractive)
  THnSparseF *fGenDDEventMatrix; //-> mcZv:multTrueMC (single diffractive)
  THnSparseF *fGenNDEventMatrix; //-> mcZv:multTrueMC (non diffractive)
  THnSparseF *fGenNSDEventMatrix; //-> mcZv:multTrueMC (non single diffractive)
  THnSparseF *fGenINEL0EventMatrix; //-> mcZv:multTrueMC (non single diffractive)

  // trigger bias corrections (fTriggerEventMatrix / fGenEventMatrix)
  THnSparseF *fTriggerEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fTriggerSDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fTriggerDDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fTriggerNDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fTriggerNSDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fTriggerINEL0EventMatrix; //-> mcZv:multTrueMC

  // event vertex rec. eff correction (fRecEventMatrix / fTriggerEventMatrix)
  THnSparseF *fRecEventMatrix; //-> mcZv:multTrueMC 
  THnSparseF *fRecSDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fRecDDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fRecNDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fRecNSDEventMatrix; //-> mcZv:multTrueMC
  THnSparseF *fRecINEL0EventMatrix; //-> mcZv:multTrueMC


  //
  // track-event level correction 
  //

  THnSparseF *fGenTrackEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fGenTrackSDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fGenTrackDDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fGenTrackNDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fGenTrackNSDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fGenTrackINEL0EventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC

  // trigger bias corrections (fTriggerTrackEventMatrix / fGenTrackEventMatrix)
  THnSparseF *fTriggerTrackEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fTriggerTrackSDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fTriggerTrackDDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fTriggerTrackNDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fTriggerTrackNSDEventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC
  THnSparseF *fTriggerTrackINEL0EventMatrix; //-> mcZv:mcPt:mcEta:multTrueMC

  // event vertex rec. corrections (fRecTrackEventMatrix / fTriggerTrackEventMatrix)
  THnSparseF *fRecTrackEventMatrix; //-> mcZv:Pt:mcEta:multTrueMC
  THnSparseF *fRecTrackSDEventMatrix; //-> mcZv:Pt:mcEta:multTrueMC
  THnSparseF *fRecTrackDDEventMatrix; //-> mcZv:Pt:mcEta:multTrueMC
  THnSparseF *fRecTrackNDEventMatrix; //-> mcZv:Pt:mcEta:multTrueMC
  THnSparseF *fRecTrackNSDEventMatrix; //-> mcZv:Pt:mcEta:multTrueMC
  THnSparseF *fRecTrackINEL0EventMatrix; //-> mcZv:Pt:mcEta:multTrueMC

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
  THnSparseF *fRecEventHist2; //-> Zv:multMB:multRecMult
  THnSparseF *fRecMCEventHist1; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv
  THnSparseF *fRecMCEventHist2; //-> Xv-mcXv:Zv-mcZv:multMB
  THnSparseF *fRecMCEventHist3; //-> multRecMult:EventType (ND, DD, SD)

  THnSparseF *fRecPrimTrackHist; // -> pt:eta:phi reconstructed
  THnSparseF *fGenPrimTrackHist; // -> pt:eta:phi generated

  
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
  THnSparseF *fRecEventHist;    //-> Zv:multMB  
  THnSparseF *fRecTrackHist; //-> Zv:pT:eta:multRecMult
  THnSparseF *fEventCount; //-> trig, trig + vertex, selected event
  THnSparseF *fMCPrimTrackHist; //-> Zv:mcpT:mceta:multTrueMC  
  TH1D       *fPileUpCount; // if pileup x<0 
  TH1D       *fSPDBGCount;  // spd background 

  //
  // candle events track corrections
  //
  THnSparseF *fRecCandleEventMatrix; //-> Zv:multTrueMC

  //binning for THNsparse
  Int_t	fMultNbins;
  Int_t fPtNbins;
  Int_t fPtCorrNbins;
  Int_t fEtaNbins;
  Int_t fZvNbins;
  
  Int_t fMultNedges;   // fMultNbins+1 uses for streaming dynamic array
  Int_t fPtNedges;     // fPtNbins+1 uses for streaming dynamic array  
  Int_t fPtCorrNedges; // fCentralityNbins+1 uses for streaming dynamic array  
  Int_t fEtaNedges;    // fEtaNbins+1 uses for streaming dynamic array
  Int_t fZvNedges;     // fZvNbins+1 uses for streaming dynamic array    
  
  Double_t* fBinsMult;//[fMultNedges]
  Double_t* fBinsPt;//[fPtNedges]
  Double_t* fBinsPtCorr;//[fPtCorrNedges]
  Double_t* fBinsEta;//[fEtaNedges]
  Double_t* fBinsZv;//[fZvNedges]
  
  Bool_t fIsInit;
  
  Bool_t fTRDTriggerRequiredHQU;
  Bool_t fTRDTriggerRequiredHJT;
  Bool_t fTRDTriggerRequiredHSE;
  Bool_t triggerResult;
//   Bool_t f2015IsIncompleteDAQ;

  // Pythia hard pT bins
  Bool_t                      fIsPythia;                   ///< trigger, if it is a PYTHIA production
  AliGenPythiaEventHeader    *fPythiaHeader;               //!<!event Pythia header
  Int_t                       fNTrials;                    //!<!event trials
  Float_t                     fXsection;                   //!<!x-section from pythia header
  Float_t		      fPtHardbin;
  TH1D                       *fHistoTrials;
  TH1D                       *fHEvents;
  TProfile                   *fProfXsection;
  // generic function to change binning
  Bool_t CanChangeBins();

  AlidNdPtAnalysis(const AlidNdPtAnalysis&); // not implemented
  AlidNdPtAnalysis& operator=(const AlidNdPtAnalysis&); // not implemented

  ClassDef(AlidNdPtAnalysis,12);
};

#endif
