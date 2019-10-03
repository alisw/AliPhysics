#ifndef ALIDNDPTCORRECTION_H
#define ALIDNDPTCORRECTION_H

//------------------------------------------------------------------------------
// AlidNdPtCorrection class to correct and
// normalised dNdPt spectra. 
//
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------

class iostream;
class TFile;
class TCint;
class TFolder;
class TObjArray;
class TString;
class THnSparse;

class AliESDtrackCuts;
class AliVertexerTracks;
class AliESD;
class AliESDfriend;
class AliESDfriendTrack;

class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;
class AlidNdPtCorrection;
class AlidNdPt;
class AlidNdPtHelper;

#include "AlidNdPt.h"

class AlidNdPtCorrection : public AlidNdPt {
public :
  AlidNdPtCorrection(); 
  AlidNdPtCorrection(Char_t* name, Char_t* title, TString corrMatrixFileName="ala.root");
  ~AlidNdPtCorrection();

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
  TFolder* GetCorrectionFolder() const {return fCorrectionFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtCorrection",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj, Double_t zv, Int_t multRec, Int_t trueMult) const;
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj, Int_t multRec) const;
  void FillHistograms(AlidNdPtHelper::EventObject eventObj, Double_t zv, Int_t multMBRec) const;

  // Get correction factors
  Double_t GetCorrFactZvPtEta(THnSparse *const hist=0, Double_t zv =0, Double_t pt=0, Double_t eta=0) const;
  Double_t GetContFactZvPtEta(THnSparse *const hist=0, Double_t zv =0, Double_t pt=0, Double_t eta=0) const;
  Double_t GetCorrFactZvMult(THnSparse *const hist=0, Double_t zv =0, Int_t mult=0) const;
  Double_t GetContFactZvMult(THnSparse *const hist=0, Double_t zv =0, Int_t mult=0) const;

  // Getters
  THnSparseF *GetMCEventHist1() const { return fMCEventHist1;};
  THnSparseF *GetRecEventHist1() const { return fRecEventHist1;};
  THnSparseF *GetRecEventMultHist1() const { return fRecEventMultHist1;};

  THnSparseF *GetMCAllEventMultHist1() const {return fMCAllEventMultHist1;}; 
  THnSparseF *GetMCAllNDEventMultHist1() const {return fMCAllNDEventMultHist1;}; 
  THnSparseF *GetMCAllNSDEventMultHist1() const {return fMCAllNSDEventMultHist1;}; 
  THnSparseF *GetMCTriggerMultHist1() const {return fMCTriggerMultHist1;}; 
  THnSparseF *GetMCEventMultHist1() const {return fMCEventMultHist1;}; 

  THnSparseF *GetMCAllPrimTrackMultHist1() const {return fMCAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCNDEventAllPrimTrackMultHist1() const {return fMCNDEventAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackMultHist1() const {return fMCNSDEventAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCTriggerPrimTrackMultHist1() const {return fMCTriggerPrimTrackMultHist1;}; 
  THnSparseF *GetMCEventPrimTrackMultHist1() const {return fMCEventPrimTrackMultHist1;}; 

  THnSparseF *GetMCAllPrimTrackTrueMultHist1() const {return fMCAllPrimTrackTrueMultHist1;}; 
  THnSparseF *GetMCNDEventAllPrimTrackTrueMultHist1() const {return fMCNDEventAllPrimTrackTrueMultHist1;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackTrueMultHist1() const {return fMCNSDEventAllPrimTrackTrueMultHist1;}; 
  THnSparseF *GetMCTriggerPrimTrackTrueMultHist1() const {return fMCTriggerPrimTrackTrueMultHist1;}; 
  THnSparseF *GetMCEventPrimTrackTrueMultHist1() const {return fMCEventPrimTrackTrueMultHist1;}; 

  THnSparseF *GetMCAllPrimTrackTrueMultHist2() const {return fMCAllPrimTrackTrueMultHist2;}; 
  THnSparseF *GetMCNDEventAllPrimTrackTrueMultHist2() const {return fMCNDEventAllPrimTrackTrueMultHist2;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackTrueMultHist2() const {return fMCNSDEventAllPrimTrackTrueMultHist2;}; 
  THnSparseF *GetMCTriggerPrimTrackTrueMultHist2() const {return fMCTriggerPrimTrackTrueMultHist2;}; 
  THnSparseF *GetMCEventPrimTrackTrueMultHist2() const {return fMCEventPrimTrackTrueMultHist2;}; 

  THnSparseF *GetMCAllPrimTrackMeanPtMult1() const {return fMCAllPrimTrackMeanPtMult1;}; 
  THnSparseF *GetMCNDEventAllPrimTrackMeanPtMult1() const {return fMCNDEventAllPrimTrackMeanPtMult1;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackMeanPtMult1() const {return fMCNSDEventAllPrimTrackMeanPtMult1;}; 
  THnSparseF *GetMCTriggerPrimTrackMeanPtMult1() const {return fMCTriggerPrimTrackMeanPtMult1;}; 
  THnSparseF *GetMCEventPrimTrackMeanPtMult1() const {return fMCEventPrimTrackMeanPtMult1;}; 

  THnSparseF *GetMCAllPrimTrackMeanPtTrueMult1() const {return fMCAllPrimTrackMeanPtTrueMult1;}; 
  THnSparseF *GetMCNDEventAllPrimTrackMeanPtTrueMult1() const {return fMCNDEventAllPrimTrackMeanPtTrueMult1;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackMeanPtTrueMult1() const {return fMCNSDEventAllPrimTrackMeanPtTrueMult1;}; 
  THnSparseF *GetMCTriggerPrimTrackMeanPtTrueMult1() const {return fMCTriggerPrimTrackMeanPtTrueMult1;}; 
  THnSparseF *GetMCEventPrimTrackMeanPtTrueMult1() const {return fMCEventPrimTrackMeanPtTrueMult1;};

  THnSparseF *GetCorrRecTrackMultHist1(Int_t i) const {return fCorrRecTrackMultHist1[i];}
  THnSparseF *GetCorrRecTrackTrueMultHist1(Int_t i) const {return fCorrRecTrackTrueMultHist1[i];}
  THnSparseF *GetCorrRecTrackTrueMultHist2(Int_t i) const {return fCorrRecTrackTrueMultHist2[i];}
  THnSparseF *GetCorrRecTrackMeanPtMultHist1(Int_t i) const {return fCorrRecTrackMeanPtMultHist1[i];}
  THnSparseF *GetCorrRecTrackMeanPtTrueMultHist1(Int_t i) const {return fCorrRecTrackMeanPtTrueMultHist1[i];}
  THnSparseF *GetCorrRecTrackPt1(Int_t i) const {return fCorrRecTrackPt1[i];}

  THnSparseF *GetCorrRecEventHist1(Int_t i) const {return fCorrRecEventHist1[i];}
  THnSparseF *GetCorrRecEventHist2(Int_t i) const {return fCorrRecEventHist2[i];}

  THnSparseF *GetEventCount()   const {return fEventCount;}

  // correlation matrix
  void SetEventMultCorrelationMatrix(THnSparseF *const matrix=0) {fEventMultCorrelationMatrix = matrix;}
  THnSparseF *GetEventMultCorrelationMatrix() const {return fEventMultCorrelationMatrix;}

  //
  void SetZvNorm(TH1D *const matrix=0) {fZvNorm = matrix;}
  TH1D *GetZvNorm() const {return fZvNorm;}

  void SetZvEmptyEventsNorm(TH1D *const matrix=0) {fZvEmptyEventsNorm = matrix;}
  TH1D *GetZvEmptyEventsNorm() const {return fZvEmptyEventsNorm;}

  void SetLHCBin0Background(TH1D *const matrix=0) {fLHCBin0Background = matrix;}
  TH1D *GetLHCBin0Background() const {return fLHCBin0Background;}
  // 
  void SetCorrEventMatrix(THnSparseF *const matrix=0) {fCorrEventMatrix = matrix;}
  THnSparseF *GetCorrEventMatrix() const {return fCorrEventMatrix;}

  void SetCorrTriggerMBtoInelEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoInelEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoInelEventMatrix() const {return fCorrTriggerMBtoInelEventMatrix;}

  void SetCorrTriggerMBtoNDEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNDEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNDEventMatrix() const {return fCorrTriggerMBtoNDEventMatrix;}

  void SetCorrTriggerMBtoNSDEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNSDEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNSDEventMatrix() const {return fCorrTriggerMBtoNSDEventMatrix;}

  // 
  void SetCorrTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTrackEventMatrix() const {return fCorrTrackEventMatrix;}

  void SetCorrTriggerMBtoInelTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoInelTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoInelTrackEventMatrix() const {return fCorrTriggerMBtoInelTrackEventMatrix;}

  void SetCorrTriggerMBtoNDTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNDTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNDTrackEventMatrix() const {return fCorrTriggerMBtoNDTrackEventMatrix;}

  void SetCorrTriggerMBtoNSDTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNSDTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNSDTrackEventMatrix() const {return fCorrTriggerMBtoNSDTrackEventMatrix;}

  void SetCorrTrackMatrix(THnSparseF *const matrix=0) {fCorrTrackMatrix = matrix;}
  THnSparseF *GetCorrTrackMatrix() const {return fCorrTrackMatrix;}

  void SetCorrHighPtTrackMatrix(THnSparseF *const matrix=0) {fCorrHighPtTrackMatrix = matrix;}
  THnSparseF *GetCorrHighPtTrackMatrix() const {return fCorrHighPtTrackMatrix;}

  void SetContTrackMatrix(THnSparseF *const matrix=0) {fContTrackMatrix = matrix;}
  THnSparseF *GetContTrackMatrix() const {return fContTrackMatrix;}

  void SetContMultTrackMatrix(THnSparseF *const matrix=0) {fContMultTrackMatrix = matrix;}
  THnSparseF *GetContMultTrackMatrix() const {return fContMultTrackMatrix;}

  void SetCorrMatrixFileName(TString name="")    { fCorrMatrixFileName = name; }

  Int_t GetTrueMult(THnSparse *const hist=0, Int_t mult=0) const;

private:

  // analysis folder 
  TFolder *fCorrectionFolder; // folder for analysed histograms

  //
  // event histograms
  //
  THnSparseF *fMCEventHist1;  //-> mcXv:mcYv:mcZv
  THnSparseF *fRecEventHist1; //-> Xv:Yv:Zv
  THnSparseF *fRecEventMultHist1; //-> track multiplicity:tracklet multiplicity

  // all MC events
  THnSparseF *fMCAllEventMultHist1; //-> mcZv:multiplicity

  // all ND MC events
  THnSparseF *fMCAllNDEventMultHist1; //-> mcZv:multiplicity

  // all NSD MC events
  THnSparseF *fMCAllNSDEventMultHist1; //-> mcZv:multiplicity

  // MC triggered events
  THnSparseF *fMCTriggerMultHist1; //-> mcZv:multiplicity

  // MC triggered and event vertex
  THnSparseF *fMCEventMultHist1; //-> mcZv:multiplicity

  //
  // track histograms
  //
  
  // all mc primary tracks in acceptance (INEL)
  THnSparseF *fMCAllPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  // all mc primary tracks in acceptance (ND events)
  THnSparseF *fMCNDEventAllPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  // all mc primary tracks in acceptance (NSD events)
  THnSparseF *fMCNSDEventAllPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  // all mc primary tracks in acceptance (triggered events)
  THnSparseF *fMCTriggerPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  // mc primary tracks in acceptance (triggered and event vertex reconstructed)
  THnSparseF *fMCEventPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  // true multiplicity 
  THnSparseF *fMCAllPrimTrackTrueMultHist1; //-> mcPt:mcEta:true multiplicity
  THnSparseF *fMCNDEventAllPrimTrackTrueMultHist1; //-> mcPt:mcEta:true multiplicity
  THnSparseF *fMCNSDEventAllPrimTrackTrueMultHist1; //-> mcPt:mcEta:true multiplicity
  THnSparseF *fMCTriggerPrimTrackTrueMultHist1; //-> mcPt:mcEta:true multiplicity
  THnSparseF *fMCEventPrimTrackTrueMultHist1; //-> mcPt:mcEta:true multiplicity

  // mcPT multiplicity vs true multiplicity 
  THnSparseF *fMCAllPrimTrackTrueMultHist2; //-> mcPt:multiplicity:true multiplicity
  THnSparseF *fMCNDEventAllPrimTrackTrueMultHist2; //-> mcPt:multiplicity:true multiplicity
  THnSparseF *fMCNSDEventAllPrimTrackTrueMultHist2; //-> mcPt:multiplicity:true multiplicity
  THnSparseF *fMCTriggerPrimTrackTrueMultHist2; //-> mcPt:multiplicity:true multiplicity
  THnSparseF *fMCEventPrimTrackTrueMultHist2; //-> mcPt:multiplicity:true multiplicity


  //
  // mc <pt> from the event
  //

  THnSparseF *fMCAllPrimTrackMeanPtMult1; //-> <mcPt>:multiplicity
  THnSparseF *fMCNDEventAllPrimTrackMeanPtMult1; //-> <mcPt>:multiplicity
  THnSparseF *fMCNSDEventAllPrimTrackMeanPtMult1; //-> <mcPt>:multiplicity
  THnSparseF *fMCTriggerPrimTrackMeanPtMult1; //-> <mcPt>:multiplicity
  THnSparseF *fMCEventPrimTrackMeanPtMult1; //-> <mcPt>:multiplicity

  // true multiplicity 

  THnSparseF *fMCAllPrimTrackMeanPtTrueMult1; //-> <mcPt>:true multiplicity
  THnSparseF *fMCNDEventAllPrimTrackMeanPtTrueMult1; //-> <mcPt>:true multiplicity
  THnSparseF *fMCNSDEventAllPrimTrackMeanPtTrueMult1; //-> <mcPt>:true multiplicity
  THnSparseF *fMCTriggerPrimTrackMeanPtTrueMult1; //-> <mcPt>:true multiplicity
  THnSparseF *fMCEventPrimTrackMeanPtTrueMult1; //-> <mcPt>:true multiplicity


  // track histograms 
  // [0]=all charged tracks, 
  // [1]=[0]+fiducual volume, 
  // [2]=[1]+after esd track cuts 
  THnSparseF *fRecTrackHist1[AlidNdPtHelper::kCutSteps]; //-> Pt:Eta:Phi

  // corrected track histograms 
  // ([0]-not corrected,
  // [1]=track cont.,
  // [2]=[1]+track eff.,
  // [3]=[2]+multple track,
  // [4]=[3]+event vertex,
  // [5]=[4]+trigger MBtoInel,
  // [6]=[4]+trigger MBtoND, 
  // [7]=[4]+trigger MBToNSD)
  THnSparseF *fCorrRecTrackMultHist1[8]; //-> Pt:Eta:mult corrected histograms 
  THnSparseF *fCorrRecTrackTrueMultHist1[8]; //-> Pt:Eta:trueMult corrected histograms
  THnSparseF *fCorrRecTrackTrueMultHist2[8]; //-> Pt:mult:trueMult corrected histograms

  // <pt> vs multiplicity from the event
  THnSparseF *fCorrRecTrackMeanPtMultHist1[8]; //-> <Pt>:mult corrected histograms
  THnSparseF *fCorrRecTrackMeanPtTrueMultHist1[8]; //-> <Pt>:trueMult corrected histograms

  THnSparseF *fCorrRecTrackPt1[8]; //-> Pt corrected histograms in the event (helper histograms)

  // corrected event histograms
  // [0]-not corrected, 
  // [1]=event vertex, 
  // [2]=[1]+trigger MBtoInel, 
  // [3]=[1]+trigger MBtoND, 
  // [4]=[1]+trigger MBtoNSD 
  THnSparseF *fCorrRecEventHist1[5]; //-> mcZv:multMB 


  // corrected event histograms (empty events)
  // [0]=not corrected,
  // [1]=trigger/trigger+vertex correction,
  // [2]=[1]+trigger MBtoInel,
  // [3]=[1]+trigger MBtoND,
  // [4]=[1]+trigger MBtoNSD
  THnSparseF *fCorrRecEventHist2[5]; //-> mcZv:multMB


  // 
  // correction matrices (histograms)
  // must be loaded
  //

  // event rec. track vs true track multiplicity correlation matrix 
  THnSparseF *fEventMultCorrelationMatrix; //-> mult:mult_true_tracks

  // histograms needed for empty events corrections
  TH1D *fZvNorm; //-> normalised reconstructed zv distribution
  TH1D *fZvEmptyEventsNorm; //-> trigger/trigger+vertex empty event correction

  // LHC background correction in 0-bin
  TH1D *fLHCBin0Background; //-> good / accepted in the 0-bin

  // trigger bias corrections
  THnSparseF *fCorrTriggerMBtoInelEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerEventMatrix / fGenEventMatrix)
  THnSparseF *fCorrTriggerMBtoNDEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerEventMatrix / fGenNDEventMatrix)
  THnSparseF *fCorrTriggerMBtoNSDEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerEventMatrix / fGenNSDEventMatrix)

  // event vertex corrections
  THnSparseF *fCorrEventMatrix;  //-> mcVz:mcPt:mcEta   (fRecEventMatrix/fTriggerEventMatrix)

  // track-event trigger bias corrections
  THnSparseF *fCorrTriggerMBtoInelTrackEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerTrackEventMatrix / fGenTrackEventMatrix)
  THnSparseF *fCorrTriggerMBtoNDTrackEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerTrackEventMatrix / fGenTrackEventNDMatrix)
  THnSparseF *fCorrTriggerMBtoNSDTrackEventMatrix; //-> mcVz:mcPt:mcEta (fTriggerTrackEventMatrix / fGenTrackEventNSDMatrix)

  // track-event vertex corrections
  THnSparseF *fCorrTrackEventMatrix;  //-> mcVz:mcPt:mcEta   (fRecTrackEventMatrix/fTriggerTrackEventMatrix)

  // track-level corrections
  THnSparseF *fCorrTrackMatrix;  //-> mcVz:mcPt:mcEta (fRecPrimTrackMatrix / fGenPrimTrackMatrix)
  THnSparseF *fCorrHighPtTrackMatrix;  //-> mcVz:mcPt:mcEta (fRecPrimTrackMatrix / fGenPrimTrackMatrix high pt tracks)
  THnSparseF *fContTrackMatrix;  //-> mcVz:mcPt:mcEta (fRecSecTrackMatrix / fRecTrackMatrix)
  THnSparseF *fContMultTrackMatrix; //-> mcVz:mcPt:mcEta (fRecMultTrackMatrix / fRecTrackMatrix)
  
  TString fCorrMatrixFileName; // file with efficiency correction matrices

  //  deta, dphi, pt1 for cosmics
  THnSparseF *fCosmicsHisto; //-> deta:dphi:pt
  
  THnSparseF *fEventCount; //-> trig, trig + vertex, selected event

  AlidNdPtCorrection(const AlidNdPtCorrection&); // not implemented
  AlidNdPtCorrection& operator=(const AlidNdPtCorrection&); // not implemented

  ClassDef(AlidNdPtCorrection,3);
};

#endif
