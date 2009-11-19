#ifndef ALIDEDPTCORRECTION_H
#define ALIDEDPTCORRECTION_H

//------------------------------------------------------------------------------
// AlidNdPtCorrection class:
//
// a. functionality:
// - applies corrections on dNdPt spectra
// - fills corrected dNdPt histograms
// - fills correction control histograms 
//
// b. data members:
// - dNdPt spectra before and after correction procedure
// - control histograms
// - correction matrices (must be loaded)
// 
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------


class TFolder;
class TObjArray;
class TString;

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

#include "THnSparse.h"
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
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms 
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * array=0);

  // Get analysis folder
  TFolder* GetCorrectionFolder() {return fCorrectionFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtCorrection",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj, Double_t zv, Int_t multRec);
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj, Int_t multRec);
  void FillHistograms(AlidNdPtHelper::EventObject eventObj, Double_t zv, Int_t multMBRec);

  // Get correction factors
  Double_t GetCorrFactZvPtEta(THnSparse *const hist=0, Double_t zv =0, Double_t pt=0, Double_t eta=0) const;
  Double_t GetContFactZvPtEta(THnSparse *const hist=0, Double_t zv =0, Double_t pt=0, Double_t eta=0) const;
  Double_t GetCorrFactZvMult(THnSparse *const hist=0, Double_t zv =0, Int_t mult=0) const;
  Double_t GetContFactZvMult(THnSparse *const hist=0, Double_t zv =0, Int_t mult=0) const;

  // Getters
  THnSparseF *GetMCAllEventMultHist1() {return fMCAllEventMultHist1;}; 
  THnSparseF *GetMCAllNDEventMultHist1() {return fMCAllNDEventMultHist1;}; 
  THnSparseF *GetMCAllNSDEventMultHist1() {return fMCAllNSDEventMultHist1;}; 
  THnSparseF *GetMCTriggerMultHist1() {return fMCTriggerMultHist1;}; 
  THnSparseF *GetMCEventMultHist1() {return fMCEventMultHist1;}; 

  THnSparseF *GetMCAllPrimTrackMultHist1() {return fMCAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCNDEventAllPrimTrackMultHist1() {return fMCNDEventAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCNSDEventAllPrimTrackMultHist1() {return fMCNSDEventAllPrimTrackMultHist1;}; 
  THnSparseF *GetMCTriggerPrimTrackMultHist1() {return fMCTriggerPrimTrackMultHist1;}; 
  THnSparseF *GetMCEventPrimTrackMultHist1() {return fMCEventPrimTrackMultHist1;}; 

  THnSparseF *GetCorrRecTrackMultHist1(Int_t i) {return fCorrRecTrackMultHist1[i];}
  THnSparseF *GetCorrRecEventHist1(Int_t i) {return fCorrRecEventHist1[i];}
  THnSparseF *GetCorrRecEventHist2(Int_t i) {return fCorrRecEventHist2[i];}

  THnSparseF *GetPtvsPt(Int_t i) {return fPtvsPt[i];}

  // correlation matrix
  void SetEventMultCorrelationMatrix(THnSparseF *const matrix=0) {fEventMultCorrelationMatrix = matrix;}
  THnSparseF *GetEventMultCorrelationMatrix() {return fEventMultCorrelationMatrix;}

  //
  void SetZvNorm(TH1D *const matrix=0) {fZvNorm = matrix;}
  TH1D *GetZvNorm() {return fZvNorm;}

  void SetZvEmptyEventsNorm(TH1D *const matrix=0) {fZvEmptyEventsNorm = matrix;}
  TH1D *GetZvEmptyEventsNorm() {return fZvEmptyEventsNorm;}

  // 
  void SetCorrEventMatrix(THnSparseF *const matrix=0) {fCorrEventMatrix = matrix;}
  THnSparseF *GetCorrEventMatrix() {return fCorrEventMatrix;}

  void SetCorrTriggerMBtoInelEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoInelEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoInelEventMatrix() {return fCorrTriggerMBtoInelEventMatrix;}

  void SetCorrTriggerMBtoNDEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNDEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNDEventMatrix() {return fCorrTriggerMBtoNDEventMatrix;}

  void SetCorrTriggerMBtoNSDEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNSDEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNSDEventMatrix() {return fCorrTriggerMBtoNSDEventMatrix;}

  // 
  void SetCorrTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTrackEventMatrix() {return fCorrTrackEventMatrix;}

  void SetCorrTriggerMBtoInelTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoInelTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoInelTrackEventMatrix() {return fCorrTriggerMBtoInelTrackEventMatrix;}

  void SetCorrTriggerMBtoNDTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNDTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNDTrackEventMatrix() {return fCorrTriggerMBtoNDTrackEventMatrix;}

  void SetCorrTriggerMBtoNSDTrackEventMatrix(THnSparseF *const matrix=0) {fCorrTriggerMBtoNSDTrackEventMatrix = matrix;}
  THnSparseF *GetCorrTriggerMBtoNSDTrackEventMatrix() {return fCorrTriggerMBtoNSDTrackEventMatrix;}

  void SetCorrTrackMatrix(THnSparseF *const matrix=0) {fCorrTrackMatrix = matrix;}
  THnSparseF *GetCorrTrackMatrix() {return fCorrTrackMatrix;}

  void SetContTrackMatrix(THnSparseF *const matrix=0) {fContTrackMatrix = matrix;}
  THnSparseF *GetContTrackMatrix() {return fContTrackMatrix;}

  void SetContMultTrackMatrix(THnSparseF *const matrix=0) {fContMultTrackMatrix = matrix;}
  THnSparseF *GetContMultTrackMatrix() {return fContMultTrackMatrix;}

  void SetCorrMatrixFileName(TString name="")    { fCorrMatrixFileName = name; }

  Int_t GetTrueMult(THnSparse *const hist=0, Int_t mult=0);

private:

  // analysis folder 
  TFolder *fCorrectionFolder; // folder for analysed histograms

  //
  // event histograms
  //
  THnSparseF *fMCEventHist1;  //-> mcXv:mcYv:mcZv
  THnSparseF *fRecEventHist1; //-> Xv:Yv:Zv

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

  //all mc primary tracks in acceptance (triggered events)
  THnSparseF *fMCTriggerPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity

  //mc primary tracks in acceptance (triggered and event vertex reconstructed)
  THnSparseF *fMCEventPrimTrackMultHist1; //-> mcPt:mcEta:multiplicity


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

  // pt vs pt to get proper pt bin (center of gravity)
  THnSparseF *fPtvsPt[8]; //-> pt:pt 

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
  THnSparseF *fContTrackMatrix;  //-> mcVz:mcPt:mcEta (fRecSecTrackMatrix / fRecTrackMatrix)
  THnSparseF *fContMultTrackMatrix; //-> mcVz:mcPt:mcEta (fRecMultTrackMatrix / fRecTrackMatrix)
  
  TString fCorrMatrixFileName; // file with efficiency correction matrices
  
  AlidNdPtCorrection(const AlidNdPtCorrection&); // not implemented
  AlidNdPtCorrection& operator=(const AlidNdPtCorrection&); // not implemented

  ClassDef(AlidNdPtCorrection,1);
};

#endif
