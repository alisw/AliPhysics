#ifndef ALIDNDPTCUTANALYSISPBPB_H
#define ALIDNDPTCUTANALYSISPBPB_H

//------------------------------------------------------------------------------
// AlidNdPtCutAnalysisPbPb class to determine 
// cuts to be used for dNdPt analysis. 
//
// Author: J.Otwinowski 12/11/2010 
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

#include "AlidNdPt.h"

class AlidNdPtCutAnalysisPbPb : public AlidNdPt {
public :
  AlidNdPtCutAnalysisPbPb(); 
  AlidNdPtCutAnalysisPbPb(Char_t* name, Char_t* title);
  ~AlidNdPtCutAnalysisPbPb();

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

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, Float_t centralityF) const;

  // Getters
  THnSparseF *GetEventCount()   const {return fEventCount;}
  THnSparseF *GetRecEventHist() const   {return fRecEventHist;}
  THnSparseF *GetMCEventHist()  const   {return fMCEventHist;}
  THnSparseF *GetRecMCEventHist() const {return fRecMCEventHist;}

  //
  THnSparseF *GetRecMCTrackHist() const {return fRecMCTrackHist;}
  

  TString GetCentralityEstimator() const {return fCentralityEstimator; }
  void SetCentralityEstimator(TString centEst="V0M") { fCentralityEstimator = centEst; }

private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  //
  // THnSparse event histograms
  //
  THnSparseF *fEventCount; //-> trig, trig + vertex

  THnSparseF *fRecEventHist;   //-> Xv:Yv:Zv:ResZv:Mult
  THnSparseF *fMCEventHist;    //-> mcXv:mcYv:mcZv
  THnSparseF *fRecMCEventHist; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult

  //
  // THnSparse track histograms
  //
  //THnSparseF *fRecMCTrackHist; //-> nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:hasStrangeMother:isFromMaterial:isPrim:charge:centr
  THnSparseF *fRecMCTrackHist; //-> nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:hasStrangeMother:isFromMaterial:isPrim:charge:centr
  TString fCentralityEstimator;     // use centrality can be "VOM" (default), "FMD", "TRK", "TKL", "CL0", "CL1", "V0MvsFMD", "TKLvsV0M", "ZEMvsZDC"

  AlidNdPtCutAnalysisPbPb(const AlidNdPtCutAnalysisPbPb&); // not implemented
  AlidNdPtCutAnalysisPbPb& operator=(const AlidNdPtCutAnalysisPbPb&); // not implemented

  ClassDef(AlidNdPtCutAnalysisPbPb,2);
};

#endif
