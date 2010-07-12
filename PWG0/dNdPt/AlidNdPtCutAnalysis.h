#ifndef ALIDNDPTCUTANALYSIS_H
#define ALIDNDPTCUTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtCutAnalysis class to determine 
// cuts to be used for dNdPt analysis. 
//
// Author: J.Otwinowski 04/11/2008 
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

class AlidNdPtCutAnalysis : public AlidNdPt {
public :
  AlidNdPtCutAnalysis(); 
  AlidNdPtCutAnalysis(Char_t* name, Char_t* title);
  ~AlidNdPtCutAnalysis();

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
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack) const;

  // Getters
  THnSparseF *GetEventCount()   const {return fEventCount;}
  THnSparseF *GetRecEventHist() const   {return fRecEventHist;}
  THnSparseF *GetMCEventHist()  const   {return fMCEventHist;}
  THnSparseF *GetRecMCEventHist() const {return fRecMCEventHist;}

  //
  THnSparseF *GetRecMCTrackHist() const {return fRecMCTrackHist;}
  
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
  THnSparseF *fRecMCTrackHist; //-> nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity

  AlidNdPtCutAnalysis(const AlidNdPtCutAnalysis&); // not implemented
  AlidNdPtCutAnalysis& operator=(const AlidNdPtCutAnalysis&); // not implemented

  ClassDef(AlidNdPtCutAnalysis,1);
};

#endif
