#ifndef ALIDEDPTCUTANALYSIS_H
#define ALIDEDPTCUTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtCutAnalysis class. 
//
// a. functionality:
// - fills generic cut histograms
// - generates cuts (selection criteria)
//
// b. data members:
// - generic cut histograms
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
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms 
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * array=0);

  // Get analysis folder
  TFolder* GetAnalysisFolder() {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack);

  // Getters
  THnSparseF *GetEventCount()   {return fEventCount;}
  THnSparseF *GetRecEventHist()   {return fRecEventHist;}
  THnSparseF *GetMCEventHist()    {return fMCEventHist;}
  THnSparseF *GetRecMCEventHist() {return fRecMCEventHist;}

  //
  THnSparseF *GetRecMCTrackHist() {return fRecMCTrackHist;}
  
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
  THnSparseF *fRecMCTrackHist; //-> nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:isNotKink:isPrim

  AlidNdPtCutAnalysis(const AlidNdPtCutAnalysis&); // not implemented
  AlidNdPtCutAnalysis& operator=(const AlidNdPtCutAnalysis&); // not implemented

  ClassDef(AlidNdPtCutAnalysis,1);
};

#endif
