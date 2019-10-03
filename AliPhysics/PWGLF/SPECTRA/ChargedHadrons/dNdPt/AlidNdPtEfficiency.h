#ifndef ALIDNDPTEFFICIENCY_H
#define ALIDNDPTEFFICIENCY_H

//------------------------------------------------------------------------------
// AlidNdPtEfficiency class to determine 
// efficiency TPC->ITS, ITS->TPC for dNdPt analysis. 
//
// Author: J.Otwinowski 18/11/2010 
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

class AlidNdPtEfficiency : public AlidNdPt {
public :
  AlidNdPtEfficiency(); 
  AlidNdPtEfficiency(Char_t* name, Char_t* title);
  ~AlidNdPtEfficiency();

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
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Bool_t isMatch, const Bool_t isTPC,const Bool_t isITSTPC) const;

  // Getters
  THnSparseF *GetRecMCTrackHistTPCITS() const {return fRecMCTrackHistTPCITS;}
  THnSparseF *GetRecMCTrackHistITSTPC() const {return fRecMCTrackHistITSTPC;}
  
private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  //
  // THnSparse event histograms
  //

  //TPC -> ITS matching efficiency
  THnSparseF *fRecMCTrackHistTPCITS; //-> eta:phi:pt:isPrim:charge:isMatch:isTPC

  //ITS -> TPC matching efficiency
  THnSparseF *fRecMCTrackHistITSTPC; //-> eta:phi:pt:isPrim:charge:isMatch

  AlidNdPtEfficiency(const AlidNdPtEfficiency&); // not implemented
  AlidNdPtEfficiency& operator=(const AlidNdPtEfficiency&); // not implemented

  ClassDef(AlidNdPtEfficiency,2);
};

#endif
