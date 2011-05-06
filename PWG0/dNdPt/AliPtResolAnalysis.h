#ifndef ALIPTRESOLANALYSIS_H
#define ALIPTRESOLANALYSIS_H

//------------------------------------------------------------------------------
// AliPtResolAnalysis class used for dNdPt analysis. 
// 
// Author: J.Otwinowski 05/05/2011 
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

class AliPtResolAnalysis : public AlidNdPt {
public :
  AliPtResolAnalysis(); 
  AliPtResolAnalysis(Char_t* name, Char_t* title);
  ~AliPtResolAnalysis();

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
  TFolder* CreateFolder(TString name,TString title);

  // Get analysis folder
  TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}
  THnSparseF *GetTrackParamHist() const {return fTrackParamHist;} 

private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  THnSparseF *fTrackParamHist;  //-> sigma(1/pT):1/pT

  AliPtResolAnalysis(const AliPtResolAnalysis&); // not implemented
  AliPtResolAnalysis& operator=(const AliPtResolAnalysis&); // not implemented

  ClassDef(AliPtResolAnalysis,1);
};

#endif
