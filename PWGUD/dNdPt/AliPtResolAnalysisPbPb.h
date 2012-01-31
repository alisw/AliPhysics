#ifndef ALIPTRESOLANALYSISPBPB_H
#define ALIPTRESOLANALYSISPBPB_H

//------------------------------------------------------------------------------
// AliPtResolAnalysisPbPb class used for dNdPt analysis. 
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

class AliPtResolAnalysisPbPb : public AlidNdPt {
public :
  AliPtResolAnalysisPbPb(); 
  AliPtResolAnalysisPbPb(Char_t* name, Char_t* title);
  ~AliPtResolAnalysisPbPb();

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
  THnSparseF *GetTrackParamHist2() const {return fTrackParamHist2;} 

  void SetCentralityEstimator(TString centEst="V0M") { fCentralityEstimator = centEst; }
  TString GetCentralityEstimator() const {return fCentralityEstimator; }

private:

  // analysis folder 
  TFolder *fAnalysisFolder;     // folder for analysed histograms
  THnSparseF *fTrackParamHist;  //-> sigma(1/pT):1/pT:centr
  THnSparseF *fTrackParamHist2; //-> sigma(1/pT)*pT:pT:centr

  TString fCentralityEstimator;     // use centrality can be "VOM" (default), "FMD", "TRK", "TKL", "CL0", "CL1", "V0MvsFMD", "TKLvsV0M", "ZEMvsZDC"

  AliPtResolAnalysisPbPb(const AliPtResolAnalysisPbPb&); // not implemented
  AliPtResolAnalysisPbPb& operator=(const AliPtResolAnalysisPbPb&); // not implemented

  ClassDef(AliPtResolAnalysisPbPb,2);
};

#endif
