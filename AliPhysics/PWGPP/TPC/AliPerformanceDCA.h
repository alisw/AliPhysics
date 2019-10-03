#ifndef ALIPERFORMANCEDCA_H
#define ALIPERFORMANCEDCA_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (DCA - Distance of Closest Approach 
// to the vertex).   
// 
// Author: J.Otwinowski 04/02/2008
// Changes by J.Salzwedel 22/10/2014
//------------------------------------------------------------------------------

class AliVEvent; 
class AliVfriendEvent; 
class AliESDVertex;
class AliVTrack;
class TH3;
class TH2;
class TH1;
class TString;
class TNamed;
class TRootIOCtor;

#include "THnSparse.h"
#include "AliPerformanceObject.h"

class AliPerformanceDCA : public AliPerformanceObject {
public :
  AliPerformanceDCA(TRootIOCtor*);
  AliPerformanceDCA(const Char_t* name="AliPerformanceDCA", const Char_t* title="AliPerformanceDCA",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE);

  virtual ~AliPerformanceDCA();

  // Init data members
  virtual void Init();

  // Execute analysis
  virtual void Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderDCA",TString title = "Analysed DCA histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  void ProcessConstrained(AliMCEvent* const mcev, AliVTrack *const vTrack);
  void ProcessTPC(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent);
  void ProcessTPCITS(AliMCEvent* const mcev, AliVTrack *const vTrack, AliVEvent* const vEvent);

  // getters
  THnSparse* GetDCAHisto() const {return fDCAHisto;}

  // Make stat histograms
  TH1F* MakeStat1D(TH2 *hist, Int_t delta1, Int_t type);
  TH2F* MakeStat2D(TH3 *hist, Int_t delta0, Int_t delta1, Int_t type);

private:

  // DCA histograms
  THnSparseF *fDCAHisto; //-> dca_r:dca_z:eta:pt:phi 
 
  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliPerformanceDCA(const AliPerformanceDCA&); // not implemented
  AliPerformanceDCA& operator=(const AliPerformanceDCA&); // not implemented

  ClassDef(AliPerformanceDCA,2);
};

#endif
