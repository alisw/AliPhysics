#ifndef ALIJETEMBEDDINGFROMAODTASK_H
#define ALIJETEMBEDDINGFROMAODTASK_H

// $Id: AliJetEmbeddingFromAODTask.h  $

class TFile;
class TObjArray;
class TClonesArray;
class TString;
class AliVCaloCells;
class AliVHeader;
class TH2;

#include "AliJetModelBaseTask.h"

class AliJetEmbeddingFromAODTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingFromAODTask();
  AliJetEmbeddingFromAODTask(const char *name, Bool_t drawqa=kFALSE); 
  virtual ~AliJetEmbeddingFromAODTask();

  void           UserCreateOutputObjects();
  Bool_t         UserNotify();

  void           SetFileList(TObjArray *list)                      { fFileList           = list  ; }
  void           SetAODTreeName(const char *t)                     { fAODTreeName        = t     ; }
  void           SetAODHeaderName(const char *t)                   { fAODHeaderName      = t     ; }
  void           SetAODTracksName(const char *n)                   { fAODTrackName       = n     ; }
  void           SetAODClusName(const char *n)                     { fAODClusName        = n     ; }
  void           SetAODCellsName(const char *n)                    { fAODCellsName       = n     ; }
  void           SetCentralityRange(Double_t min, Double_t max)    { fMinCentrality      = min   ; fMaxCentrality    = max ; }
  void           SetTriggerMask(UInt_t mask)                       { fTriggerMask        = mask  ; }
  void           SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)      { fAODfilterBits[0]   = b0    ; fAODfilterBits[1] = b1  ; }
  void           SetIncludeNoITS(Bool_t f)                         { fIncludeNoITS       = f     ; }
  void           SetTotalFiles(Int_t n)                            { fTotalFiles         = n     ; }

 protected:
  Bool_t         ExecOnce()            ;// intialize task
  void           Run()                 ;// do jet model action
  Bool_t         OpenNextFile()        ;// open next file in fFileList
  Bool_t         GetNextEntry()        ;// get next entry in current tree
  Bool_t         IsAODEventSelected()  ;// AOD event trigger/centrality selection

  TObjArray     *fFileList         ;//  List of AOD files
  TString        fAODTreeName      ;//  Name of the tree in the AOD file
  TString        fAODHeaderName    ;//  Name of the header in the AOD tree
  TString        fAODVertexName    ;//  Name of the vertex branch in the AOD tree
  TString        fAODTrackName     ;//  Name of the track collection branch in the AOD tree
  TString        fAODClusName      ;//  Name of the cluster collection branch in the AOD tree
  TString        fAODCellsName     ;//  Name of the cell collection branch in the AOD tree
  Double_t       fMinCentrality    ;//  Minimum centrality
  Double_t       fMaxCentrality    ;//  Maximum centrality
  UInt_t         fTriggerMask      ;//  Trigger selection mask
  Double_t       fZVertexCut       ;//  Z vertex cut
  Int_t          fAODfilterBits[2] ;//  AOD track filter bit map
  Bool_t         fIncludeNoITS     ;//  True = includes tracks with failed ITS refit
  Int_t          fTotalFiles       ;//  Total number of files per pt hard bin
  Bool_t         fEsdTreeMode      ;//! True = embed from ESD (must be a skimmed ESD!)
  Int_t          fCurrentFileID    ;//! Current file being processed (trough the event handler)
  Int_t          fCurrentAODFileID ;//! Current file ID
  TFile         *fCurrentAODFile   ;//! Current open file
  AliVHeader    *fAODHeader        ;//! AOD header
  TClonesArray  *fAODVertex        ;//! AOD vertex
  TClonesArray  *fAODTracks        ;//! AOD track collection
  TClonesArray  *fAODClusters      ;//! AOD cluster collection
  AliVCaloCells *fAODCaloCells     ;//! AOD cell collection
  TH2           *fHistFileIDs      ;//! Current file ID vs. AOD file ID (to be embedded)

 private:
  AliJetEmbeddingFromAODTask(const AliJetEmbeddingFromAODTask&);            // not implemented
  AliJetEmbeddingFromAODTask &operator=(const AliJetEmbeddingFromAODTask&); // not implemented

  ClassDef(AliJetEmbeddingFromAODTask, 1) // Jet embedding from AOD task
};
#endif
