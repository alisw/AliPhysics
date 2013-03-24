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
class TH1;
class TLorentzVector;

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
  void           SetAODMCParticlesName(const char *n)              { fAODMCParticlesName = n     ; }
  void           SetCentralityRange(Double_t min, Double_t max)    { fMinCentrality      = min   ; fMaxCentrality    = max ; }
  void           SetTriggerMask(UInt_t mask)                       { fTriggerMask        = mask  ; }
  void           SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)      { fAODfilterBits[0]   = b0    ; fAODfilterBits[1] = b1  ; }
  void           SetIncludeNoITS(Bool_t f)                         { fIncludeNoITS       = f     ; }
  void           SetUseNegativeLabels(Bool_t f)                    { fUseNegativeLabels  = f     ; }
  void           SetTrackEfficiency(Double_t eff = 0.95)           { fTrackEfficiency    = eff   ; }
  void           SetTotalFiles(Int_t n)                            { fTotalFiles         = n     ; }
  void           SetAttempts(Int_t n)                              { fAttempts           = n     ; }
  void           SetRandomAccess(Bool_t r=kTRUE)                   { fRandomAccess       = r     ; }
  void           SetAODMC(Bool_t a)                                { fIsAODMC            = a     ; }
  void           SetJetMinPt(Double_t pt)                          { fJetMinPt           = pt    ; }
  void           SetJetEtaRange(Double_t emi, Double_t ema)        { fJetMinEta = emi; fJetMaxEta = ema; }
  void           SetJetPhiRange(Double_t pmi, Double_t pma)        { fJetMinPhi = pmi; fJetMaxPhi = pma; }
  void           SetJetConstituentMinPt(Double_t pt)               { fJetConstituentMinPt= pt    ; }
  void           SetJetType(Byte_t t)                              { fJetType            = t     ; }
  void           SetJetAlgo(Byte_t t)                              { fJetAlgo            = t     ; }
  void           SetZVertexCut(Double_t z)                         { fZVertexCut         = z     ; }

 protected:
  Bool_t          ExecOnce()            ;// intialize task
  void            Run()                 ;// do jet model action
  virtual TFile  *GetNextFile()         ;// get next file from fFileList
  virtual Bool_t  OpenNextFile()        ;// open next file
  virtual Bool_t  GetNextEntry()        ;// get next entry in current tree
  virtual Bool_t  IsAODEventSelected()  ;// AOD event trigger/centrality selection
  TLorentzVector  GetLeadingJet(TClonesArray *tracks, TClonesArray *clusters=0);  // get the leading jet

  TObjArray     *fFileList            ;//  List of AOD files 
  Bool_t         fRandomAccess        ;//  Random access to file number and event
  TString        fAODTreeName         ;//  Name of the tree in the AOD file
  TString        fAODHeaderName       ;//  Name of the header in the AOD tree
  TString        fAODVertexName       ;//  Name of the vertex branch in the AOD tree
  TString        fAODTrackName        ;//  Name of the track collection branch in the AOD tree
  TString        fAODClusName         ;//  Name of the cluster collection branch in the AOD tree
  TString        fAODCellsName        ;//  Name of the cell collection branch in the AOD tree
  TString        fAODMCParticlesName  ;//  Name of the cell collection branch in the AOD tree
  Double_t       fMinCentrality       ;//  Minimum centrality
  Double_t       fMaxCentrality       ;//  Maximum centrality
  UInt_t         fTriggerMask         ;//  Trigger selection mask
  Double_t       fZVertexCut          ;//  Z vertex cut
  Double_t       fJetMinPt            ;//  Select events with a minimum jet pt
  Double_t       fJetMinEta           ;//  Min eta for jets
  Double_t       fJetMaxEta           ;//  Max eta for jets
  Double_t       fJetMinPhi           ;//  Min phi for jets
  Double_t       fJetMaxPhi           ;//  Max phi for jets
  Double_t       fJetConstituentMinPt ;//  Jet constituent min pt
  Double_t       fJetRadius           ;//  Jet radius
  Byte_t         fJetType             ;//  Jet type (0=full, 1=charged, 2=neutral)
  Byte_t         fJetAlgo             ;//  Jet algorithm (0=kT, 1=anti-kT)
  Bool_t         fJetParticleLevel    ;//  Trigger, look at particle level jets
  Int_t          fAODfilterBits[2]    ;//  AOD track filter bit map
  Bool_t         fIncludeNoITS        ;//  True = includes tracks with failed ITS refit
  Bool_t         fUseNegativeLabels   ;//  Whether or not should use negative MC labels
  Double_t       fTrackEfficiency     ;//  Track efficiency
  Bool_t         fIsAODMC             ;//  Whether the embedding AOD is MC or not
  Int_t          fTotalFiles          ;//  Total number of files per pt hard bin
  Int_t          fAttempts            ;//  Attempts to be tried before giving up in opening the next file
  Bool_t         fEsdTreeMode         ;//! True = embed from ESD (must be a skimmed ESD!)
  Int_t          fCurrentFileID       ;//! Current file being processed (via the event handler)
  Int_t          fCurrentAODFileID    ;//! Current file ID
  TFile         *fCurrentAODFile      ;//! Current open file
  Int_t          fPicoTrackVersion    ;//! Version of the PicoTrack class (if any) in fCurrentAODFile
  TTree         *fCurrentAODTree      ;//! Current open tree
  AliVHeader    *fAODHeader           ;//! AOD header
  TClonesArray  *fAODVertex           ;//! AOD vertex
  TClonesArray  *fAODTracks           ;//! AOD track collection
  TClonesArray  *fAODClusters         ;//! AOD cluster collection
  AliVCaloCells *fAODCaloCells        ;//! AOD cell collection
  TClonesArray  *fAODMCParticles      ;//! AOD MC particles collection
  Int_t          fCurrentAODEntry     ;//! Current entry in the AOD tree
  TH2           *fHistFileMatching    ;//! Current file ID vs. AOD file ID (to be embedded)
  TH1           *fHistAODFileError    ;//! AOD file ID (to be embedded) error
  TH1           *fHistNotEmbedded     ;//! File ID not embedded
  TH1           *fHistEmbeddingQA     ;//! Embedding QA

 private:
  AliJetEmbeddingFromAODTask(const AliJetEmbeddingFromAODTask&);            // not implemented
  AliJetEmbeddingFromAODTask &operator=(const AliJetEmbeddingFromAODTask&); // not implemented

  ClassDef(AliJetEmbeddingFromAODTask, 8) // Jet embedding from AOD task
};
#endif
