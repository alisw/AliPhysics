#ifndef AliEmbeddingEventForHFTask_H
#define AliEmbeddingEventForHFTask_H

// $Id$

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
//Embedding task modified to check D-meson candidates in the base sample
//Antônio Silva (antonio.silva@cern.ch)
//-----------------------------------------------------------------------

#include "AliAnalysisTaskEmcal.h"

class TFile;
class TObjArray;
class TClonesArray;
class AliRDHFCuts;
class TString;
class AliVCaloCells;
class AliVHeader;
class TH2;
class TH1;
class TLorentzVector;
class AliNamedString;
class AliAnalysisDataContainer;
class AliStack;
class AliAODRecoCascadeHF;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecay;
class TRandom3;


class AliEmbeddingEventForHFTask : public AliAnalysisTaskEmcal
{

 public:
    
    enum ACandidateType{ kD0toKpi, kDstartoKpipi };
    enum AParticleOrigin { kQuarkNotFound, kFromCharm, kFromBottom };

    AliEmbeddingEventForHFTask();
    AliEmbeddingEventForHFTask(const Char_t* name,AliRDHFCuts* cuts,ACandidateType candtype);
    virtual ~AliEmbeddingEventForHFTask();

    void     UserCreateOutputObjects();
    
    static Int_t CheckOrigin(AliAODRecoDecay* cand, TClonesArray* mcArray); // AOD
    static Int_t CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray); // AOD
    static Int_t CheckOrigin(AliAODRecoDecay* cand, AliStack* stack); // ESD
    static Int_t CheckOrigin(Int_t ipart, AliStack* stack); // ESD
  
    void SetAODTreeName(const char *t) { fAODTreeName = t; }
    void SetAODHeaderName(const char *t) { fAODHeaderName = t; }
    void SetAODTracksName(const char *n) { fAODTrackName = n; }
    void SetCentralityRange(Double_t min, Double_t max) { fMinCentrality = min; fMaxCentrality = max ;}
    void SetTriggerMask(UInt_t mask) { fTriggerMask  = mask; }
    void SetZVertexCut(Double_t z) { fZVertexCut = z; }
    void SetMaxVertexDist(Double_t d) { fMaxVertexDist = d ; }
    void SetRandomEventAccess(Bool_t c) {fRandomAccess = c;}
  
    
    void SetCheckDmeson(Bool_t c) { fCheckDmeson = c; }
    Bool_t GetCheckDmeson() const { return fCheckDmeson; }
    
    void SetRejectQuarkNotFound(Bool_t c) { fRejectQuarkNotFound = c; }
    Bool_t GetRejectQuarkNotFound() const { return fRejectQuarkNotFound; }
    
    void SetRejectDfromB(Bool_t c) { fRejectDfromB = c; }
    Bool_t GetRejectDfromB() const { return fRejectDfromB; }
    
    void SetKeepOnlyDfromB(Bool_t c) { fKeepOnlyDfromB = c; }
    Bool_t GetKeepOnlyDfromB() const { return fKeepOnlyDfromB; }
    
    void SetDSignalOnly(Bool_t c) { fDSignalOnly = c; }
    Bool_t GetDSignalOnly() const { return fDSignalOnly; }
    
    void SetGridDir(const char *t) { fGridDir = t; }
    void SetPassAndDataType(const char *t) { fPassAndDataType = t; }
    void SetDataFile(const char *t) { fDataFile = t; }
    void SetGridDir(Bool_t c) { fIsSimulation = c; }
    
    void SetEmbeddedTracksName(const char *t) {fEmbeddedTracksName = t;}
    
    void SetMaxFileNumber(Int_t c) { fMaxFileNumber = c; }
    
    void SetEmbTrackEff(Bool_t c)  { fEmbTrackEff = c; }
    void SetEmbTrackEffPath(const char *t) { fEmbTrackEffPath = t; }
    void SetEmbTrackEffHistName(const char *t) { fEmbTrackEffHistName = t; }
    
    void SetBaseTrackEff(Bool_t c) { fBaseTrackEff = c; }
    void SetBaseTrackEffPath(const char *t) { fBaseTrackEffPath = t; }
    void SetBaseTrackEffHistName(const char *t) { fBaseTrackEffHistName = t; }
    
    // inizializations

 protected:
  void ExecOnce();
  Bool_t   Run();
    virtual TFile  *GetNextFile()         ;// get next file 
    virtual Bool_t  OpenNextFile()        ;// open next file
    virtual Bool_t  GetNextEntry()        ;// get next entry in current tree
    virtual Bool_t  IsAODEventSelected()  ;// AOD event trigger/centrality selection
    void AddObjectToEvent(TObject *obj, Bool_t attempt); //Add object to the event
    Bool_t CheckForDmesons(); //Check if there is a D-meson candidate in the event
    
    Int_t          fAttempts            ;//  Attempts to be tried before giving up in opening the next file
    Int_t          fCurrentAODFileID    ;//! Current file ID
    TFile         *fCurrentAODFile      ;//! Current open file
    Int_t          fCurrentAODEntry     ;//! Current entry in the AOD tree
    Int_t          fFirstAODEntry       ;//! First entry in the AOD tree
    Int_t          fLastAODEntry        ;//! Last entry in the AOD tree
    Int_t          fEmbeddingCount      ;//! Number of embedded events from the current file
    TTree         *fCurrentAODTree      ;//! Current open tree
    TString        fAODTreeName         ;// Name of the tree in the AOD file
    TString        fAODHeaderName       ;// Name of the header in the AOD tree
    TString        fAODVertexName       ;// Name of the vertex branch in the AOD tree
    TString        fAODTrackName        ;// Name of the track collection branch in the AOD tree
    TString        fAODVZEROName        ;// Name of the VZERO collection branch in the AOD tree
    TString        fGridDir             ;// Name of the directory to be embedded (Ex.: alien:///alice/data/2015/LHC15o)
    TString        fPassAndDataType     ;// Pass and data type (Ex.: pass1/AOD)
    TString        fDataFile            ;// Name of the data file to be embedded (Ex.: AliAOD.root)
    Int_t          fMaxFileNumber       ;// Maximum number of the file to be accessed in the directory (from 1 to fMaxFileNumber)
    AliVHeader    *fAODHeader           ;//! AOD header
    AliNamedString *fAODFilePath        ;//! Current AOD file path being embedded
    Double_t       fZVertexCut          ;//  Z vertex cut
    Double_t       fMaxVertexDist       ;//  Maximum distance allowed between the vertices of the current and the embedded events
    Double_t       fMinCentrality       ;//  Minimum centrality
    Double_t       fMaxCentrality       ;//  Maximum centrality
    UInt_t         fTriggerMask         ;//  Trigger selection mask

    TClonesArray  *fAODVertex           ;//! AOD vertex
    TClonesArray  *fAODTracks           ;//! AOD track collection
    AliAODVZERO  *fAODVZERO           ;//! AOD VZERO
    TClonesArray   *fEmbeddedTracks;     //! contains tracks for embedding
    TString        fEmbeddedTracksName;  // Embedded tracks name (output)
    AliRDHFCuts    *fCuts;                   //  cuts
    AliAODEvent    *fAodEvent;               //!
    TClonesArray   *fArrayDStartoD0pi;       //!
    Bool_t          fInhibitTask;            //
    UInt_t          fCandidateType;          //  Dstar or D0
    TString         fBranchName;             //  AOD branch name
    TClonesArray   *fMCarray;                //!
    Int_t           fPDGmother;              //  PDG code of D meson
    Int_t           fNProngs;                //  number of prong of the decay channel
    Int_t           fPDGdaughters[4];        //  PDG codes of daughters
    Bool_t          fRejectQuarkNotFound;    //  reject D mesons for which the original charm or bottom quark could not be found (MC)
    Bool_t          fRejectDfromB;           //  reject D mesons coming from a B meson decay (MC)
    Bool_t          fKeepOnlyDfromB;         //  only accept D mesons coming from a B meson decay (MC)
    Bool_t          fCheckDmeson;           // Check if there are D-meson candidates in the event
    Bool_t          fDSignalOnly;           // Embed only if D-meson candidate matches MC signal
    Bool_t          fRandomAccess;          // Access Events randomically
    Bool_t          fIsSimulation;          // Is simulation (kTRUE) or data (kFALSE)
    Bool_t          fEmbTrackEff;         // Apply track efficiency to embedded tracks
    TString         fEmbTrackEffPath;
    TString         fEmbTrackEffHistName;
    Bool_t          fBaseTrackEff;         // Apply track efficiency to embedded tracks
    TString         fBaseTrackEffPath;
    TString         fBaseTrackEffHistName;
    TRandom3        *fRandomGetter;         // Random number 0-1
    TH1F            *fTrackEffHist;         //! Track efficiency
    TH1F            *fTrackEffON;           //! With track eff
    TH1F            *fTrackEffOFF;          //! Without track eff
    TH1F            *fEmbTrackEffHist;      //! Histogram with the efficiency to be applied to embedded tracks
    TH1F            *fEmbTrackEffAccepted;  //! Embedded tracks that were accepted after applied track efficiency
    TH1F            *fEmbTrackEffAll;       //! All tracks (accepted and rejected) of the embedded tracks
    TH1F            *fBaseTrackEffHist;     //! Histogram with the efficiency to be applied to base tracks
    TH1F            *fBaseTrackEffAccepted; //! Base tracks that were accepted after applied track efficiency
    TH1F            *fBaseTrackEffAll;      //! All tracks (accepted and rejected) of the base tracks
  
    
 private:
  
  AliEmbeddingEventForHFTask(const AliEmbeddingEventForHFTask &source);
  AliEmbeddingEventForHFTask& operator=(const AliEmbeddingEventForHFTask& source); 

  ClassDef(AliEmbeddingEventForHFTask, 1);
};

#endif
