#ifndef ALIMCTRUTHTRACKMAKER_H
#define ALIMCTRUTHTRACKMAKER_H

class TClonesArray;
class AliESDEvent;
class AliESDMuonTrack;
class AliMCEvent;
class AliAODTrack;
class AliVParticle;
class AliStack;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliMCTruthTrackMaker : public AliAnalysisTaskSE {
 public:
  AliMCTruthTrackMaker();
  AliMCTruthTrackMaker(const char *name);
  virtual ~AliMCTruthTrackMaker();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  void SetChargedMC(Bool_t c = kTRUE)                { fChargedMC        = c    ; }
  void SetEtaMax(Double_t e)                         { fEtaMax           = e    ; }
  void SetFillMuMothers(Bool_t c = kTRUE)            { fFillMuMothers    = c    ; }
  void SetTriggerMatch(Bool_t b)                     { fTriggerMatch     = b    ; }
  void SetTracksOutName(const char *name)            { fTracksOutName    = name ; }

 protected:
  Int_t              GetNumberOfTracks() const;
  AliVParticle      *GetTrack(Int_t i);
  void               AddTrack(AliVParticle *track, Int_t nacc);
  Bool_t             IsGoodMUONtrack(AliESDMuonTrack &track);
  Bool_t             IsGoodMUONtrack(AliAODTrack &track);
  Int_t              GetFirstPrimaryMother(Int_t muonlabel);
  
  TString            fTracksOutName;        // name of output track array
  Bool_t             fChargedMC;            // true = only charged particles
  Bool_t             fFillMuMothers;        // true = fill primary mother of reconstructed muon tracksmuon cuts
  Bool_t             fTriggerMatch;         // require trigger match as muon cut?
  Double_t           fEtaMax;               // maximum eta to accept tracks
  Bool_t             fInit;                 //!true = task initialized
  Bool_t             fEsdMode;              //!switch for ESD/AOD mode
  TClonesArray      *fTracksIn;             //!track array in (AOD only)
  TClonesArray      *fTracksOut;            //!track array out
  AliESDEvent       *fESD;                  //! ESD object
  AliMCEvent        *fMC;                   //! MC object
  AliStack          *fStack;                //! MC stack
  
 private:
  AliMCTruthTrackMaker(const AliMCTruthTrackMaker&);            // not implemented
  AliMCTruthTrackMaker &operator=(const AliMCTruthTrackMaker&); // not implemented

  ClassDef(AliMCTruthTrackMaker, 1); // Task to select tracks in MC events
};
#endif
