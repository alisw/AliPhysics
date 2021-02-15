#ifndef ALIANALYSISTASKTRACKSINJET_H
#define ALIANALYSISTASKTRACKSINJET_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliAnalysisTaskSE.h>

class TTree;
class AliAnalysisUtils;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliGenPythiaEventHeader;
class AliMCEvent;
class AliVParticle;
class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {


/**
 * @class AliAnalysisTaskTracksInJet
 * @brief Stores p-vector of jet, leading track and subleading track
 */
class AliAnalysisTaskTracksInJet: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskTracksInJet();
  AliAnalysisTaskTracksInJet(const char *taskname);
  virtual ~AliAnalysisTaskTracksInJet();

  virtual void UserCreateOutputObjects();
  virtual Bool_t UserNotify();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {}

  void SetMC(Bool_t isMC) { fIsMC = isMC; }
  void SetOutlierCut(double fracpthard = 1.2) { fFracPtHard = fracpthard; }


protected:
  struct JetData{
    Double_t        fPvecJet[3];
    Double_t        fPvecLead[3];
    Double_t        fPvecSubLead[3];
    Int_t           fIsData;

    JetData()
    {
      Reset();
    }

    void Reset(){
      memset(fPvecJet, 0, sizeof(Double_t) * 3);
      memset(fPvecLead, 0, sizeof(Double_t) * 3);
      memset(fPvecSubLead, 0, sizeof(Double_t) * 3);
      fIsData = 0;
    }
  };

  Bool_t PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const;
  Bool_t IsPhysicalPrimary(const AliVParticle* const part, AliMCEvent* const mcevent) const;
  AliGenPythiaEventHeader *GetPythiaHeader() const;
  Bool_t IsOutlier(AliGenPythiaEventHeader * const header) const;
  Bool_t TrackSelectionESDHybrid(AliESDtrack* track) const;
  Bool_t TrackSelectionESDDefault(AliESDtrack* track) const;
  Bool_t TrackSelectionAODHybrid(AliAODTrack* track) const;
  Bool_t TrackSelectionAODDefault(AliAODTrack* track) const;

  JetData                   fJetStructure;
  TTree                     *fJetTree;
  AliAnalysisUtils          *fAnalysisUtils;
  AliESDtrackCuts           *fTrackCutsDefault;
  AliESDtrackCuts           *fHybridCutsCat1;
  AliESDtrackCuts           *fHybridCutsCat2;

  Bool_t                    fIsMC;
  Double_t                  fFracPtHard;

  // Histos for MC
  THistManager    *fHistosMC;

private:

  AliAnalysisTaskTracksInJet(AliAnalysisTaskTracksInJet &ref);
  AliAnalysisTaskTracksInJet &operator=(const AliAnalysisTaskTracksInJet &ref);

  ClassDef(AliAnalysisTaskTracksInJet, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKTRACKSINJET_H */
