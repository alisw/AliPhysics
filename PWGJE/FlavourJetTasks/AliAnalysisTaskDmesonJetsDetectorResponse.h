/// \class AliAnalysisTaskDmesonJetsDetectorResponse
/// \brief Analysis task used to build a detector response for D meson jets
///
/// This task derives from AliAnalysisTaskDmesonJets.
/// Most of the analysis is performed there. This task only
/// takes care of matching detector level D meson jets with generator level.
/// The matching is done using the method MatchToMC of AliAODRecoDecayHF2Prong.
///
/// The main output is stored in a THnSparse histogram or in a TTree.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date June 9, 2016

#ifndef ALIANALYSISTASKDMESONJETSDETECTORRESPONSE_H
#define ALIANALYSISTASKDMESONJETSDETECTORRESPONSE_H

/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskDmesonJets.h"

class AliAnalysisTaskDmesonJetsDetectorResponse : public AliAnalysisTaskDmesonJets
{
 public:

  /// \class AliDmesonMatchInfoSummary
  /// \brief Lightweight class that encapsulates matching between reconstructed and generated D mesons
  ///
  /// This virtual class provides a common interface to AliD0MatchInfoSummary and AliDStarMatchInfoSummary
  class AliDmesonMatchInfoSummary {
  public:
    AliDmesonMatchInfoSummary() {;}
    virtual ~AliDmesonMatchInfoSummary() {}

    virtual void Reset() = 0;
    virtual void SetReconstructed(const AliDmesonJetInfo& reco) = 0;
    virtual void SetGenerated(const AliDmesonJetInfo& truth) = 0;
    virtual AliDmesonInfoSummary* GetReconstructed() = 0;
    virtual AliDmesonInfoSummary* GetGenerated() = 0;

    /// \cond CLASSIMP
    ClassDef(AliDmesonMatchInfoSummary, 1);
    /// \endcond
  };

  /// \class AliD0MatchInfoSummary
  /// \brief Lightweight class that encapsulates matching between reconstructed and generated D0 mesons
  ///
  /// This class encapsulates matching information between reconstructed and generated D0 mesons
  /// in a compact data structure
  class AliD0MatchInfoSummary : public AliDmesonMatchInfoSummary {
  public:
    AliD0MatchInfoSummary() : AliDmesonMatchInfoSummary(), fReconstructed(), fGenerated() {}

    virtual void Reset();
    virtual void SetReconstructed(const AliDmesonJetInfo& reco);
    virtual void SetGenerated(const AliDmesonJetInfo& truth);

    virtual AliDmesonInfoSummary* GetReconstructed() { return &fReconstructed; }
    virtual AliDmesonInfoSummary* GetGenerated()     { return &fGenerated    ; }

    AliD0InfoSummary     fReconstructed ; ///<  Reconstructed D meson
    AliDmesonInfoSummary fGenerated     ; ///<  Generated D meson

    /// \cond CLASSIMP
    ClassDef(AliD0MatchInfoSummary, 1);
    /// \endcond
  };

  /// \class AliDStarMatchInfoSummary
  /// \brief Lightweight class that encapsulates matching between reconstructed and generated D* mesons
  ///
  /// This class encapsulates matching information between reconstructed and generated D* mesons
  /// in a compact data structure
  class AliDStarMatchInfoSummary : public AliDmesonMatchInfoSummary {
  public:
    AliDStarMatchInfoSummary() : AliDmesonMatchInfoSummary(), fReconstructed(), fGenerated() {}

    virtual void Reset();
    virtual void SetReconstructed(const AliDmesonJetInfo& reco);
    virtual void SetGenerated(const AliDmesonJetInfo& truth);

    virtual AliDmesonInfoSummary* GetReconstructed() { return &fReconstructed; }
    virtual AliDmesonInfoSummary* GetGenerated()     { return &fGenerated    ; }

    AliDStarInfoSummary     fReconstructed ; ///<  Reconstructed D meson
    AliDmesonInfoSummary    fGenerated     ; ///<  Generated D meson

    /// \cond CLASSIMP
    ClassDef(AliDStarMatchInfoSummary, 1);
    /// \endcond
  };

  /// \class ResponseEngine
  /// \brief Analysis engine to produce detector response matrix in the D meson jet analysis
  ///
  /// Analysis engine to produce detector response matrix in the D meson jet analysis
  class ResponseEngine : public TObject {
  public:
    ResponseEngine();
    ResponseEngine(ECandidateType_t type);
    ResponseEngine(const ResponseEngine &source);
    ResponseEngine& operator=(const ResponseEngine& source);
    virtual ~ResponseEngine() {;}

    void SetMaxJetDmesonDistance(Double_t d) { fMaxJetDmesonDistance = d; }

    void SetReconstructedAnalysisEngine(AnalysisEngine* reco) { fRecontructed = reco; }
    void SetGeneratedAnalysisEngine(AnalysisEngine* truth) { fGenerated = truth; }

    Bool_t CheckInit();
    Bool_t IsInhibit() const { return fInhibit; }

    const char* GetName() const { return fName.Data(); }

    void RunAnalysis();

    TTree* BuildTree(const char* taskName);
    TTree* GetTree() const { return fTree; }
    Bool_t FillTree(Bool_t applyKinCuts);

    void AssignDataSlot(Int_t n) { fDataSlotNumber = n; }
    Int_t GetDataSlotNumber() const { return fDataSlotNumber; }

    friend bool        operator==(const ResponseEngine& lhs, const ResponseEngine& rhs) { return (lhs.fCandidateType == rhs.fCandidateType);}
    friend inline bool operator!=(const ResponseEngine& lhs, const ResponseEngine& rhs) { return !(lhs == rhs); }

  protected:

    ECandidateType_t              fCandidateType       ; ///<  D meson candidate type
    Bool_t                        fInhibit             ; ///<  Inhibit the task
    Double_t                      fMaxJetDmesonDistance; ///<  Maximum distance between a generated D meson and a reconstructed jet, used for geometrical matching (in units of R)

    TString                       fName                ; //!<! D meson candidate name
    TTree                        *fTree                ; //!<! Output tree
    AliDmesonMatchInfoSummary    *fCurrentDmeson       ; //!<! Tree branch
    AliJetInfoSummary           **fCurrentJetInfoReco  ; //!<! Tree branch
    AliJetInfoSummary           **fCurrentJetInfoTruth ; //!<! Tree branch
    Int_t                         fDataSlotNumber      ; //!<! Data slot where the tree output is posted

    AnalysisEngine               *fRecontructed        ; //!<! Reconstructed level analysis engine
    AnalysisEngine               *fGenerated           ; //!<! Generated level analysis engine

    friend class AliAnalysisTaskDmesonJetsDetectorResponse;

  private:

    /// \cond CLASSIMP
    ClassDef(ResponseEngine, 2);
    /// \endcond
  };

  AliAnalysisTaskDmesonJetsDetectorResponse();
  AliAnalysisTaskDmesonJetsDetectorResponse(const char* name, Int_t nOutputTrees=2);
  virtual ~AliAnalysisTaskDmesonJetsDetectorResponse() {;}

  virtual void         UserCreateOutputObjects();
  virtual void         ExecOnce();
  virtual Bool_t       Run();
  virtual Bool_t       FillHistograms();

 protected:

  virtual void SetOutputTypeInternal(EOutputType_t b);

  Int_t PostDataFromResponseEngine(const ResponseEngine& eng);

  std::map<ECandidateType_t, ResponseEngine>  fResponseEngines  ; //!<! Response engines

 private:

  AliAnalysisTaskDmesonJetsDetectorResponse(const AliAnalysisTaskDmesonJetsDetectorResponse &source);
  AliAnalysisTaskDmesonJetsDetectorResponse& operator=(const AliAnalysisTaskDmesonJetsDetectorResponse& source);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDmesonJetsDetectorResponse, 1);
  /// \endcond
};

#endif
