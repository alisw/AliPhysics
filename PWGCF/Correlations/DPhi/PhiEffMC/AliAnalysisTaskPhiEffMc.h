#ifndef ALIANALYSISTASKPHIEFFMC_H
#define ALIANALYSISTASKPHIEFFMC_H

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"
#include "AliHelperPID.h"
#include "AliAnalyseLeadingTrackUE.h"
class AliHelperPID;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;
class AliAnalyseLeadingTrackUE;

class AliAnalysisTaskPhiEffMc : public AliAnalysisTaskSE
{
 public:
  
  // constructors
 AliAnalysisTaskPhiEffMc() : AliAnalysisTaskSE(), 
    fAOD(0x0),
    fIsMC(0),
    fOutput(0x0),
    fHelperPID(0x0),
    fTrackCuts(0x0),
    fEventCuts(0x0)
      {}
  AliAnalysisTaskPhiEffMc(const char *name);
  virtual ~AliAnalysisTaskPhiEffMc() {}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; }
  Bool_t GetIsMC()           const           { return fIsMC;}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHelperPID(AliHelperPID* pid)                     { fHelperPID = pid; }
  AliHelperPID                   * GetHelperPID()          { return fHelperPID; }

  AliSpectraAODTrackCuts      * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts      * GetEventCuts()         {  return fEventCuts; }
  void SetTrackCuts(AliSpectraAODTrackCuts * tc)       { fTrackCuts = tc; }
  void SetEventCuts(AliSpectraAODEventCuts * vc)       { fEventCuts = vc; }

  void UnlikeSign(TObjArray* kaonsPos, TObjArray* kaonsNeg, THnSparseF* h, Double_t cent = -1);
  void LikeSign(TObjArray* kaonsPos, TObjArray* kaonsNeg, THnSparseF* h, Double_t cent = -1);
  TLorentzVector* makePhi(AliVParticle* p1, AliVParticle* p2);

  void SetPtCut(Double_t ptcut){fPtCut = ptcut;}
  Double_t GetPtCut(){return fPtCut;}

 private:
  
  AliAODEvent           *fAOD;         //! AOD object
  Bool_t          fIsMC;// true if processing MC
  TList *fOutput; //! tlist with output  
  AliHelperPID                * fHelperPID;     // points to class for PID
  AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts

  Double_t fPtCut;  // min pt cut on tracks (pt cut in AliSpectraAODTrackCuts is max pt)

  AliAnalysisTaskPhiEffMc(const AliAnalysisTaskPhiEffMc&);
  AliAnalysisTaskPhiEffMc& operator=(const AliAnalysisTaskPhiEffMc&);
  
  ClassDef(AliAnalysisTaskPhiEffMc, 2);
};

#endif

