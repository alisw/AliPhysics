#ifndef ALIANALYSISTASKJETSPECTRAAOD_H
#define ALIANALYSISTASKJETSPECTRAAOD_H

class TH1F;
class TH2F;
class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;
class AliHelperPID;
class AliAODJet;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskJetSpectraAOD : public AliAnalysisTaskSE
{
 public:
  
  // constructors
  AliAnalysisTaskJetSpectraAOD() : AliAnalysisTaskSE(),
     fAOD(0),
     fIsMC(0),
     fEventCuts(0), 
     fTrackCuts(0),
     fVZEROside(0),
     fOutput(0),
     fAODJets(0),
     fJetBranchName(""),
     fListJets(0),
     fBackgroundBranch(""),
     fOfflineTrgMask(AliVEvent::kMB),
     fFilterMask(0),
     fJetPtMin(0),
     fJetEtaMin(0x0),
     fJetEtaMax(0x0),
     fLeadPtMin(0x0),
     fnCentBins(20),
     fnQvecBins(20),
     fnptLeadBins(4),
     fIsQvecCalibMode(0),
     fQvecUpperLim(100),
     fIsQvecCut(0),
     fQvecMin(0),
     fQvecMax(100),
     fHistEvtSelection(0x0),
     fDebug(0),
     fMinNcontributors(0),
     fRejectPileup(0),
     fR(0.4),
     fZvertexDiff(1),
     fZvertex(10.)
       {}
  AliAnalysisTaskJetSpectraAOD(const char *name);
  virtual ~AliAnalysisTaskJetSpectraAOD();
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; }
  Bool_t GetIsMC()           const           { return fIsMC;}
  
  AliSpectraAODTrackCuts * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts * GetEventCuts()         {  return fEventCuts; }
  TList * GetOutputList()                         {  return fOutput;    }
  
  void SetTrackCuts(AliSpectraAODTrackCuts * tc)   {   fTrackCuts = tc; }
  void SetEventCuts(AliSpectraAODEventCuts * vc)   {   fEventCuts = vc; }
  void SetnCentBins(Int_t val)                     {   fnCentBins = val; }
  void SetnQvecBins(Int_t val)                     {   fnQvecBins = val; }
  void SetnPtLeadBins(Int_t val)                   {   fnptLeadBins = val; }
  void SetQvecCalibMode(Bool_t mode)               {   fIsQvecCalibMode = mode; }
  void SetQvecUpperLimit(Double_t val)             {   fQvecUpperLim = val; }
  void SetDebugMode (Bool_t val) { fDebug = val; }
  
  //jet getter
  void     GetBranchNames(TString &branch) const { branch = fJetBranchName; }
  void     GetBackgroundBranchNames(TString &branch) const { branch = fBackgroundBranch; }
  Float_t  GetJetPtMin() const { return fJetPtMin; }
  Float_t  GetJetEtaMin() const { return fJetEtaMin; }
  Float_t  GetJetEtaMax() const { return fJetEtaMax; }
  //jet setter
  void     SetBranchNames(const TString &branch);
  void     SetRecBackgroundBranch(const TString &bckbranch);
  void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; } 
  void     SetFilterMask(UInt_t i){fFilterMask = i;}
  void     SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
  void     SetEtaJet(Float_t etamin,Float_t etamax)   { fJetEtaMin = etamin; fJetEtaMax = etamax; }
  
  void SetVZEROside(Int_t side = 0)    {fVZEROside = side; }
  Int_t GetVZEROside()           const           { return fVZEROside;}
  
  void SetQvecCut(Bool_t qcut) { fIsQvecCut = qcut; }
  void SetQvecCutLimits(Float_t qmin,Float_t qmax)   { fQvecMin = qmin; fQvecMax = qmax; }
  
  AliVParticle *LeadingTrackFromJetRefs(AliAODJet* jet);
  void SetleadingTrackPtMin(Double_t val) { fLeadPtMin = val;}
  
  
  void    SetMinNcontributors(Int_t val)  {fMinNcontributors = val;}
  void    SetPileupRejection(Bool_t val)  {fRejectPileup = val;}
  
  void    SetZvertexDiff (Bool_t val) { fZvertexDiff = val; }
  
  void SetJetParameterR(Double_t val) {fR = val;}

  void   UserCreateOutputObjects();
  void   UserExec(Option_t *option);
  void   Terminate(Option_t *);
  
 private:
  
  AliAODEvent           *fAOD;         //! AOD object
  Bool_t          fIsMC;// true if processing MC
  AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts
  AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
  
  Int_t                      fVZEROside;                  // 0: VZERO-A 1: VZERO-C
  
  TList                       * fOutput;        // output list
  
  //jet
  AliAODEvent                 * fAODJets;         //! AOD jet object
  TString                       fJetBranchName;   //  name of jet branches to compare
  TList                       * fListJets;        //! jet lists
  TString                       fBackgroundBranch;
  
  AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline trigs 
   
  UInt_t  fFilterMask;       // filter bit for slecected tracks
  Float_t fJetPtMin;         // minimum jet pT
  Float_t fJetEtaMin;        // lower bound on eta for found jets
  Float_t fJetEtaMax;        // upper bound on eta for found jets
  Float_t fLeadPtMin;        // leading track pt min.
  
  Int_t                            fnCentBins;            // number of bins for the centrality axis
  Int_t                            fnQvecBins;            // number of bins for the q vector axis
  Int_t                            fnptLeadBins;          // number of bins for the pt of leading particle axis
  Bool_t                           fIsQvecCalibMode;      // calib mode for Qvector percentile
  Double_t                         fQvecUpperLim;         // Upper limit for Qvector
  
  Bool_t                           fIsQvecCut;            // Q-vec cut switch
  Double_t                         fQvecMin;              // lower bound for Qvec
  Double_t                         fQvecMax;              // upper bound for Qvec
 
  TH1I  *fHistEvtSelection;                  //! event selection statistic 
  Bool_t fDebug;
  Int_t   fMinNcontributors;
  Bool_t  fRejectPileup;
  
  Double_t fR; // cone radius
  
  Bool_t fZvertexDiff;
  Double_t fZvertex;
  
  AliAnalysisTaskJetSpectraAOD(const AliAnalysisTaskJetSpectraAOD&);
  AliAnalysisTaskJetSpectraAOD& operator=(const AliAnalysisTaskJetSpectraAOD&);
  
  ClassDef(AliAnalysisTaskJetSpectraAOD, 4);
};

#endif
