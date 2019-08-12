/// \class AliRsnMiniAnalysisTask
/// \brief Analysis task for resonance analysis ('mini' package)
/// Analysis task for 'mini' sub-package
/// Contains all definitions needed for running an analysis:
/// -- global event cut
/// -- a list of track cuts (any number)
/// -- definitions of output histograms
/// -- values to be computed.
///
/// \author: Alberto Pulvirenti
/// \author: Francesca Bellini <fbellini@cern.ch>
/// 


#ifndef ALIRSNMINIANALYSISTASK_H
#define ALIRSNMINIANALYSISTASK_H

#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"

#include "AliRsnEvent.h"
#include "AliRsnMiniValue.h"
#include "AliRsnMiniOutput.h"
#include "AliRsnCutEventUtils.h"
#include "AliRsnCutPrimaryVertex.h"
#include "AliRsnMiniResonanceFinder.h"

#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"

class TList;

class AliTriggerAnalysis;
class AliRsnMiniEvent;
class AliRsnCutSet;
class AliQnCorrectionsManager;
class AliQnCorrectionsQnVector;

class AliRsnMiniAnalysisTask : public AliAnalysisTaskSE {

public:

   AliRsnMiniAnalysisTask();
   AliRsnMiniAnalysisTask(const char *name, Bool_t isMC = kFALSE,Bool_t saveRsnTreeInFile=kFALSE);
   AliRsnMiniAnalysisTask(const AliRsnMiniAnalysisTask &copy);
   AliRsnMiniAnalysisTask &operator=(const AliRsnMiniAnalysisTask &copy);
   virtual ~AliRsnMiniAnalysisTask();

   void                UseMC(Bool_t yn = kTRUE)           {fUseMC = yn;}
   void                UseESDTriggerMask(UInt_t trgMask)     {fTriggerMask = trgMask;}
   void                SkipTriggerMask(UInt_t skip)    {fSkipTriggerMask = skip;}
   void                UseCentrality(const char *type)    {fUseCentrality = kTRUE; fCentralityType = type; fCentralityType.ToUpper();}
   void                UseReferenceMultiplicity(const char *type)    {fRefMultiType = type; fRefMultiType.ToUpper();}
   void                SetUseCentralityPatch(Bool_t isAOD049) {fUseAOD049CentralityPatch = isAOD049;}
   void                SetUseCentralityPatchPbPb2011(Int_t centralityPatchPbPb2011) {fUseCentralityPatchPbPb2011 = centralityPatchPbPb2011;}
   void                SetFlowQnVectorSubDet(const char *s) { fFlowQnVectorSubDet = s;}
   void                SetFlowQnVectorExpStep(const char *s) { fFlowQnVectorExpStep = s;}
   void                UseMultiplicity(const char *type)  {fUseCentrality = kFALSE; fCentralityType = type; if(!fCentralityType.Contains("AliMultSelection")) fCentralityType.ToUpper();}
   void                UseContinuousMix()                 {fContinuousMix = kTRUE;}
   void                UseBinnedMix()                     {fContinuousMix = kFALSE;}
   void                SetNMix(Int_t nmix)                {fNMix = nmix;}
   void                SetMaxDiffMult (Double_t val)      {fMaxDiffMult  = val;}
   void                SetMaxDiffVz   (Double_t val)      {fMaxDiffVz    = val;}
   void                SetMaxDiffAngle(Double_t val)      {fMaxDiffAngle = val;}
   void                SetEventCuts(AliRsnCutSet *cuts)   {fEventCuts    = cuts;}
   void                SetMixPrintRefresh(Int_t n)        {fMixPrintRefresh = n;}
   void                SetCheckDecay(Bool_t checkDecay = kTRUE) {fCheckDecay = checkDecay;}
   void                SetMaxNDaughters(Short_t n)        {fMaxNDaughters = n;}
   void                SetCheckMomentumConservation(Bool_t checkP) {fCheckP = checkP;}
   void                SetCheckFeedDown(Bool_t checkFeedDown)      {fCheckFeedDown = checkFeedDown;}
   void                SetDselection(UShort_t originDselection);
   void 	              SetRejectCandidateIfNotFromQuark(Bool_t opt){fRejectIfNoQuark=opt;}
   void                SetMotherAcceptanceCutMinPt(Float_t minPt)  {fMotherAcceptanceCutMinPt = minPt;}
   void                SetMotherAcceptanceCutMaxEta(Float_t maxEta){fMotherAcceptanceCutMaxEta = maxEta;}
   void                KeepMotherInAcceptance(Bool_t keepMotherInAcceptance) {fKeepMotherInAcceptance = keepMotherInAcceptance;}
   void                SaveRsnTreeInFile(Bool_t saveInFile=kTRUE) {fRsnTreeInFile = saveInFile;}
   void                SetComputeSpherocity(Bool_t doit=kTRUE) {fComputeSpherocity = doit;}
   void                SetTrackCuts(AliAnalysisFilter* fTrackFilter);

   Int_t               AddTrackCuts(AliRsnCutSet *cuts);
   TClonesArray       *Outputs()                          {return &fHistograms;}
   TClonesArray       *Values()                           {return &fValues;}
   Short_t             GetMaxNDaughters()                 {return fMaxNDaughters;}
   void                SetEventQAHist(TString type,TH1 *histo);
   void                UseBigOutput(Bool_t b=kTRUE) { fBigOutput = b; }
   Int_t               GetNumberOfTrackCuts() { return fTrackCuts.GetEntries(); }

   virtual void        UserCreateOutputObjects();
   virtual void        UserExec(Option_t *);
   virtual void        Terminate(Option_t *);
   virtual void        FinishTaskOutput();

   Int_t               ValueID(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   Int_t               CreateValue(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   AliRsnMiniOutput   *CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src);
   AliRsnMiniOutput   *CreateOutput(const char *name, const char *outType, const char *compType);

   Int_t               AddResonanceFinder(AliRsnMiniResonanceFinder* f);
   Int_t               GetNResonanceFinders() {return fResonanceFinders.GetEntries();}

private:
   Char_t   CheckCurrentEvent();
   void     FillMiniEvent(Char_t evType);
   Double_t ComputeAngle();
   Double_t ComputeCentrality(Bool_t isESD);
   Double_t ComputeMultiplicity(Bool_t isESD,TString type);
   Double_t ComputeReferenceMultiplicity(Bool_t isESD,TString type);
   Double_t ComputeTracklets();
   Double_t ComputeSpherocity();
   Double_t ApplyCentralityPatchAOD049();
   Double_t ApplyCentralityPatchPbPb2011();
   void     FillTrueMotherESD(AliRsnMiniEvent *event);
   void     FillTrueMotherAOD(AliRsnMiniEvent *event);
   void     StoreTrueMother(AliRsnMiniPair *pair, AliRsnMiniEvent *event);
   Bool_t   EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2);
   AliQnCorrectionsQnVector * GetQnVectorFromList(const TList *list, const char *subdetector, const char *expectedstep) const;

   Bool_t               fUseMC;           ///<  use or not MC info
   Int_t                fEvNum;           ///< absolute event counter
   UInt_t               fTriggerMask;   ///< trigger mask
   UInt_t               fSkipTriggerMask; ///< skip events with this trigger mask, even if they are consistent with fTriggerMask
   Bool_t               fUseCentrality;   ///<  if true, use centrality for event, otherwise use multiplicity
   TString              fCentralityType;  ///<  definition used to choose what centrality or multiplicity to use
   TString              fRefMultiType;    ///< reference multiplicity to use, TRACKLETS (SPD only) or GLOBAL (ITS+TPC)
   Bool_t               fUseAOD049CentralityPatch; ///< flag to enable AOD049 centrality patch
   Int_t                fUseCentralityPatchPbPb2011; ///< for PbPb 2011 centrality flattening
   AliQnCorrectionsManager *fFlowQnVectorMgr; ///< Qn vector manager
   TString              fFlowQnVectorSubDet; ///< Sub Detector used to extract Qn vector, e.g. VZEROA
   TString              fFlowQnVectorExpStep; ///< Desired step in Qn vector calculation
   Bool_t               fContinuousMix;   ///<  mixing --> technique chosen (continuous or binned)
   Int_t                fNMix;            ///<  mixing --> required number of mixes
   Double_t             fMaxDiffMult;     ///<  mixing --> max difference in multiplicity
   Double_t             fMaxDiffVz;       ///<  mixing --> max difference in Vz of prim vert
   Double_t             fMaxDiffAngle;    ///<  mixing --> max difference in reaction plane angle

   TList               *fOutput;          ///< output list
   TClonesArray         fHistograms;      ///< list of histogram definitions
   TClonesArray         fValues;          ///< list of values to be computed
   TH1F                *fHEventStat;      //!<! histogram of event statistics
   TH1F                *fHAEventsVsMulti; //!<! histogram of event statistics
   TH1F                *fHAEventsVsTracklets; //!<! histogram of event statistics
   TH2F                *fHAEventVzCent;       ///< histogram of vertex-z vs. multiplicity/centrality
   TH2F                *fHAEventSpherocityCent; ///< histogram of spherocity vs. multiplicity/centrality
   TH2F                *fHAEventMultiCent;    ///< histogram of multiplicity vs. centrality
   TH2F                *fHAEventRefMultiCent; //!<! histogram of reference multiplicity vs. centrality
   TH2F                *fHAEventPlane;        //!<! histogram of event plane vs. multiplicity/centrality

   AliRsnCutSet        *fEventCuts;       ///< cuts on events
   TObjArray            fTrackCuts;       ///< list of single track cuts
   AliRsnEvent          fRsnEvent;        ///< interface object to the event
   TTree               *fEvBuffer;        //!<! mini-event buffer
   AliTriggerAnalysis  *fTriggerAna;      //!<! trigger analysis
   AliESDtrackCuts     *fESDtrackCuts;    //!<! quality cut for ESD tracks
   AliRsnMiniEvent     *fMiniEvent;       ///< mini-event cursor
   Bool_t               fBigOutput;       ///< flag if open file for output list
   Int_t                fMixPrintRefresh; ///< how often info in mixing part is printed
   Bool_t               fCheckDecay;      ///< check if the mother decayed via the requested channel
   Short_t              fMaxNDaughters;   ///< maximum number of allowed mother's daughter
   Bool_t               fCheckP;          ///< flag to set in order to check the momentum conservation for mothers
   
   Bool_t               fCheckFeedDown;     ///< flag to set in order to check the particle feed down (specific for D meson analysis)
   UShort_t 		      fOriginDselection;  ///< flag to select D0 origins. 0 Only from charm 1 only from beauty 2 both from charm and beauty (specific for D meson analysis)
   Bool_t   		      fKeepDfromB;  	     ///< flag for the feed down from b quark decay (specific for D meson analysis)
   Bool_t   		      fKeepDfromBOnly;    ///< flag to keep only the charm particles that comes from beauty decays (specific for D meson analysis)
   Bool_t 		         fRejectIfNoQuark;   ///< flag to remove events not generated with PYTHIA
   Float_t              fMotherAcceptanceCutMinPt;  ///< cut value to apply when selecting the mothers inside a defined acceptance
   Float_t              fMotherAcceptanceCutMaxEta; ///< cut value to apply when selecting the mothers inside a defined acceptance
   Bool_t               fKeepMotherInAcceptance;    ///< flag to keep also mothers in acceptance
   Bool_t               fRsnTreeInFile;     ///< flag rsn tree should be saved in file instead of memory
   Bool_t               fComputeSpherocity; ///< compute spherocity, false by default since spherocity calculation is time-consuming
   AliAnalysisFilter   *fTrackFilter;       //!<! track filter for spherocity estimator 
   Double_t             fSpherocity;        ///< stores value of spherocity
   TObjArray            fResonanceFinders;  ///< list of AliRsnMiniResonanceFinder objects

/// \cond CLASSIMP
   ClassDef(AliRsnMiniAnalysisTask, 20);     
/// \endcond
};


#endif
