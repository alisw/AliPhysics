#ifndef AliRsnMiniTaskPhiVn_H
#define AliRsnMiniTaskPhiVn_H

//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
//

#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliRsnEvent.h"
#include "AliRsnMiniValue.h"
#include "AliRsnMiniOutput.h"

class TList;

class AliTriggerAnalysis;
class AliRsnMiniEvent;
class AliRsnCutSet;

class AliRsnMiniTaskPhiVn : public AliAnalysisTaskSE {

public:

   AliRsnMiniTaskPhiVn();
   AliRsnMiniTaskPhiVn(const char *name, Bool_t isMC = kFALSE);
   AliRsnMiniTaskPhiVn(const AliRsnMiniTaskPhiVn &copy);
   AliRsnMiniTaskPhiVn &operator=(const AliRsnMiniTaskPhiVn &copy);
   virtual ~AliRsnMiniTaskPhiVn();

   void                UseMC(Bool_t yn = kTRUE)           {fUseMC = yn;}
   void                UseESDTriggerMask(UInt_t trgMask)     {fTriggerMask = trgMask;}
   void                UseCentrality(const char *type)    {fUseCentrality = kTRUE; fCentralityType = type; fCentralityType.ToUpper();}
   void                SetUseCentralityPatch(Bool_t isAOD049) {fUseAOD049CentralityPatch = isAOD049;}
   void                SetUseCentralityPatchPbPb2011(Int_t centralityPatchPbPb2011) {fUseCentralityPatchPbPb2011 = centralityPatchPbPb2011;}
   void                UseMultiplicity(const char *type)  {fUseCentrality = kFALSE; fCentralityType = type; fCentralityType.ToUpper();}
   void                UseContinuousMix()                 {fContinuousMix = kTRUE;}
   void                UseBinnedMix()                     {fContinuousMix = kFALSE;}
   void                SetNMix(Int_t nmix)                {fNMix = nmix;}
   void                SetMaxDiffMult (Double_t val)      {fMaxDiffMult  = val;}
   void                SetMaxDiffVz   (Double_t val)      {fMaxDiffVz    = val;}
   void                SetMaxDiffAngle(Double_t val)      {fMaxDiffAngle = val;}
   void                SetEventCuts(AliRsnCutSet *cuts)   {fEventCuts    = cuts;}
   void                SetMixPrintRefresh(Int_t n)        {fMixPrintRefresh = n;}
   void                SetMaxNDaughters(Short_t n)        {fMaxNDaughters = n;}
   void                SetCheckMomentumConservation(Bool_t checkP) {fCheckP = checkP;}
   void                SetCheckFeedDown(Bool_t checkFeedDown)      {fCheckFeedDown = checkFeedDown;}
   void                SetDselection(UShort_t originDselection);
   void                SetRejectCandidateIfNotFromQuark(Bool_t opt){fRejectIfNoQuark=opt;}
   void                SetMotherAcceptanceCutMinPt(Float_t minPt)  {fMotherAcceptanceCutMinPt = minPt;}
   void                SetMotherAcceptanceCutMaxEta(Float_t maxEta){fMotherAcceptanceCutMaxEta = maxEta;}
   void                KeepMotherInAcceptance(Bool_t keepMotherInAcceptance) {fKeepMotherInAcceptance = keepMotherInAcceptance;}
   void                SetNharmToProcess (Int_t value){fnHarmToProcess = value; return;}//new
   void                SetSelCharge (Char_t charge1, Char_t charge2){fCharge[0] = charge1; fCharge[1] = charge2;}
   void                SetSelPid (RSNPID pid1, RSNPID pid2){fDaughter[0] = pid1; fDaughter[1] = pid2;}
   Int_t               AddTrackCuts(AliRsnCutSet *cuts);
   TClonesArray       *Outputs()                          {return &fHistograms;}
   TClonesArray       *Values()                           {return &fValues;}
   Short_t             GetMaxNDaughters()                 {return fMaxNDaughters;}
   void                SetEventQAHist(TString type,TH2F *histo);
   void                UseBigOutput(Bool_t b=kTRUE) { fBigOutput = b; }

   virtual void        UserCreateOutputObjects();

   virtual void        UserExec(Option_t *);
   virtual void        Terminate(Option_t *);
   virtual void        FinishTaskOutput();

   RSNPID              GetDaughter(Int_t i) const {if (i <= 0) return fDaughter[0]; else return fDaughter[1];}//new
   Double_t            GetMass(Int_t i) const {return AliRsnDaughter::SpeciesMass(GetDaughter(i));}//new


   Int_t               ValueID(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   Int_t               PtIndex(Double_t pt, Double_t ptmax, Double_t ptmin, Double_t totalBinsPt);
   Int_t               InvMassIndex(Double_t invmass, Double_t invmassmax, Double_t invmassmin, Double_t totalBinsInvMass);
   Int_t               CentralityIndex(Double_t centrality, Double_t centmax, Double_t centmin, Double_t totalBinsCentrality);

   Int_t               CreateValue(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   AliRsnMiniOutput   *CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src);
   AliRsnMiniOutput   *CreateOutput(const char *name, const char *outType, const char *compType);

private:

   Char_t   CheckCurrentEvent();
   void     FillMiniEvent(Char_t evType);
   Double_t ComputeAngle();
   Double_t ComputeCentrality(Bool_t isESD);
   Double_t ComputeMultiplicity(Bool_t isESD,TString type);
   Double_t ComputeTracklets();
   Double_t ApplyCentralityPatchAOD049();
   Double_t ApplyCentralityPatchPbPb2011();
   void     FillTrueMotherESD(AliRsnMiniEvent *event);
   void     FillTrueMotherAOD(AliRsnMiniEvent *event);
   void     StoreTrueMother(AliRsnMiniPair *pair, AliRsnMiniEvent *event);
   Bool_t   EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2);


   Bool_t               fUseMC;           //  use or not MC info
   Int_t                fEvNum;           //! absolute event counter
   UInt_t               fTriggerMask;   //trigger mask
   Bool_t               fUseCentrality;   //  if true, use centrality for event, otherwise use multiplicity
   TString              fCentralityType;  //  definition used to choose what centrality or multiplicity to use
   Bool_t               fUseAOD049CentralityPatch; //flag to enable AOD049 centrality patch
   Int_t                fUseCentralityPatchPbPb2011; //for PbPb 2011 centrality flattening

   Bool_t               fContinuousMix;   //  mixing --> technique chosen (continuous or binned)
   Int_t                fNMix;            //  mixing --> required number of mixes
   Double_t             fMaxDiffMult;     //  mixing --> max difference in multiplicity
   Double_t             fMaxDiffVz;       //  mixing --> max difference in Vz of prim vert
   Double_t             fMaxDiffAngle;    //  mixing --> max difference in reaction plane angle

   RSNPID fDaughter[2];                   //  species of daughters, used to assign mass
   Char_t fCharge[2];                     //  required track charge  
   Bool_t fUseStoredMass[2];


   TList               *fOutput;          //  output list
   TClonesArray         fHistograms;      //  list of histogram definitions
   TClonesArray         fValues;          //  list of values to be computed
   TH1F                *fHEventStat;      //  histogram of event statistics
   TH1F                *fHAEventsVsMulti; //  histogram of event statistics
   TH1F                *fHAEventsVsTracklets; //  histogram of event statistics
   TH1F                *fHQVectorPosReTest;//new test histo
   TH1F                *fHQVectorPosImTest;  
   TH1F                *fHQVectorNegReTest;
   TH1F                *fHQVectorNegImTest; 
   TH1F                *fHCentrality;
   TH1F                *fHMass1;
   TH1F                *fHMass2;
   TH1F                *fHMassPhi;
   TH3F                *fHTestPtMassCentrality;
   TH1F                *fHQQVectorDenominator[10];//new: array of histograms of qvector with dimension 10 as maximal number of harmonic, which will never be exceeded. 
   TProfile            *fHQQVectorCentralityNorm[10];
   TProfile            *fHQQVectorCentrality[10];
   TProfile            *fHPposQnegCharged[10][10];
   TProfile            *fHPposQnegChargedWeight[10][10];
   TProfile            *fHPnegQposCharged[10][10];
   TProfile            *fHPnegQposChargedWeight[10][10];
 //  TProfile            *fHPVectorTest[10][10];
   TH1F                *fHPVectorPosTest[10];
   TH1F                *fHPVectorNegTest[10];
  // TProfile2D          *fHPposQnegNoNorm[10][10];// 10 centrality bins
   //TProfile2D          *fHPnegQposNoNorm[10][10];// 10 centrality bins
   TProfile2D          *fHPposQneg[10][10];// 10 centrality bins
   TProfile2D          *fHPnegQpos[10][10];// 10 centrality bins   
   TProfile2D          *fHPposQnegWeight[10][10];// 10 centrality bins
   TProfile2D          *fHPnegQposWeight[10][10];
   TH2F                *fHAEventVz;       //  histogram of vertex-z vs. multiplicity/centrality
   TH2F                *fHAEventMultiCent;//  histogram of multiplicity vs. centrality
   TH2F                *fHAEventPlane;    //  histogram of event plane vs. multiplicity/centrality

   TArrayI              fSel1;//! list of selected particles for definition 1    
   TArrayI              fSel2;//! list of selected particles for definition 2
   AliRsnCutSet        *fEventCuts;       //  cuts on events
   TObjArray            fTrackCuts;       //  list of single track cuts
   AliRsnEvent          fRsnEvent;        //! interface object to the event
   TTree               *fEvBuffer;        //! mini-event buffer
   AliTriggerAnalysis  *fTriggerAna;      //! trigger analysis
   AliESDtrackCuts     *fESDtrackCuts;    //! quality cut for ESD tracks
   AliRsnMiniEvent     *fMiniEvent;       //! mini-event cursor
   Double_t             fMotherMass; //  nominal resonance mass
   Bool_t               fBigOutput;       // flag if open file for output list
   Int_t                fMixPrintRefresh; // how often info in mixing part is printed
   Short_t              fMaxNDaughters;   // maximum number of allowed mother's daughter
   Bool_t               fCheckP;          // flag to set in order to check the momentum conservation for mothers
   
   Bool_t               fCheckFeedDown;      // flag to set in order to check the particle feed down (specific for D meson analysis)
   UShort_t             fOriginDselection;   // flag to select D0 origins. 0 Only from charm 1 only from beauty 2 both from charm and beauty (specific for D meson analysis)
   Bool_t               fKeepDfromB;        // flag for the feed down from b quark decay (specific for D meson analysis)         
   Bool_t               fKeepDfromBOnly;     // flag to keep only the charm particles that comes from beauty decays (specific for D meson analysis)
   Bool_t               fRejectIfNoQuark;    // flag to remove events not generated with PYTHIA
   Float_t              fMotherAcceptanceCutMinPt;              // cut value to apply when selecting the mothers inside a defined acceptance
   Float_t              fMotherAcceptanceCutMaxEta;             // cut value to apply when selecting the mothers inside a defined acceptance
   Bool_t               fKeepMotherInAcceptance;                // flag to keep also mothers in acceptance
   Int_t                fnHarmToProcess;//new

   ClassDef(AliRsnMiniTaskPhiVn, 1);   // AliRsnMiniTaskPhiVn
};


#endif
