#ifndef ALIDNDPT_H
#define ALIDNDPT_H

//------------------------------------------------------------------------------
// Abstract class for dNdPt analysis. All dNdPt components should derive from it.   
// 
// Author: J.Otwinowski 03/11/2008 
// last change: 2011-06-13 by M.Knichel
//------------------------------------------------------------------------------

class AliESDEvent; 
class AliMCEvent; 
class AliESDtrackCuts; 
class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;
class AliPhysicsSelection;
class AlidNdPtBackgroundCuts;

#include "TNamed.h"
#include "TFolder.h"
#include "AliTriggerAnalysis.h"
#include "AlidNdPtHelper.h"

class AlidNdPt : public TNamed {
public:
  AlidNdPt(); 
  AlidNdPt(Char_t* name, Char_t* title);
  ~AlidNdPt();

  // Init data members
  virtual void Init() = 0;

  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list=0) = 0;

  // Analyse output histograms 
  virtual void Analyse() = 0;

  // Export analysed output objects to folder
  virtual TFolder *ExportToFolder(TObjArray * const array=0) = 0;

  //

  //
  void SetEventCuts(AlidNdPtEventCuts* const cuts)              { fdNdPtEventCuts = cuts; }
  void SetAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts)    { fdNdPtAcceptanceCuts = cuts; }  
  void SetRecAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts) { fdNdPtRecAcceptanceCuts = cuts; }
  void SetMultAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts){ fMultAcceptanceCuts = cuts; }  
  void SetTrackCuts(AliESDtrackCuts* const cuts)                { fEsdTrackCuts = cuts; }
  void SetMultTrackCuts(AliESDtrackCuts* const cuts)            { fMultTrackCuts = cuts; }
  void SetUseMCInfo(const Bool_t info)                          { fUseMCInfo = info; }
  void SetUsePileUpRejection(Bool_t pileup)			{ fUsePileUpRejection = pileup;}  
  void SetUseSPDClusterVsTrackletRejection(Bool_t spdclutrkls)  { fUseSPDClusterVsTrackletRejection = spdclutrkls; }
  void SetUseTOFBunchCrossing(Bool_t tofbunchcross)             { fUseTOFBunch = tofbunchcross; }
  void SetUseTOFExpectedTimeDiff(Bool_t tofexptimediff)         { fUseTOFExpectedTimeDiff = tofexptimediff; }
  void SetUseKinkMotherRejection(Bool_t kinkmother)             { fUseKinkMother = kinkmother; }               
  void SetAnalysisMode(const AlidNdPtHelper::AnalysisMode mode) { fAnalysisMode = mode; }
  void SetTrigger(const AliTriggerAnalysis::Trigger trigger)    { fTrigger = trigger; }
  void SetTriggerClass(const Char_t *triggerClass)              { fTriggerClass = triggerClass; }
  void SetParticleMode(const AlidNdPtHelper::ParticleMode mode) { fParticleMode = mode; }
  void SetPhysicsTriggerSelection(AliPhysicsSelection* const selection)  { fPhysicsSelection = selection; }
  void SetBackgroundCuts(AlidNdPtBackgroundCuts* const cuts)    { fdNdPtBackgroundCuts = cuts; }
  void SetRequireCompleteDAQ(const Bool_t req = kFALSE){ f2015IsIncompleteDAQ = req;     }

  
  AlidNdPtEventCuts* GetEventCuts() const                       { return fdNdPtEventCuts; }
  AlidNdPtAcceptanceCuts* GetAcceptanceCuts() const             { return fdNdPtAcceptanceCuts; }
  AlidNdPtAcceptanceCuts* GetMultAcceptanceCuts() const         { return (fMultAcceptanceCuts) ? fMultAcceptanceCuts : fdNdPtAcceptanceCuts; }
  AlidNdPtAcceptanceCuts* GetRecAcceptanceCuts() const          { return fdNdPtRecAcceptanceCuts; }    
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  AliESDtrackCuts* GetMultTrackCuts() const                     { return (fMultTrackCuts) ? fMultTrackCuts: fEsdTrackCuts; }  
  Bool_t IsUseMCInfo() const                                    { return fUseMCInfo; }
  Bool_t IsUsePileUpRejection() const				{ return fUsePileUpRejection;}
  Bool_t IsUseSPDClusterVsTrackletRejection() const             { return fUseSPDClusterVsTrackletRejection;}
  Bool_t IsUseTOFBunchCrossing() const                          { return fUseTOFBunch;}
  Bool_t IsUseTOFExpectedTimeDiff() const                       { return fUseTOFExpectedTimeDiff;}
  Bool_t IsRequiredCompleteDAQ() const  { return f2015IsIncompleteDAQ;}
  Bool_t IsUseKinkMotherReject() const                          { return fUseKinkMother; }
  AlidNdPtHelper::AnalysisMode GetAnalysisMode() const          { return fAnalysisMode; }
  AliTriggerAnalysis::Trigger GetTrigger() const                { return fTrigger; }
  const Char_t* GetTriggerClass() const                         { return fTriggerClass; }
  AlidNdPtHelper::ParticleMode GetParticleMode() const          { return fParticleMode; }
  AliPhysicsSelection* GetPhysicsTriggerSelection() const       { return fPhysicsSelection; }
  AlidNdPtBackgroundCuts* GetBackgroundCuts() const             { return fdNdPtBackgroundCuts; }
  Double_t* CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax);

  void SetAnalyseOutput(Bool_t analyseoutput)  { fAnalyseOutput = analyseoutput; }
  Bool_t GetAnalyseOutput() const { return fAnalyseOutput; }

  void SetMergeTHnSparse(Bool_t mergethnsparse)  { fMergeTHnSparse = mergethnsparse; }
  Bool_t GetMergeTHnSparse() const { return fMergeTHnSparse; }

  void SetTriggerMask(UInt_t triggermask)  { fTriggerMask = triggermask; }
  UInt_t GetTriggerMask()  { return fTriggerMask; }

protected:
   static Double_t* CloneArray(Int_t n, Double_t* source);

private:

  AlidNdPt(const AlidNdPt&); // not implemented
  AlidNdPt& operator=(const AlidNdPt&); // not implemented

  AlidNdPtEventCuts      *fdNdPtEventCuts;      // event cuts
  AlidNdPtAcceptanceCuts *fdNdPtAcceptanceCuts; // acceptance cuts  
  AlidNdPtAcceptanceCuts *fdNdPtRecAcceptanceCuts; // additional recontruction acceptance cuts (not used for MC truth)
  AlidNdPtAcceptanceCuts *fMultAcceptanceCuts; // acceptance cuts for multiplicity estimator
  AliESDtrackCuts *fEsdTrackCuts;               // esd track cuts
  AliESDtrackCuts *fMultTrackCuts;               // esd track cuts for multiplicity estimator

  Bool_t fUseMCInfo;                            // use MC information
  Bool_t fUsePileUpRejection;			// use Pile up rejection (pp)
  Bool_t fUseSPDClusterVsTrackletRejection;     // use SPD Clusters vs Tracklet
  Bool_t fUseTOFBunch;                          // use TOF bunch crossing
  Bool_t fUseTOFExpectedTimeDiff;		// use TOF expected time difference
  Bool_t fUseKinkMother;                        // use Rejection of Kink Mothers
  Bool_t f2015IsIncompleteDAQ;                  // Needed for Run2-2015 data analysis, where a problem with incomplete events from the daq-side exists -- kFALSE will reject incomplete events
  AlidNdPtHelper::AnalysisMode fAnalysisMode;   // analysis mode TPC only, TPC + ITS
  AliTriggerAnalysis::Trigger fTrigger;         // trigger definition MB1, MB2 ...
  const Char_t * fTriggerClass;                 // trigger class
  AlidNdPtHelper::ParticleMode fParticleMode;   // selected particle (pion, kaon, ...)

  AliPhysicsSelection* fPhysicsSelection; // physics trigger selection class
  AlidNdPtBackgroundCuts *fdNdPtBackgroundCuts; // background cuts (cosmics and splitted tracks)
  
 
  Bool_t fAnalyseOutput;  // call Analyse() function in the FinishTaskOutput
  Bool_t fMergeTHnSparse; // merge THnSparse histograms in Merge() function

  UInt_t fTriggerMask;    // trigger mask

  ClassDef(AlidNdPt,8);
};

#endif
