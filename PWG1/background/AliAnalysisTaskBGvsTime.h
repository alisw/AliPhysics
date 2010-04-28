#ifndef ALIANALYSISTASKBGVSTIME_H
#define ALIANALYSISTASKBGVSTIME_H



#include "AliAnalysisTaskSE.h"
#include "TH2F.h"


//-------------------------------------------------------------------------
//                      AliAnalysisTaskBGvsTime
// 
// This Task produces histograms of rate of events vs time and
// distributions (pt,vz,multiplicity) for all trigger classes and for
// different selections. Those are used for luminosity and background
// studies. The class does not need to know in advance the trigger
// classes: histos are created on the fly. This complicates the
// merging a bit, which is taken care of by AliHistoListWrapper
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------


class AliESDEvent;
class TF1;
class AliMCParticle;
class AliVParticle;
class AliESDVertex;
class AliHistoListWrapper;
class AliPhysicsSelection;
class AliESDtrack;
class AliAnalysisTaskBGvsTime : public AliAnalysisTaskSE {

public:

  enum {kDistSPDMult, kDistTPCMult, kDistPtLoose, kDistPt, kDistVertex, kDistVertexZ, kDistVertex3D, kDistDCATPC, kDistClsITSLayer}; // Used to get and book histograms of "distributions" (not vs time")

  enum {kEffPt05, kEffPt1, kNEff}; // used to bin histo for efficiency calculation.

  enum {kEffStepGen, kEffStepTrig, kEffStepRec, kEffNSteps}; // used to bin histo for efficiency calculation.

  AliAnalysisTaskBGvsTime();
  AliAnalysisTaskBGvsTime(const char * name);
  AliAnalysisTaskBGvsTime(const AliAnalysisTaskBGvsTime& obj) ;
  ~AliAnalysisTaskBGvsTime();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  const char * GetVsTimeHistoName(const char * triggerClass, Int_t nclusters, Int_t bx, const char * v0flag, const char * prefix = "");
  const char * GetVsTimeHistoNameAll(const char * triggerClass);

  const char * GetVsTimeHistoForLuminosityName(const char * triggerClass, Float_t ptmin);

  TH1F * GetVsTimeHisto(const char * triggerClass, Int_t nclusters, Int_t bx, const char * v0flag, const char * prefix = "") 
  {return GetVsTimeHisto(GetVsTimeHistoName(triggerClass,nclusters,bx,v0flag,prefix));}
  TH1F * GetVsTimeHistoAll(const char * triggerClass)
  {return GetVsTimeHisto(GetVsTimeHistoNameAll(triggerClass));}

  TH1F * GetVsTimeHisto(const char * name) ;
  
  TH1F * BookVsTimeHisto(const char * name, const char * title);  

//   TH1F * GetMultHisto(const char * triggerClass) ;
//   TH1F * BookMultHisto(const char * name, const char * title);  

  TH1F * GetDistributionHisto (const char * triggerClass, Int_t dist, const char * suffix = 0) ;
  TH1F * BookDistributionHisto(const char * name, const char * title,const char * xtitle, Int_t nbin, Float_t min, Float_t max);  

  TH1F * GetDeadTimeHisto(const char * triggerClass); //dead time vs time stamp
//   TH2F * BookDeadTimeHisto(const char * name, const char * title);   

  void SetTimes(UInt_t start, UInt_t end)      { fStartTime=start;   fEndTime =end;   }
  void SetMultBins(Int_t nbins, Int_t * bins) { fNMultBins = nbins; fMultBins = bins;}
  void SetBinWidth(Int_t deltatime) { fBinWidth = deltatime;} // set bin width in seconds

  void SetNoCuts(Bool_t cuts = kTRUE) { fNoCuts = cuts;}
  void SetUsePhysicsSelection(Bool_t sel = kTRUE) { fUsePhysicsSelection = sel;}
  void SetUseZeroBin(Bool_t sel = kTRUE) { fUseZeroBin = sel;} 
  void SetSkipV0(Bool_t sel = kTRUE) { fSkipV0 = sel;} 
  void SetSkipZeroBin(Bool_t sel = kTRUE) { fSkipZeroBin = sel;} // skip zero bin events
  void SetMC(Bool_t sel = kTRUE) { fIsMC = sel;} // analyze MC
  void SetUseBI(Bool_t useBunchInt = kTRUE){fUseBunchInt=useBunchInt;}

  Bool_t IsEventInBinZero(); // returns true if the event has to be put in the bin0.
  Bool_t SelectOnImpPar(AliESDtrack* t) ;

  TH1F * GetEfficiencyHisto(Int_t step);  // if the first argument is 0, returns the histo at generation level, 1 returns the histo after the HW trigger, 2 returns the histo at rec level. Fill this histo with elements of the enum kEffPt1, kEffPt05, ...


private:

  //
  AliESDEvent *  fESD;    //! ESD object  AliVEvent*     fEvent;
  TList * fListHisto;     // list of output object
  AliHistoListWrapper * fListWrapper; // wrapper for the list, takes care of merging

  UInt_t fStartTime;    // run start time, used to fill and book histos 
  UInt_t fEndTime;      // run end time, used to book histos
  
  Int_t    fNMultBins;  // number of multiplicity bins
  Int_t * fMultBins;    //[fNMultBins] edges of multiplicity bins
  
  Long64_t fFirstTimeStamp; // time stamp of first event found (works only locally)
  Long64_t fLastTimeStamp;  // time stamp of first event found (works only locally)

  Int_t fBinWidth; // Bin width of the histos vs time (in seconds)

  Bool_t fNoCuts; // if true, no selection is applied

  AliPhysicsSelection * fPhysicsSelection; // physics selection
  Bool_t fUsePhysicsSelection;   // use physics selection to actually select events
  Bool_t fUseZeroBin;   // use only the zero bin
  Bool_t fIsMC; // true if running on MC
  Bool_t fSkipV0; // Ignore V0 information
  Bool_t fSkipZeroBin;  // skip the events in the zero bin
  Bool_t fUseBunchInt; // use the BI information to compute Background

  TH2F * fHistoTimeStampVsUTC; // correlation between timestamp computed with bx, period, orbit and GDC time
  TH1F * fHistoTimeStampDiffUTC; // difference  between timestamp computed with bx, period, orbit and GDC time

  AliAnalysisTaskBGvsTime& operator=(const AliAnalysisTaskBGvsTime& task);
  
  ClassDef(AliAnalysisTaskBGvsTime, 2)


};

#endif
