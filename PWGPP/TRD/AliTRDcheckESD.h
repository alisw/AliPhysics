#ifndef ALITRDCHECKESD_H
#define ALITRDCHECKESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDcheckESD.h 27496 2008-07-22 08:35:45Z cblume $ */

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//   Ionut Arsene <iarsene@cern.ch>
//
//////////////////////////////////////////////////////

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliMCEvent;
class AliESDpid;
class AliCFContainer;
class AliAnalysisCuts;
class AliESDtrack;
class TH1;
class TH2;
class TH1F;
class TH1D;
class TH2F;
class TH3F;
class TH3;
class TObjArray;
class TGraph;
class TGraphErrors;
class TAxis;

class AliTRDcheckESD : public AliAnalysisTaskSE {
public:
  enum ETRDcheckESDstatus {
     kMC        = BIT(0)  // use MC info
    ,kCollision = BIT(1)  // 
  };
    
  enum ETrdCfVariables {
    // Event wise variables
    kEventVtxZ=0,                // event vtx. Z        ---
    kEventMult,                  // event multiplicity  ---
    kEventTrigger,               // trigger class
    kEventBC,                    // event BC        --- 
    // Track wise variables
    kTrackTOFBC,                 // track TOF BC    ---
    kTrackDCAxy,                 // dca xy          ---
    kTrackDCAz,                  // dca z           ---
    kTrackCharge,                // charge          ---
    kTrackOuterParamRadius,      // outer param radius
    kTrackPhi,                   // phi at the vertex  ---
    kTrackPhiTRD,                // phi at the entrance of TRD  ---
    kTrackEta,                   // eta at the vertex  ---
    kTrackEtaTRD,                // eta at the entrance of TRD  ---
    kTrackPt,                    // pt at the vertex     ---
    kTrackPtTRD,                 // pt at the entrance of TRD ---
    kTrackP,                     // p at the vertex      ---
    kTrackPTRD,                  // p at the entrance of TRD  ---
    kTrackTrdChi2,               // TRD track chi2 ---
    kTrackTrdTracklets,          // number of TRD tracklets  ---
    kTrackTrdClusters,           // number of TRD clusters   ---
    kTrackTrdQuality,            // TRD quality for TOF
    kTrackTRDBudget,             // TRD material budget
    kTrackTOFchi2,               // TOF chi2
    // Tracklets wise variables
    kTrackletQtot,               // tracklet qtot in each layer
    kTrackletClustersVsRows,     // clusters / crossed rows
    kTrackletClusters,           // number of clusters per tracklet
    kTrackletP,                  // tracklet p in each layer
    kTrackPlossTRDlayer,         // p loss at each layer
    kTrackletLayer,              // layer of the current tracklet
    // Tracklet slice variables
    kTrackletSlice,              // tracklet slice number
    kTrackletPHslice,            // charge per tracklet slice
    kNTrdCfVariables,
    kNMaxAssignedTriggers = 200
  };
  
  enum ETrdCfSteps {
    kTPCreference=0,
    kTRD,
    kTOF,
    kTOFin,
    kTOFout,
    kNSteps
  };
     
  AliTRDcheckESD();
  AliTRDcheckESD(char* name);
  virtual ~AliTRDcheckESD();
  
  void          UserCreateOutputObjects();
  void          UserExec(Option_t *);

  void          SetRefTrackFilter(AliAnalysisCuts* const filter) {fReferenceTrackFilter = filter;}
  
  Bool_t        HasMC() const { return TESTBIT(fStatus, kMC);}
  Bool_t        IsCollision() const {return TESTBIT(fStatus, kCollision);}
  void          SetCollision(Bool_t set=kTRUE) {set ? SETBIT(fStatus, kCollision) : CLRBIT(fStatus, kCollision);}
  TObjArray*    Histos();
  Int_t         GetTriggerCounter(const Char_t* triggerName) const;
  void          PrintTriggers() const;
  Bool_t        Load(const Char_t *fn="AnalysisResults.root", const Char_t *dir="TRD_Performance", const Char_t *name=NULL);
  void          SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  Bool_t        PutTrendValue(const Char_t *name, Double_t val);
  void          Terminate(Option_t *);
  void          MakeSummaryFromCF(Double_t* trendValues=0x0, const Char_t* triggerName="", Bool_t useIsolatedBC=kFALSE, Bool_t cutTOFbc=kFALSE);
  Int_t         GetNAssignedTriggers();
  void          AddUserTrigger(const Char_t* name) {fUserEnabledTriggers += name; fUserEnabledTriggers += ";";}
  
  void          AddCFContainer(const Char_t* name, const Char_t* title, Int_t nSteps, Int_t* steps, 
		               Int_t nVars, UInt_t* vars, TArrayD* binLimits);
  
private:
  static const Float_t fgkxTPC; // end radial position of TPC
  static const Float_t fgkxTOF; // start radial position of TOF
  static const Char_t* fgkVarNames[kNTrdCfVariables];
  static const Char_t* fgkStepNames[kNSteps];
  
  void PlotTrackingSummaryFromCF(Double_t* trendValues=0x0,
				 const Char_t* triggerName="",
                                 Bool_t useIsolatedBC=kFALSE, Bool_t cutTOFbc=kFALSE);   // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotPidSummaryFromCF(Double_t* trendValues=0x0,
			    const Char_t* triggerName="",
                            Bool_t useIsolatedBC=kFALSE, Bool_t cutTOFbc=kFALSE);        // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotCentSummaryFromCF(Double_t* trendValues=0x0, const Char_t* triggerName="",
                             Bool_t useIsolatedBC=kFALSE, Bool_t cutTOFbc=kFALSE);       // centrality dependent plots
  void PlotOtherSummaryFromCF(Double_t* trendValues);
  
  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  Int_t         Pdg2Idx(Int_t pdg) const;
  void          Process(TH1 **h, TGraphErrors *g);
  void          Process2D(TH2 * const h, TGraphErrors **g);
  void          PrintStatus(ULong_t s);
  TH2F*         Proj3D(TH3* hist, TH2* accMap, Int_t binLow, Int_t binHigh, Float_t &entries);
  TH1D*         Proj2D(TH2* hist, TH1* mpvErr=0x0, TH1* widthErr=0x0, TH1* chi2=0x0);
  TH1F*         EfficiencyTRD(TH3* tpc3D, TH3* trd3D, Bool_t useAcceptance=kTRUE);
  TH1F*         EfficiencyFromPhiPt(AliCFContainer* cf, Int_t minNtrkl, Int_t maxNtrkl, Int_t stepNom, Int_t stepDenom, Int_t var=kTrackPt);
  void          DrawTRDGrid();
  void          SetStyle(TH1* hist, Int_t lineStyle, Int_t lineColor, Int_t lineWidth, 
                         Int_t markerStyle, Int_t markerColor, Int_t markerSize);
  void          SetStyle(TAxis* axis, const Char_t* title, Float_t titleSize, Float_t titleOffset, Bool_t centerTitle, 
                         Float_t labelSize);
  void          CheckActiveSM(TH1D* phiProj, Bool_t activeSM[18]);
  void          FindIsolatedBCs(TH1D* bcHist, Bool_t isIsolated[3500]);
  void          InitializeCFContainers();
  Int_t         GetTriggerIndex(const Char_t* name, Bool_t createNew=kTRUE);
  void          FillEventInfo(Double_t* values);
  void          FillTrackInfo(Double_t* values, AliESDtrack* esdTrack);
  void          FillTrackletInfo(Double_t* values, AliESDtrack* esdTrack, Int_t iPlane,
                                 Double_t* localSagitaPhi, Double_t localMom[][3], Bool_t* localMomGood);
  void          FillTrackletSliceInfo(Double_t* values, AliESDtrack* esdTrack, Int_t iSlice);
  //void          FillCFContainer(AliCFContainer* cf, Double_t* values, Int_t step);
  void          FillCFContainer(AliCFContainer* cf, Double_t* values, Bool_t* stepSelections);
  Bool_t        IsTrackSelected(AliESDtrack* track, Double_t* values, Int_t step);
  void          FillGlobalTrackContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig);
  void          FillTrdTrackletContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig);
  void          FillTrdSliceContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig);
  
  Int_t            fStatus;            // bit mask for controlling the task
  Int_t            fNRefFigures;       // number of current ref plots
  AliESDEvent      *fESD;              //! ESD event
  AliMCEvent       *fMC;               //! MC event
  AliESDpid        *fESDpid;           //  ESD pid object 
  static FILE      *fgFile;            //! trend file streamer
  TObjArray*        fHistos;           //! QA histograms
    
  AliAnalysisCuts* fReferenceTrackFilter;     // reference track filter
  Bool_t           fPhysSelTriggersEnabled;   // flag wheter physics selection triggers were enabled
  TString          fUserEnabledTriggers;      // list of user enabled triggers
  Int_t            fNAssignedTriggers;        // number of assigned triggers
    
  ClassDef(AliTRDcheckESD, 10)          // user oriented TRD analysis based on ESD-MC data
};
#endif
