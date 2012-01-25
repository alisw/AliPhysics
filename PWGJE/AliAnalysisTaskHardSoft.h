#ifndef AliAnalysisTaskHardSoft_cxx
#define AliAnalysisTaskHardSoft_cxx
 

class TH1F;
class TH2F;
class TH3F;
class TList;
class TNtuple;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#include "TFile.h"
#include "TNtuple.h"

class AliAnalysisTaskHardSoft : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHardSoft(const char *name = "AliAnalysisTaskHardSoft");
  virtual ~AliAnalysisTaskHardSoft() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  virtual void   SetUseMc(Bool_t useMc)                    {fUseMc = useMc;}  
  virtual void   SetUseNeutralMC(Bool_t  useNeutral)       {fUseNeutral = useNeutral;}    

  virtual void   SetRadiusCut(Float_t radiusCut)           {fRadiusCut = radiusCut;}  
  virtual void   SetTriggerPtCut(Float_t triggerPtCut)     {fTriggerPtCut = triggerPtCut;}  
  virtual void   SetAssociatePtCut(Float_t associatePtCut) {fAssociatePtCut = associatePtCut;}  

  virtual void   SetMap(TH2F* map)                         {fMap = map;}  
  virtual void   SetMapLeading(TH2F* mapLeading)           {fMapLeading = mapLeading;}  
  
  
  virtual void   SetCuts(AliESDtrackCuts* cuts)
  {fCuts = cuts;}

  virtual void   SetFieldOn(Bool_t b = kTRUE){fFieldOn = b;} 

  
 private:

  Bool_t       fUseMc;                     // for simulated data: calculate the same for MCpartciles aswell
  Bool_t       fUseNeutral;                // use also neutral particle in MC case
  Float_t      fRadiusCut;                 // radius cut for hard/soft event determination (CDF approach)
  Float_t      fTriggerPtCut;              // first pt cut for hard/soft event determination (CDF approach)
  Float_t      fAssociatePtCut;            // first pt cut for hard/soft event determination (CDF approach)
  
  TH2F       * fMap;                       // map for eta-phi acceptance (used for random particles)
  TH2F       * fMapLeading;                // map for eta-phi acceptance (used for random particles)

  AliESDtrackCuts* fCuts;                  // List of cuts for ESDs
  Bool_t      fFieldOn;


  TList       * fHists;                    // List of histos

  //properties of particles(0)/esdtracks(1)
  TH1F       * fPt[2];                     // pt 
  TH1F       * fEta[2];                    // eta
  TH1F       * fPhi[2];                    // phi
  TH2F       * fEtaPhi[2];                 // eta-phi (needed as input for random position -> will be fMap)
  TH2F       * fEtaPhiLeading[2];          // eta-phi (needed as input for random position -> will be fMapLeading)
  TH1F       * fNch[2];                    // all accepted tracks/particles
  TH1F       * fPtLeading[2];              // pt of leading track/particle

  TProfile   * fPtLeadingNch[2];           // pt of leading track/particle vs Nch
  TProfile   * fPtSumNch[2];               // pt sum track/particle vs Nch
  TProfile   * fPtAvNch[2];                // average pt track/particle vs Nch

  TH1F       * fDPhiLeading[2];            // delta phi of associate tracks to leading track
  TH1F       * fRadiusLeading[2];          // radius of associate tracks to leading track

  TH1F       * fDPhiLeadingR[2];           // delta phi of associate tracks to leading track for random pos.
  TH1F       * fRadiusLeadingR[2];         // radius of associate tracks to leading track for random pos.

  TH1F       * fDPhiLeadingRS[2];          // delta phi of associate tracks to leading track for random seed pos.
  TH1F       * fRadiusLeadingRS[2];        // radius of associate tracks to leading track for random seed pos

  TProfile   * fNchAssInR[2];              // number of tracks within R around leading track vs Nch
  TH1F       * fTrigger[2];                // number of triggers with at accepted-track number

  //per Nch bin of all accepted tracks
  TH1F       * fDPhiLeadingNchBin[2][100]; // delta phi of associate tracks to leading track per Nch bin

  TH1F       * fNchHardSoft[2][2];         // Nch for hard and soft events (classified with CDF algorithm)
  TH1F       * fPtHardSoft[2][2];          // Pt for hard and soft events (classified with CDF algorithm)
  TProfile   * fPtAvHardSoftNch[2][2];     // <Pt> for hard and soft events (classified with CDF algorithm)


  AliAnalysisTaskHardSoft(const AliAnalysisTaskHardSoft&); // not implemented
  AliAnalysisTaskHardSoft& operator=(const AliAnalysisTaskHardSoft&); // not implemented
  
  ClassDef(AliAnalysisTaskHardSoft, 1);    // Hard and Soft event characteristics
};

#endif
