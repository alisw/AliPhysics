#ifndef ALITRDRECOTASK_H
#define ALITRDRECOTASK_H

///////////////////////////////////////////////////////
//
// Basic class for Performance/Calibration TRD tasks
//
// Author: Alexandru Bercuci, 10/09/2008
//
//////////////////////////////////////////////////

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ALITRDTRACKINFO_H
#include "info/AliTRDtrackInfo.h"
#endif

#ifndef ALITRDEVENTINFO_H
#include "info/AliTRDeventInfo.h"
#endif

class TH1;
class TF1;
class TList;
class TObjArray;
class TTreeSRedirector;
class AliTRDtrackV1;

class AliTRDrecoTask : public AliAnalysisTaskSE 
{
public:
  enum AliTRDrecoSteeringBits{
    kMCdata       = BIT(18)
    ,kFriends     = BIT(19)
    ,kPostProcess = BIT(20)
    ,kHeavyIon     = BIT(21)
  };
  
  AliTRDrecoTask();
  AliTRDrecoTask(const char *name, const char *title);
  virtual ~AliTRDrecoTask();
  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual void   SetDebugLevel(Int_t level);
  
    
  Int_t          GetNRefFigures() const; 
  const Char_t*  GetNameId() const       { return fNameId;}
  TList*         GetPlotFunctors() const { return fPlotFuncList;}
  virtual Bool_t GetRefFigure(Int_t ifig);
  virtual void   MakeSummary();
  void           MakeDetectorPlot(Int_t ly=0);
  Bool_t         IsHeavyIon() const      { return TestBit(kHeavyIon);};
  Bool_t         IsPP() const            { return !TestBit(kHeavyIon);};
  Bool_t         HasFriends() const      { return TestBit(kFriends);};
  Bool_t         HasMCdata() const       { return TestBit(kMCdata);};
  Bool_t         HasPostProcess() const  { return TestBit(kPostProcess);};
  Bool_t         HasRunTerminate() const { return fRunTerminate; }
  virtual TObjArray* Histos()            { return fContainer;}

  virtual Bool_t Load(const Char_t *file = "AnalysisResults.root", const Char_t *dir = "TRD_Performance");
  virtual Bool_t LoadDetectorMap(const Char_t *file = "AnalysisResults.root", const Char_t *dir = "TRD_Performance");
  virtual Bool_t Save(TObjArray * const res);
  virtual Bool_t PostProcess();
  virtual Bool_t PutTrendValue(const Char_t *name, Double_t val);
  virtual void   SetFriends(Bool_t fr = kTRUE) {SetBit(kFriends, fr);}
  virtual void   SetMCdata(Bool_t mc = kTRUE) {SetBit(kMCdata, mc);}
  virtual void   SetNameId(const Char_t *nid) {snprintf(fNameId, 10, "%s", nid);}
  virtual void   SetPostProcess(Bool_t pp = kTRUE) {SetBit(kPostProcess, pp);}
  void SetRunTerminate(Bool_t runTerminate = kTRUE) { fRunTerminate = runTerminate; }
  virtual void   Terminate(Option_t *);

protected:
  static TTreeSRedirector* DebugStream() { return fgDebugStream;}
  void           InitFunctorList();
  void           Adjust(TF1 *f, TH1 * const h);
  Bool_t         HasFunctorList() const { return fPlotFuncList != NULL; }
  Char_t                fNameId[10];       // unique identifier of task particularity
  UChar_t               fNRefFigures;      // no of reference figures reported by task
  TObjArray             *fDets;            //! container to store detector position and status
  TObjArray             *fContainer;       //! container to store results
  AliTRDeventInfo       *fEvent;           //! Event Info
  TObjArray             *fTracks;          //! Array of tracks
  const AliTRDtrackV1   *fkTrack;          //! current track
  const AliTRDtrackInfo::AliMCinfo  *fkMC; //! MC info
  const AliTRDtrackInfo::AliESDinfo *fkESD;//! ESD info
  Char_t                 fSpecies;         //! species index +1 with charge sign
  Float_t                fPt;              //! p_t of the track being analyzed
  Float_t                fPhi;             //! phi of the track being analyzed
  Float_t                fEta;             //! eta of the track being analyzed

private:
  AliTRDrecoTask(const AliTRDrecoTask&);
  AliTRDrecoTask& operator=(const AliTRDrecoTask&);

  TList             *fPlotFuncList;//! plot functors list
  Bool_t            fRunTerminate;  // Switch for Terminate Function
  static TList      *fgTrendPoint;          //! trend point
  static TTreeSRedirector *fgDebugStream;  //! Debug stream 

  ClassDef(AliTRDrecoTask, 5) // base TRD reconstruction task
};

#endif

