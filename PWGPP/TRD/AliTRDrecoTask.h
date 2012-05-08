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

class TAxis;
class TH1;
class TH2;
class TH3;
class TF1;
class TList;
class TObjArray;
class TTreeSRedirector;
class AliTRDtrackV1;
template <typename Value> class TVectorT;
typedef class TVectorT<Float_t> TVector;
class AliTRDrecoTask : public AliAnalysisTaskSE 
{
friend class AliEveTRDTrackList;
public:
  enum AliTRDrecoSteeringBits{
     kMCdata      = BIT(18)
    ,kFriends     = BIT(19)
    ,kPostProcess = BIT(20)
    ,kHeavyIon    = BIT(21)
  };
  
  class AliTRDrecoProjection : public TNamed
  {
  public:
    AliTRDrecoProjection();
    virtual ~AliTRDrecoProjection();
    AliTRDrecoProjection& operator+=(const AliTRDrecoProjection& other);
    AliTRDrecoProjection& operator=(const AliTRDrecoProjection& other);
    void  Build(const Char_t *n, const Char_t *t, Int_t ix, Int_t iy, Int_t iz, TAxis *aa[]);
    void  Increment(Int_t bin[], Double_t v);
    TH3*  H() const { return fH;}
    TH2*  Projection2D(const Int_t nstat, const Int_t ncol, const Int_t mid=0, Bool_t del=kTRUE);
    void  SetRebinStrategy(Int_t n, Int_t rebx[], Int_t reby[]);
    void  SetShowRange(Float_t zm, Float_t zM, Float_t em=0., Float_t eM=0.) {fRange[0] = zm; fRange[1] = zM; fRange[2] = em; fRange[3] = eM;}
  private:
    AliTRDrecoProjection(const AliTRDrecoProjection&);
  protected:
    TH3  *fH;          // data container
    Int_t fAx[3];      // projection axes
    Int_t fNrebin;     // no. of rebinning steps
    Int_t *fRebinX;    //[fNrebin] rebinning of the X axis
    Int_t *fRebinY;    //[fNrebin] rebinning of the Y axis
    Float_t fRange[4]; //show range of the z processed

    ClassDef(AliTRDrecoProjection, 2)  // wrapper for a projection container THnSparse -> TH3
  };

  AliTRDrecoTask();
  AliTRDrecoTask(const char *name, const char *title);
  virtual ~AliTRDrecoTask();
  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual void   SetDebugLevel(Int_t level);
  
    
  static Float_t GetMeanStat(TH1 *h, Float_t cut=0., Option_t *opt="");
  Int_t          GetNRefFigures() const; 
  const Char_t*  GetNameId() const       { return fNameId;}
  TList*         GetPlotFunctors() const { return fPlotFuncList;}
  static Int_t   GetPtBinSignificant(Float_t pt);
  virtual Bool_t GetRefFigure(Int_t ifig);
  virtual void   MakeSummary();
  void           MakeDetectorPlot(Int_t ly=0, const Option_t *opt="eta");
  void           MakeDetectorPlotNEW(Int_t ly=0, const Option_t *opt="eta");
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
  static void    SetNormZ(TH2 *h2, Int_t bxmin=1, Int_t bxmax=-1, Int_t bymin=1, Int_t bymax=-1, Float_t thr=0.);
  static void    SetRangeZ(TH2 *h2, Float_t m, Float_t M, Float_t thr=0.);
  void SetRunTerminate(Bool_t runTerminate = kTRUE) { fRunTerminate = runTerminate; }
  virtual void   Terminate(Option_t *);

protected:
  static TTreeSRedirector* DebugStream() { return fgDebugStream;}
  virtual void   InitFunctorList();
  void           Adjust(TF1 *f, TH1 * const h);
  Bool_t         HasFunctorList() const { return fPlotFuncList != NULL; }

  Char_t                fNameId[10];       // unique identifier of task particularity
  UChar_t               fNRefFigures;      // no of reference figures reported by task
  TObjArray             *fDets;            //! OLD container to store detector position and status support should be discontinued 
  TVector               *fDetsV;           //! NEW container to store detector position and status
  TObjArray             *fContainer;       //! container to store results
  AliTRDeventInfo       *fEvent;           //! Event Info
  TObjArray             *fTracks;          //! Array of tracks
  TObjArray             *fClusters;        //! Array of clusters
  const TObjArray       *fkClusters;       //! current detector clusters array
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

  TList             *fPlotFuncList;        //! track functors list
  TList             *fDetFuncList;         //! detector functors list
  Bool_t            fRunTerminate;         // Switch for Terminate Function
  static TList      *fgTrendPoint;         //! trend point
  static TTreeSRedirector *fgDebugStream;  //! Debug stream
  static const Int_t fgNPt0 = 4;           // No of significant pt bins 
  static Float_t fgPt0[fgNPt0];            // Array with limits for significant pt bins 

  ClassDef(AliTRDrecoTask, 5) // base TRD reconstruction task
};

#endif

