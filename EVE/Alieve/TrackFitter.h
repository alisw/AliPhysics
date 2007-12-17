// $Header$

#ifndef ALIEVE_TrackFitter_H
#define ALIEVE_TrackFitter_H

#include <Reve/PointSet.h>
#include <TQObject.h>
#include <map>

class TGraphErrors;
class TGraph;
class AliRieman;

namespace Reve 
{
class TrackList;
}

namespace Alieve {

class TrackFitter : public Reve::PointSet
{
private:
  TrackFitter(const TrackFitter&);            // Not implemented
  TrackFitter& operator=(const TrackFitter&); // Not implemented

  TGraph       *fGraphSelected;  // graph of selected points
  TGraphErrors *fGraphFitted;    // graph of fitted points

protected:
  struct Point_t
  {
    // inner structure to check duplicates
    Reve::PointSet* fPS;   // selected pointset
    Int_t           fIdx;  // location in the point set array
    Point_t(Reve::PointSet* ps, Int_t i): fPS(ps), fIdx(i){} 
    bool operator<(const Point_t& o) const
    { if (fPS != o.fPS) return fPS < o.fPS; return fIdx < o.fIdx; }
  };

  Float_t    fAlpha;          // transformation agle to local system (where x>>y)
  AliRieman* fRieman;         // rieman fitter

  Bool_t     fConnected;      // object connected to pointset Ctrl-shift signal 
  
  Reve::TrackList* fTrackList; // track list created with rieman fit 

  std::map<Point_t, Int_t> fMapPS; // map of selected points from different PointSet
public:
  TrackFitter(const Text_t* name, Int_t n_points=0, TreeVarType_e tv_type=TVT_XYZ);
  virtual ~TrackFitter();
  
  void AddFitPoint(Reve::PointSet*,Int_t);  // slot for PointCtrlClicked() signal

  virtual void  DestroyElements(); // *MENU*

  void       Start();
  void       Stop();
  void       FitTrack();
  virtual    void  Reset(Int_t n_points=0, Int_t n_int_ids=0);

  Bool_t     GetConnected(){ return fConnected; }
  AliRieman* GetRieman(){ return fRieman; }

  void       DrawRiemanGraph();
  TGraph*    GetGraphSelected(){ return fGraphSelected; }
  TGraphErrors* GetGraphFitted(){ return fGraphFitted; }

  ClassDef(TrackFitter, 0); // Interface to AliRieman fit.
}; // endclass TrackFitter

}

#endif
