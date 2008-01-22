// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TrackFitter_H
#define ALIEVE_TrackFitter_H

#include <TEvePointSet.h>
#include <TQObject.h>
#include <map>

class TGraphErrors;
class TGraph;
class AliRieman;

class TEveTrackList;


class AliEveTrackFitter : public TEvePointSet
{
private:
  AliEveTrackFitter(const AliEveTrackFitter&);            // Not implemented
  AliEveTrackFitter& operator=(const AliEveTrackFitter&); // Not implemented

  TGraph       *fGraphSelected;  // graph of selected points
  TGraphErrors *fGraphFitted;    // graph of fitted points

protected:
  struct Point_t
  {
    // inner structure to check duplicates
    TEvePointSet   *fPS;   // selected pointset
    Int_t           fIdx;  // location in the point set array

    Point_t(TEvePointSet* ps, Int_t i) : fPS(ps), fIdx(i) {}
    Point_t(const Point_t& p) : fPS(p.fPS), fIdx(p.fIdx)  {}
    Point_t& operator=(const Point_t& p)
    { fPS = p.fPS; fIdx = p.fIdx; return *this; }

    bool operator<(const Point_t& o) const
    { if (fPS != o.fPS) return fPS < o.fPS; return fIdx < o.fIdx; }
  };

  Float_t    fAlpha;                // transformation agle to local system (where x>>y)
  AliRieman* fRieman;               // rieman fitter

  Bool_t     fConnected;            // object connected to pointset Ctrl-shift signal

  TEveTrackList* fTrackList;        // track list created with rieman fit

  std::map<Point_t, Int_t> fMapPS;  // map of selected points from different TEvePointSet

public:
  AliEveTrackFitter(const Text_t* name, Int_t n_points=0);
  virtual ~AliEveTrackFitter();

  virtual void AddFitPoint(TEvePointSet*,Int_t);  // slot for PointCtrlClicked() signal

  virtual void DestroyElements();   // *MENU*

  virtual void Start();
  virtual void Stop();
  virtual void FitTrack();
  virtual void Reset(Int_t n_points=0, Int_t n_int_ids=0);

  Bool_t        GetConnected(){ return fConnected; }
  AliRieman*    GetRieman(){ return fRieman; }

  void          DrawRiemanGraph();

  TGraph*       GetGraphSelected() { return fGraphSelected; }
  TGraphErrors* GetGraphFitted()   { return fGraphFitted; }

  ClassDef(AliEveTrackFitter, 0); // Interface to AliRieman fit.
}; // endclass AliEveTrackFitter

#endif
