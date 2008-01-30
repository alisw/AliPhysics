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

protected:

  struct Point_t
  {
    TEvePointSet   *fPS;   // point set
    Int_t           fIdx;  // id in point set

    Point_t(TEvePointSet* ps=0, Int_t i=0) : fPS(ps), fIdx(i) {}
    Point_t(const Point_t& p) : fPS(p.fPS), fIdx(p.fIdx)  {}

    Point_t& operator=(const Point_t& p) {
      fPS = p.fPS; fIdx = p.fIdx; return *this;
    }

    bool operator<(const Point_t& o) const {
      if (fPS != o.fPS) return fPS < o.fPS;
      return fIdx < o.fIdx;
    }
  };

  typedef std::map<Point_t, Int_t>          PointMap_t;

  Float_t            fAlpha;                // transformation angle to AliRieman local system (where x>>y)
  AliRieman*         fRieman;               // rieman fitter

  Bool_t             fConnected;            // connection to the TEvePointSet signal

  PointMap_t         fSPMap;                // map of selected points
  TEveTrackList*     fTrackList;            // list of tracks removed in the destructor

  TGraph            *fGraphPicked;          // graph of selected points debug info
  TGraphErrors      *fGraphHelix;           // graph of fitted points for debug info

public:
  AliEveTrackFitter(const Text_t* name = "TrackFitter", Int_t n_points=0);
  virtual ~AliEveTrackFitter();

  virtual void  AddFitPoint(TEvePointSet*,Int_t);  // slot for PointCtrlClicked() signal

  virtual void  Start();
  virtual void  Stop();
  virtual void  FitTrack();
  virtual void  Reset(Int_t n_points=0, Int_t n_int_ids=0);

  Bool_t        GetConnected(){ return fConnected; }
  AliRieman*    GetRieman(){ return fRieman; }

  TGraph*       GetGraphPicked()   { return fGraphPicked; }
  TGraphErrors* GetGraphHelix()    { return fGraphHelix; }
  void          DrawDebugGraph();

  virtual void  DestroyElements();   // *MENU*

  ClassDef(AliEveTrackFitter, 0); // Interface of TEvePointSet allowing helix fit.
};

#endif
