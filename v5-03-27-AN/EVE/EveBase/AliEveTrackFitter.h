// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTrackFitter_H
#define AliEveTrackFitter_H

#include <TEvePointSet.h>
#include <map>

class TGraphErrors;
class TGraph;
class AliRieman;

class TEveTrackList;

class AliEveTrackFitter : public TEvePointSet
{
public:
  AliEveTrackFitter(const Text_t* name = "TrackFitter", Int_t nPoints=0);
  virtual ~AliEveTrackFitter();

  virtual void DestroyElements();

  virtual void  AddFitPoint(Int_t pointId);  // slot for TEvePointSet::PointSelected() signal

  virtual void  Start();
  virtual void  Stop();
  virtual void  FitTrack();
  virtual void  Reset(Int_t nPoints=0, Int_t nIntIds=0);

  Bool_t        GetConnected() const { return fConnected; }
  AliRieman*    GetRieman()    const { return fRieman; }

  TGraph*       GetGraphPicked() const { return fGraphPicked; }
  TGraphErrors* GetGraphHelix()  const { return fGraphHelix; }
  void          DrawDebugGraph();


protected:

  struct Point_t
  {
    TEvePointSet   *fPS;   // point set
    Int_t           fIdx;  // id in point set

    Point_t(TEvePointSet* ps=0, Int_t i=0) : fPS(ps), fIdx(i) {}
    Point_t(const Point_t& p) : fPS(p.fPS), fIdx(p.fIdx)  {}

    Point_t& operator=(const Point_t& p)
    {
      fPS = p.fPS; fIdx = p.fIdx; return *this;
    }

    bool operator<(const Point_t& o) const
    {
      if (fPS != o.fPS) return fPS < o.fPS;
      return fIdx < o.fIdx;
    }
  };

  typedef std::map<Point_t, Int_t> PointMap_t;

  Float_t            fAlpha;           // transformation angle to AliRieman local system (where x>>y)
  AliRieman*         fRieman;          // rieman fitter

  Bool_t             fConnected;       // connection to the TEvePointSet signal

  PointMap_t         fSPMap;           // map of selected points
  TEveTrackList*     fTrackList;       // list of tracks removed in the destructor

  TGraph            *fGraphPicked;     // graph of selected points debug info
  TGraphErrors      *fGraphHelix;      // graph of fitted points for debug info

private:
  AliEveTrackFitter(const AliEveTrackFitter&);            // Not implemented
  AliEveTrackFitter& operator=(const AliEveTrackFitter&); // Not implemented

  ClassDef(AliEveTrackFitter, 0); // Interface of TEvePointSet allowing helix fit.
};

#endif
