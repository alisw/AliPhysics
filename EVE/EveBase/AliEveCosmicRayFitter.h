// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_CosmicRayFitter_H
#define ALIEVE_CosmicRayFitter_H

#include <TEvePointSet.h>
#include <map>

class TGraphErrors;
class TGraph;
class TGraph2DErrors;
class TGraph2D;
class TLinearFitter;
class TVirtualFitter;

class AliEveCosmicRayFitter : public TEvePointSet
{
public:
  AliEveCosmicRayFitter(const Text_t* name = "CosmicRayFitter", Int_t n_points=0);
  virtual ~AliEveCosmicRayFitter();
  
  void AddFitPoint(Int_t); // slot for PointCtrlClicked() signal

  virtual void       Start();
  virtual void       Stop();
  virtual void       FitTrack();
  virtual void  Reset(Int_t n_points=0, Int_t n_int_ids=0);

  Bool_t       GetConnected() { return fConnected; }
  TVirtualFitter* GetLine3DFitter(){ return fLine3DFitter; }

  TGraph*         GetGraphSelected1() { return fGraphPicked1; }
  TGraph*         GetGraphSelected2() { return fGraphPicked2; }
  TGraph2D*       GetGraphSelected3() { return fGraphPicked3; }
  TGraphErrors*   GetGraphFitted1() { return fGraphLinear1; }
  TGraphErrors*   GetGraphFitted2() { return fGraphLinear2; }
  TGraph2DErrors* GetGraphFitted3() { return fGraphLinear3; }
  void          DrawDebugGraph();

  virtual void  DestroyElements(); // *MENU*

protected:
  struct Point_t
  {
    // inner structure to check duplicates
    TEvePointSet   *fPS;   // selected pointset
    Int_t           fIdx;  // location in the point set array

    Point_t(TEvePointSet* ps=0, Int_t i=0): fPS(ps), fIdx(i){} 
    Point_t(const Point_t& p) : fPS(p.fPS), fIdx(p.fIdx)  {}

    Point_t& operator=(const Point_t& p) {
      fPS = p.fPS; fIdx = p.fIdx; return *this;
    }

    bool operator<(const Point_t& o) const {
      if (fPS != o.fPS) return fPS < o.fPS;
      return fIdx < o.fIdx;
    }
  };

  typedef std::map<Point_t, Int_t>          PointMap_t; // Map of registered points.

  TVirtualFitter* fLine3DFitter; // 3D straight line fitter

  Bool_t     fConnected;         // object connected to pointset Ctrl-shift signal 
  
  PointMap_t fSPMap;             // map of selected points from different PointSet

  TGraph            *fGraphPicked1; // graph of selected points debug info
  TGraphErrors      *fGraphLinear1; // graph of fitted points for debug info
  TGraph            *fGraphPicked2; // graph of selected points debug info
  TGraphErrors      *fGraphLinear2; // graph of fitted points for debug info
  TGraph2D          *fGraphPicked3; // graph of selected points debug info
  TGraph2DErrors    *fGraphLinear3; // graph of fitted points for debug info

private:
  AliEveCosmicRayFitter(const AliEveCosmicRayFitter&);            // Not implemented
  AliEveCosmicRayFitter& operator=(const AliEveCosmicRayFitter&); // Not implemented

  ClassDef(AliEveCosmicRayFitter, 0); // Interface to TEvePointSet allowing 3D straight linear fit.
};

#endif
