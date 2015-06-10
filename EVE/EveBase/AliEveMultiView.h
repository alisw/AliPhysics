// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMultiView_H
#define AliEveMultiView_H

#include <TEveManager.h>

#include <TEveViewer.h>
#include <TGLViewer.h>

#include <TEveScene.h>
#include <TEveGeoShape.h>

#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>

#include <TEveBrowser.h>
#include <TEveWindow.h>

//______________________________________________________________________________
// Short description of AliEveMultiView
//

class AliEveMultiView
{
public:
  AliEveMultiView(Bool_t setMuonView = kFALSE);
  virtual ~AliEveMultiView();

  static AliEveMultiView* Instance();

  void InitGeomGentle(TEveGeoShape* g3d, TEveGeoShape* grphi, TEveGeoShape* grhoz, TEveGeoShape* gmuon);
  void InitGeomGentleTrd(TEveGeoShape* gtrd);
  void InitGeomGentleMuon(TEveGeoShape* gmuon, Bool_t showRPhi, Bool_t showRhoZ, Bool_t showMuon);

  //-------------------------------------------------------------------------

  void SetDepth(Float_t d);

  //-------------------------------------------------------------------------

  void ImportGeomRPhi(TEveElement* el);
  void ImportGeomRhoZ(TEveElement* el);
  void ImportGeomMuon(TEveElement* el);
  void ImportEventRPhi(TEveElement* el);
  void ImportEventRhoZ(TEveElement* el);
  void ImportEventMuon(TEveElement* el);

  void DestroyEventRPhi();
  void DestroyEventRhoZ();
  void DestroyEventMuon();

  void SetCenterRPhi(Double_t x, Double_t y, Double_t z);
  void SetCenterRhoZ(Double_t x, Double_t y, Double_t z);
  void SetCenterMuon(Double_t x, Double_t y, Double_t z);

  void DestroyAllGeometries();

  //-------------------------------------------------------------------------

  TEveViewer* Get3DView()   { return f3DView; }
  TEveViewer* GetRPhiView() { return fRPhiView; }
  TEveViewer* GetRhoZView() { return fRhoZView; }
  TEveViewer* GetMuonView() { return fMuonView; }

    TEveProjectionManager* GetRhoZMgr() { return fRhoZMgr;}
    
  void SetMuonView(Bool_t set) { fIsMuonView = set; }
  Bool_t IsMuonView() { return fIsMuonView; }


protected:
  TEveProjectionManager *fRPhiMgr; // Obvious meaning.
  TEveProjectionManager *fRhoZMgr; // Obvious meaning.
  TEveProjectionManager *fMuonMgr; // Obvious meaning.

  TEveViewer            *f3DView;   // Obvious meaning.
  TEveViewer            *fRPhiView; // Obvious meaning.
  TEveViewer            *fRhoZView; // Obvious meaning.
  TEveViewer            *fMuonView; // Obvious meaning.
    
    TEveWindowPack *fPack;

  TEveScene             *fRPhiGeomScene;  // Obvious meaning.
  TEveScene             *fRhoZGeomScene;  // Obvious meaning.
  TEveScene             *fMuonGeomScene;  // Obvious meaning.
  TEveScene             *fRPhiEventScene; // Obvious meaning.
  TEveScene             *fRhoZEventScene; // Obvious meaning.
  TEveScene             *fMuonEventScene; // Obvious meaning.

  TEveGeoShape          *fGeomGentle;     // Obvious meaning.
  TEveGeoShape          *fGeomGentleRPhi; // Obvious meaning.
  TEveGeoShape          *fGeomGentleRhoZ; // Obvious meaning.
  TEveGeoShape          *fGeomGentleTrd;  // Obvious meaning.
  TEveGeoShape          *fGeomGentleMuon; // Obvious meaning.

  Bool_t		 fIsMuonView;

  static AliEveMultiView* fgInstance;     // Obvious meaning.

private:
  AliEveMultiView(const AliEveMultiView&);            // Not implemented
  AliEveMultiView& operator=(const AliEveMultiView&); // Not implemented

  ClassDef(AliEveMultiView, 0); // Multiple-views.
};

#endif
