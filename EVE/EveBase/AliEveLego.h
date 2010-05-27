// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveLego_H
#define AliEveLego_H

#include "TEveElement.h"

class AliESDEvent;
class AliEveEventSelector;
class AliEveMultiView;
class AliPhysicsSelection;

class TEveCaloDataHist;
class TEveCaloLego;
class TEveCaloLegoOverlay;
class TEveCalo3D;
class TEveScene;
class TEveViewer;
class TEveWindowSlot;
class TGLOverlayButton;
class TGLViewer;
class TH2F;

//______________________________________________________________________________
// Creates a 2D lego histogram on a new tab.
// Allow 3D visualization and projection of the lego histogram on the 3D Viewer and the MultiView.
//

class AliEveLego : public TEveElementList
{
public:
  AliEveLego(const char* name="AliEveLego");
  virtual ~AliEveLego();

  // Get Methods
  TEveCaloDataHist* GetData()           {   return fData;             }
  TEveCaloDataHist* GetDataAllEvents()  {   return fData_all_events;  }
  AliESDEvent*      GetESD()            {   return fEsd;              }
  TEveCaloLego*     GetLego()           {   return fLego;             }
  TEveCaloLego*     GetLegoAllEvents()  {   return fLego_all_events;  }
  TEveCalo3D*       GetCalo3D()         {   return fCalo3d;           }
  TEveCalo3D*       GetCalo3DAllEvents(){   return fCalo3d_all_events;}
  AliEveMultiView*  GetMultiView()      {   return fAl;               }
  Float_t           GetPtMax();
  Float_t           GetPtMaxAE();

  // Set Methods
  void SetCharge(Int_t id) {  fChargeId = id;  Update();  }
  void SetAllEventsCharge(Int_t id) {  fChargeIdAE = id; FilterAllData();  }
  void SetTracks(Int_t id) {  fTracksId = id;  Update();  }
  void SetTracksAE(Int_t id) {  fTracksIdAE = id;  FilterAllData();  }
  void SetEvents(Int_t id) {  fEventsId = id;  Update();  }
  void SetMaxPt(Double_t val);
  void SetMaxPtAE(Double_t val);
  void SetThreshold(Double_t val);
  void SetThresholdAE(Double_t val);
  void SetEventSelection();
  void ShowPrevEvent();
  void ShowNextEvent();

  // Functions
  void Update();
  TEveCaloDataHist* LoadData();
  TEveCaloDataHist* LoadAllData();
  TEveCaloDataHist* FilterData();
  TEveCaloDataHist* FilterAllData();
  TEveCaloLego*     CreateHistoLego();
  TEveCaloLego*     CreateHistoLego(TEveWindowSlot* slot);
  TEveCalo3D*       Create3DView();
  TEveCalo3D*       Create3DView(TEveWindowSlot* slot);
  void              CreateProjections(TEveWindowSlot* slot1,
                                      TEveWindowSlot* slot2);
  TEveCaloDataHist* LoadAllEvents();
  void              ShowEventSeletion(Bool_t show, Bool_t updateonly = kFALSE);
  void              SelectEventSelection(Int_t id);

private:
  Int_t fChargeId;
  Int_t fTracksId;
  Int_t fEventsId;
  Float_t fMaxPt;

  Int_t fChargeIdAE;
  Int_t fTracksIdAE;
  Float_t fMaxPtAE;

  AliESDEvent *fEsd;
  AliPhysicsSelection *fPhysicsSelection;
  TH2F *fHistopos;
  TH2F *fHistoposclone;
  TH2F *fHistopos_all_events;
  TH2F *fHistoneg;
  TH2F *fHistonegclone;
  TH2F *fHistoneg_all_events;

  TEveCaloDataHist *fData;
  TEveCaloDataHist *fData_all_events;
  TEveCaloLego *fLego;
  TEveCaloLego *fLego_all_events;
  TEveCalo3D *fCalo3d;
  TEveCalo3D *fCalo3d_all_events;
  TGLViewer  *fGlv;

  TEveViewer *fHisto2d_v;
  TEveScene  *fHisto2d_s;
  TEveScene  *fHisto2d_s2;
  TEveViewer *fHisto2d_all_events_v0;
  TEveViewer *fHisto2d_all_events_v1;
  TEveViewer *fHisto2d_all_events_v2;
  TEveViewer *fHisto2d_all_events_v3;
  TEveScene  *fHisto2d_all_events_s0;
  TEveScene  *fHisto2d_all_events_s1;
  TEveScene  *fHisto2d_all_events_s2;
  TEveScene  *fHisto2d_all_events_s3;

  AliEveMultiView *fAl;
  TEveCaloLegoOverlay *fHisto2d_lego_overlay;
  TEveCaloLegoOverlay* fHisto2d_all_events_lego_overlay;
  TEveWindowSlot* fHisto2d_all_events_slot;

  AliEveEventSelector *fEventSelector;
  Bool_t fShowEventsInfo;
  TGLOverlayButton *fGButton;
  TGLOverlayButton *fB1;
  TGLOverlayButton *fB2;

  AliEveLego(const AliEveLego&);            // Not implemented
  AliEveLego& operator=(const AliEveLego&); // Not implemented


  ClassDef(AliEveLego, 0); // Short description.
};

#endif
