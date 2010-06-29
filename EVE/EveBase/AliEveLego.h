// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVELEGO_H
#define ALIEVELEGO_H

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
  TEveCaloDataHist* GetDataAllEvents()  {   return fDataAllEvents;  }
  AliESDEvent*      GetESD()            {   return fEsd;              }
  TEveCaloLego*     GetLego()           {   return fLego;             }
  TEveCaloLego*     GetLegoAllEvents()  {   return fLegoAllEvents;  }
  TEveCalo3D*       GetCalo3D()         {   return fCalo3d;           }
  TEveCalo3D*       GetCalo3DAllEvents(){   return fCalo3dAllEvents;}
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
  Int_t fChargeId;  // determine if the charge is +/-
  Int_t fTracksId;  // determine tracks selection
  Int_t fEventsId;  // determine events selection
  Float_t fMaxPt;   // set maximum pT

  Int_t fChargeIdAE;  // determine if the charge is +/- for all events
  Int_t fTracksIdAE;  // determine tracks selection for all events
  Float_t fMaxPtAE;   // determine maximum pT for all events

  AliESDEvent *fEsd;                      // ESD tree
  AliPhysicsSelection *fPhysicsSelection; // physics selection object
  TH2F *fHistopos;                        // positive charge histogram
  TH2F *fHistoposAllEvents;               // positive charge histogram for all events
  TH2F *fHistoneg;                        // negative charge histogram
  TH2F *fHistonegAllEvents;               // negative charge histogram for all events

  TEveCaloDataHist *fData;          // calo data for 2D, 3D histograms
  TEveCaloDataHist *fDataAllEvents; // calo data for all events
  TEveCaloLego *fLego;              // calo lego for histograms
  TEveCaloLego *fLegoAllEvents;     // calo lego for all events histograms
  TEveCalo3D *fCalo3d;              // 3D histogram for single event
  TEveCalo3D *fCalo3dAllEvents;     // 3D histogram for all events
  TGLViewer  *fGlv;                 // viewer object

  TEveViewer *fHisto2dv;           // viewer for histograms
  TEveScene  *fHisto2ds;           // scene for 3d histogram
  TEveScene  *fHisto2ds2;          // scene for 3d histogram new tab
  TEveViewer *fHisto2dAllEventsv0; // viewer 1 for all events tab
  TEveViewer *fHisto2dAllEventsv1; // viewer 2 for all events tab
  TEveViewer *fHisto2dAllEventsv2; // viewer 3 for all events tab
  TEveViewer *fHisto2dAllEventsv3; // viewer 4 for all events tab
  TEveScene  *fHisto2dAllEventss0; // scene for all events tab
  TEveScene  *fHisto2dAllEventss1; // scene for all events tab
  TEveScene  *fHisto2dAllEventss2; // scene for all events tab
  TEveScene  *fHisto2dAllEventss3; // scene for all events tab

  AliEveMultiView *fAl;            // AliEveMultiView object to create 2D projections
  TEveCaloLegoOverlay *fHisto2dLegoOverlay; // Overlay for calo lego
  TEveCaloLegoOverlay* fHisto2dAllEventsLegoOverlay; // Overlay for calo lego all events
  TEveWindowSlot* fHisto2dAllEventsSlot;  // window slot for 2d all events histogram

  AliEveEventSelector *fEventSelector;  // event selector
  Bool_t fShowEventsInfo;               // determine if show events info
  TGLOverlayButton *fGButton;           // button box with collision candidate information
  TGLOverlayButton *fB1;                // button box with beam 1 information
  TGLOverlayButton *fB2;                // button box with beam 2 information

  AliEveLego(const AliEveLego&);            // Not implemented
  AliEveLego& operator=(const AliEveLego&); // Not implemented


  ClassDef(AliEveLego, 0); // Short description.
};

#endif
