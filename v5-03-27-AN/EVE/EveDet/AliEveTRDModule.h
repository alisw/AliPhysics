// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTRDModule_H
#define AliEveTRDModule_H

/////////////////////////////////////////////////////////////////////////
//
// - AliEVE implementation -
// The common structure of a TRD module (SM, Stack or Chamber)
//    - AliEveTRDModule - structure of TRD module for visualisation
//    - AliEveTRDModuleEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
///////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TGedFrame.h>

class TObject;
class TGWindow;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;


class AliEveTRDModule : public TNamed
{
  friend class AliEveTRDModuleEditor;
  friend class AliEveTRDNode;
  friend class AliEveTRDChamber;

public:
  AliEveTRDModule(const char *typ="XXX", Int_t id=0);
  virtual ~AliEveTRDModule() {}

  virtual Bool_t GetDigitsBox() const {return fDigitsBox;}
  virtual Bool_t GetDigitsLog() const {return fDigitsLog;}
  virtual UShort_t GetDigitsThreshold() const {return fDigitsThreshold;}
  virtual Int_t	GetID() const {return fDet;}

protected:
  // UI section
  Bool_t	fLoadHits, fRnrHits, fLoadDigits;   // What to load.
  Bool_t	fRnrDigits, fDigitsLog, fDigitsBox; // What to show.
  Bool_t	fDigitsNeedRecompute;               // Need to recompute digits.

  Bool_t	fLoadRecPoints, fRnrRecPoints; // What to do with recpoints.
  Bool_t	fLoadTracklets, fRnrTracklets; // What to do with tracklets.

  Int_t         fDet;             // detector number
  UShort_t	fDigitsThreshold; // digits threshold

private:
  AliEveTRDModule(const AliEveTRDModule&);            // Not implemented
  AliEveTRDModule& operator=(const AliEveTRDModule&); // Not implemented

  ClassDef(AliEveTRDModule, 0); // Structure holder for TRD chamber.
};


class AliEveTRDModuleEditor : public TGedFrame
{
public:
  AliEveTRDModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDModuleEditor() {}

  virtual void SetModel(TObject* obj);

  void	ModifyDigitsView();
  void	SetThreshold(Long_t thres);
  void	UpdateChamber();
  void	UpdateClusters(Pixel_t);
  void	UpdateHits(Pixel_t);

protected:
  AliEveTRDModule* fM; // Model object.

private:
  TGCheckButton *fDisplayHits; // Hit control.
  TGColorSelect *fHitsColor;   // Hit color.
  TGCheckButton *fDisplayDigits, *fToggleLog, *fToggleBox, *fThreshold; // Display toggles.
  TGNumberEntry	*fThresValue;       // Threshold weed.
  TGCheckButton *fDisplayClusters;  // Cluster control.
  TGColorSelect *fClustersColor;    // Cluster color.
  TGCheckButton *fDisplayTracks;    // Track control.

  AliEveTRDModuleEditor(const AliEveTRDModuleEditor&);            // Not implemented
  AliEveTRDModuleEditor& operator=(const AliEveTRDModuleEditor&); // Not implemented

  ClassDef(AliEveTRDModuleEditor, 0); // Editor for AliEveTRDModule.
};

#endif
