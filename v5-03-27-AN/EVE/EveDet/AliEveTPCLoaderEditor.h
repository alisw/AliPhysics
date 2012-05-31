// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCLoaderEditor_H
#define AliEveTPCLoaderEditor_H

#include <TGedFrame.h>

class TGTextButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;

class TEveGValuator;
class TEveGDoubleValuator;

class AliEveTPCLoader;

//------------------------------------------------------------------------------
// AliEveTPCLoaderEditor
//
// GUI editor for AliEveTPCLoader.
//

class AliEveTPCLoaderEditor : public TGedFrame
{
  AliEveTPCLoaderEditor(const AliEveTPCLoaderEditor&);            // Not implemented
  AliEveTPCLoaderEditor& operator=(const AliEveTPCLoaderEditor&); // Not implemented

public:
  AliEveTPCLoaderEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                        UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTPCLoaderEditor() {}

  virtual void SetModel(TObject* obj);

  void FileSelect();
  void FileChanged();
  void DoOpen();

  void DoEvent();
  void DoDoubleSR();

  void DoDataLoadThreshold();
  void DoDataLoadPedestal();
  void DoDataAutoPedestal();

  void DoUpdateSectors();
  void DoReloadSectors();
  void DoCreateSectors3D();
  void DoDeleteSectors3D();
  void DoShowSectors2D();
  void DoHideSectors2D();

protected:
  AliEveTPCLoader *fM;                  // Model object.

  TGTextEntry     *fFile;               // Text entry for file-name.
  TGTextButton    *fOpenFile;           // Button to open the file.

  TEveGValuator   *fEvent;              // Valueator for event number.
  TGCheckButton   *fDoubleSR;           // Check-box for double sampling-rate.

  // AliEveTPCData loading settings
  TEveGValuator   *fDataLoadThreshold;  // Valuator for threshold.
  TEveGValuator   *fDataLoadPedestal;   // Valuator for pedestal.
  TGCheckButton   *fDataAutoPedestal;   // Check-box for auto pedestal.

  TGTextButton    *fUpdateSectors;      // Button to update sectors.
  TGTextButton    *fReloadSectors;      // Button to reload sectors.
  TGTextButton    *fCreateSectors3D;    // Button to create 3D sectors.
  TGTextButton    *fDeleteSectors3D;    // Button to delete 3D sectors.

  TEveGDoubleValuator *gEtaRange;       // Slider to set eta cuts
  TGCheckButton *gCutOnEta;             // Checkbutton to apply eta cuts

  ClassDef(AliEveTPCLoaderEditor, 0); // Editor for AliEveTPCLoader.
};

#endif
