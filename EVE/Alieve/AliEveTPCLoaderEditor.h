// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_TPCLoaderEditor_H
#define ALIEVE_TPCLoaderEditor_H

#include <TGedFrame.h>

class TGTextButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;

class TEveGValuator;


class AliEveTPCLoader;

class AliEveTPCLoaderEditor : public TGedFrame
{
  AliEveTPCLoaderEditor(const AliEveTPCLoaderEditor&);            // Not implemented
  AliEveTPCLoaderEditor& operator=(const AliEveTPCLoaderEditor&); // Not implemented

protected:
  AliEveTPCLoader* fM; // fModel dynamic-casted to AliEveTPCLoaderEditor

  TGTextEntry*  fFile;
  TGTextButton* fOpenFile;

  TEveGValuator* fEvent;
  TGCheckButton*    fDoubleSR;

  // AliEveTPCData loading settings
  TEveGValuator* fDataLoadThreshold;
  TEveGValuator* fDataLoadPedestal;
  TGCheckButton*    fDataAutoPedestal;

  TGTextButton* fUpdateSectors;
  TGTextButton* fReloadSectors;
  TGTextButton* fCreateSectors3D;
  TGTextButton* fDeleteSectors3D;

public:
  AliEveTPCLoaderEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTPCLoaderEditor();

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

  ClassDef(AliEveTPCLoaderEditor, 0); // Editor for AliEveTPCLoader
}; // endclass AliEveTPCLoaderEditor

#endif
