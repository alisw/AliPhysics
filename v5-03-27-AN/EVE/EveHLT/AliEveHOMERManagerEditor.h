// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveAliEVEHOMERManagerEditor_H
#define AliEveAliEVEHOMERManagerEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGTextButton;
class TGNumberEntry;
class TGColorSelect;
class TGPictureButton;
class TGComboBox;
class AliEveHOMERManager;

class AliEveHOMERManagerEditor : public TGedFrame
{
public:
  AliEveHOMERManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, 
			   UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveHOMERManagerEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void ConnectToHLT();
  void NextEvent();

  void SetTriggerString(int id);

protected:

  AliEveHOMERManager  *fM; // Model object.
  
  TGTextButton     *fButtonConnect; // Button to connect to HOMER.
  TGTextButton     *fButtonWriteToFile; // Button to write block list to file
  TGTextButton     *fButtonNextEvent; // Button to call next Even
  TGTextButton     *fButtonPrintScreens;  // Button to print viewers
  TGComboBox       *fBoxTriggerSelector; // Drop down menu to select trigger criteria.
  

private:
  AliEveHOMERManagerEditor(const AliEveHOMERManagerEditor&);            // Not implemented
  AliEveHOMERManagerEditor& operator=(const AliEveHOMERManagerEditor&); // Not implemented

  ClassDef(AliEveHOMERManagerEditor, 0); // Editor for AliEveHOMERManager
};

#endif
