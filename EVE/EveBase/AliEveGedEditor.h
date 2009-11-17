// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveGedEditor_H
#define AliEveGedEditor_H

#include <TGedFrame.h>
#include <TEveGedEditor.h>


//==============================================================================
// AliEveGedFrame
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveGedFrame
//

class AliEveGedFrame  : public TGedFrame
{
public:
  AliEveGedFrame(const TGWindow *p=0);
  virtual ~AliEveGedFrame() {}

  virtual void SetModel(TObject* obj);

protected:
  TGTextButton *fB; // Info button.

private:
  AliEveGedFrame(const AliEveGedFrame&);            // Not implemented
  AliEveGedFrame& operator=(const AliEveGedFrame&); // Not implemented

  ClassDef(AliEveGedFrame, 0); // Short description.
};


//==============================================================================
// AliEveGedEditor
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveGedEditor
//

class AliEveGedEditor : public TEveGedEditor
{
public:
  AliEveGedEditor();
  virtual ~AliEveGedEditor() {}

protected:
  virtual TGedFrame* CreateNameFrame(const TGWindow* parent, const char* tab_name);

private:
  AliEveGedEditor(const AliEveGedEditor&);            // Not implemented
  AliEveGedEditor& operator=(const AliEveGedEditor&); // Not implemented

  ClassDef(AliEveGedEditor, 0); // Short description.
};

#endif
