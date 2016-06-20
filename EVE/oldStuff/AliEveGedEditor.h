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
// AliEveGedNameFrame
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveGedNameFrame
//

class AliEveGedNameFrame  : public TGedFrame
{
public:
  AliEveGedNameFrame(const TGWindow *p=0);
  virtual ~AliEveGedNameFrame() {}

  virtual void SetModel(TObject* obj);

protected:
  TGTextButton *fB; // Info button.

private:
  AliEveGedNameFrame(const AliEveGedNameFrame&);            // Not implemented
  AliEveGedNameFrame& operator=(const AliEveGedNameFrame&); // Not implemented

  ClassDef(AliEveGedNameFrame, 0); // Specialization of GED top name-frame for AliEve.
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

  ClassDef(AliEveGedEditor, 0); // // Specialization of GED editor for AliEve.
};

#endif
