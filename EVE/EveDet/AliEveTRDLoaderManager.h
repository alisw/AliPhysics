// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TRDLoaderManager_H
#define ALIEVE_TRDLoaderManager_H

////////////////////////////////////////////////////////////////////////
//                                                                      
// - ALIEVE implementation -
// Loader manager for the TRD detector
//    - AliEveTRDLoaderManager - manager of TRD data loaders (simulation + measured)
//    - AliEveTRDLoaderManagerEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Mon Feb 26 2007
////////////////////////////////////////////////////////////////////////

#include <TEveElement.h>

#include <TNamed.h>
#include <TGedFrame.h>

class TGComboBox;
class TGTextButton;
class TClonesArray;


class AliEveTRDLoaderManager : public TEveElementList
{
  friend class AliEveTRDLoaderManagerEditor;
private:
  AliEveTRDLoaderManager(const AliEveTRDLoaderManager&);            // Not implemented
  AliEveTRDLoaderManager& operator=(const AliEveTRDLoaderManager&); // Not implemented
public:
  AliEveTRDLoaderManager(const Text_t* name="AliEveTRDLoader", const Text_t* title=0x0);
  ~AliEveTRDLoaderManager();
  void 	Paint(Option_t *option);

protected:
  void	Add(Int_t type, const Text_t *name, const Text_t *title=0x0);
  void	Remove(Int_t entry);

  ClassDef(AliEveTRDLoaderManager, 1); // Alieve loaders manager for TRD
};


class AliEveTRDLoaderManagerEditor : public TGedFrame
{
private:
  AliEveTRDLoaderManagerEditor(const AliEveTRDLoaderManagerEditor&);            // Not implemented
  AliEveTRDLoaderManagerEditor& operator=(const AliEveTRDLoaderManagerEditor&); // Not implemented
public:
  AliEveTRDLoaderManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			       UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTRDLoaderManagerEditor();

  virtual void	Add();
  virtual void	Remove(Int_t entry);
  virtual void	SetModel(TObject* obj);

protected:
  AliEveTRDLoaderManager* fM;

private:
  TGComboBox	*fSelector;
  TGTextButton	*fAdd, *fRemoveButton;
  TGGroupFrame 	*fGroupFrame;
  TClonesArray	*fRemove;              

  ClassDef(AliEveTRDLoaderManagerEditor, 1); // Editor for AliEveTRDLoaderManager
};

#endif

