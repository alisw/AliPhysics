// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveTRDLoaderImp_H
#define AliEveTRDLoaderImp_H

////////////////////////////////////////////////////////////////////////
//
// Single event loader for the TRD detector
//    - AliEveTRDLoaderSim - loader for simulations based on gAlice
//    - AliEveTRDLoaderRaw - loader for raw data
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#include "AliEveTRDLoader.h"

class AliRunLoader;
class AliTRDrawData;
class AliRawReaderDate;
class AliRawReaderRoot;

class TGCheckButton;

class AliEveTRDLoaderSim : public AliEveTRDLoader
{
  friend class AliEveTRDLoaderSimEditor;

public:
  AliEveTRDLoaderSim(const Text_t* n="AliEveTRDLoaderSim", const Text_t* t=0);
  virtual ~AliEveTRDLoaderSim() {}

  Bool_t	GoToEvent(int ev);
  Bool_t	Open(const char *file, const char *dir=".");

protected:
  AliRunLoader *fRunLoader; // Run Loader

private:
  AliEveTRDLoaderSim(const AliEveTRDLoaderSim&);            // Not implemented
  AliEveTRDLoaderSim& operator=(const AliEveTRDLoaderSim&); // Not implemented

  ClassDef(AliEveTRDLoaderSim, 0); // Alieve loader for the TRD detector (gAlice)
};


class AliEveTRDLoaderRaw : public AliEveTRDLoader
{
public:
  AliEveTRDLoaderRaw(const Text_t* n="AliEveTRDLoaderRaw", const Text_t* t=0);
  virtual ~AliEveTRDLoaderRaw() {}

  Bool_t	GoToEvent(int ev);
  //Bool_t  NextEvent(Bool_t rewindOnEnd=kTRUE);
  Bool_t	Open(const char *file, const char *dir=".");

private:
  Bool_t 	LoadEvent();

  AliRawReaderDate	*fRawDateReader; // raw data reader
  AliRawReaderRoot	*fRawRootReader; // raw root reader
  AliTRDrawData		  *fRaw;           // raw data
  Int_t             fEventCnt;       // event contor

  AliEveTRDLoaderRaw(const AliEveTRDLoaderRaw&);            // Not implemented
  AliEveTRDLoaderRaw& operator=(const AliEveTRDLoaderRaw&); // Not implemented

  ClassDef(AliEveTRDLoaderRaw, 0); // Alieve loader for the TRD detector (raw)
};


class AliEveTRDLoaderSimEditor : public TGedFrame
{
public:
  AliEveTRDLoaderSimEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			   UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDLoaderSimEditor() {}

  virtual void	SetModel(TObject* obj);
  virtual void	Toggle(Int_t id);

protected:
  AliEveTRDLoaderSim  *fM; // Model object.
  TGCheckButton       *fCheckedHits, *fCheckedDigits, *fCheckedClusters, *fCheckedTracklets; // What data to load.

private:
  AliEveTRDLoaderSimEditor(const AliEveTRDLoaderSimEditor&);            // Not implemented
  AliEveTRDLoaderSimEditor& operator=(const AliEveTRDLoaderSimEditor&); // Not implemented

  ClassDef(AliEveTRDLoaderSimEditor, 0); // Editor for AliEveTRDLoaderSim
};

#endif
