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
private:
  AliEveTRDLoaderSim(const AliEveTRDLoaderSim&);            // Not implemented
  AliEveTRDLoaderSim& operator=(const AliEveTRDLoaderSim&); // Not implemented
public:
  AliEveTRDLoaderSim(const Text_t* n="AliEveTRDLoaderSim", const Text_t* t=0);
  virtual ~AliEveTRDLoaderSim() {}

  Bool_t	GoToEvent(int ev);
  Bool_t	LoadHits(TTree *tH);
  Bool_t	Open(const char *file, const char *dir=".");

private:
  AliRunLoader *fRunLoader; // Run Loader

  ClassDef(AliEveTRDLoaderSim, 1); // Alieve loader for the TRD detector (gAlice)
};


class AliEveTRDLoaderRaw : public AliEveTRDLoader
{
private:
  AliEveTRDLoaderRaw(const AliEveTRDLoaderRaw&);            // Not implemented
  AliEveTRDLoaderRaw& operator=(const AliEveTRDLoaderRaw&); // Not implemented

public:
  AliEveTRDLoaderRaw(const Text_t* n="AliEveTRDLoaderRaw", const Text_t* t=0);
  virtual ~AliEveTRDLoaderRaw() {}

  Bool_t	GoToEvent(int ev);
  Bool_t 	LoadEvent();
  Bool_t	Open(const char *file, const char *dir=".");
  void		SetDataType(TRDDataTypes type);

private:
  void NextEvent(Bool_t rewindOnEnd=kTRUE);

private:
  AliRawReaderDate	*fRawDateReader; // raw data reader
  AliRawReaderRoot	*fRawRootReader; // raw root reader
  AliTRDrawData		*fRaw;           // raw data
  Bool_t		 fDataRoot;      // data in root format
  Int_t			 fEventOld;      // old event

  ClassDef(AliEveTRDLoaderRaw, 1); // Alieve loader for the TRD detector (raw)
};


class AliEveTRDLoaderSimEditor : public TGedFrame
{
private:
  AliEveTRDLoaderSimEditor(const AliEveTRDLoaderSimEditor&);            // Not implemented
  AliEveTRDLoaderSimEditor& operator=(const AliEveTRDLoaderSimEditor&); // Not implemented

public:
  AliEveTRDLoaderSimEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			   UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDLoaderSimEditor() {}

  virtual void	SetModel(TObject* obj);
  virtual void	Toggle(Int_t id);

protected:
  AliEveTRDLoaderSim  *fM; // Model object.
  TGCheckButton       *fLoadHits, *fLoadDigits, *fLoadClusters, *fLoadTracks; // What data to load.

  ClassDef(AliEveTRDLoaderSimEditor,1); // Editor for AliEveTRDLoaderSim
};

#endif
