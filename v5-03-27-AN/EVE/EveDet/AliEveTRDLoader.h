// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveTRDLoader_H
#define AliEveTRDLoader_H

////////////////////////////////////////////////////////////////////////
// - ALIEVE implementation -
// Loader for the TRD detector - base class
//    - AliEveTRDLoader - loader of TRD data (simulation + measured)
//    - AliEveTRDLoaderEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#include <TEveElement.h>
#include <TGedFrame.h>
#include <TString.h>

class AliTRDv1;
class AliTRDgeometry;

class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;
class TTree;

class TEveGValuator;

class AliEveTRDChamber;
class AliEveTRDLoaderManager;

class AliEveTRDLoader : public TEveElementList
{
  friend class AliEveTRDLoaderEditor;

public:
  enum TRDDataTypes {
    kTRDHits     = BIT(0),
    kTRDDigits   = BIT(1),
    kTRDClusters = BIT(2),
    kTRDTracklets= BIT(3),
    kTRDRawRoot  = BIT(4),
    kTRDRawDate  = BIT(5)
  };

  AliEveTRDLoader(const Text_t* n="AliEveTRDLoader", const Text_t* t=0);
  virtual ~AliEveTRDLoader() {}

  virtual void		AddChambers(int sm=-1, int stk=-1, int ly=-1);
  virtual AliEveTRDChamber*	GetChamber(int d);
  virtual Int_t   GetDataType() const {return fDataType;}
  virtual Int_t   GetEvent() const {return fEvent;}
  virtual Bool_t	GoToEvent(int ev);
          Bool_t  IsDataLinked() const {return TestBit(1);}

  virtual Bool_t	Open(const char *file, const char *dir = ".");
  inline virtual  Bool_t  NextEvent(Bool_t rewind = kTRUE);
  virtual void 		Paint(Option_t *option="");

  virtual void		SetDataType(UChar_t type = 0){fDataType = type;}


protected:
  virtual Bool_t	LoadHits(TTree *tH);
  virtual Bool_t	LoadClusters(TTree *tC);
  virtual Bool_t	LoadDigits(TTree *tD);
  virtual Bool_t	LoadTracklets(TTree *trklTree);
          void    SetDataLinked(Bool_t linked = kTRUE) {SetBit(1, linked);}

  virtual void		Unload();

  UChar_t fDataType;        // data type
  //Bool_t	fLoadHits, fLoadDigits, fLoadClusters, fLoadTracks; // flags for data-loading
  Char_t	fSM, fStack, fLy; // supermodule, stack, layer
  Short_t		fEvent;         // current event to be displayed
  AliTRDgeometry		*fGeo;  // the TRD geometry
  TString	fFilename;        // name of data file
  TString	fDir;             // data directory

private:
  AliEveTRDLoader(const AliEveTRDLoader&);            // Not implemented
  AliEveTRDLoader& operator=(const AliEveTRDLoader&); // Not implemented

  ClassDef(AliEveTRDLoader, 0); // Alieve Loader class for the TRD detector.
};

Bool_t AliEveTRDLoader::NextEvent(Bool_t)
{
  fEvent++;
  return GoToEvent(fEvent);
}

class AliEveTRDLoaderEditor : public TGedFrame
{
public:
  AliEveTRDLoaderEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDLoaderEditor() {}

  virtual void	AddChambers();
  virtual void	FileOpen();
  virtual void	GoTo();
  virtual void	Next();
  virtual void	SetEvent(Double_t ev){fM->fEvent = (Int_t)ev;}
  virtual void	SetModel(TObject* obj);

protected:
  AliEveTRDLoader	*fM;     // Model object.
  TGTextEntry		*fFile;    // File name weed.
  TGTextButton  *fBrowse;  // browse button
  TEveGValuator		*fEvent; // Event no weed.
  TEveGValuator		*fSMNumber, *fStackNumber, *fPlaneNumber; // Detector id weeds.

private:
  AliEveTRDLoaderEditor(const AliEveTRDLoaderEditor&);            // Not implemented
  AliEveTRDLoaderEditor& operator=(const AliEveTRDLoaderEditor&); // Not implemented

  ClassDef(AliEveTRDLoaderEditor, 0); // Editor for AliEveTRDLoader.
};

#endif
