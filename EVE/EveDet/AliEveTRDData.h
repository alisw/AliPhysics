// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
/////////////////////////////////////////////////////////////////////////
//
// - AliEVE implementation -
// Containers for visualisation of TRD data structures
//    - AliEveTRDHits - visualisation of MC Hits, Clusters (RecPoints)
//    - AliEveTRDDigits - visualisation of TRD digits
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
///////////////////////////////////////////////////////////////////////

#ifndef AliEveTRDData_H
#define AliEveTRDData_H

#include <TEveQuadSet.h>
#include <TEveBoxSet.h>
#include <TEvePointSet.h>

#include <TGedFrame.h>

#include "AliTRDdataArrayI.h"

class AliTRDdigitsManager;
class AliEveTRDChamber;

class AliEveTRDHits : public TEvePointSet
{
public:
  AliEveTRDHits(AliEveTRDChamber *p);

  void PointSelected(Int_t n);

protected:
  AliEveTRDChamber *fParent; // Chaber holding the hits.

private:
  AliEveTRDHits(const AliEveTRDHits&);            // Not implemented
  AliEveTRDHits& operator=(const AliEveTRDHits&); // Not implemented

  ClassDef(AliEveTRDHits, 0); // Base class for TRD hits visualisation
};


class AliEveTRDHitsEditor : public TGedFrame
{
public:
  AliEveTRDHitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		      UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDHitsEditor() {}

  virtual void SetModel(TObject* obj);

protected:
  AliEveTRDHits* fM; // Model object.

private:
  AliEveTRDHitsEditor(const AliEveTRDHitsEditor&);            // Not implemented
  AliEveTRDHitsEditor& operator=(const AliEveTRDHitsEditor&); // Not implemented

  ClassDef(AliEveTRDHitsEditor, 0); // Editor for AliEveTRDHits.
};


class AliEveTRDDigits : public TEveQuadSet
{
  friend class AliEveTRDDigitsEditor;

public:
  AliEveTRDDigits(AliEveTRDChamber *p);

  void			ComputeRepresentation();
  void			Paint(Option_t *opt="");
  void			Reset();
  void			SetData(AliTRDdigitsManager *digits);

protected:
  AliEveTRDChamber *fParent;

private:
  TEveBoxSet		fBoxes; // Boxset for didigit representation.
  AliTRDdataArrayI	fData;  // Raw-data array.

  AliEveTRDDigits(const AliEveTRDDigits&);            // Not implemented
  AliEveTRDDigits& operator=(const AliEveTRDDigits&); // Not implemented

  ClassDef(AliEveTRDDigits, 0); // Digits visualisation for TRD
};


class AliEveTRDDigitsEditor : public TGedFrame
{
public:
  AliEveTRDDigitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTRDDigitsEditor() {}

  virtual void SetModel(TObject* obj);

protected:
  AliEveTRDDigits* fM; // Model object.

private:
  AliEveTRDDigitsEditor(const AliEveTRDDigitsEditor&);            // Not implemented
  AliEveTRDDigitsEditor& operator=(const AliEveTRDDigitsEditor&); // Not implemented

  ClassDef(AliEveTRDDigitsEditor, 0); // Editor for AliEveTRDDigits
};


class AliEveTRDClusters : public AliEveTRDHits
{
public:
  AliEveTRDClusters(AliEveTRDChamber *p);

  void PointSelected(Int_t n);

private:
  AliEveTRDClusters(const AliEveTRDClusters&);            // Not implemented
  AliEveTRDClusters& operator=(const AliEveTRDClusters&); // Not implemented

  ClassDef(AliEveTRDClusters, 0); // Base class for TRD clusters visualisation
};

#endif
