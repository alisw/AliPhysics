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

#ifndef ALIEVE_TRDData_H
#define ALIEVE_TRDData_H

#include <TEveQuadSet.h>
#include <TEveBoxSet.h>
#include <TEvePointSet.h>

#include <TGedFrame.h>


#include "AliTRDdataArrayI.h"

class AliTRDdigitsManager;
class AliEveTRDChamber;


class AliEveTRDHits : public TEvePointSet
{
private:
  AliEveTRDHits(const AliEveTRDHits&);            // Not implemented
  AliEveTRDHits& operator=(const AliEveTRDHits&); // Not implemented

public:
  AliEveTRDHits(AliEveTRDChamber *p);

  void PointSelected(Int_t n);

protected:
  AliEveTRDChamber *fParent;

  ClassDef(AliEveTRDHits,1); // Base class for TRD hits visualisation
};


class AliEveTRDHitsEditor : public TGedFrame
{
private:
  AliEveTRDHitsEditor(const AliEveTRDHitsEditor&);            // Not implemented
  AliEveTRDHitsEditor& operator=(const AliEveTRDHitsEditor&); // Not implemented

public:
  AliEveTRDHitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		      UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTRDHitsEditor();

  virtual void SetModel(TObject* obj);

protected:
  AliEveTRDHits* fM;

  ClassDef(AliEveTRDHitsEditor,1); // Editor for AliEveTRDHits
};


class AliEveTRDDigits : public TEveQuadSet
{
  friend class AliEveTRDDigitsEditor;

private:
  AliEveTRDDigits(const AliEveTRDDigits&);            // Not implemented
  AliEveTRDDigits& operator=(const AliEveTRDDigits&); // Not implemented

public:
  AliEveTRDDigits(AliEveTRDChamber *p);

  void			ComputeRepresentation();
  void			Paint(Option_t *opt="");
  void			Reset();
  void			SetData(AliTRDdigitsManager *digits);

protected:
  AliEveTRDChamber *fParent;

private:
  TEveBoxSet		fBoxes;
  AliTRDdataArrayI	fData;

  ClassDef(AliEveTRDDigits,1); // Digits visualisation for TRD
};


class AliEveTRDDigitsEditor : public TGedFrame
{
private:
  AliEveTRDDigitsEditor(const AliEveTRDDigitsEditor&);            // Not implemented
  AliEveTRDDigitsEditor& operator=(const AliEveTRDDigitsEditor&); // Not implemented

public:
  AliEveTRDDigitsEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~AliEveTRDDigitsEditor();

  virtual void SetModel(TObject* obj);

protected:
  AliEveTRDDigits* fM;

  ClassDef(AliEveTRDDigitsEditor,1); // Editor for AliEveTRDDigits
};


class AliEveTRDClusters : public AliEveTRDHits
{
private:
  AliEveTRDClusters(const AliEveTRDClusters&);            // Not implemented
  AliEveTRDClusters& operator=(const AliEveTRDClusters&); // Not implemented

public:
  AliEveTRDClusters(AliEveTRDChamber *p);

  void PointSelected(Int_t n);

  ClassDef(AliEveTRDClusters,1); // Base class for TRD clusters visualisation
};

#endif
