// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTRDModuleImp_H
#define AliEveTRDModuleImp_H

/////////////////////////////////////////////////////////////////////////
//
// Implementation of AliEveTRDModule:
//    - AliEveTRDChamber - Data holder
//    - AliEveTRDNode    - Node structure
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
/////////////////////////////////////////////////////////////////////////

#include <vector>

#include <TEveElement.h>

#include "AliEveTRDModule.h"

class AliTRDpadPlane;
class AliTRDgeometry;
class AliTRDhit;
class AliTRDdataArrayI;
class AliTRDdigitsManager;
class TObjArray;

class TEveTrack;

class AliEveTRDHits;
class AliEveTRDDigits;

class AliEveTRDChamber : public TEveElement, public AliEveTRDModule
{
  friend class AliEveTRDDigits;

public:

  AliEveTRDChamber(Int_t det=0);
  virtual ~AliEveTRDChamber() {}

  void  AddHit(AliTRDhit *hit);
  Int_t GetRowMax()  const {return fRowMax;}
  Int_t GetColMax()  const {return fColMax;}
  Int_t GetTimeMax() const {return fTimeMax;}
  Int_t GetSM() const;
  Int_t GetSTK() const;
  Int_t GetPlane() const {return fPla;}

  void  LoadClusters(TObjArray *cs);
  void  LoadDigits(AliTRDdigitsManager *digits);
  void  LoadTracklets(TObjArray *ts);
  void  Paint(Option_t* option="");
  void  Reset();
  void  SetGeometry(AliTRDgeometry *geo);

protected:
  AliEveTRDDigits         *fDigits;    // digits representation
  AliEveTRDHits           *fHits;      // hits representation
  AliEveTRDHits           *fRecPoints; // cluster representation
  std::vector<TEveTrack*> *fTracklets; // mcm tracklets

  // data representation section
  Int_t           fRowMax;   // number of rows for this pad plane
  Int_t           fColMax;   // number of columns for this pad plane
  Int_t           fTimeMax;  // number of timebins
  Float_t         fSamplingFrequency; // sampling frequency
  Float_t         fX0;       // radial distance from vertex to the chamber
  Int_t           fPla;      // detector plane
  AliTRDpadPlane *fPadPlane; // pad plane object
  AliTRDgeometry *fGeo;      // TRD geometry

private:
  AliEveTRDChamber(const AliEveTRDChamber&);            // Not implemented.
  AliEveTRDChamber& operator=(const AliEveTRDChamber&); // Not implemented.

  ClassDef(AliEveTRDChamber, 0); // Holder for TRD chamber data
};


class AliEveTRDNode : public TEveElement, public AliEveTRDModule
{
public:
  AliEveTRDNode(const char *typ, Int_t det=0);

  void Paint(Option_t* option="");
  void Reset();

  void Collapse();            // *MENU*
  void Expand();              // *MENU*
  void EnableListElements();  // *MENU*
  void DisableListElements(); // *MENU*
  void UpdateLeaves();
  void UpdateNode();

  List_i begin() { return fChildren.begin(); }
  List_i end()   { return fChildren.end(); }

  ClassDef(AliEveTRDNode, 0);
};

#endif
