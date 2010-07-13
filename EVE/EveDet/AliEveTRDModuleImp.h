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

class TTree;
class TObjArray;
class TClonesArray;

class TEveTrack;
class TEveGeoTopNode;

class AliTRDgeometry;
class AliTRDdataArrayI;
class AliTRDdigitsManager;

class AliEveTRDHits;
class AliEveTRDDigits;
class AliTRDtrackingChamber;
class AliEveTRDChamber : public TEveElement, public AliEveTRDModule
{
  friend class AliEveTRDDigits;

public:

  AliEveTRDChamber(Int_t det=-1);
  virtual ~AliEveTRDChamber() {}

  Int_t GetNrows()  const {return fNrows;}
  Int_t GetNcols()  const {return fNcols;}
  Int_t GetNtime() const {return fNtime;}

  void  LoadHits(TClonesArray *hits, Int_t &idx);
  void  LoadClusters(TObjArray *cs);
  void  LoadClusters(AliTRDtrackingChamber *tc);
  void  LoadDigits(AliTRDdigitsManager *digits);
  void  LoadTracklets(TTree *trklTree);

  void  SetGeometry(AliTRDgeometry *geo);
  void  SetNtime(Int_t nt) {fNtime = nt;}

protected:
  AliEveTRDDigits  *fDigits;    // digits representation
  AliEveTRDHits    *fHits;      // hits representation
  AliEveTRDHits    *fRecPoints; // cluster representation
  TClonesArray     *fTracklets; // mcm tracklets
  AliTRDgeometry   *fGeo;      // TRD geometry
  TEveGeoTopNode   *fShape;    // rendarable geometry of the chamber 
  // data representation section
  Int_t           fNrows;   // number of rows for this pad plane
  Int_t           fNcols;   // number of columns for this pad plane
  Int_t           fNtime;  // number of timebins

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
  void EnableListElements(Bool_t rnr_self = kTRUE, Bool_t rnr_children = kTRUE);  // *MENU*
  void DisableListElements(Bool_t rnr_self = kTRUE, Bool_t rnr_children = kTRUE); // *MENU*
  void UpdateLeaves();
  void UpdateNode();

  List_i begin() { return fChildren.begin(); }
  List_i end()   { return fChildren.end(); }

  ClassDef(AliEveTRDNode, 0);
};

#endif
