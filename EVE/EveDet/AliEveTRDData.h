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
#include <TEveLine.h>


#ifndef ALITRDDATAARRAYI_H
#include "AliTRDdataArrayI.h"
#endif

class AliEveTRDChamber;
class AliEveTRDHits : public TEvePointSet
{
public:
  AliEveTRDHits();
  ~AliEveTRDHits();

  void PointSelected(Int_t n);


private:
  AliEveTRDHits(const AliEveTRDHits&);            // Not implemented
  AliEveTRDHits& operator=(const AliEveTRDHits&); // Not implemented

  ClassDef(AliEveTRDHits, 0); // Base class for TRD hits visualisation
};


class AliTRDdigitsManager;
class AliEveTRDDigits : public TEveQuadSet
{
  friend class AliEveTRDDigitsEditor;

public:
  AliEveTRDDigits(AliEveTRDChamber *p);
  ~AliEveTRDDigits();

  void			ComputeRepresentation();
  const AliTRDdataArrayI*	GetData() const {return fData.GetNelems() ? &fData : 0x0;}
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



class AliEveTRDClusters : public AliEveTRDHits
{
public:
  AliEveTRDClusters();

  void Load(Char_t *what="all", Bool_t stkwise=kTRUE) const; // *MENU*
  void PointSelected(Int_t n);
  void Print(Option_t *o = "") const; // *MENU*

private:
  AliEveTRDClusters(const AliEveTRDClusters&);            // Not implemented
  AliEveTRDClusters& operator=(const AliEveTRDClusters&); // Not implemented

  ClassDef(AliEveTRDClusters, 0); // Base class for TRD clusters visualisation
};



class AliTRDseedV1;
class AliEveTRDTracklet : public TEveLine
{
public:
  AliEveTRDTracklet(AliTRDseedV1 *trklt);
//  ~AliEveTRDTracklet();
  AliEveTRDClusters* GetClusters() const {return fClusters;}
  void               Print(Option_t *o="") const; // *MENU*
private:
  AliEveTRDClusters *fClusters;  // clusters

  AliEveTRDTracklet(const AliEveTRDTracklet&);            // Not implemented
  AliEveTRDTracklet& operator=(const AliEveTRDTracklet&); // Not implemented

  ClassDef(AliEveTRDTracklet, 0); // TRD tracklet visualisation
};


class AliTrackPoint;
class AliTRDtrackV1;
class AliEveTRDTrack : public TEveLine
{
public:
  enum AliEveTRDTrackState{
    kSource        = 0 
    ,kPID          = 1 
    ,kTrackCosmics = 2
    ,kTrackModel   = 3
  };
  enum AliEveTRDTrackModel{
    kRieman  = 0
    ,kKalman = 1
  };


  AliEveTRDTrack(AliTRDtrackV1 *trk);
  virtual ~AliEveTRDTrack();
  //AliEveTRDTracklet*  GetTracklet(Int_t plane) const {return plane <6 && plane >= 0 ? fTracklet[plane] : 0x0;}
  void    Print(Option_t *opt="a") const; // *MENU*
  void    SetStatus(UChar_t s);
  void    SetESDstatus(ULong_t stat) {fESDStatus = stat;} 
private:
  AliEveTRDTrack(const AliEveTRDTrack&);            // Not implemented
  AliEveTRDTrack& operator=(const AliEveTRDTrack&); // Not implemented

  UChar_t        fTrackState;   // bit map for the track drawing state
  ULong_t        fESDStatus;    // ESD status bit for this track
  Float_t        fAlpha;        // sector angle for this track  
  AliTrackPoint* fPoints;       // track crossing points

  ClassDef(AliEveTRDTrack, 0);  // TRD track visualisation
};


#endif
