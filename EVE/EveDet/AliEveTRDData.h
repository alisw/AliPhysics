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


#include "AliTRDdataArrayI.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"

class AliEveTRDChamber;
class AliEveTRDHits : public TEvePointSet
{
public:
  AliEveTRDHits(AliEveTRDChamber *p);
  ~AliEveTRDHits();

  void PointSelected(Int_t n);

protected:
  AliEveTRDChamber *fParent; // Chamber holding the hits.

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
  AliEveTRDClusters(AliEveTRDChamber *p);

  void PointSelected(Int_t n);

private:
  AliEveTRDClusters(const AliEveTRDClusters&);            // Not implemented
  AliEveTRDClusters& operator=(const AliEveTRDClusters&); // Not implemented

  ClassDef(AliEveTRDClusters, 0); // Base class for TRD clusters visualisation
};




class AliEveTRDTracklet : public TEveLine, public AliTRDseedV1 
{
public:
  AliEveTRDTracklet();
//  ~AliEveTRDTracklet();
  

private:
  AliEveTRDTracklet(const AliEveTRDTracklet&);            // Not implemented
  AliEveTRDTracklet& operator=(const AliEveTRDTracklet&); // Not implemented

  ClassDef(AliEveTRDTracklet, 0); // TRD tracklet visualisation
};




class AliEveTRDTrack : public TEveLine, public AliTRDtrackV1
{
public:
  AliEveTRDTrack();
//  ~AliEveTRDTrack();
  
private:
  AliEveTRDTrack(const AliEveTRDTrack&);            // Not implemented
  AliEveTRDTrack& operator=(const AliEveTRDTrack&); // Not implemented

  ClassDef(AliEveTRDTrack, 0); // TRD track visualisation
};


#endif
