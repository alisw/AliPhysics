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

#ifndef ALITRDARRAYADC_H
#include "AliTRDarrayADC.h"
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
  //friend class AliEveTRDDigitsEditor;

public:
  AliEveTRDDigits(AliEveTRDChamber *p);
  ~AliEveTRDDigits();
  void			SetData(AliTRDdigitsManager *digits);

protected:
  AliEveTRDChamber *fParent;

private:
  AliEveTRDDigits(const AliEveTRDDigits&);            // Not implemented
  AliEveTRDDigits& operator=(const AliEveTRDDigits&); // Not implemented

  ClassDef(AliEveTRDDigits, 0); // Digits visualisation for TRD
};



class AliEveTRDClusters : public AliEveTRDHits
{
public:
  AliEveTRDClusters();

  void Load(const Char_t *what="all") const; // *MENU*
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
  ~AliEveTRDTracklet();
  AliEveTRDClusters* GetClusters() const {return fClusters;}
  inline void        Load(const Char_t *what="all") const; // *MENU*
  void               Print(Option_t *o="") const; // *MENU*

private:
  AliEveTRDClusters *fClusters;  // clusters

  AliEveTRDTracklet(const AliEveTRDTracklet&);            // Not implemented
  AliEveTRDTracklet& operator=(const AliEveTRDTracklet&); // Not implemented

  ClassDef(AliEveTRDTracklet, 0); // TRD tracklet visualisation
};
void AliEveTRDTracklet::Load(const Char_t *what) const
{
  if(fClusters) fClusters->Load(what);
}


class AliTrackPoint;
class AliTRDtrackV1;
class AliRieman;
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
  //AliEveTRDTracklet*  GetTracklet(Int_t plane) const {return plane <6 && plane >= 0 ? fTracklet[plane] : NULL;}
  void    Print(Option_t *opt="a") const; // *MENU*
  void    Load(const Char_t *what="all") const; // *MENU*
  void    SetStatus(UChar_t s=0);         // *MENU*
  void    SetESDstatus(ULong_t stat) {fESDStatus = stat;} 
private:
  AliEveTRDTrack(const AliEveTRDTrack&);            // Not implemented
  AliEveTRDTrack& operator=(const AliEveTRDTrack&); // Not implemented

  UChar_t        fTrackState;   // bit map for the track drawing state
  ULong_t        fESDStatus;    // ESD status bit for this track
  Float_t        fAlpha;        // sector angle for this track  
  AliTrackPoint *fPoints;       // track crossing points
  AliRieman     *fRim;          // rieman fitter
  ClassDef(AliEveTRDTrack, 0);  // TRD track visualisation
};


#include "TEveElement.h"

class AliTRDtrackletMCM;
class AliTRDtrackletWord;
class AliEveTRDTrackletOnline : public TEveLine
{
public:
  AliEveTRDTrackletOnline(AliTRDtrackletMCM *tracklet);
  AliEveTRDTrackletOnline(AliTRDtrackletWord *tracklet);
  ~AliEveTRDTrackletOnline();

  void               ShowMCM(Option_t *opt = "RHT") const; // *MENU*
private:
  Int_t fDetector;
  Int_t fROB;
  Int_t fMCM;

  AliEveTRDTrackletOnline(const AliEveTRDTrackletOnline&);            // Not implemented
  AliEveTRDTrackletOnline& operator=(const AliEveTRDTrackletOnline&); // Not implemented

  ClassDef(AliEveTRDTrackletOnline, 0); // TRD tracklet visualisation
};


class AliTRDmcmSim;
class AliEveTRDmcm : public TEveElement, public TNamed
{
 public:
  AliEveTRDmcm();
  ~AliEveTRDmcm();

  Bool_t Init(Int_t det, Int_t rob, Int_t mcm); // *MENU*
  Bool_t LoadDigits(); // *MENU*
  Bool_t Filter(); // *MENU*
  Bool_t Tracklet(); // *MENU*
  void Draw(Option_t* option = "FHT"); // *MENU*
  Bool_t AssignPointer(const char* ptrname = "mcmtest"); // *MENU*

 protected:
  AliTRDmcmSim *fMCM;

  AliEveTRDmcm(const AliEveTRDmcm&); // not implemented
  AliEveTRDmcm& operator=(const AliEveTRDmcm&); // not implemented

  ClassDef(AliEveTRDmcm, 0);
};

#endif
