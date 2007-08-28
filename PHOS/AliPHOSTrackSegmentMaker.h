#ifndef ALIPHOSTRACKSEGMENTMAKER_H
#define ALIPHOSTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.42  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.41  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.40  2006/08/29 11:41:19  kharlov
 * Missing implementation of ctors and = operator are added
 *
 * Revision 1.39  2006/08/25 16:00:53  kharlov
 * Compliance with Effective C++AliPHOSHit.cxx
 *
 * Revision 1.38  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Algorithm Base class to construct PHOS track segments
// Associates EMC and CPV clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Kurchatov Institute  & SUBATECH)

// --- ROOT system ---
#include <TObject.h>
class TTree;

// --- AliRoot header files ---
class AliPHOSGeometry ;
class AliESDEvent ;
class AliPHOSQualAssDataMaker ; 

class  AliPHOSTrackSegmentMaker : public TObject {

public:

  AliPHOSTrackSegmentMaker();
  AliPHOSTrackSegmentMaker(AliPHOSGeometry *geom);
  AliPHOSTrackSegmentMaker(const AliPHOSTrackSegmentMaker & tsmaker) ;
  virtual ~ AliPHOSTrackSegmentMaker() ;
  AliPHOSTrackSegmentMaker & operator = (const AliPHOSTrackSegmentMaker & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

  virtual void   Clusters2TrackSegments(Option_t *option) = 0;

  void    SetInput(TTree *clustersTree);

  virtual void    Print(const Option_t * = "")const {Warning("Print", "Not Defined" ) ; }

  void SetESD(AliESDEvent *esd) { fESD = esd; }

  AliESDEvent *GetESD()             const {return fESD;            }

  AliPHOSQualAssDataMaker * GetQualAssDataMaker() const { return fQADM ; } 

  virtual TClonesArray * GetTrackSegments() const = 0;

protected:

  AliESDEvent * fESD;              //! ESD object
  AliPHOSQualAssDataMaker * fQADM ; //!Quality Assurance Data Maker
  AliPHOSGeometry *fGeom;           //! Pointer to the PHOS geometry
  TObjArray *fEMCRecPoints;         //  Array with the EMC clusters
  TObjArray *fCPVRecPoints;         //  Array with the CPV clusters

  ClassDef( AliPHOSTrackSegmentMaker,6)  // Algorithm class to make PHOS track segments (Base Class)
};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
