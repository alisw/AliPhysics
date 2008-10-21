#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.41  2007/08/28 12:55:08  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.40  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.39  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.38  2007/04/01 15:40:15  kharlov
 * Correction for actual vertex position implemented
 *
 * Revision 1.37  2006/08/29 11:41:19  kharlov
 * Missing implementation of ctors and = operator are added
 *
 * Revision 1.36  2006/08/25 16:00:53  kharlov
 * Compliance with Effective C++AliPHOSHit.cxx
 *
 * Revision 1.35  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class                             
//  of identified particles                
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TObject.h"
class TTree;

// --- Standard library ---

// --- AliRoot header files ---
class AliESDEvent ;
class AliPHOSGeometry ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;

class AliPHOSPID : public TObject {

 public:

  AliPHOSPID() ;          // ctor            
  AliPHOSPID (AliPHOSGeometry *geom);
  AliPHOSPID(const AliPHOSPID & pid) ;
  virtual ~AliPHOSPID() ; // dtor
  AliPHOSPID & operator = (const AliPHOSPID & /*rvalue*/)  {
    Fatal("operator =", "not implemented") ; return *this ; }

  virtual void TrackSegments2RecParticles(Option_t * option) = 0;

  void SetInput(TTree *clustersTree, TClonesArray *trackSegments);
  TClonesArray* GetRecParticles() const { return fRecParticles; }

  virtual void Print(const Option_t * = "") const = 0;

  void SetESD(AliESDEvent *esd) { fESD = esd; }

  void SetEnergyCorrectionOn(Bool_t on=kTRUE) {fEnergyCorrectionOn = on;}
  Bool_t GetEnergyCorrectionOn() const  {return fEnergyCorrectionOn;}

  virtual const char * Version() const = 0;

protected:

  AliPHOSGeometry * fGeom;    //! Pointer to PHOS Geometry
  AliESDEvent * fESD;         //! ESD object

  TObjArray *fEMCRecPoints;      //!Array with EMC clusters
  TObjArray *fCPVRecPoints;      //!Array with CPV clusters

  TClonesArray *fTrackSegments;     //!Array with found track segments
  TClonesArray *fRecParticles;      //!Array with reconstructed particles (PID)
  
  Bool_t   fEnergyCorrectionOn;     // Do energy correction in GetCalibratedEnergy()
  
private: 

  ClassDef(AliPHOSPID,7)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
