#ifndef ALIPHOSRECONSTRUCTOR_H
#define ALIPHOSRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.15  2007/10/01 20:24:08  kharlov
 * Memory leaks fixed
 *
 * Revision 1.14  2007/09/26 14:22:18  cvetan
 * Important changes to the reconstructor classes. Complete elimination of the run-loaders, which are now steered only from AliReconstruction. Removal of the corresponding Reconstruct() and FillESD() methods.
 *
 * Revision 1.13  2007/08/30 10:40:27  cvetan
 * Minor
 *
 * Revision 1.12  2007/08/28 12:55:08  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.11  2007/07/24 17:20:35  policheh
 * Usage of RecoParam objects instead of hardcoded parameters in reconstruction.
 * (See $ALICE_ROOT/PHOS/macros/BeamTest2006/RawReconstruction.C).
 *
 * Revision 1.10  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.9  2006/11/15 16:05:03  kharlov
 * New FillESD() for raw data is added
 *
 * Revision 1.8  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//--
//-- Author: Yves Schutz (SUBATECH) 
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
//_________________________________________________________________________

// --- ROOT system ---

#include "AliReconstructor.h" 
#include "AliPHOSRecoParam.h"
class AliPHOSDigitizer ;
class AliPHOSClusterizer ;
class AliPHOSClusterizerv1 ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSPID ;
class AliPHOSSDigitizer ;
class AliESDEvent ;
class AliRawReader; 
class AliPHOSRecoParam;
class AliPHOSGeometry;
class AliPHOSCalibData ;
class AliPHOSTriggerParameters;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructor : public AliReconstructor {

public:

  AliPHOSReconstructor() ; //ctor            
  virtual ~AliPHOSReconstructor() ; //dtor            

  static void                SetDebug()   { fgDebug = kTRUE ; }
  static void                ResetDebug() { fgDebug = kFALSE ; }
  static Bool_t              Debug() { return fgDebug ; }
  AliTracker *CreateTracker() const;
  using AliReconstructor::FillESD;
  virtual void               FillESD(TTree* digitsTree, TTree* clustersTree, 
				     AliESDEvent* esd) const;
  using AliReconstructor::Reconstruct;
  virtual void               Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual Bool_t             HasDigitConversion() const {return kTRUE;};
  virtual void               ConvertDigits   (AliRawReader* rawReader, TTree* digitsTree) const;
  virtual void               ConvertDigitsEMC(AliRawReader* rawReader, TClonesArray* digits) const;
  virtual void               ConvertDigitsCPV(AliRawReader* rawReader, TClonesArray* digits) const;
  virtual Float_t            Calibrate(Float_t amp, Int_t absId) const ;
  virtual Float_t            CalibrateT(Float_t time, Int_t absId, Bool_t isLG) const ;

  void FillMisalMatrixes(AliESDEvent* esd)const ;
  
  static const AliPHOSRecoParam* GetRecoParam() {
    return dynamic_cast<const AliPHOSRecoParam*>(AliReconstructor::GetRecoParam(4)); }
  static Float_t CorrectNonlinearity(Float_t oldEnergy) ;
  
  void readTRUParameters(AliPHOSTriggerParameters* parameters) const;
  
private:
  AliPHOSReconstructor(const AliPHOSReconstructor & rec); // Not implemented
  AliPHOSReconstructor & operator = (const AliPHOSReconstructor &); // Not implemented
  
  static Bool_t fgDebug ; //! verbosity controller
  AliPHOSGeometry          *fGeom;           // pointer to the PHOS geometry
  AliPHOSClusterizerv1     *fClusterizer;    //! PHOS clusterizer
  AliPHOSTrackSegmentMaker *fTSM;            //! PHOS TrackSegmentMaker
  AliPHOSPID               *fPID;            //! PHOS PID maker
  TClonesArray             *fTmpDigLG;       //! Temporary array of LG digits
  static TClonesArray      *fgDigitsArray;   //! Array of PHOS digits
  static TObjArray         *fgEMCRecPoints;  //! Array of EMC rec.points
  static TObjArray         *fgCPVRecPoints;  //! Array of CPV rec.points
  static AliPHOSCalibData * fgCalibData ;    //! Calibration database if aval.
  static TClonesArray      *fgTriggerDigits; //! Array of PHOS trigger digits

  ClassDef(AliPHOSReconstructor,11)  // PHOS Reconstruction class

}; 

#endif // ALIPHOSRECONSTRUCTOR_H
