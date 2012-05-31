#ifndef ALIPHOSPIDV0_H
#define ALIPHOSPIDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.12  2007/03/06 06:57:05  kharlov
 * DP:calculation of distance to CPV done in TSM
 *
 * Revision 1.11  2006/09/07 18:31:08  kharlov
 * Effective c++ corrections (T.Pocheptsov)
 *
 * Revision 1.10  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version v0 of the PHOS particle identifier 
// Identification is based on information from PPSD and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
class TFormula ;
class TVector3 ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;

#include "AliPHOSPID.h"

class  AliPHOSPIDv0 : public AliPHOSPID {

public:

  AliPHOSPIDv0() ;          // ctor            
  AliPHOSPIDv0(AliPHOSGeometry *geom);
  virtual ~AliPHOSPIDv0() ; // dtor, empty, memory leak in fFormula member

  //Compiler generated should be ok, because destructor is empty.
  AliPHOSPIDv0(const AliPHOSPIDv0 & rhs);
  AliPHOSPIDv0 & operator = (const AliPHOSPIDv0 & rhs);

  virtual void TrackSegments2RecParticles(Option_t * option);

  virtual void PlotDispersionCuts()const ;
  virtual void Print(const Option_t * = "")const ; 
  virtual void SetIdentificationMethod(const char * option = "CPV DISP" ){fIDOptions = option ;} 
  virtual void SetShowerProfileCut(const char * formula = "0.35*0.35 - (x-1.386)*(x-1.386) - 1.707*1.707*(y-1.008)*(y-1.008)") ;
  virtual void SetDispersionCut(Float_t cut){fDispersion = cut ; } 
  virtual void SetCpvtoEmcDistanceCut(Float_t cut )      {fCpvEmcDistance = cut ;}
  virtual void SetTimeGate(Float_t gate)                 {fTimeGate = gate ;}

  virtual const char * Version() const { return "pid-v0" ; }  
                     
 private:
  
  void     MakeRecParticles(void ) ;
//  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; // Relative Distance CPV-EMC
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;

 private:

  TString                fIDOptions ;         // PID option

  AliPHOSClusterizer   * fClusterizer ;       // !
  AliPHOSTrackSegmentMaker * fTSMaker ;       // !

  TFormula             * fFormula ;           // formula to define cut on the shower elips axis
  Float_t                fDispersion ;        // dispersion cut
  Float_t                fCpvEmcDistance ;    // Max EMC-CPV distance
  Float_t                fTimeGate ;          // Time of the latest EmcRecPoint accepted as EM

  ClassDef( AliPHOSPIDv0,2)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV0_H
