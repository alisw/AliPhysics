#ifndef ALIFMDV1_H
#define ALIFMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDv1.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:48:51 2006
    @brief   Concrete implementation of FMD detector driver - detailed
    version 
*/
//____________________________________________________________________
//
//  Manager class for the FMD - Detailed version. 
//  Implements the full geometry, 
//  And does stepping
// 
#ifndef ALIFMD_H 
# include "AliFMD.h"
#endif
#ifndef ROOT_TLorentzVector
# include <TLorentzVector.h>
#endif
 
//____________________________________________________________________
/** @brief Forward Multiplicity Detector based on Silicon wafers. 

    This class contains the base procedures for the Forward
    Multiplicity detector Detector consists of 3 sub-detectors FMD1,
    FMD2, and FMD3, each of which has 1 or 2 rings of silicon sensors.
                                                           
    This class contains the detailed version of the FMD - that is,
    hits are produced during simulation. 
    @ingroup FMD_sim
*/
class AliFMDv1 : public AliFMD 
{
public:
  /** CTOR */
  AliFMDv1()
    : AliFMD(),
      fCurrentDeltaE(0),
      fCurrentV(),
      fCurrentP(),
      fCurrentPdg(0) { fDetailed = kTRUE; }
  /** CTOR 
      @param name Name 
      @param title Title */
  AliFMDv1(const char *name, const char *title="Detailed geometry") 
    : AliFMD(name, title),
      fCurrentDeltaE(0),
      fCurrentV(),
      fCurrentP(),
      fCurrentPdg(0) { fDetailed = kTRUE; }
  /** DTOR */
  virtual ~AliFMDv1() {}

  // Required member functions 
  /** Get version number 
      @return always 1 */
  virtual Int_t  IsVersion() const {return 1;}
  /** Member function that is executed each time a hit is made in the 
      FMD.  None-charged particles are ignored.   Dead tracks  are
      ignored.  
      
      The procedure is as follows: 
      - IF NOT track is alive THEN RETURN ENDIF
      - IF NOT particle is charged THEN RETURN ENDIF
      - IF NOT volume name is "STRI" or "STRO" THEN RETURN ENDIF 
      - Get strip number (volume copy # minus 1)
      - Get phi division number (mother volume copy #)
      - Get module number (grand-mother volume copy #)
      - section # = 2 * module # + phi division # - 1
      - Get ring Id from volume name 
      - Get detector # from grand-grand-grand-mother volume name 
      - Get pointer to sub-detector object. 
      - Get track position 
      - IF track is entering volume AND track is inside real shape THEN
      -   Reset energy deposited 
      -   Get track momentum 
      -   Get particle ID # 
      - ENDIF
      - IF track is inside volume AND inside real shape THEN 
      -   Update energy deposited 
      - ENDIF 
      - IF track is inside real shape AND (track is leaving volume,
           or it died, or it is stopped  THEN
      -   Create a hit 
      - ENDIF 
  */
  virtual void   StepManager();
protected:
  /** Translate VMC coordinates to detector coordinates
      @param v        On output, Current position
      @param detector On output, detector #
      @param ring     On output, ring id
      @param sector   On output, sector #
      @param strip    On output, strip #
      @return @c true on success */
  Bool_t VMC2FMD(TLorentzVector& v, UShort_t& detector,
		 Char_t& ring, UShort_t& sector, UShort_t& strip) const;
  /** Translate VMC coordinates to detector coordinates
      @param copy     Volume copy number 
      @param v        On output, Current position
      @param detector On output, detector #
      @param ring     On output, ring id
      @param sector   On output, sector #
      @param strip    On output, strip #
      @return @c true on success */
  Bool_t VMC2FMD(Int_t copy, TLorentzVector& v,
		 UShort_t& detector, Char_t& ring,
		 UShort_t& sector, UShort_t& strip) const;
  /** Check if hit is bad.  A hit is bad if 
      @f[
      \Delta E > |Q|^2 p / m > 1 
      @f]
      holds, where @f$ \Delta E@f$ is the energy loss in this step, 
      @f$ Q@f$ is the particle charge, @f$ p@f$ is the track momentum,
      and @f$ m@f$ is the particle mass.   If a track is marked as
      bad, it's kept in a cache, and can be printed at the end of the
      event. 
      @param trackno Track number 
      @param pdg     PDG particle type ID
      @param absQ    Absolute value of particle charge
      @param p       Track momentum
      @param edep    Energy loss in this step.
      @return @c true if hit is `bad' */
  Bool_t CheckHit(Int_t trackno, Int_t pdg, Float_t absQ, 
		  const TLorentzVector& p, Float_t edep) const;

  Double_t       fCurrentDeltaE;    // The current accumulated energy loss
  TLorentzVector fCurrentV;         // Current production vertex 
  TLorentzVector fCurrentP;         // Current momentum vector 
  Int_t          fCurrentPdg;       // Current PDG code 
  
  ClassDef(AliFMDv1,5)  // Detailed FMD geometry
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
