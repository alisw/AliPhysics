#ifndef ALIFMDV0_H
#define ALIFMDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Manager class for the FMD - Coarse version. 
//
#ifndef ALIFMD_H 
# include "AliFMD.h"
#endif

//____________________________________________________________________
/** Forward Multiplicity Detector based on Silicon wafers. This class 
    contains the base procedures for the Forward Multiplicity detector
    Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of 
    which has 1 or 2 rings of silicon sensors.  
                                                           
    This contains the coarse version of the FMD - that is, the
    simulation produces no hits in the FMD volumes, and the sensors
    are not divided into strips and sectors.   Useful for material
    budget calculations
    @ingroup FMD_sim
 */
class AliFMDv0 : public AliFMD 
{
public:
  /** CTOR */
  AliFMDv0() {}
  /** CTOR
      @param name Name
      @param title Title  */
  AliFMDv0(const char *name, const char *title="Coarse geometry") 
    : AliFMD(name, title)
  {}
  /** DTOR */
  virtual ~AliFMDv0() {}

  // Required member functions 
  /** @return Version number - always 0 */
  virtual Int_t  IsVersion() const {return 0;}
  /** Function called at each hit.  Empty */
  virtual void   StepManager() {}

  ClassDef(AliFMDv0,1) // Coarse FMD geometry 
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
