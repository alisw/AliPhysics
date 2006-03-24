//
// $Id$
//
#ifndef ALIFMD2_H
#define ALIFMD2_H

#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif

//____________________________________________________________________
/** @class AliFMD2 AliFMD2.h <FMD/AliFMD2.h>
    Geometry parameters of the FMD2 detector. This has two rings. 
    @ingroup FMD_base
*/
class AliFMD2 : public AliFMDDetector 
{
protected: 
public: 
  /** Constructor 
      @param inner Pointer to inner ring description 
      @param outer Pointer to outer ring description */
  AliFMD2(AliFMDRing* inner, AliFMDRing* outer);
  /** Initialize the geometry */
  virtual void Init();
  ClassDef(AliFMD2, 1);
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
//
// EOF
//
