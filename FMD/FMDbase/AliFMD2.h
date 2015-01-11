//
// $Id$
//
#ifndef ALIFMD2_H
#define ALIFMD2_H
/** @file    AliFMD2.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:25:28 2006
    @brief   Geometry parameters of the FMD2 detector. 
*/
// Geometry parameters of the FMD2 detector. This has two rings.
// Other than that it's fairly straight forward.   Needs to make the
// full support stuff. 
//
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif

//____________________________________________________________________
/** @class AliFMD2 AliFMD2.h <FMD/AliFMD2.h>
    @brief Geometry parameters of the FMD2 detector. 
    This has two rings. 
    @todo Flesh out support once it's defined 
    @image html FMD2.png 
    @ingroup FMD_base
*/
class AliFMD2 : public AliFMDDetector 
{
public: 
  /** Constructor 
      @param inner Pointer to inner ring description 
      @param outer Pointer to outer ring description */
  AliFMD2(AliFMDRing* inner, AliFMDRing* outer);
  /** Destructor */
  virtual ~AliFMD2() {}
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
