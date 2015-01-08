#ifndef ALIFMD1_H
#define ALIFMD1_H
/* $Id$ */
/** @file    AliFMD1.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:00:56 2006
    @brief   Declaration of FMD1 declaration     
*/
// Geometry description and parameters of the FMD1 detector. The FMD1
// has only one ring.     Need to flesh out the support once it's
// defined 
//
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif
class AliFMDRing;

//__________________________________________________________________
/** @class AliFMD1 AliFMD1.h <FMD/AliFMD1.h>
    @brief Geometry description and parameters of the FMD1 detector. 
    The FMD1 has only one ring.     
    @todo Need to flesh out the support once it's defined 
    @image html FMD1.png 
    @ingroup FMD_base
*/
class AliFMD1 : public AliFMDDetector 
{
public:
  /** Constructor 
      @param inner Pointer to inner ring description  */
  AliFMD1(AliFMDRing* inner);
  /** Destructor */
  virtual ~AliFMD1() {}
  /** Initialize */
  virtual void Init();
  ClassDef(AliFMD1,1)
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
