#ifndef ALIFMD1_H
#define ALIFMD1_H
//
// $Id$
//
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif
class AliFMDRing;

//__________________________________________________________________
/** @class AliFMD1 AliFMD1.h <FMD/AliFMD1.h>
    Geometry description and parameters of the FMD1 detector. 
    The FMD1 has only one ring.     
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
  virtual void Init() { AliFMDDetector::Init(); }
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
