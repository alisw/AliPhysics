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
/** Geometry description and parameters of the FMD1
    detector. 
    
    The FMD1 only has one ring.     
*/
class AliFMD1 : public AliFMDDetector 
{
public:
  AliFMD1(AliFMDRing* inner);
  virtual ~AliFMD1() {}
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
