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
#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif

//____________________________________________________________________
/** @class AliFMD2 AliFMD2.h <FMD/AliFMD2.h>
    @brief Geometry parameters of the FMD2 detector. 
    This has two rings. 
    @image html FMD2.png 
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
