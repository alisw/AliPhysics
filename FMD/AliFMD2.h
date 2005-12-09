//
// $Id$
//
#ifndef ALIFMD2_H
#define ALIFMD2_H

#ifndef ALIFMDDETECTOR_H
# include "AliFMDDetector.h"
#endif

// Geometry description and parameters of the FMD2
// detector. 
// This has two rings. 
class AliFMD2 : public AliFMDDetector 
{
protected: 
public: 
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
