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
class AliFMDv0 : public AliFMD 
{
public:
  AliFMDv0() {}
  AliFMDv0(const char *name, const char *title="Coarse geometry") 
    : AliFMD(name, title, false)
  {}
  virtual ~AliFMDv0() 
  {}

  // Required member functions 
  virtual Int_t  IsVersion() const {return 0;}
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
