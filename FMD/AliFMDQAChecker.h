#ifndef ALIFMDQACHECKER_H
#define ALIFMDQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */

class TFile; 
class TH1F; 
class TH1I; 

#include "AliQACheckerBase.h"

class AliFMDQAChecker : public AliQACheckerBase 
{
public:
  AliFMDQAChecker() 
    : AliQACheckerBase("FMD","FMD Quality Assurance Checker") 
  {}          
 
  virtual ~AliFMDQAChecker() {}

private:
  ClassDef(AliFMDQAChecker,0)  // Yves? what to do? 
};

#endif // AliFMDQAChecker_H
// Local Variables:
//  mode: c++
// End:
