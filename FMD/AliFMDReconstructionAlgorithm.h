#ifndef ALIFMDRECONSTRUCTIONALGORITHM_H
#define ALIFMDRECONSTRUCTIONALGORITHM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/* $Id$ */
//____________________________________________________________________
// 
// Base class for FMD reconstruction algorithms. 
//
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 

#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliFMDDigit;

//____________________________________________________________________
class AliFMDReconstructionAlgorithm : public TNamed
{
public:
  AliFMDReconstructionAlgorithm(const char* name, const char* title);
  virtual ~AliFMDReconstructionAlgorithm() {}
  
  virtual void PreEvent() {};
  virtual void ProcessDigit(AliFMDDigit* digit, 
			    Float_t eta, 
			    Float_t phi, 
			    UShort_t counts) = 0;
  virtual void PostEvent() {};
protected:
  ClassDef(AliFMDReconstructionAlgorithm, 0) // Base class for algorithms
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
