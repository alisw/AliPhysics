#ifndef ALIFMDPOISSONALGORITHM_H
#define ALIFMDPOISSONALGORITHM_H
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
// Class to do multiplicity reconstruction using the Poisson method.
// That is, we count the number of empty strips in a region, and
// derive the charge particle multiplicity from that number. 
// 

#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliFMDDigit;

//____________________________________________________________________
class AliFMDPoissonAlgorithm : public TNamed
{
public:
  AliFMDPoissonAlgorithm();
  virtual ~AliFMDPoissonAlgorithm() {}
  
  virtual void Reset();
  virtual void ProcessDigit(AliFMDDigit* digit, Float_t ipZ);
  
protected:
  ClassDef(AliFMDPoissonAlgorithm, 0) // Poisson algorithm
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
