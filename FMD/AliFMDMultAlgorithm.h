#ifndef ALIFMDMULTALGORITHM_H
#define ALIFMDMULTALGORITHM_H
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
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 

#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliFMDDigit;
class TTree;
class AliFMD;
class TClonesArray;

//____________________________________________________________________
class AliFMDMultAlgorithm : public TNamed
{
public:
  AliFMDMultAlgorithm(const char* name, const char* title);
  virtual ~AliFMDMultAlgorithm();
  
  virtual void PreRun(AliFMD* fmd) { fFMD = fmd; }
  virtual void PreEvent(TTree* treeR, Float_t ipZ);
  virtual void ProcessDigit(AliFMDDigit* digit, 
			    Float_t eta, 
			    Float_t phi, 
			    UShort_t counts) = 0;
  virtual void PostEvent() {}
  virtual void PostRun() {}
protected:
  TTree*        fTreeR;     //! Reconstruction tree  
  TClonesArray* fMult;      //! Reconstructed multiplicities
  Int_t         fNMult;     //! Number of reconstructed multiplicities
  AliFMD*       fFMD;       //! Detector information 
  
  ClassDef(AliFMDMultAlgorithm, 0) // Base class for algorithms
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
