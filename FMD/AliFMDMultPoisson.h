#ifndef ALIFMDMULTPOISSON_H
#define ALIFMDMULTPOISSON_H
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

#ifndef ALIFMDMULTALGORITHM_H
# include "AliFMDMultAlgorithm.h"
#endif
#ifndef ALIFMDBOOLMAP_H
# include "AliFMDBoolMap.h"
#endif

//____________________________________________________________________
class AliFMDMultPoisson : public AliFMDMultAlgorithm
{
public:
  AliFMDMultPoisson();
  virtual ~AliFMDMultPoisson() {}
  
  virtual void PreEvent(TTree* treeR, Float_t ipZ);
  virtual void ProcessDigit(AliFMDDigit* digit, 
			    Float_t eta, 
			    Float_t phi, 
			    UShort_t counts);
  virtual void PostEvent();
  
  void         SetDeltaEta(Float_t deta=.1)  { fDeltaEta = deta;  }
  void         SetDeltaPhi(Float_t dphi=360) { fDeltaPhi = dphi;  } 
  void         SetThreshold(UShort_t t=6)    { fThreshold = t; }
protected:
  AliFMDBoolMap fEmpty;          //! Map of empty channels
  Float_t       fCurrentVertexZ; //! Current IP's Z-coordinate
  Float_t       fDeltaEta;       // Bin size in eta
  Float_t       fDeltaPhi;       // Bin size in phi
  UShort_t      fThreshold;      // Threshold for Poisson recon.
  
  ClassDef(AliFMDMultPoisson, 0) // Poisson algorithm
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
