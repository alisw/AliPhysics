#ifndef ALIFMDMULTNAIIVE_H
#define ALIFMDMULTNAIIVE_H
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
// Class to do multiplicity reconstruction using the Naiive method.
// That is, we count the number of empty strips in a region, and
// derive the charge particle multiplicity from that number. 
// 

#ifndef ALIFMDMULTALGORITHM_H
# include "AliFMDMultAlgorithm.h"
#endif

//____________________________________________________________________
class AliFMDMultNaiive : public AliFMDMultAlgorithm
{
public:
  AliFMDMultNaiive();
  virtual ~AliFMDMultNaiive() {}
  
  virtual void PreRun(AliFMD* fmd);
  virtual void PreEvent(TTree* treeR, Float_t ipZ);
  virtual void ProcessDigit(AliFMDDigit* digit, 
			    Float_t eta, 
			    Float_t phi, 
			    UShort_t counts);
  void    SetGain(Float_t g) { fGain = g; }
  void    SetEdepMip(Float_t e) { fEdepMip = e; }
  Float_t GetGain() const { return fGain; }
  Float_t GetEdepMip() const { return fEdepMip; }
protected:
  Float_t Adc2Energy(AliFMDDigit* digit, Float_t eta, UShort_t count);
  Float_t Energy2Multiplicity(AliFMDDigit* digit, Float_t edep);
  Float_t fGain;           //  GeV per ADC count
  Float_t fEdepMip;        //  Energy deposited per MIP
  
  ClassDef(AliFMDMultNaiive, 0) // Naiive algorithm
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
