#ifndef ALIFMDV1_H
#define ALIFMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Manager class for the FMD - Detailed version. 
//
#ifndef ALIFMD_H 
# include "AliFMD.h"
#endif
#ifndef ROOT_TLorentzVector
# include <TLorentzVector.h>
#endif
 
//____________________________________________________________________
class AliFMDv1 : public AliFMD 
{
public:
  AliFMDv1() { fDetailed = kTRUE; }
  AliFMDv1(const char *name, const char *title="Detailed geometry") 
    : AliFMD(name, title) 
  { fDetailed = kTRUE; }
  virtual ~AliFMDv1() {}

  // Required member functions 
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   StepManager();
protected:
  Double_t   fCurrentDeltaE;        // The current accumelated energy loss
  TLorentzVector fCurrentV;         // Current production vertex 
  TLorentzVector fCurrentP;         // Current momentum vector 
  Int_t          fCurrentPdg;       // Current PDG code 
  
  ClassDef(AliFMDv1,4)  // Detailed FMD geometry
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
