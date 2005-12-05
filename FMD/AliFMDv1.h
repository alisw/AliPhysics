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
  AliFMDv1()
    : AliFMD(),
      fCurrentDeltaE(0),
      fCurrentV(),
      fCurrentP(),
      fCurrentPdg(0) { fDetailed = kTRUE; }
  AliFMDv1(const char *name, const char *title="Detailed geometry") 
    : AliFMD(name, title),
      fCurrentDeltaE(0),
      fCurrentV(),
      fCurrentP(),
      fCurrentPdg(0) { fDetailed = kTRUE; }
  virtual ~AliFMDv1() {}

  // Required member functions 
  virtual Int_t  IsVersion() const {return 1;}
  virtual void   StepManager();
protected:
#ifndef USE_PRE_MOVE
  Bool_t VMC2FMD(TLorentzVector& v, UShort_t& detector,
		 Char_t& ring, UShort_t& sector, UShort_t& strip) const;
  Bool_t VMC2FMD(Int_t copy, TLorentzVector& v,
		 UShort_t& detector, Char_t& ring,
		 UShort_t& sector, UShort_t& strip) const;
#endif  

  Double_t       fCurrentDeltaE;    // The current accumulated energy loss
  TLorentzVector fCurrentV;         // Current production vertex 
  TLorentzVector fCurrentP;         // Current momentum vector 
  Int_t          fCurrentPdg;       // Current PDG code 
  
  ClassDef(AliFMDv1,5)  // Detailed FMD geometry
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
