/*
  Copyright   : The FASTMC and SPHMC Collaboration
  Author      : Ionut Cristian Arsene 
  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
  e-mail      : i.c.arsene@fys.uio.no
  Date        : 2007/05/30

  This class is using the particle and decays lists provided by the 
  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
*/

#ifndef DECAY_CHANNEL
#define DECAY_CHANNEL

#include "Rtypes.h"

const Int_t kMaxDaughters = 3;
const Int_t kNonsensePDG = 1000000000;

class DecayChannel {
 private:
  Int_t    fMotherPDG;
  Double_t fBranchingRatio;
  Int_t    fNDaughters;
  Int_t    fDaughtersPDG[kMaxDaughters];
  
 public:
  DecayChannel();                                                                           // default constructor
  DecayChannel(const DecayChannel &copy);                                                   // copy constructor
  DecayChannel(Int_t mother, Double_t branching, Int_t nDaughters, Int_t *daughters);       // explicit constructor
  ~DecayChannel() {};                                                                       // destructor
  
  void     SetMotherPDG(Int_t value)              {fMotherPDG = value;}
  void     SetBranching(Double_t value)           {fBranchingRatio = value;}
  void     SetDaughters(Int_t *values, Int_t n);
  void     AddDaughter(Int_t pdg);
  Int_t    GetMotherPDG()                         {return fMotherPDG;}
  Double_t GetBranching()                         {return fBranchingRatio;}
  Int_t    GetNDaughters()                        {return fNDaughters;}
  Int_t*   GetDaughters()                         {return fDaughtersPDG;}
  Int_t    GetDaughterPDG(Int_t i);                                                         // i --> must be the zero-based index of daughter
};

#endif
