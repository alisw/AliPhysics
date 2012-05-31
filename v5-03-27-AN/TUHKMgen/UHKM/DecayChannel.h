////////////////////////////////////////////////////////////////////////////////////////////////
//  Copyright   : The FASTMC and SPHMC Collaboration                                          //
//  Author      : Ionut Cristian Arsene                                                       //
//  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania  //
//  e-mail      : i.c.arsene@fys.uio.no                                                       //
//  Date        : 2007/05/30                                                                  //
//                                                                                            //
//  This class is using the particle and decays lists provided by the                         //
//  THERMINATOR (Computer Physics Communications 174 669 (2006)) and                          //
//  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.                    //
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DECAYCHANNEL_H
#define DECAYCHANNEL_H

#include <Rtypes.h>

const Int_t kMaxDaughters = 3;
const Int_t kNonsensePDG = 1000000000;

class DecayChannel {
 public:
  DecayChannel();                                                                           // default constructor
  DecayChannel(const DecayChannel &copy);                                                   // copy constructor
  DecayChannel(Int_t mother, Double_t branching, Int_t nDaughters, const Int_t *daughters); // explicit constructor
  ~DecayChannel() {};                                                                       // destructor
  
  void     SetMotherPDG(Int_t value)              {fMotherPDG = value;}
  void     SetBranching(Double_t value)           {fBranchingRatio = value;}
  void     SetDaughters(const Int_t *values, Int_t n);
  void     AddDaughter(Int_t pdg);
  Int_t    GetMotherPDG() const                   {return fMotherPDG;}
  Double_t GetBranching() const                   {return fBranchingRatio;}
  Int_t    GetNDaughters() const                  {return fNDaughters;}
  const Int_t*   GetDaughters() const                   {return fDaughtersPDG;}
  Int_t    GetDaughterPDG(Int_t i);                                                         // i --> must be the zero-based index of daughter

 private:
  Int_t    fMotherPDG;                          // PDG code of the mother particle
  Double_t fBranchingRatio;                     // branching ratio
  Int_t    fNDaughters;                         // number of daughter particles
  Int_t    fDaughtersPDG[kMaxDaughters];        // array with daughters PDG
};

#endif
