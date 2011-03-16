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
#include <iostream>
using std::cout;
using std::endl;

#include "DecayChannel.h"

//_______________________________________________________________________________
DecayChannel::DecayChannel() :
  fMotherPDG(kNonsensePDG),
  fBranchingRatio(0.0),
  fNDaughters(0)
{
  //
  //  default constructor
  //
  for(Int_t i=0; i<kMaxDaughters; i++)
    fDaughtersPDG[i] = kNonsensePDG;
}

//_______________________________________________________________________________
DecayChannel::DecayChannel(const DecayChannel &copy):
  fMotherPDG(copy.fMotherPDG),
  fBranchingRatio(copy.fBranchingRatio),
  fNDaughters(copy.fNDaughters)
{
  //
  // copy constructor
  //
  for(Int_t i=0; i<kMaxDaughters; i++)
    fDaughtersPDG[i] = copy.fDaughtersPDG[i];
}

//_______________________________________________________________________________
DecayChannel::DecayChannel(Int_t mother, Double_t branching, Int_t nDaughters, const Int_t *daughters):
  fMotherPDG(mother),
  fBranchingRatio(branching),
  fNDaughters(0)
{
  //
  // constructor
  //
  for(Int_t i=0; i<kMaxDaughters; i++) {
    if(i<nDaughters)
      fDaughtersPDG[fNDaughters++] = *(daughters+i);
    else
      fDaughtersPDG[i] = kNonsensePDG;
  }
}

//_______________________________________________________________________________
void DecayChannel::SetDaughters(const Int_t *daughters, Int_t n) {
  //
  // set the daughter PDGs
  //
  for(Int_t i=0; i<n; i++) {
    if(i >= kMaxDaughters) {
      cout << "ERROR in DecayChannel::SetDaughters() :" << endl;
      cout << "Number of daughters bigger than the maximum allowed one (" << kMaxDaughters << ") !!" << endl;
    }
    fDaughtersPDG[fNDaughters++] = *(daughters+i);
  }
}

//_______________________________________________________________________________
void DecayChannel::AddDaughter(Int_t pdg) {
  //
  // add a daughter to this channel
  //
  if(fNDaughters >= kMaxDaughters) {
    cout << "ERROR in DecayChannel::AddDaughter() :" << endl;
    cout << "Number of daughters is already >= than the maximum allowed one (" << kMaxDaughters << ") !!" << endl;
  }
  fDaughtersPDG[fNDaughters++] = pdg;
}


//_______________________________________________________________________________
Int_t DecayChannel::GetDaughterPDG(Int_t i) {
  //
  // get the PDG of the i'th daughter
  //
  if((i >= fNDaughters) || (i<0)) {
    cout << "ERROR in DecayChannel::GetDaughterPDG() :" << endl;
    cout << "Daughter index required is too big or less than zero!! There are only " << fNDaughters << " secondaries in this channel !!" << endl;
    return kNonsensePDG;
  }
  return fDaughtersPDG[i];
}

