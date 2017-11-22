/*
 * **********************************************************
 * Implementation of AliSignalMC class.
 * Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
 * 2017/08/26
 *********************************************************
 */

#ifndef ALISIGNALMC_H
#include "AliSignalMC.h"
#endif

#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliSignalMC)


//________________________________________________________________________________________________________________
AliSignalMC::AliSignalMC(Int_t nProngs /*=1*/, Int_t nGenerations /*=kNMaxGenerations*/) :
  TNamed(),
  fNProngs(1),
  fNGenerations(kNMaxGenerations),
  fPDGcodes(),
  fCheckBothCharges(),
  fSourceBits(),
  fCommonAncestorIdx(0)
{
   if(nProngs<1 || nProngs>kNMaxProngs) {
      cout << "WARNING in AliSignalMC(nProngs, nGenerations):  Attempt to build an AliSignalMC object with " 
      << nProngs << " prongs!" << endl; 
      cout << "WARNING in AliSignalMC(nProngs, nGenerations):  Allowed number of prongs: 1 - " << kNMaxProngs << endl;
      fNProngs = 1;
   }
   else 
      fNProngs = nProngs;
   if(nGenerations<1 || nGenerations>kNMaxGenerations) {
      cout << "WARNING in AliSignalMC(nProngs, nGenerations):  Attempt to build an AliReducedSignal object with " 
      << nGenerations << " generations!" << endl; 
      cout << "WARNING in AliSignalMC(nProngs, nGenerations):  Allowed number of generations: 1 - " << kNMaxGenerations << endl;
      fNGenerations = kNMaxGenerations;
   }
   else 
      fNGenerations = nGenerations;
   
   for(UInt_t p=0;p<fNProngs;++p) {
      for(UInt_t g=0;g<fNGenerations;++g) {
         fPDGcodes[p][g] = kPDGnotAssigned;
         fCheckBothCharges[p][g] = kFALSE;
         fSourceBits[p][g] = 0;
      }
   }
}


//________________________________________________________________________________________________________________
AliSignalMC::AliSignalMC(const Char_t* name, const Char_t* title, Int_t nProngs /*=1*/, Int_t nGenerations /*=kNMaxGenerations*/) :
TNamed(name, title),
fNProngs(1),
fNGenerations(kNMaxGenerations),
fPDGcodes(),
fCheckBothCharges(),
fSourceBits(),
fCommonAncestorIdx(0)
{
   if(nProngs<1 || nProngs>kNMaxProngs) {
      cout << "WARNING in AliSignalMC(name, title, nProngs, nGenerations):  Attempt to build an AliSignalMC object with " 
      << nProngs << " prongs!" << endl; 
      cout << "WARNING in AliSignalMC(name, title, nProngs, nGenerations):  Allowed number of prongs: 1 - " << kNMaxProngs << endl;
      fNProngs = 1;
   }
   else
      fNProngs = nProngs;
   if(nGenerations<1 || nGenerations>kNMaxGenerations) {
      cout << "WARNING in AliSignalMC(name, title, nProngs, nGenerations):  Attempt to build an AliSignalMC object with " 
      << nGenerations << " generations!" << endl; 
      cout << "WARNING in AliSignalMC(name, title, nProngs, nGenerations):  Allowed number of generations: 1 - " << kNMaxGenerations << endl;
      fNGenerations = kNMaxGenerations;
   }
   else 
      fNGenerations = nGenerations;
   
   for(UInt_t p=0;p<fNProngs;++p) {
      for(UInt_t g=0;g<fNGenerations;++g) {
         fPDGcodes[p][g] = kPDGnotAssigned;
         fCheckBothCharges[p][g] = kFALSE;
         fSourceBits[p][g] = 0;
      }
   }
}


//________________________________________________________________________________________________________________
AliSignalMC::AliSignalMC(const AliSignalMC &c) :
fNProngs(c.GetNProngs()),
fNGenerations(c.GetNGenerations()),
fPDGcodes(),
fCheckBothCharges(),
fSourceBits(),
fCommonAncestorIdx(c.GetCommonAncestorIdx())
{
   //
   // copy constructor
   //
   for(UInt_t iprong=0;iprong<fNProngs;++iprong) {
      for(UInt_t igen=0; igen<fNGenerations;++igen) {
         fPDGcodes[iprong][igen] = c.GetPDGcode(iprong,igen);
         fCheckBothCharges[iprong][igen] = c.GetCheckBothCharges(iprong,igen);
         fSourceBits[iprong][igen] = c.GetSources(iprong,igen);
      }
   }
}

//________________________________________________________________________________________________________________
AliSignalMC::~AliSignalMC() {
   //
   // destructor
   //
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetCommonAncestorIdx(UInt_t idx) {
   //
   // set the common ancestor
   //
   fCommonAncestorIdx = idx;
   if(idx>=fNGenerations) fCommonAncestorIdx = fNGenerations-1;
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetProngHistory(UInt_t prong, UInt_t pdgCodes[], Bool_t checkBothCharges[], UInt_t sourceBits[]) {
   //
   // set the entire prong history
   //
   if(prong>=fNProngs) {
      cout << "WARNING in AliSignalMC::SetProngHistory:  Invalid prong index, should be a number between 0 and " << fNProngs << endl;
      return;
   }

   for(UInt_t i=0; i<fNGenerations; ++i) {
      fPDGcodes[prong][i] = pdgCodes[i];
      fCheckBothCharges[prong][i] = checkBothCharges[i];
      fSourceBits[prong][i] = sourceBits[i];
   }
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetPDGcode(UInt_t prong, UInt_t generation, Int_t pdgCode, Bool_t checkBothCharges /*= kFALSE*/) {
   //
   // set PDG code for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   fPDGcodes[prong][generation] = pdgCode;
   fCheckBothCharges[prong][generation] = checkBothCharges;
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetSources(UInt_t prong, UInt_t generation, UInt_t bits) {
   //
   // set the entire sources bit map for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   fSourceBits[prong][generation] = bits;
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit) {
   //
   // set a specific source bit for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   if(sourceBit>=kNSources) return;
   fSourceBits[prong][generation] |= (UInt_t(1)<<sourceBit);
}

