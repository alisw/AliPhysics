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
  fExcludePDG(),
  fSourceBits(),
  fExcludeSource(),
  fUseANDonSourceBitMap(),
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
         fExcludePDG[p][g] = 0;
         fSourceBits[p][g] = 0;
         fExcludeSource[p][g] = 0;
         fUseANDonSourceBitMap[p][g] = kTRUE;
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
fExcludePDG(),
fSourceBits(),
fExcludeSource(),
fUseANDonSourceBitMap(),
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
         fExcludePDG[p][g] = kFALSE;
         fSourceBits[p][g] = 0;
         fExcludeSource[p][g] = 0;
         fUseANDonSourceBitMap[p][g] = kTRUE;
      }
   }
}


//________________________________________________________________________________________________________________
AliSignalMC::AliSignalMC(const AliSignalMC &c) :
fNProngs(c.GetNProngs()),
fNGenerations(c.GetNGenerations()),
fPDGcodes(),
fCheckBothCharges(),
fExcludePDG(),
fSourceBits(),
fExcludeSource(),
fUseANDonSourceBitMap(),
fCommonAncestorIdx(c.GetCommonAncestorIdx())
{
   //
   // copy constructor
   //
   for(UInt_t iprong=0;iprong<fNProngs;++iprong) {
      for(UInt_t igen=0; igen<fNGenerations;++igen) {
         fPDGcodes[iprong][igen] = c.GetPDGcode(iprong,igen);
         fCheckBothCharges[iprong][igen] = c.GetCheckBothCharges(iprong,igen);
         fExcludePDG[iprong][igen] = c.GetPDGExclude(iprong, igen);
         fSourceBits[iprong][igen] = c.GetSources(iprong,igen);
         fExcludeSource[iprong][igen] = c.GetSourcesExclude(iprong, igen);
         fUseANDonSourceBitMap[iprong][igen] = c.GetUseANDonSourceBits(iprong, igen);
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
void AliSignalMC::SetProngHistory(UInt_t prong, UInt_t pdgCodes[], Bool_t checkBothCharges[], UInt_t sourceBits[], 
                                  Bool_t excludePDG[] /* = 0x0*/, UInt_t excludeSources[] /* = 0x0*/, Bool_t useANDonSourceBits[] /*=0x0*/) {
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
      if(excludePDG)
         fExcludePDG[prong][i] = excludePDG[i];
      if(excludeSources)
         fExcludeSource[prong][i] = excludeSources[i];
      if(useANDonSourceBits)
         fUseANDonSourceBitMap[prong][i] = useANDonSourceBits[i];
   }
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetPDGcode(UInt_t prong, UInt_t generation, Int_t pdgCode, Bool_t checkBothCharges /*= kFALSE*/, Bool_t exclude /*=kFALSE*/) {
   //
   // set PDG code for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   fPDGcodes[prong][generation] = pdgCode;
   fCheckBothCharges[prong][generation] = checkBothCharges;
   fExcludePDG[prong][generation] = exclude;
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetSources(UInt_t prong, UInt_t generation, UInt_t bits, UInt_t exclude /*=0*/, Bool_t useANDonSourceBits /*=kTRUE*/) {
   //
   // set the entire sources bit map for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   fSourceBits[prong][generation] = bits;
   fExcludeSource[prong][generation] = exclude;
   fUseANDonSourceBitMap[prong][generation] = useANDonSourceBits;
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit, Bool_t exclude /*=kFALSE*/) {
   //
   // set a specific source bit for a prong at a given generation
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   if(sourceBit>=kNSources) return;
   fSourceBits[prong][generation] |= (UInt_t(1)<<sourceBit);
   // NOTE: The exclude factor is set here via the inidividual source bit, but actually overrides the flag for all sources (in case multiple sources have been defined)
   // TODO: implement a bit map corresponding to each source, to hold an exclude flag for each source
   if(exclude) fExcludeSource[prong][generation] |= (UInt_t(1)<<sourceBit);
}

//________________________________________________________________________________________________________________
void AliSignalMC::SetUseANDonSourceBits(UInt_t prong, UInt_t generation, Bool_t option /*=kTRUE*/) {
   //
   // set the option to be used when evaluating the source bits
   //
   if(prong<0 || prong>=fNProngs) return;
   if(generation<0 || generation>=fNGenerations) return;
   fUseANDonSourceBitMap[prong][generation] = option;
}

//________________________________________________________________________________________________________________
Bool_t AliSignalMC::TestPDG(UInt_t prong, UInt_t generation, Int_t pdg) {
   //
   // test if the code pdg matches this signal
   //
   if(prong<0 || prong>=fNProngs) return kFALSE;
   if(generation<0 || generation>=fNGenerations) return kFALSE;
   
   Bool_t decision = kTRUE;
   Int_t absPDG = TMath::Abs(pdg);
   
   switch(TMath::Abs(fPDGcodes[prong][generation])) {
      case kPDGnotAssigned:
         // PDG not required (any code will do fine)
         break;
      case 100:     // light flavoured mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=100 && absPDG<=199;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=100 && pdg<=199;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-199 && pdg<=-100;
         }
         break;
      case 1000:     // light flavoured baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=1000 && absPDG<=1999;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=1000 && pdg<=1999;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-1999 && pdg<=-1000;
         }
         break;
      case 200:     // light flavoured mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=200 && absPDG<=299;
         else {
            if(fPDGcodes[prong][generation]>0)decision = pdg>=200 && pdg<=299;
            if(fPDGcodes[prong][generation]<0)decision = pdg>=-299 && pdg<=-200;
         }
         break;
      case 2000:     // light flavoured baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=2000 && absPDG<=2999;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=2000 && pdg<=2999;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-2999 && pdg<=-2000;
         }
         break;
      case 300:     // all strange mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=300 && absPDG<=399;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=300 && pdg<=399;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-399 && pdg<=-300;
         }
         break;
      case 3000:     // all strange baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=3000 && absPDG<=3999;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=3000 && pdg<=3999;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-3999 && pdg<=-3000;
         }
         break;
      case 400:     // all charmed mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=400 && absPDG<=499;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=400 && pdg<=499;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-499 && pdg<=-400;
         }
         break;
      case 401:     // open charm mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=400 && absPDG<=439;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=400 && pdg<=439;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-439 && pdg<=-400;
         }
         break;
      case 402:     // open charm mesons and baryons together
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=400 && absPDG<=439) ||
            (absPDG>=4000 && absPDG<=4399);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=400 && pdg<=439) ||
               (pdg>=4000 && pdg<=4399);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-439 && pdg<=-400) ||
               (pdg>=-4399 && pdg<=-4000);
         }
         break;
      case 403:     // all charm hadrons
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=400 && absPDG<=499) ||
            (absPDG>=4000 && absPDG<=4999);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=400 && pdg<=499) ||
               (pdg>=4000 && pdg<=4999);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-499 && pdg<=-400) ||
               (pdg>=-4999 && pdg<=-4000);
         }
         break;
      case 404:     // charged open charmed mesons NO s quark
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=410 && absPDG<=419);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=410 && pdg<=419);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-419 && pdg<=-410);
         }
         break;
      case 405:     // neutral open charmed mesons
         if(fCheckBothCharges[prong][generation])
            decision =absPDG>=420 && absPDG<=429;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=420 && pdg<=429 ;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-429 && pdg<=-420;
         }
         break;
         
      case 406:     // charged open charmed mesons with s quark
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=430 && absPDG<=439);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=430 && pdg<=439);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-439 && pdg<=-430);
         }
         break;
         
      case 4000:     // all charmed baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=4000 && absPDG<=4999;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=4000 && pdg<=4999;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-4999 && pdg<=-4000;
         }
         break;
      case 4001:     // open charm baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=4000 && absPDG<=4399;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=4000 && pdg<=4399;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-4399 && pdg<=-4000;
         }
         break;
      case 500:      // all beauty mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=500 && absPDG<=599;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=500 && pdg<=599;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-599 && pdg<=-500;
         }
         break;
      case 501:      // open beauty mesons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=500 && absPDG<=549;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=500 && pdg<=549;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-549 && pdg<=-500;
         }
         break;
      case 502:      // open beauty mesons and baryons
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=500 && absPDG<=549) ||
            (absPDG>=5000 && absPDG<=5499);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=500 && pdg<=549) ||
               (pdg>=5000 && pdg<=5499);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-549 && pdg<=-500) ||
               (pdg>=-5499 && pdg<=-5000);
         }
         break;
      case 503:      // all beauty hadrons
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=500 && absPDG<=599) ||
            (absPDG>=5000 && absPDG<=5999);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=500 && pdg<=599) ||
               (pdg>=5000 && pdg<=5999);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-599 && pdg<=-500) ||
               (pdg>=-5999 && pdg<=-5000);
         }
         break;
         
      case 504:     // neutral open beauty mesons NO s quark
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=510 && absPDG<=519);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=510 && pdg<=519);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-519 && pdg<=-510);
         }
         break;
      case 505:     // charged open beauty mesons
         if(fCheckBothCharges[prong][generation])
            decision =absPDG>=520 && absPDG<=529;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=520 && pdg<=529 ;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-529 && pdg<=-520;
         }
         break;
         
      case 506:     // charged open beauty mesons with s quark
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=530 && absPDG<=539);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=530 && pdg<=539);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-539 && pdg<=-530);
         }
         break;
         
      case 5000:      // all beauty baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=5000 && absPDG<=5999;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=5000 && pdg<=5999;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-5999 && pdg<=-5000;
         }
         break;
      case 5001:      // open beauty baryons
         if(fCheckBothCharges[prong][generation])
            decision = absPDG>=5000 && absPDG<=5499;
         else {
            if(fPDGcodes[prong][generation]>0) decision = pdg>=5000 && pdg<=5499;
            if(fPDGcodes[prong][generation]<0) decision = pdg>=-5499 && pdg<=-5000;
         }
         break;
      case 902:      // // open charm,beauty  mesons and baryons together
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=400 && absPDG<=439) ||
            (absPDG>=4000 && absPDG<=4399) ||
            (absPDG>=500 && absPDG<=549) ||
            (absPDG>=5000 && absPDG<=5499);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=400 && pdg<=439) ||
               (pdg>=4000 && pdg<=4399)      ||
               (pdg>=500 && pdg<=549)        ||
               (pdg>=5000 && pdg<=5499);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-439 && pdg<=-400) ||
               (pdg>=-4399 && pdg<=-4000)      ||
               (pdg>=-549 && pdg<=-500)        ||
               (pdg>=-5499 && pdg<=-5000);
         }
         break;
      case 903:      // // all hadrons in the code range 100-599, 1000-5999
         if(fCheckBothCharges[prong][generation])
            decision = (absPDG>=100 && absPDG<=599) ||
            (absPDG>=1000 && absPDG<=5999);
         else {
            if(fPDGcodes[prong][generation]>0) decision = (pdg>=100 && pdg<=599) ||
               (pdg>=1000 && pdg<=5999);
            if(fPDGcodes[prong][generation]<0) decision = (pdg>=-599 && pdg<=-100) ||
               (pdg>=-5999 && pdg<=-1000);
         }
         break;
      default:          // all specific cases
         if(fCheckBothCharges[prong][generation])
            decision = (TMath::Abs(fPDGcodes[prong][generation])==absPDG);
         else
            decision = (fPDGcodes[prong][generation]==pdg);
   }
   
   if(fPDGcodes[prong][generation] && fExcludePDG[prong][generation]) decision = !decision;
   return decision;
}
