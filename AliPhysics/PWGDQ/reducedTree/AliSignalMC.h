//
// Creation date: 2017/08/26
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
/*
Monte Carlo signal definition:

Prong(1,0) <-- Prong(1,1) <-- ... Prong(1,X) ... <-- Prong(1,kNMaxGenerations)
                                                         |
Prong(2,0) <-- Prong(2,1) <-- ... Prong(2,X) ... <-- Prong(2,kNMaxGenerations)
                                                         |
Prong(3,0) <-- Prong(3,1) <-- ... Prong(3,X) ... <-- Prong(3,kNMaxGenerations)

The MC signal model foresees the usage of up to kNMaxProngs (currently 3) particles
which are followed back in history for maximum kNMaxGenerations (currently 10).
The most recent common ancestor is at generation X. Older ancestors then X are all assumed to be common. 

For every (prong,generation) doublet, the user can specify a PDG code (see below), a source bit map and an exclusion flag.
If fCheckBothCharges is set to true for a given doublet, then both signs of the PDG code will be used as selection criteria.

See below a few more details.

1.) For the PDG codes, the PYTHIA standard is used. 
A few non-existent PYTHIA codes are used to select more than one PYTHIA code. 
This is a convention, and the way this will be implemented in physics analyses is under the responsability of the user.  

0 - default, accepts all PYTHIA codes
100 - light unflavoured mesons in the code range 100-199
200 -        --"--                               200-299
300 - strange mesons in the code range           300-399
400 - charmed mesons in the code range           400-499
401 - open charm mesons (all D and D* mesons)    400-439
402 - open charm mesons and baryons together     400-439, 4000-4399
403 - all (open- or hidden-) charm hadrons (mesons and baryons) in the range  400-499, 4000-4999    (! no psi' here)
500 - beauty mesons in the code range            500-599
501 - open beauty mesons                         500-549
502 - open beauty mesons and baryons             500-549, 5000-5499
503 - all beauty hadrons                         500-599, 5000-5999
902 - all open charm open beauty mesons+baryons  400-439, 500-549, 4000-4399, 5000-5499
903 - all hadrons in the code range              100-599, 1000-5999
1000 - light unflavoured baryons in the code range 1000-1999
2000 -        --"--                                2000-2999
3000 - strange baryons in the code range           3000-3999
4000 - charmed baryons in the code range           4000-4999
4001 - open charm baryons                          4000-4399
5000 - beauty baryons in the code range            5000-5999
5001 - open beauty baryons                         5000-5499

2.) If the exclusion flags are turned ON then the requested criteria (PDG code and sources) for a given doublet 
are used to exclude the selected cases 

3.) If the selection of both charges is switched ON then the PDG codes act on both particles and anti-particles.

4.) Particles sources implemented:
     See the items defined in Source. Most of the sources are defined using the AliMCEvent methods, but
     other sources can be defined in a customized way. Currently there is a limit of 32 source types which can be defined.
     
*/

#ifndef ALISIGNALMC_H
#define ALISIGNALMC_H

#include <TNamed.h>


class AliSignalMC : public TNamed {
   
public:
   
   enum Source {
      kPhysicalPrimary=0,                         // AliMCEvent::IsPhysicalPrimary()
      kFromBGEvent,                                 // AliMCEvent::IsFromBGEvent()
      kSecondaryFromWeakDecay,           // AliMCEvent::IsSecondaryFromWeakDecay()
      kSecondaryFromMaterial,                 // AliMCEvent::IsSecondaryFromMaterial()
      kFromSubsidiaryEvent,                     // AliMCEvent::IsFromSubsidiaryEvent()
      kRadiativeDecay,                              // particle decayed in QED radiative process (e.g. J/psi -> e+ e-  + photons)
      kFirstInStack,                                    // first particle in stack
      kSecondInStack,                               // second particle in stack
      kFirstTenInStack,                              // one of the first ten particles in stack 
      kNSources
   };
   
   enum Constants {
      kNMaxProngs=3,                       // maximum numbers of prongs allowed for the signal
      kNMaxGenerations=10,            // maximum number of generations to look back
      kPDGnotAssigned=0         
   };
   
   AliSignalMC(Int_t nProngs = 1, Int_t nGenerations = kNMaxGenerations);
   AliSignalMC(const Char_t* name, const Char_t* title, Int_t nProngs = 1, Int_t nGenerations = kNMaxGenerations);
   AliSignalMC(const AliSignalMC &c);
   virtual ~AliSignalMC();
   
   void SetCommonAncestorIdx(UInt_t idx);
   void SetProngHistory(UInt_t prong, UInt_t pdgCodes[], Bool_t checkBothCharges[], UInt_t sourceBits[], Bool_t excludePDG[]=0x0, UInt_t excludeSources[]=0x0, Bool_t useANDonSourceBits[]=0x0);
   void SetPDGcode(UInt_t prong, UInt_t generation, Int_t pdgCode, Bool_t checkBothCharges = kFALSE, Bool_t exclude = kFALSE);
   void SetSources(UInt_t prong, UInt_t generation, UInt_t bits, UInt_t exclude=0, Bool_t useANDonSourceBits=kTRUE);
   void SetSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit, Bool_t exclude=kFALSE);
   void SetUseANDonSourceBits(UInt_t prong, UInt_t generation, Bool_t option=kTRUE);
   
   UInt_t GetCommonAncestorIdx() const {return fCommonAncestorIdx;}
   UInt_t GetNProngs() const {return fNProngs;}
   UInt_t GetNGenerations() const {return fNGenerations;}
   Int_t GetPDGcode(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fPDGcodes[prong][generation] : 0);}
   Bool_t GetCheckBothCharges(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fCheckBothCharges[prong][generation] : 0);}
   Bool_t GetPDGExclude(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fExcludePDG[prong][generation] : 0);}
   Bool_t CheckSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit) const {return (prong<fNProngs && generation<fNGenerations && sourceBit<kNSources ? fSourceBits[prong][generation] & (UInt_t(1)<<sourceBit) : kFALSE);}
   UInt_t GetSources(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fSourceBits[prong][generation] : 0);}
   UInt_t GetSourcesExclude(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fExcludeSource[prong][generation] : 0);}
   Bool_t GetSourceExclude(UInt_t prong, UInt_t generation, UInt_t sourceBit) const {return (prong<fNProngs && generation<fNGenerations && sourceBit<kNSources ? fExcludeSource[prong][generation] & (UInt_t(1)<<sourceBit) : kFALSE);}
   Bool_t GetUseANDonSourceBits(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fUseANDonSourceBitMap[prong][generation] : kFALSE);}
   
   Bool_t TestPDG(UInt_t prong, UInt_t generation, Int_t code);
   
private:

   UInt_t     fNProngs;                                                                 // number of prongs
   UInt_t     fNGenerations;                                                         // number of generations to look back in history
   
   Int_t     fPDGcodes[kNMaxProngs][kNMaxGenerations];      // PDG codes for all particles in the defined signal
   Bool_t  fCheckBothCharges[kNMaxProngs][kNMaxGenerations];   // include both charge signs of the specified PDG code
   Bool_t  fExcludePDG[kNMaxProngs][kNMaxGenerations];     // if TRUE, the specified PDG criteria are used to exclude the particle
   
   UInt_t   fSourceBits[kNMaxProngs][kNMaxGenerations];     // bit maps encoding physical sources/processes of the particles (see ESource) 
   UInt_t   fExcludeSource[kNMaxProngs][kNMaxGenerations];     // if TRUE, the specified source criteria are used to exclude the particle
   Bool_t   fUseANDonSourceBitMap[kNMaxProngs][kNMaxGenerations]; // if TRUE request all enabled source bits (AND); if FALSE request at least one of the enabled source bits (OR)   
   
   UInt_t   fCommonAncestorIdx;                                           // index of first common ancestor for all prongs; defaults to -1
                                                                                               // all older ancestors are considered to be common
                                                                                               // NOTE: this signal model is not suited for situations where for a given prong there is more than one ancestor
                                                                                               // These situations could be encoded by defining more Source bits
   
   AliSignalMC& operator= (const AliSignalMC &c);
                                                                                               
   ClassDef(AliSignalMC, 2)
};

#endif
