//
// Creation date: 2017/08/26
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALISIGNALMC_H
#define ALISIGNALMC_H

#include <TNamed.h>


class AliSignalMC : public TNamed {
   
public:
   
   enum Source {
      kPhysicalPrimary=0,                         // AliMCEvent::IsPhysicalPrimary()
      kSecondaryFromWeakDecay,           // AliMCEvent::IsSecondaryFromWeakDecay()
      kSecondaryFromMaterial,                 // AliMCEvent::IsSecondaryFromMaterial()
      kFromSubsidiaryEvent,                     // AliMCEvent::IsFromSubsidiaryEvent()
      kRadiativeDecay,                              // particle decayed in QED radiative process (e.g. J/psi -> e+ e-  + photons)
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
   void SetProngHistory(UInt_t prong, UInt_t pdgCodes[], Bool_t checkBothCharges[], UInt_t sourceBits[]);
   void SetPDGcode(UInt_t prong, UInt_t generation, Int_t pdgCode, Bool_t checkBothCharges = kFALSE);
   void SetSources(UInt_t prong, UInt_t generation, UInt_t bits);
   void SetSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit);
   
   UInt_t GetCommonAncestorIdx() const {return fCommonAncestorIdx;}
   UInt_t GetNProngs() const {return fNProngs;}
   UInt_t GetNGenerations() const {return fNGenerations;}
   Int_t GetPDGcode(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fPDGcodes[prong][generation] : 0);}
   Int_t GetCheckBothCharges(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fCheckBothCharges[prong][generation] : 0);}
   Bool_t CheckSourceBit(UInt_t prong, UInt_t generation, UInt_t sourceBit) const {return (prong<fNProngs && generation<fNGenerations && sourceBit<kNSources ? fSourceBits[prong][generation] & (UInt_t(1)<<sourceBit) : kFALSE);}
   UInt_t GetSources(UInt_t prong, UInt_t generation) const {return (prong<fNProngs && generation<fNGenerations ? fSourceBits[prong][generation] : 0);}
   
private:

   UInt_t     fNProngs;                                                                 // number of prongs
   UInt_t     fNGenerations;                                                         // number of generations to look back in history
   Int_t     fPDGcodes[kNMaxProngs][kNMaxGenerations];      // PDG codes for all particles in the defined signal
   Bool_t  fCheckBothCharges[kNMaxProngs][kNMaxGenerations];   // include both charge signs of the specified PDG code
   UInt_t   fSourceBits[kNMaxProngs][kNMaxGenerations];     // bit maps encoding physical sources/processes of the particles (see ESource) 
   
   UInt_t     fCommonAncestorIdx;                                             // index of first common ancestor for all prongs; defaults to -1
                                                                                               // all older ancestors are considered to be common
                                                                                               // NOTE: this signal model is not suited for situations where for a given prong there is more than one ancestor
                                                                                               // These situations could be encoded by defining more Source bits
   
   AliSignalMC& operator= (const AliSignalMC &c);
                                                                                               
   ClassDef(AliSignalMC, 1)
};

#endif
