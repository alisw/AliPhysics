#ifndef ALIHBTPAIR_H
#define ALIHBTPAIR_H
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliHBTPair
//
// class implements pair of particles and taking care of caluclation (almost)
// all of pair properties (Qinv, InvMass,...)
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////

#include <AliAODPair.h>

class AliVAODParticle;

class AliHBTPair: public AliAODPair
{
 public:
   AliHBTPair(Bool_t rev = kFALSE); //contructor
   AliHBTPair(AliVAODParticle* part1, AliVAODParticle* part2, Bool_t rev = kFALSE); //contructor
   AliHBTPair(const AliHBTPair& in);
   
   virtual ~AliHBTPair(){}
   
   AliHBTPair& operator=(const AliHBTPair& in);

   void     Changed();
   Double_t GetWeight();
   
 protected:

   Double_t fWeight;//Value of the weight
   Bool_t   fWeightNotCalc;//flag indicating if fWeight is calculated for current pair
   
 private:
  ClassDef(AliHBTPair,2)
};

inline
void AliHBTPair::Changed()
{
 // Resel all calculations (flags)
 
 AliAODPair::Changed();
 fWeightNotCalc = kTRUE;
}
#endif
