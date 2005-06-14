#include "AliHBTPair.h"
//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
// class AliHBTPair
//
// class implements pair of particles and taking care of caluclation (almost)
// all of pair properties (Qinv, InvMass,...)
// 
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////

#include "AliVAODParticle.h"
#include "AliHBTWeights.h"

ClassImp(AliHBTPair)

/************************************************************************/
AliHBTPair::AliHBTPair(Bool_t rev):
 AliAODPair(rev),
 fWeight(0.0),
 fWeightNotCalc(kTRUE)
 {
//Constructor
  
 }
/************************************************************************/

AliHBTPair::AliHBTPair(AliVAODParticle* part1, AliVAODParticle* part2, Bool_t rev):
 AliAODPair(part1,part2,rev),
 fWeight(0.0),
 fWeightNotCalc(kTRUE)
 {
//constructor
  
 }
/************************************************************************/
AliHBTPair::AliHBTPair(const AliHBTPair& in):
 AliAODPair(in),
 fWeight(in.fWeight),
 fWeightNotCalc(in.fWeightNotCalc)
{
 //cpy constructor
}
/************************************************************************/

AliHBTPair& AliHBTPair::operator=(const AliHBTPair& in)
{
 //Assigment operator
 in.Copy(*this);
 return *this;
}
/************************************************************************/

Double_t AliHBTPair::GetWeight()
{
  //returns and buffers weight for this pair
  if (fWeightNotCalc)
   {
      fWeight = AliHBTWeights::Weight(this);
      fWeightNotCalc = kFALSE;  
   }
  return fWeight; 
}
/************************************************************************/
