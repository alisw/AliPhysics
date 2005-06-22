// $Id$

#ifndef ALIJFMCJETFINDERH
#define ALIJFMCJETFINDERH

#include "AliJFJetFinder.h"
#include <TObject.h>
class TClonesArray;
class TParticle;
class AliJFMCJet;

class AliJFMCJetFinder : public AliJFJetFinder
{
 public:
  AliJFMCJetFinder(Int_t n=25);
  virtual ~AliJFMCJetFinder();

  Int_t Init(TClonesArray *particles);
  Int_t Run();

  void Debug();
  //void Clean();

 protected:
  void FollowDaughters(Int_t first, Int_t last);

  TClonesArray *fParticles; //!
  AliJFMCJet *fJet; //!
  
  ClassDef(AliJFMCJetFinder,1) //AliJFMCJetFinder class
};

#endif /*ALIJFMCJETFINDERH*/
