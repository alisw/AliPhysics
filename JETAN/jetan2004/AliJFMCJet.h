// $Id$

#ifndef ALIJFMCJETH
#define ALIJFMCJETH

#include "AliJFJet.h"

class AliJFMCJet : public AliJFJet 
{
 public:
  AliJFMCJet(Int_t n=250);
  virtual ~AliJFMCJet();

  void AddParticle(TParticle *p);
  void Update();
  
  void Debug();
  void Clean();

 protected:

  ClassDef(AliJFMCJet,1) //AliJFMCJet class
};

#endif /*ALIJFMCJETH*/
