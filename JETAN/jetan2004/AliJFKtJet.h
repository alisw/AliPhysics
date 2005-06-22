// $Id$

#ifndef ALIJFKTJETH
#define ALIJFKTJETH

#include "AliJFJet.h"

class AliJFCluster;

class AliJFKtJet : public AliJFJet 
{
 public:
  AliJFKtJet(Int_t n=250);
  virtual ~AliJFKtJet();

  void AddJet(AliJFCluster *c);
  void Update();
  
  //void Debug();
  //void Clean();

 protected:

  ClassDef(AliJFKtJet,1) //AliJFKtJet class
};

#endif /*ALIJFKTJETH*/
