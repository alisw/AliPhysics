#ifndef ALIJFKTJETFINDERH
#define ALIJFKTJETFINDERH

#include <TObject.h>
#include <vector>
#include <set>

#include "AliJFJetFinder.h"

class TClonesArray;
class TParticle;
class AliJFPreCluster;
class AliJFCluster;
class AliJFClusterDifference;
class AliJFKtJet;
//class multiset<AliJFClusterDifference>;
class AliJFKtJetFinder : public AliJFJetFinder
{
 public:
  AliJFKtJetFinder(Int_t n=50);
  virtual ~AliJFKtJetFinder();

  Int_t Init(TClonesArray *particles);
  Int_t Run();

  Bool_t IsAcceptedParticle(TParticle *p);
  //Bool_t IsAcceptedTower(JFTower*);

  void Debug();
  void Clean();

 protected:
  vector<AliJFPreCluster*> fPreClusterList; //!
  vector<AliJFCluster*>    fClusterList; //!
  multiset<AliJFClusterDifference> fClusterDiffSet; //!

  AliJFKtJet *fJet; //!

  ClassDef(AliJFKtJetFinder,1) //AliJFKtJetFinder class
};

#endif /*ALIJFKTJETFINDERH*/
