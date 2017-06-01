#ifndef AliEventClassifierMult_cxx
#define AliEventClassifierMult_cxx

#include "AliEventClassifierBase.h"

class AliEventClassifierMult : public AliEventClassifierBase {
 public:
  AliEventClassifierMult()
    : AliEventClassifierBase() {}
  AliEventClassifierMult(const char* name, const char* title,
			  Float_t regions[],
			  Int_t lengthRegions,   //how many points are in "regions"
			  Bool_t regionsAreInclusive,
			  Bool_t countCharged,
			  TList *taskOutputList,
                          Int_t collisionSystem);
  virtual ~AliEventClassifierMult() {}

 private:
  // each region is defined by a start and endpoint in eta in a vector.
  // fRegions is a vector of these (region) vectors and can have arbitrary length
  std::vector< std::vector<Float_t> > fRegions;
  Bool_t fRegionsAreInclusive;
  Bool_t fCountCharged;
  void CalculateClassifierValue(AliMCEvent *event, AliStack *stack);

  ClassDef(AliEventClassifierMult, 1);
};

#endif
  
