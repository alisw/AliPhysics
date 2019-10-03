#ifndef AliEventClassifierSpherocity_cxx
#define AliEventClassifierSpherocity_cxx

#include "AliEventClassifierBase.h"

class AliEventClassifierSpherocity : public AliEventClassifierBase {
 public:
  AliEventClassifierSpherocity()
    : AliEventClassifierBase() {}
  AliEventClassifierSpherocity(const char* name, const char* title,
			TList *taskOutputList);
  virtual ~AliEventClassifierSpherocity() {}

 private:
  Bool_t TrackPassesSelection(AliMCParticle* track, AliStack *stack, Int_t iTrack);
  void CalculateClassifierValue(AliMCEvent *event, AliStack *stack);
  
  ClassDef(AliEventClassifierSpherocity, 1);
};

#endif
  
