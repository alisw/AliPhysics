////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelFreezeOutGenerator - abstract base class for freeze-out     ///
/// coordinates generator                                                    ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelFreezeOutGenerator_hh
#define AliFemtoModelFreezeOutGenerator_hh

#include "TRandom2.h"
#include "AliFemtoPair.h"

class AliFemtoModelFreezeOutGenerator 
{
 public:
  AliFemtoModelFreezeOutGenerator();
  AliFemtoModelFreezeOutGenerator(AliFemtoModelFreezeOutGenerator &aModel);
  
  AliFemtoModelFreezeOutGenerator& operator=(const AliFemtoModelFreezeOutGenerator& aGen);
  
  virtual ~AliFemtoModelFreezeOutGenerator();
  virtual void GenerateFreezeOut(AliFemtoPair *aPair) = 0;
  
  virtual AliFemtoModelFreezeOutGenerator* Clone() const;
  
 protected:
  TRandom2 *fRandom;
  
 private:
  
#ifdef __ROOT__
  ClassDef(AliFemtoModelFreezeOutGenerator, 1)
#endif
    
};

inline AliFemtoModelFreezeOutGenerator::AliFemtoModelFreezeOutGenerator(): fRandom(0) { /* no-op */ };
inline AliFemtoModelFreezeOutGenerator::AliFemtoModelFreezeOutGenerator(AliFemtoModelFreezeOutGenerator &aModel): fRandom(0){/* no-op */};
inline AliFemtoModelFreezeOutGenerator& AliFemtoModelFreezeOutGenerator::operator=(const AliFemtoModelFreezeOutGenerator& aGen) 
{ 
  if (this == &aGen) return *this; 
  if (aGen.fRandom) 
    fRandom = new TRandom2(*aGen.fRandom);
  else 
    fRandom=0; 
  return *this; 
};

inline AliFemtoModelFreezeOutGenerator::AliFemtoModelFreezeOutGenerator* Clone() { return 0; };


#endif


