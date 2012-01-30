////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelFreezeOutGenerator - abstract base class for freeze-out     ///
/// coordinates generator                                                    ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelFreezeOutGenerator, 1)
#endif

#include "AliFemtoModelFreezeOutGenerator.h"

//____________________________
AliFemtoModelFreezeOutGenerator::AliFemtoModelFreezeOutGenerator(): 
  fRandom(0) 
{ /* no-op */ }
//____________________________
AliFemtoModelFreezeOutGenerator::AliFemtoModelFreezeOutGenerator(const AliFemtoModelFreezeOutGenerator &/* aModel */): 
  fRandom(0)
{/* no-op */}
//____________________________
AliFemtoModelFreezeOutGenerator::~AliFemtoModelFreezeOutGenerator()
{
  if (fRandom) delete fRandom;
}
//____________________________
AliFemtoModelFreezeOutGenerator& AliFemtoModelFreezeOutGenerator::operator=(const AliFemtoModelFreezeOutGenerator& aGen) 
{ 
  if (this == &aGen) return *this; 
  if (aGen.fRandom) 
    fRandom = new TRandom2(*aGen.fRandom);
  else 
    fRandom=0; 
  return *this; 
}
//____________________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelFreezeOutGenerator::Clone() const 
{ 
  return 0; 
}



