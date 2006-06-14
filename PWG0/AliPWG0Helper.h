/* $Id$ */

#ifndef ALIPWG0HELPER_H
#define ALIPWG0HELPER_H

#include <TObject.h>

// static helper functions

class AliESD;
class TParticle;
class TH3F;

class AliPWG0Helper : public TObject
{
  public:
    static Bool_t IsEventTriggered(AliESD* aEsd);
    static Bool_t IsVertexReconstructed(AliESD* aEsd);
    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries);

    static void CreateProjections(TH3F* hist);
    static void CreateDividedProjections(TH3F* hist, TH3F* hist2, const char* axis = 0);

  protected:
    ClassDef(AliPWG0Helper, 0)
};

#endif

