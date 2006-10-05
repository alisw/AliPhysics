/* $Id$ */

#ifndef ALIPWG0HELPER_H
#define ALIPWG0HELPER_H

#include <TObject.h>

#include <AliHeader.h>

// static helper functions

class AliESD;
class TParticle;
class TH3;

class AliPWG0Helper : public TObject
{
  public:
    static Bool_t IsEventTriggered(AliESD* aEsd);
    static Bool_t IsVertexReconstructed(AliESD* aEsd);
    static Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries, Bool_t adebug = kFALSE);

    static void CreateProjections(TH3* hist);
    static void CreateDividedProjections(TH3* hist, TH3* hist2, const char* axis = 0, Bool_t putErrors = kFALSE);
    static const char* GetAxisTitle(TH3* hist, const char axis);
    
  protected:
    ClassDef(AliPWG0Helper, 0)

  private:
    AliPWG0Helper(const AliPWG0Helper&);
    AliPWG0Helper& operator=(const AliPWG0Helper&);
};

#endif

