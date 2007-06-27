/* $Id$ */

#ifndef ALIPWG0DEPHELPER_H
#define ALIPWG0DEPHELPER_H

#include <TObject.h>

// static helper functions that depend on more than ESD

class AliHeader;
class TParticle;
class AliStack;

class AliPWG0depHelper : public TObject
{
  public:
    static Int_t GetPythiaEventProcessType(AliHeader* aHeader, Bool_t adebug = kFALSE);
    static TParticle* FindPrimaryMother(AliStack* stack, Int_t label);
    static Int_t FindPrimaryMotherLabel(AliStack* stack, Int_t label);

  protected:
    ClassDef(AliPWG0depHelper, 0)

  private:
    AliPWG0depHelper(const AliPWG0depHelper&);
    AliPWG0depHelper& operator=(const AliPWG0depHelper&);
};

#endif

