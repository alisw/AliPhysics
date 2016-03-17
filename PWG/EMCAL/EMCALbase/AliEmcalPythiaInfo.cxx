#include "AliEmcalPythiaInfo.h"

/// \cond CLASSIMP
ClassImp(AliEmcalPythiaInfo)
/// \endcond

//_______________________________________________
AliEmcalPythiaInfo::AliEmcalPythiaInfo() :
  TNamed("AliEmcalPythiaInfo", "AliEmcalPythiaInfo"),
  fPartonFlag6(0),
  fParton6(),
  fPartonFlag7(0),
  fParton7(),
  fPythiaEventWeight(1)
{
}

//_______________________________________________
AliEmcalPythiaInfo::AliEmcalPythiaInfo(const char* name) :
  TNamed(name, name),
  fPartonFlag6(0),
  fParton6(),
  fPartonFlag7(0),
  fParton7(),
  fPythiaEventWeight(1)
{
}

