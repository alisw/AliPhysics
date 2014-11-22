#include "AliPythiaInfo.h" 

ClassImp(AliPythiaInfo)

//_______________________________________________
AliPythiaInfo::AliPythiaInfo() :
TNamed("AliPythiaInfo","AliPythiaInfo"),
  fPartonFlag6(0),
  fPartonPt6(0),
  fPartonEta6(0),
  fPartonPhi6(0),
  fPartonFlag7(0),
  fPartonPt7(0),
  fPartonEta7(0),
  fPartonPhi7(0),
  fPythiaEventWeight(1)
{
  
}

//_______________________________________________
AliPythiaInfo::AliPythiaInfo(const char* name) :
  TNamed(name,name),
  fPartonFlag6(0),
  fPartonPt6(0),
  fPartonEta6(0),
  fPartonPhi6(0),
  fPartonFlag7(0),
  fPartonPt7(0),
  fPartonEta7(0),
  fPartonPhi7(0),
  fPythiaEventWeight(1)

{
  
}

