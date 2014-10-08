#include "AliStackPartonInfo.h" 

ClassImp(AliStackPartonInfo)

AliStackPartonInfo::AliStackPartonInfo() :
  TNamed(),
  fPartonFlag6(0),
  fPartonPt6(0),
  fPartonEta6(0),
  fPartonPhi6(0),
  fPartonFlag7(0),
  fPartonPt7(0),
  fPartonEta7(0),
  fPartonPhi7(0)
{
  
}

//_______________________________________________

AliStackPartonInfo::AliStackPartonInfo(const char* name) :
  TNamed(),
  fPartonFlag6(0),
  fPartonPt6(0),
  fPartonEta6(0),
  fPartonPhi6(0),
  fPartonFlag7(0),
  fPartonPt7(0),
  fPartonEta7(0),
  fPartonPhi7(0)
{
  
}

