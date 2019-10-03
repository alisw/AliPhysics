///
/// \file AliFemtoModelGlobalHiddenInfo.cxx
///

#include "AliFemtoModelGlobalHiddenInfo.h"

//_____________________________________________
AliFemtoModelGlobalHiddenInfo::AliFemtoModelGlobalHiddenInfo():
  AliFemtoModelHiddenInfo(),
  fGlobalEmissionPoint(NULL)
{ // Default constructor
}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo::AliFemtoModelGlobalHiddenInfo(const AliFemtoModelGlobalHiddenInfo &aInfo):
  AliFemtoModelHiddenInfo(aInfo),
  fGlobalEmissionPoint(NULL)
{ // Copy constructor

  if (aInfo.GetGlobalEmissionPoint() != NULL) {
    fGlobalEmissionPoint = new AliFemtoThreeVector(*aInfo.GetGlobalEmissionPoint());
  }

}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo::~AliFemtoModelGlobalHiddenInfo()
{
  // Destructor
  delete fGlobalEmissionPoint;
}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo& AliFemtoModelGlobalHiddenInfo::operator=(const AliFemtoModelGlobalHiddenInfo& aInfo)
{
  // assignment operator
  if (this == &aInfo) {
    return *this;
  }

  AliFemtoModelHiddenInfo::operator=(aInfo);

  if (aInfo.GetGlobalEmissionPoint() != NULL) {
    SetGlobalEmissionPoint(*aInfo.GetGlobalEmissionPoint());
  } else {
    delete fGlobalEmissionPoint;
    fGlobalEmissionPoint = NULL;
  }

  return *this;
}
//_____________________________________________
AliFemtoThreeVector *AliFemtoModelGlobalHiddenInfo::GetGlobalEmissionPoint() const
{
  return fGlobalEmissionPoint;
}
//_____________________________________________
void AliFemtoModelGlobalHiddenInfo::SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos)
{
  // set position from vector
  if (fGlobalEmissionPoint) {
    *fGlobalEmissionPoint = aPos;
  } else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aPos);
  }
}
//_____________________________________________
void AliFemtoModelGlobalHiddenInfo::SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz)
{
  // Set position from components
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->SetX(aRx);
    fGlobalEmissionPoint->SetY(aRy);
    fGlobalEmissionPoint->SetZ(aRz);
  } else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aRx, aRy, aRz);
  }
}
//_____________________________________________
 AliFemtoHiddenInfo* AliFemtoModelGlobalHiddenInfo::GetParticleHiddenInfo() const
{
  // return copy of this hidden info
  AliFemtoModelGlobalHiddenInfo* tBuf = new AliFemtoModelGlobalHiddenInfo(*this);
  return tBuf;
}
