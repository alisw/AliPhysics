////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGlobalHiddenInfo - the hidden info for model calculations    //
/// Stores information needed for the weight generation - the true           ///
/// simulated momenta, freeze-out coordinates from model and particle PID    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelGlobalHiddenInfo.h"

//_____________________________________________
AliFemtoModelGlobalHiddenInfo::AliFemtoModelGlobalHiddenInfo() :
  AliFemtoModelHiddenInfo(),
  fGlobalEmissionPoint(0)
{
  // Default constructor
}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo::AliFemtoModelGlobalHiddenInfo(const AliFemtoModelGlobalHiddenInfo &aInfo) :
  AliFemtoModelHiddenInfo(aInfo),
  fGlobalEmissionPoint(0)
{
  // Copy constructor
  if (aInfo.GetGlobalEmissionPoint())
    SetGlobalEmissionPoint((*aInfo.GetGlobalEmissionPoint()));
}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo::~AliFemtoModelGlobalHiddenInfo()
{
  // Destructor
  if (fTrueMomentum) delete fTrueMomentum;
  if (fEmissionPoint) delete fEmissionPoint;
  if (fGlobalEmissionPoint) delete fGlobalEmissionPoint;
}
//_____________________________________________
AliFemtoModelGlobalHiddenInfo& AliFemtoModelGlobalHiddenInfo::operator=(const AliFemtoModelGlobalHiddenInfo& aInfo)
{
  // assignment operator
  if (this == &aInfo)
    return *this;

  if (aInfo.GetTrueMomentum())
    SetTrueMomentum(aInfo.GetTrueMomentum());
  else SetTrueMomentum(0);
  if (aInfo.GetEmissionPoint())
    SetEmissionPoint(aInfo.GetEmissionPoint());
  if (aInfo.GetGlobalEmissionPoint())
    SetGlobalEmissionPoint(*aInfo.GetGlobalEmissionPoint());
  else SetEmissionPoint(0);
  fPDGPid = aInfo.GetPDGPid();
  fMass = aInfo.GetMass();

  return *this;
}
//_____________________________________________
AliFemtoThreeVector *AliFemtoModelGlobalHiddenInfo::GetGlobalEmissionPoint() const
{
  return fGlobalEmissionPoint;
}
//_____________________________________________
void                   AliFemtoModelGlobalHiddenInfo::SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos)
{
  // set position from vector
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->setX(aPos.x());
    fGlobalEmissionPoint->setY(aPos.y());
    fGlobalEmissionPoint->setZ(aPos.z());
  }
  else {
    fGlobalEmissionPoint = new AliFemtoThreeVector();
    *fGlobalEmissionPoint = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelGlobalHiddenInfo::SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz)
{
  // set position from components
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->setX(aRx);
    fGlobalEmissionPoint->setY(aRy);
    fGlobalEmissionPoint->setZ(aRz);
  }
  else {
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
