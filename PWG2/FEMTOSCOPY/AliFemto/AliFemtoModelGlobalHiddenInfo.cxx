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
//   if (fTrueMomentum) delete fTrueMomentum;
//   if (fEmissionPoint) delete fEmissionPoint;
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
  else SetEmissionPoint(0);
  if (aInfo.GetGlobalEmissionPoint())
    SetGlobalEmissionPoint(*aInfo.GetGlobalEmissionPoint());
  else fGlobalEmissionPoint = 0;
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
    fGlobalEmissionPoint->SetX(aPos.x());
    fGlobalEmissionPoint->SetY(aPos.y());
    fGlobalEmissionPoint->SetZ(aPos.z());
  }
  else {
    fGlobalEmissionPoint = new AliFemtoThreeVector();
    *fGlobalEmissionPoint = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelGlobalHiddenInfo::SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz)
{
  // Set position from components
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->SetX(aRx);
    fGlobalEmissionPoint->SetY(aRy);
    fGlobalEmissionPoint->SetZ(aRz);
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
