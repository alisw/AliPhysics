////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelHiddenInfo - the hidden info for model calculations         ///
/// Stores information needed for the weight generation - the true           ///
/// simulated momenta, freeze-out coordinates from model and particle PID    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoModelHiddenInfo.h"

//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo() :
  fTrueMomentum(0),
  fEmissionPoint(0),
  fPDGPid(0),
  fMass(0)
{
  // Default constructor
}
//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo) :
  AliFemtoHiddenInfo(aInfo),
  fTrueMomentum(0),
  fEmissionPoint(0),
  fPDGPid(0),
  fMass(0)
{
  // Copy constructor
  if (aInfo.GetTrueMomentum())
    SetTrueMomentum(aInfo.GetTrueMomentum());
  if (aInfo.GetEmissionPoint())
    SetEmissionPoint(aInfo.GetEmissionPoint());
  fPDGPid = aInfo.GetPDGPid();
  fMass = aInfo.GetMass();
}
//_____________________________________________
AliFemtoModelHiddenInfo::~AliFemtoModelHiddenInfo()
{
  // Destructor
  if (fTrueMomentum) delete fTrueMomentum;
  if (fEmissionPoint) delete fEmissionPoint;
}
//_____________________________________________
AliFemtoModelHiddenInfo& AliFemtoModelHiddenInfo::operator=(const AliFemtoModelHiddenInfo& aInfo)
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
  fPDGPid = aInfo.GetPDGPid();
  fMass = aInfo.GetMass();

  return *this;
}
//_____________________________________________
AliFemtoThreeVector   *AliFemtoModelHiddenInfo::GetTrueMomentum() const
{
  return fTrueMomentum;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelHiddenInfo::GetEmissionPoint() const
{
  return fEmissionPoint;
}
//_____________________________________________
Int_t                  AliFemtoModelHiddenInfo::GetPDGPid() const
{
  return fPDGPid;
}
//_____________________________________________
Double_t                  AliFemtoModelHiddenInfo::GetMass() const
{
  return fMass;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(AliFemtoThreeVector *aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    fTrueMomentum->SetX(aMom->x());
    fTrueMomentum->SetY(aMom->y());
    fTrueMomentum->SetZ(aMom->z());
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    fTrueMomentum->SetX(aMom.x());
    fTrueMomentum->SetY(aMom.y());
    fTrueMomentum->SetZ(aMom.z());
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector();
    *fTrueMomentum = aMom;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (!fTrueMomentum) fTrueMomentum = new AliFemtoThreeVector();
    fTrueMomentum->SetX(aPx);
    fTrueMomentum->SetY(aPy);
    fTrueMomentum->SetZ(aPz);
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aPos->px());
    fEmissionPoint->SetY(aPos->py());
    fEmissionPoint->SetZ(aPos->pz());
    fEmissionPoint->SetT(aPos->e());
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aPos.px());
    fEmissionPoint->SetY(aPos.py());
    fEmissionPoint->SetZ(aPos.pz());
    fEmissionPoint->SetT(aPos.e());
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector();
    *fEmissionPoint = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetPDGPid(Int_t aPid)
{
  fPDGPid = aPid;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMass(Double_t aMass)
{
  fMass = aMass;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // Set position from components
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aRx);
    fEmissionPoint->SetY(aRy);
    fEmissionPoint->SetZ(aRz);
    fEmissionPoint->SetT(aT);
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aRx, aRy, aRz, aT); 
  }
}
//_____________________________________________
 AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::GetParticleHiddenInfo() const
{
  // return copy of this hidden info
  AliFemtoModelHiddenInfo* tBuf = new AliFemtoModelHiddenInfo(*this);
  return tBuf;
}
