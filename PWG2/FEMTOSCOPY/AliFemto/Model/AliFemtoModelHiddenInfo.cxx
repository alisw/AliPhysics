#include "AliFemtoModelHiddenInfo.h"

//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo() :
  fTrueMomentum(0),
  fEmissionPoint(0),
  fPDGPid(0),
  fMass(0)
{
};
//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo) :
  fTrueMomentum(0),
  fEmissionPoint(0),
  fPDGPid(0),
  fMass(0)
{
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
  if (fTrueMomentum) delete fTrueMomentum;
  if (fEmissionPoint) delete fEmissionPoint;
}
//_____________________________________________
AliFemtoModelHiddenInfo& AliFemtoModelHiddenInfo::operator=(const AliFemtoModelHiddenInfo& aInfo)
{
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
  if (fTrueMomentum) {
    fTrueMomentum->setX(aMom->x());
    fTrueMomentum->setY(aMom->y());
    fTrueMomentum->setZ(aMom->z());
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(const AliFemtoThreeVector& aMom)
{
  if (fTrueMomentum) {
    fTrueMomentum->setX(aMom.x());
    fTrueMomentum->setY(aMom.y());
    fTrueMomentum->setZ(aMom.z());
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector();
    *fTrueMomentum = aMom;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  if (!fTrueMomentum) fTrueMomentum = new AliFemtoThreeVector();
    fTrueMomentum->setX(aPx);
    fTrueMomentum->setY(aPy);
    fTrueMomentum->setZ(aPz);
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  if (fEmissionPoint) {
    fEmissionPoint->setX(aPos->px());
    fEmissionPoint->setY(aPos->py());
    fEmissionPoint->setZ(aPos->pz());
    fEmissionPoint->setT(aPos->e());
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(const AliFemtoLorentzVector& aPos)
{
  if (fEmissionPoint) {
    fEmissionPoint->setX(aPos.px());
    fEmissionPoint->setY(aPos.py());
    fEmissionPoint->setZ(aPos.pz());
    fEmissionPoint->setT(aPos.e());
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
  fEmissionPoint->setX(aRx);
  fEmissionPoint->setY(aRy);
  fEmissionPoint->setZ(aRz);
  fEmissionPoint->setT(aT);
}
//_____________________________________________
 AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::getParticleHiddenInfo() const
{
  AliFemtoModelHiddenInfo* tBuf = new AliFemtoModelHiddenInfo(*this);
  return tBuf;
}
