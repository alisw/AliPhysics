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
  fMotherPdg(0),
  fMass(0),
  fTrueMomentumPos(0),
  fEmissionPointPos(0),
  fPDGPidPos(0),
  fMassPos(0),
  fTrueMomentumNeg(0),
  fEmissionPointNeg(0),
  fPDGPidNeg(0),
  fMassNeg(0)
{
  // Default constructor
}
//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo) :
  AliFemtoHiddenInfo(aInfo),
  fTrueMomentum(0),
  fEmissionPoint(0),
  fPDGPid(0),
  fMotherPdg(0),
  fMass(0),
  fTrueMomentumPos(0),
  fEmissionPointPos(0),
  fPDGPidPos(0),
  fMassPos(0),
  fTrueMomentumNeg(0),
  fEmissionPointNeg(0),
  fPDGPidNeg(0),
  fMassNeg(0)
{
  // Copy constructor
  if (aInfo.GetTrueMomentum())
    SetTrueMomentum(aInfo.GetTrueMomentum());
  if (aInfo.GetEmissionPoint())
    SetEmissionPoint(aInfo.GetEmissionPoint());
  fPDGPid = aInfo.GetPDGPid();
  fMotherPdg = aInfo.GetMotherPdgCode();
  fMass = aInfo.GetMass();

  if (aInfo.GetTrueMomentumPos())
    SetTrueMomentumPos(aInfo.GetTrueMomentumPos());
  if (aInfo.GetEmissionPointPos())
    SetEmissionPointPos(aInfo.GetEmissionPointPos());
  fPDGPidPos = aInfo.GetPDGPidPos();
  fMassPos = aInfo.GetMassPos();

  if (aInfo.GetTrueMomentumNeg())
    SetTrueMomentumNeg(aInfo.GetTrueMomentumNeg());
  if (aInfo.GetEmissionPointNeg())
    SetEmissionPointNeg(aInfo.GetEmissionPointNeg());
  fPDGPidNeg = aInfo.GetPDGPidNeg();
  fMassNeg = aInfo.GetMassNeg();
}
//_____________________________________________
AliFemtoModelHiddenInfo::~AliFemtoModelHiddenInfo()
{
  // Destructor
  if (fTrueMomentum) delete fTrueMomentum;
  if (fEmissionPoint) delete fEmissionPoint;
  if (fTrueMomentumPos) delete fTrueMomentumPos;
  if (fEmissionPointPos) delete fEmissionPointPos;
  if (fTrueMomentumNeg) delete fTrueMomentumNeg;
  if (fEmissionPointNeg) delete fEmissionPointNeg;
}
//_____________________________________________
AliFemtoModelHiddenInfo& AliFemtoModelHiddenInfo::operator=(const AliFemtoModelHiddenInfo& aInfo)
{
  // assignment operator
  if (this == &aInfo)
    return *this;

  if (fTrueMomentum) delete fTrueMomentum;
  if (aInfo.GetTrueMomentum())
    SetTrueMomentum(aInfo.GetTrueMomentum());
  else SetTrueMomentum(0);
  if (fEmissionPoint) delete fEmissionPoint;
  if (aInfo.GetEmissionPoint())
    SetEmissionPoint(aInfo.GetEmissionPoint());
  else SetEmissionPoint(0);
  fPDGPid = aInfo.GetPDGPid();
  fMotherPdg = aInfo.GetMotherPdgCode();
  fMass = aInfo.GetMass();

  if (fTrueMomentumPos) delete fTrueMomentumPos;
  if (aInfo.GetTrueMomentumPos())
    SetTrueMomentumPos(aInfo.GetTrueMomentumPos());
  else SetTrueMomentumPos(0);
  if (fEmissionPointPos) delete fEmissionPointPos;
  if (aInfo.GetEmissionPointPos())
    SetEmissionPointPos(aInfo.GetEmissionPointPos());
  else SetEmissionPointPos(0);
  fPDGPidPos = aInfo.GetPDGPidPos();
  fMassPos = aInfo.GetMassPos();

  if (fTrueMomentumNeg) delete fTrueMomentumNeg;
  if (aInfo.GetTrueMomentumNeg())
    SetTrueMomentumNeg(aInfo.GetTrueMomentumNeg());
  else SetTrueMomentumNeg(0);
  if (fEmissionPointNeg) delete fEmissionPointNeg;
  if (aInfo.GetEmissionPointNeg())
    SetEmissionPointNeg(aInfo.GetEmissionPointNeg());
  else SetEmissionPointNeg(0);
  fPDGPidNeg = aInfo.GetPDGPidNeg();
  fMassNeg = aInfo.GetMassNeg();

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
Int_t                  AliFemtoModelHiddenInfo::GetMotherPdgCode() const
{
  return fMotherPdg;
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
void                   AliFemtoModelHiddenInfo::SetMotherPdgCode(Int_t aMotherPdg)
{
  fMotherPdg = aMotherPdg;
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
AliFemtoThreeVector   *AliFemtoModelHiddenInfo::GetTrueMomentumPos() const
{
  return fTrueMomentumPos;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelHiddenInfo::GetEmissionPointPos() const
{
  return fEmissionPointPos;
}
//_____________________________________________
Int_t                  AliFemtoModelHiddenInfo::GetPDGPidPos() const
{
  return fPDGPidPos;
}
//_____________________________________________
Double_t                  AliFemtoModelHiddenInfo::GetMassPos() const
{
  return fMassPos;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumPos(AliFemtoThreeVector *aMom)
{
  // Set momentum from vector
  if (fTrueMomentumPos) {
    fTrueMomentumPos->SetX(aMom->x());
    fTrueMomentumPos->SetY(aMom->y());
    fTrueMomentumPos->SetZ(aMom->z());
  }
  else {
    fTrueMomentumPos = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumPos(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentumPos) {
    fTrueMomentumPos->SetX(aMom.x());
    fTrueMomentumPos->SetY(aMom.y());
    fTrueMomentumPos->SetZ(aMom.z());
  }
  else {
    fTrueMomentumPos = new AliFemtoThreeVector();
    *fTrueMomentumPos = aMom;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumPos(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (!fTrueMomentumPos) fTrueMomentumPos = new AliFemtoThreeVector();
    fTrueMomentumPos->SetX(aPx);
    fTrueMomentumPos->SetY(aPy);
    fTrueMomentumPos->SetZ(aPz);
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointPos(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPointPos) {
    fEmissionPointPos->SetX(aPos->px());
    fEmissionPointPos->SetY(aPos->py());
    fEmissionPointPos->SetZ(aPos->pz());
    fEmissionPointPos->SetT(aPos->e());
  }
  else {
    fEmissionPointPos = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointPos(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPointPos) {
    fEmissionPointPos->SetX(aPos.px());
    fEmissionPointPos->SetY(aPos.py());
    fEmissionPointPos->SetZ(aPos.pz());
    fEmissionPointPos->SetT(aPos.e());
  }
  else {
    fEmissionPointPos = new AliFemtoLorentzVector();
    *fEmissionPointPos = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetPDGPidPos(Int_t aPid)
{
  fPDGPidPos = aPid;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMassPos(Double_t aMass)
{
  fMassPos = aMass;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointPos(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // Set position from components
  if (fEmissionPointPos) {
    fEmissionPointPos->SetX(aRx);
    fEmissionPointPos->SetY(aRy);
    fEmissionPointPos->SetZ(aRz);
    fEmissionPointPos->SetT(aT);
  }
  else {
    fEmissionPointPos = new AliFemtoLorentzVector(aRx, aRy, aRz, aT);
  }
}

//_____________________________________________
AliFemtoThreeVector   *AliFemtoModelHiddenInfo::GetTrueMomentumNeg() const
{
  return fTrueMomentumNeg;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelHiddenInfo::GetEmissionPointNeg() const
{
  return fEmissionPointNeg;
}
//_____________________________________________
Int_t                  AliFemtoModelHiddenInfo::GetPDGPidNeg() const
{
  return fPDGPidNeg;
}
//_____________________________________________
Double_t                  AliFemtoModelHiddenInfo::GetMassNeg() const
{
  return fMassNeg;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumNeg(AliFemtoThreeVector *aMom)
{
  // Set momentum from vector
  if (fTrueMomentumNeg) {
    fTrueMomentumNeg->SetX(aMom->x());
    fTrueMomentumNeg->SetY(aMom->y());
    fTrueMomentumNeg->SetZ(aMom->z());
  }
  else {
    fTrueMomentumNeg = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumNeg(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentumNeg) {
    fTrueMomentumNeg->SetX(aMom.x());
    fTrueMomentumNeg->SetY(aMom.y());
    fTrueMomentumNeg->SetZ(aMom.z());
  }
  else {
    fTrueMomentumNeg = new AliFemtoThreeVector();
    *fTrueMomentumNeg = aMom;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumNeg(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (!fTrueMomentumNeg) fTrueMomentumNeg = new AliFemtoThreeVector();
    fTrueMomentumNeg->SetX(aPx);
    fTrueMomentumNeg->SetY(aPy);
    fTrueMomentumNeg->SetZ(aPz);
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointNeg(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPointNeg) {
    fEmissionPointNeg->SetX(aPos->px());
    fEmissionPointNeg->SetY(aPos->py());
    fEmissionPointNeg->SetZ(aPos->pz());
    fEmissionPointNeg->SetT(aPos->e());
  }
  else {
    fEmissionPointNeg = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointNeg(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPointNeg) {
    fEmissionPointNeg->SetX(aPos.px());
    fEmissionPointNeg->SetY(aPos.py());
    fEmissionPointNeg->SetZ(aPos.pz());
    fEmissionPointNeg->SetT(aPos.e());
  }
  else {
    fEmissionPointNeg = new AliFemtoLorentzVector();
    *fEmissionPointNeg = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetPDGPidNeg(Int_t aPid)
{
  fPDGPidNeg = aPid;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMassNeg(Double_t aMass)
{
  fMassNeg = aMass;
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointNeg(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // Set position from components
  if (fEmissionPointNeg) {
    fEmissionPointNeg->SetX(aRx);
    fEmissionPointNeg->SetY(aRy);
    fEmissionPointNeg->SetZ(aRz);
    fEmissionPointNeg->SetT(aT);
  }
  else {
    fEmissionPointNeg = new AliFemtoLorentzVector(aRx, aRy, aRz, aT);
  }
}


//_____________________________________________
 AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::GetParticleHiddenInfo() const
{
  // return copy of this hidden info
  AliFemtoModelHiddenInfo* tBuf = new AliFemtoModelHiddenInfo(*this);
  return tBuf;
}
