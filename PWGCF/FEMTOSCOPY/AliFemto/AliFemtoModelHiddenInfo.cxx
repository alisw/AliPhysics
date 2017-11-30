///
/// \file AliFemtoModelHiddenInfo.cxx
///

#include "AliFemtoModelHiddenInfo.h"

//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo():
  fTrueMomentum(NULL),
  fMotherMomentum(NULL),
  fEmissionPoint(NULL),
  fPDGPid(0),
  fMotherPdg(0),
  fMass(0.0),
  fTrueMomentumPos(NULL),
  fEmissionPointPos(NULL),
  fPDGPidPos(0),
  fMassPos(0.0),
  fTrueMomentumNeg(NULL),
  fEmissionPointNeg(NULL),
  fPDGPidNeg(0),
  fMassNeg(0.0),
  fOrigin(-900)
{
  // Default constructor
}
//_____________________________________________
AliFemtoModelHiddenInfo::AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo):
  AliFemtoHiddenInfo(aInfo),
  fTrueMomentum(NULL),
  fMotherMomentum(NULL),
  fEmissionPoint(NULL),
  fPDGPid(aInfo.fPDGPid),
  fMotherPdg(aInfo.fMotherPdg),
  fMass(aInfo.fMass),
  fTrueMomentumPos(NULL),
  fEmissionPointPos(NULL),
  fPDGPidPos(aInfo.fPDGPidPos),
  fMassPos(aInfo.fMassPos),
  fTrueMomentumNeg(NULL),
  fEmissionPointNeg(NULL),
  fPDGPidNeg(aInfo.fPDGPidNeg),
  fMassNeg(aInfo.fMassNeg),
  fOrigin(-900)
{
  // Copy constructor
  SetTrueMomentum(aInfo.GetTrueMomentum());
  SetMotherMomentum(aInfo.GetMotherMomentum());
  SetEmissionPoint(aInfo.GetEmissionPoint());
  SetTrueMomentumPos(aInfo.GetTrueMomentumPos());
  SetEmissionPointPos(aInfo.GetEmissionPointPos());
  SetTrueMomentumNeg(aInfo.GetTrueMomentumNeg());
  SetEmissionPointNeg(aInfo.GetEmissionPointNeg());
  SetOrigin(aInfo.GetOrigin());
}
//_____________________________________________
AliFemtoModelHiddenInfo::~AliFemtoModelHiddenInfo()
{
  // Destructor
  if (fTrueMomentum) delete fTrueMomentum;
  if (fMotherMomentum) delete fMotherMomentum;
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
  if (this == &aInfo) {
    return *this;
  }

  SetTrueMomentum(aInfo.GetTrueMomentum());
  SetMotherMomentum(aInfo.GetMotherMomentum());
  SetEmissionPoint(aInfo.GetEmissionPoint());
  fPDGPid = aInfo.GetPDGPid();
  fMotherPdg = aInfo.GetMotherPdgCode();
  fMass = aInfo.GetMass();

  SetTrueMomentumPos(aInfo.GetTrueMomentumPos());
  SetEmissionPointPos(aInfo.GetEmissionPointPos());
  fPDGPidPos = aInfo.GetPDGPidPos();
  fMassPos = aInfo.GetMassPos();

  SetTrueMomentumNeg(aInfo.GetTrueMomentumNeg());
  SetEmissionPointNeg(aInfo.GetEmissionPointNeg());
  fPDGPidNeg = aInfo.GetPDGPidNeg();
  fMassNeg = aInfo.GetMassNeg();

  fOrigin = aInfo.GetOrigin();

  return *this;
}
//_____________________________________________
AliFemtoThreeVector   *AliFemtoModelHiddenInfo::GetTrueMomentum() const
{
  return fTrueMomentum;
}
//_____________________________________________
AliFemtoThreeVector   *AliFemtoModelHiddenInfo::GetMotherMomentum() const
{
    return fMotherMomentum;
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

Int_t AliFemtoModelHiddenInfo::GetOrigin() const
{
  return fOrigin;
}

//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(AliFemtoThreeVector *aMom)
{
  // neither are null - simple set
  if (aMom != NULL && fTrueMomentum != NULL) {
    *fTrueMomentum = *aMom;
  }
  // parameter is NULL - delete & null member
  else if (fTrueMomentum != NULL) {
    delete fTrueMomentum;
    fTrueMomentum = NULL;
  }
  // member is NULL - allocate vector
  else if (aMom != NULL) {
    fTrueMomentum = new AliFemtoThreeVector(*aMom);
  }

}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    *fTrueMomentum = aMom;
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (fTrueMomentum == NULL) {
     fTrueMomentum = new AliFemtoThreeVector(aPx, aPy, aPz);
  }
  else {
    fTrueMomentum->SetX(aPx);
    fTrueMomentum->SetY(aPy);
    fTrueMomentum->SetZ(aPz);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMotherMomentum(AliFemtoThreeVector *aMom)
{
    // neither are null - simple set
    if (aMom != NULL && fMotherMomentum != NULL) {
        *fMotherMomentum = *aMom;
    }
    // parameter is NULL - delete & null member
    else if (fMotherMomentum != NULL) {
        delete fMotherMomentum;
        fMotherMomentum = NULL;
    }
    // member is NULL - allocate vector
    else if (aMom != NULL) {
        fMotherMomentum = new AliFemtoThreeVector(*aMom);
    }
    
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMotherMomentum(const AliFemtoThreeVector& aMom)
{
    // Set momentum from vector
    if (fMotherMomentum) {
        *fMotherMomentum = aMom;
    }
    else {
        fMotherMomentum = new AliFemtoThreeVector(aMom);
    }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetMotherMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
    // Set momentum from components
    if (fMotherMomentum == NULL) {
        fMotherMomentum = new AliFemtoThreeVector(aPx, aPy, aPz);
    }
    else {
        fMotherMomentum->SetX(aPx);
        fMotherMomentum->SetY(aPy);
        fMotherMomentum->SetZ(aPz);
    }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  // neither are null - simple set
  if (aPos != NULL && fEmissionPoint != NULL) {
    *fEmissionPoint = *aPos;
  }
  // parameter is NULL - delete & null member
  else if (fEmissionPoint != NULL) {
    delete fEmissionPoint;
    fEmissionPoint = NULL;
  }
  // member is NULL - allocate vector
  else if (aPos != NULL) {
    fEmissionPoint = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPoint(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    *fEmissionPoint = aPos;
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aPos);
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
  // neither are null - simple set
  if (aMom != NULL && fTrueMomentumPos != NULL) {
    *fTrueMomentumPos = *aMom;
  }
  // parameter is NULL - delete & null member
  else if (fTrueMomentumPos != NULL) {
    delete fTrueMomentumPos;
    fTrueMomentumPos = NULL;
  }
  // member is NULL - allocate vector
  else if (aMom != NULL) {
    fTrueMomentumPos = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumPos(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentumPos) {
    *fTrueMomentumPos = aMom;
  }
  else {
    fTrueMomentumPos = new AliFemtoThreeVector(aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumPos(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (fTrueMomentumPos == NULL) {
    fTrueMomentumPos = new AliFemtoThreeVector(aPx, aPy, aPz);
  } else {
    fTrueMomentumPos->SetX(aPx);
    fTrueMomentumPos->SetY(aPy);
    fTrueMomentumPos->SetZ(aPz);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointPos(AliFemtoLorentzVector *aPos)
{
  // neither are null - simple set
  if (aPos != NULL && fEmissionPointNeg != NULL) {
    *fEmissionPointPos = *aPos;
  }
  // parameter is NULL - delete & null member
  else if (fEmissionPointPos != NULL) {
    delete fEmissionPointPos;
    fEmissionPointPos = NULL;
  }
  // member is NULL - allocate vector
  else if (aPos != NULL) {
    fEmissionPointPos = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointPos(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPointPos) {
    *fEmissionPointPos = aPos;
  }
  else {
    fEmissionPointPos = new AliFemtoLorentzVector(aPos);
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
  // neither are null - simple set
  if (aMom != NULL && fTrueMomentumNeg != NULL) {
    *fTrueMomentumNeg = *aMom;
  }
  // parameter is NULL - delete & null member
  else if (fTrueMomentumNeg != NULL) {
    delete fTrueMomentumNeg;
    fTrueMomentumNeg = NULL;
  }
  // member is NULL - allocate vector
  else if (aMom != NULL) {
    fTrueMomentumNeg = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumNeg(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentumNeg) {
    *fTrueMomentumNeg = aMom;
  }
  else {
    fTrueMomentumNeg = new AliFemtoThreeVector(aMom);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetTrueMomentumNeg(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (fTrueMomentumNeg == NULL) {
    fTrueMomentumNeg = new AliFemtoThreeVector(aPx, aPy, aPz);
  } else {
    fTrueMomentumNeg->SetX(aPx);
    fTrueMomentumNeg->SetY(aPy);
    fTrueMomentumNeg->SetZ(aPz);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointNeg(AliFemtoLorentzVector *aPos)
{
  // neither are null - simple set
  if (aPos != NULL && fEmissionPointNeg != NULL) {
    *fEmissionPointNeg = *aPos;
  }
  // parameter is NULL - delete & null member
  else if (fEmissionPointNeg != NULL) {
    delete fEmissionPointNeg;
    fEmissionPointNeg = NULL;
  }
  // member is NULL - allocate vector
  else if (aPos != NULL) {
    fEmissionPointNeg = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelHiddenInfo::SetEmissionPointNeg(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPointNeg) {
    *fEmissionPointNeg = aPos;
  }
  else {
    fEmissionPointNeg = new AliFemtoLorentzVector(aPos);
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


void AliFemtoModelHiddenInfo::SetOrigin(Int_t origin)
{
  fOrigin = origin;
}
