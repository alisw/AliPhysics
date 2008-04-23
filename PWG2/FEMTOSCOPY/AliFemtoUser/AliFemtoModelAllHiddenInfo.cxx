//K//////////////////////////////////////////////////////////////////////////M//
//K                                                                          M//
//K AliFemtoModelAllHiddenInfo -                                             M//
//K derived class inherits  the base class AliFemtoModelHiddenInfo           M//
//K the hidden info for model calculations                                   M//
//K Stores information needed for the weight generation -                    M//
//K                                                                          M//
//K in addition to  the base class AliFemtoModelHiddenInfo - the true        M//
//K simulated momenta, freeze-out coordinates from model and particle PID    M//
//K New information was added                                                M//
//K 1. Mother ID                                                             M//
//K 2. Mother 4-Momentum                                                     M//
//K 3. Mother emission point 4-vector                                        M//
//K 4. Childs IDs                                                            M//
//K 5. Childs 4-Momentum                                                     M//
//K--------------------------------------------------------------------------M//                                                                            //
//K APR2008  Konstantin Mikhailov Konstantin.Mikhailov@itep.ru               M//
//K                                                                          M//
//K//////////////////////////////////////////////////////////////////////////M//
#include "AliFemtoModelAllHiddenInfo.h"

//_____________________________________________
AliFemtoModelAllHiddenInfo::AliFemtoModelAllHiddenInfo() :
  fTrueMomentumMother(0),
  fEmissionPointMother(0),
  fPDGPidMother(0),
  fTrueMomentumChild1(0),
  fTrueMomentumChild2(0),
  fPDGPidChild1(0),
  fPDGPidChild2(0)
{
  // Default constructor
}
//_____________________________________________
AliFemtoModelAllHiddenInfo::AliFemtoModelAllHiddenInfo(const AliFemtoModelAllHiddenInfo &aInfo) :
  AliFemtoModelHiddenInfo(aInfo),
  fTrueMomentumMother(0),
  fEmissionPointMother(0),
  fPDGPidMother(0),
  fTrueMomentumChild1(0),
  fTrueMomentumChild2(0),
  fPDGPidChild1(0),
  fPDGPidChild2(0)
{
  // Copy constructor
  if (aInfo.GetTrueMomentumMother())
    SetTrueMomentumMother(aInfo.GetTrueMomentumMother());
  if (aInfo.GetEmissionPointMother())
    SetEmissionPointMother(aInfo.GetEmissionPointMother());
  fPDGPidMother = aInfo.GetPDGPidMother();
  if (aInfo.GetTrueMomentumChild1())
    SetTrueMomentumChild1(aInfo.GetTrueMomentumChild1());
  if (aInfo.GetTrueMomentumChild2())
    SetTrueMomentumChild2(aInfo.GetTrueMomentumChild2());
  fPDGPidChild1 = aInfo.GetPDGPidChild1();
  fPDGPidChild2 = aInfo.GetPDGPidChild2();
}
//_____________________________________________
AliFemtoModelAllHiddenInfo::~AliFemtoModelAllHiddenInfo()
{
  // Destructor
  if (fTrueMomentumMother)  delete fTrueMomentumMother;
  if (fEmissionPointMother) delete fEmissionPointMother;
  if (fTrueMomentumChild1) delete fTrueMomentumChild1;
  if (fTrueMomentumChild2) delete fTrueMomentumChild2;
}
//_____________________________________________
AliFemtoModelAllHiddenInfo& AliFemtoModelAllHiddenInfo::operator=(const AliFemtoModelAllHiddenInfo& aInfo)
{
  // assignment operator
  if (this == &aInfo)
    return *this;

  if (aInfo.GetTrueMomentumMother())
    SetTrueMomentumMother(aInfo.GetTrueMomentumMother());
  else SetTrueMomentumMother(0);
  if (aInfo.GetEmissionPointMother())
    SetEmissionPointMother(aInfo.GetEmissionPointMother());
  else SetEmissionPointMother(0);
  fPDGPidMother = aInfo.GetPDGPidMother();
  if (aInfo.GetTrueMomentumChild1())
    SetTrueMomentumChild1(aInfo.GetTrueMomentumChild1());
  else SetTrueMomentumChild1(0);
  if (aInfo.GetTrueMomentumChild2())
    SetTrueMomentumChild2(aInfo.GetTrueMomentumChild2());
  else SetTrueMomentumChild2(0);
  fPDGPidChild1 = aInfo.GetPDGPidChild1();
  fPDGPidChild2 = aInfo.GetPDGPidChild2();

  return *this;
}
//
//   GET
//
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelAllHiddenInfo::GetTrueMomentumMother() const
{
return fTrueMomentumMother;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelAllHiddenInfo::GetEmissionPointMother() const
{
  return fEmissionPoint;
}
//_____________________________________________
  Int_t                AliFemtoModelAllHiddenInfo::GetPDGPidMother() const
{
  return fPDGPidMother;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelAllHiddenInfo::GetTrueMomentumChild1() const
{
return fTrueMomentumChild1;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoModelAllHiddenInfo::GetTrueMomentumChild2() const
{
return fTrueMomentumChild2;
}
//_____________________________________________
  Int_t                AliFemtoModelAllHiddenInfo::GetPDGPidChild1() const
{
  return fPDGPidChild1;
}
//_____________________________________________
  Int_t                AliFemtoModelAllHiddenInfo::GetPDGPidChild2() const
{
  return fPDGPidChild2;
}
//
//   SET
//
//_____________________________________________
//  Mother momentum
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumMother(AliFemtoLorentzVector *aMomMother)
{
  // set momentum from vector
  if (fTrueMomentumMother) {
    fTrueMomentumMother->setX(aMomMother->px());
    fTrueMomentumMother->setY(aMomMother->py());
    fTrueMomentumMother->setZ(aMomMother->pz());
    fTrueMomentumMother->setT(aMomMother->e());
  }
  else {
    fTrueMomentumMother = new AliFemtoLorentzVector(*aMomMother);
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumMother(const AliFemtoLorentzVector& aMomMother)
{
  // set momentum from vector and energy
  if (fTrueMomentumMother) {
    fTrueMomentumMother->setX(aMomMother.px());
    fTrueMomentumMother->setY(aMomMother.py());
    fTrueMomentumMother->setZ(aMomMother.pz());
    fTrueMomentumMother->setT(aMomMother.e());
  }
  else {
    fTrueMomentumMother = new AliFemtoLorentzVector();
    *fTrueMomentumMother = aMomMother;
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumMother(Double_t aMotherPx, Double_t aMotherPy, Double_t aMotherPz, Double_t aMotherE)
{
  // set momentum from components and energy
  if (!fTrueMomentumMother) fTrueMomentumMother = new AliFemtoLorentzVector();
    fTrueMomentumMother->setX(aMotherPx);
    fTrueMomentumMother->setY(aMotherPy);
    fTrueMomentumMother->setZ(aMotherPz);
    fTrueMomentumMother->setT(aMotherE);
}
//_____________________________________________
//   Mother Emissin Point
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetEmissionPointMother(AliFemtoLorentzVector *aPos)
{
  // set position from vector
  if (fEmissionPointMother) {
    fEmissionPointMother->setX(aPos->px());
    fEmissionPointMother->setY(aPos->py());
    fEmissionPointMother->setZ(aPos->pz());
    fEmissionPointMother->setT(aPos->e());
  }
  else {
    fEmissionPointMother = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetEmissionPointMother(const AliFemtoLorentzVector& aPos)
{
  // set position from vector
  if (fEmissionPointMother) {
    fEmissionPointMother->setX(aPos.px());
    fEmissionPointMother->setY(aPos.py());
    fEmissionPointMother->setZ(aPos.pz());
    fEmissionPointMother->setT(aPos.e());
  }
  else {
    fEmissionPointMother = new AliFemtoLorentzVector();
    *fEmissionPointMother = aPos;
  }
}
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetEmissionPointMother(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // set position from components
  if (fEmissionPointMother) {
    fEmissionPointMother->setX(aRx);
    fEmissionPointMother->setY(aRy);
    fEmissionPointMother->setZ(aRz);
    fEmissionPointMother->setT(aT);
  }
  else {
    fEmissionPointMother = new AliFemtoLorentzVector(aRx, aRy, aRz, aT); 
  }
}
//_____________________________________________
//  Mother PID
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetPDGPidMother(Int_t aPidMother)
{
  fPDGPidMother = aPidMother;
}
//_____________________________________________
//  Child1 momentum
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild1(AliFemtoLorentzVector *aMomChild1)
{
  // set momentum from vector
  if (fTrueMomentumChild1) {
    fTrueMomentumChild1->setX(aMomChild1->px());
    fTrueMomentumChild1->setY(aMomChild1->py());
    fTrueMomentumChild1->setZ(aMomChild1->pz());
    fTrueMomentumChild1->setT(aMomChild1->e());
  }
  else {
    fTrueMomentumChild1 = new AliFemtoLorentzVector(*aMomChild1);
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild1(const AliFemtoLorentzVector& aMomChild1)
{
  // set momentum from vector and energy
  if (fTrueMomentumChild1) {
    fTrueMomentumChild1->setX(aMomChild1.px());
    fTrueMomentumChild1->setY(aMomChild1.py());
    fTrueMomentumChild1->setZ(aMomChild1.pz());
    fTrueMomentumChild1->setT(aMomChild1.e());
  }
  else {
    fTrueMomentumChild1 = new AliFemtoLorentzVector();
    *fTrueMomentumChild1 = aMomChild1;
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild1(Double_t aChild1Px, Double_t aChild1Py, Double_t aChild1Pz, Double_t aChild1E)
{
  // set momentum from components and energy
  if (!fTrueMomentumChild1) fTrueMomentumChild1 = new AliFemtoLorentzVector();
    fTrueMomentumChild1->setX(aChild1Px);
    fTrueMomentumChild1->setY(aChild1Py);
    fTrueMomentumChild1->setZ(aChild1Pz);
    fTrueMomentumChild1->setT(aChild1E);
}
//_____________________________________________
//  Child2 momentum
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild2(AliFemtoLorentzVector *aMomChild2)
{
  // set momentum from vector
  if (fTrueMomentumChild2) {
    fTrueMomentumChild2->setX(aMomChild2->px());
    fTrueMomentumChild2->setY(aMomChild2->py());
    fTrueMomentumChild2->setZ(aMomChild2->pz());
    fTrueMomentumChild2->setT(aMomChild2->e());
  }
  else {
    fTrueMomentumChild2 = new AliFemtoLorentzVector(*aMomChild2);
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild2(const AliFemtoLorentzVector& aMomChild2)
{
  // set momentum from vector and energy
  if (fTrueMomentumChild2) {
    fTrueMomentumChild2->setX(aMomChild2.px());
    fTrueMomentumChild2->setY(aMomChild2.py());
    fTrueMomentumChild2->setZ(aMomChild2.pz());
    fTrueMomentumChild2->setT(aMomChild2.e());
  }
  else {
    fTrueMomentumChild2 = new AliFemtoLorentzVector();
    *fTrueMomentumChild2 = aMomChild2;
  }
}
//_____________________________________________
void AliFemtoModelAllHiddenInfo::SetTrueMomentumChild2(Double_t aChild2Px, Double_t aChild2Py, Double_t aChild2Pz, Double_t aChild2E)
{
  // set momentum from components and energy
  if (!fTrueMomentumChild2) fTrueMomentumChild2 = new AliFemtoLorentzVector();
    fTrueMomentumChild2->setX(aChild2Px);
    fTrueMomentumChild2->setY(aChild2Py);
    fTrueMomentumChild2->setZ(aChild2Pz);
    fTrueMomentumChild2->setT(aChild2E);
}
//_____________________________________________
//  Child1 PID
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetPDGPidChild1(Int_t aPidChild1)
{
  fPDGPidMother = aPidChild1;
}
//_____________________________________________
//  Child2 PID
//_____________________________________________
void                   AliFemtoModelAllHiddenInfo::SetPDGPidChild2(Int_t aPidChild2)
{
  fPDGPidMother = aPidChild2;
}
//
//  RETURN COPY
//
//_____________________________________________
 AliFemtoModelHiddenInfo* AliFemtoModelAllHiddenInfo::GetParticleHiddenInfo() const
{
  // return copy of this hidden info
  AliFemtoModelAllHiddenInfo* tBuf = new AliFemtoModelAllHiddenInfo(*this);
  return tBuf;
}
