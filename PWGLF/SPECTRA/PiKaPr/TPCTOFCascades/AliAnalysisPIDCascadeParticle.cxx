#include "AliAnalysisPIDCascadeParticle.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"

ClassImp(AliAnalysisPIDCascadeParticle)

//___________________________________________________________

TLorentzVector AliAnalysisPIDCascadeParticle::fgLorentzVector;

//___________________________________________________________

Double_t
AliAnalysisPIDCascadeParticle::GetY() const
{

  fgLorentzVector.SetPtEtaPhiM(fPt, fEta, fPhi, GetMass());
  return fgLorentzVector.Rapidity();
}

//___________________________________________________________

Int_t
AliAnalysisPIDCascadeParticle::GetSign() const
{

  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TParticlePDG *ppdg = dbpdg->GetParticle(fPdgCode);
  if (!ppdg)
    return 0;
  return TMath::Nint(ppdg->Charge());
}

//___________________________________________________________

Int_t
AliAnalysisPIDCascadeParticle::GetPID() const
{
  /*
   * get PID
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    if (TMath::Abs(fPdgCode) == AliPID::ParticleCode(ipart))
      return ipart;
  return -1;

}

//___________________________________________________________

Double_t
AliAnalysisPIDCascadeParticle::GetMass() const
{
  /*
   * get mass
   */

  if (GetPID() == -1)
    return 0.;
  return AliPID::ParticleMass(GetPID());
}

//___________________________________________________________

AliAnalysisPIDCascadeParticle::AliAnalysisPIDCascadeParticle() :
  TObject(),
  fLabel(0),
  fPt(0.),
  fEta(0.),
  fPhi(0.),
  fPdgCode(0),
  fMotherPdgCode(0),
  fPrimary(kFALSE)
{
  /*
   * default constructor
   */
}

//___________________________________________________________

AliAnalysisPIDCascadeParticle::AliAnalysisPIDCascadeParticle(const AliAnalysisPIDCascadeParticle &source) :
  TObject(source),
  fLabel(source.fLabel),
  fPt(source.fPt),
  fEta(source.fEta),
  fPhi(source.fPhi),
  fPdgCode(source.fPdgCode),
  fMotherPdgCode(source.fMotherPdgCode),
  fPrimary(source.fPrimary)
{
  /*
   * copy constructor
   */
}

//___________________________________________________________

AliAnalysisPIDCascadeParticle &
AliAnalysisPIDCascadeParticle::operator=(const AliAnalysisPIDCascadeParticle &source)
{
  /*
   * operator=
   */

  if (&source == this) return *this;
  TObject::operator=(source);

  fLabel = source.fLabel;
  fPt = source.fPt;
  fEta = source.fEta;
  fPhi = source.fPhi;
  fPdgCode = source.fPdgCode;
  fMotherPdgCode = source.fMotherPdgCode;
  fPrimary = source.fPrimary;
  return *this;
}

//___________________________________________________________

AliAnalysisPIDCascadeParticle::~AliAnalysisPIDCascadeParticle()
{
  /*
   * default destructor
   */
}

//___________________________________________________________

void
AliAnalysisPIDCascadeParticle::Reset()
{
  /*
   * reset
   */

  fLabel = 0;
  fPt = 0.;
  fEta = 0.;
  fPhi = 0.;
  fPdgCode = 0;
  fMotherPdgCode = 0;
  fPrimary = kFALSE;
}

//___________________________________________________________

void
AliAnalysisPIDCascadeParticle::Update(AliMCParticle *particle, Int_t label, Int_t MotherPdg, Bool_t PrimCheck)
{
  /*
   * update
   */

  fLabel = label;
  fMotherPdgCode = MotherPdg;
  fPt = particle->Pt();
  fEta = particle->Eta();
  fPhi = particle->Phi();
  fPdgCode = particle->PdgCode();
  fPrimary = PrimCheck;
}
