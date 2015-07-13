#include "AliAnalysisParticle.h"
#include "TParticle.h"
#include "AliPID.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"

ClassImp(AliAnalysisParticle)

//___________________________________________________________

TLorentzVector AliAnalysisParticle::fgLorentzVector;

//___________________________________________________________

Double_t
AliAnalysisParticle::GetY() const
{
  
  fgLorentzVector.SetPtEtaPhiM(fPt, fEta, fPhi, GetMass());
  return fgLorentzVector.Rapidity();
}

//___________________________________________________________

Float_t
AliAnalysisParticle::GetSign() const
{
  
  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TParticlePDG *ppdg = dbpdg->GetParticle(fPdgCode);
  if (!ppdg)
    return 0.;
  return TMath::Sign(1., ppdg->Charge());
}

//___________________________________________________________

Int_t
AliAnalysisParticle::GetPID() const
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
AliAnalysisParticle::GetMass() const
{
  /*
   * get mass
   */
  
  if (GetPID() == -1)
    return 0.;
  return AliPID::ParticleMass(GetPID());
}

//___________________________________________________________

AliAnalysisParticle::AliAnalysisParticle() :
  TObject(),
  fLabel(0),
  fPt(0.),
  fEta(0.),
  fPhi(0.),
  fPdgCode(0)
{
  /*
   * default constructor
   */
}

//___________________________________________________________

AliAnalysisParticle::AliAnalysisParticle(const AliAnalysisParticle &source) :
  TObject(source),
  fLabel(source.fLabel),
  fPt(source.fPt),
  fEta(source.fEta),
  fPhi(source.fPhi),
  fPdgCode(source.fPdgCode)
{
  /*
   * copy constructor
   */
}

//___________________________________________________________

AliAnalysisParticle &
AliAnalysisParticle::operator=(const AliAnalysisParticle &source)
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

  return *this;
}

//___________________________________________________________

AliAnalysisParticle::~AliAnalysisParticle()
{
  /*
   * default destructor
   */
}

//___________________________________________________________

void
AliAnalysisParticle::Reset()
{
  /*
   * reset
   */

  fLabel = 0;
  fPt = 0.;
  fEta = 0.;
  fPhi = 0.;
  fPdgCode = 0;
  
}

//___________________________________________________________

void
AliAnalysisParticle::Update(TParticle *particle, Int_t label)
{
  /*
   * update
   */

  fLabel = label;
  fPt = particle->Pt();
  fEta = particle->Eta();
  fPhi = particle->Phi();
  fPdgCode = particle->GetPdgCode();
  
}

