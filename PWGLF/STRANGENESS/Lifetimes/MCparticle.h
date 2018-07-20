#ifndef Lifetimes_MCparticle_h
#define Lifetimes_MCparticle_h

#include <cmath>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace Lifetimes
{

class MCparticle
{
public:
  enum Status
  {
    kPrimary = 1,
    kSecondaryFromMaterial = 2,
    kSecondaryFromWeakDecay = 4
  };

  MCparticle();

  float
  GetDistOverP() const;
  float GetEta() const { return fEta; }
  float GetMass() const;
  int GetPDGcode() const { return fPDGcode; }
  float GetPt() const { return fPt; }
  int GetRecoIndex() const { return fRecoIndex; }
  float GetY() const;
  bool IsPrimary() const { return fStatus & kPrimary; }
  bool IsSecondaryFromMaterial() const { return fStatus & kSecondaryFromMaterial; }
  bool IsSecondaryFromWeakDecay() const { return fStatus & kSecondaryFromWeakDecay; }

  void SetDistOverP(float d) { fDistOverTotMom = d; }
  void SetEta(float eta) { fEta = eta; }
  void SetPt(float pt) { fPt = pt; }
  void SetPDGcode(int pdg) { fPDGcode = pdg; }
  void SetRecoIndex(int idx) { fRecoIndex = idx; }
  void SetStatus(Status st) { fStatus = st; }

private:
  float fPt;
  float fEta;
  float fDistOverTotMom;
  int fPDGcode;
  int fRecoIndex;
  unsigned char fStatus;
};

inline MCparticle::MCparticle() :
  fPt{0},
  fEta{-10},
  fDistOverTotMom{0},
  fPDGcode{0},
  fRecoIndex{-1},
  fStatus{0u}
{}

inline float
MCparticle::GetDistOverP() const
{
  return fDistOverTotMom;
}

inline float MCparticle::GetMass() const
{
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(fPDGcode);
  if (part)
    return part->Mass();
  else
    return -1.;
}

inline float
MCparticle::GetY() const
{
  return std::log((std::hypot(GetMass(), fPt * std::cosh(fEta)) + fPt * std::sinh(fEta)) / std::hypot(GetMass(), fPt));
}

} // namespace Lifetimes

#endif