///
/// \file AliFemtoModelHiddenInfo.h
///

#ifndef ALIFEMTOMODELHIDDENINFO_H
#define ALIFEMTOMODELHIDDENINFO_H

#include "AliFemtoHiddenInfo.h"
#include "AliFemtoTypes.h"

/// \class AliFemtoModelHiddenInfo
/// \brief The hidden info for model calculations
///
/// Stores information needed for the weight generation - the true
/// simulated momenta, freeze-out coordinates from model and particle PID
///
class AliFemtoModelHiddenInfo : public AliFemtoHiddenInfo {
public:

  AliFemtoModelHiddenInfo();
  AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo);
  virtual ~AliFemtoModelHiddenInfo();

  AliFemtoModelHiddenInfo& operator=(const AliFemtoModelHiddenInfo& aInfo);
  AliFemtoThreeVector   *GetTrueMomentum() const;
  AliFemtoThreeVector   *GetMotherMomentum() const;
  AliFemtoLorentzVector *GetEmissionPoint() const;
  Int_t                  GetPDGPid() const;
  Int_t                  GetMotherPdgCode() const;
  Double_t               GetMass() const;

  AliFemtoThreeVector   *GetTrueMomentumPos() const;
  AliFemtoLorentzVector *GetEmissionPointPos() const;
  Int_t                  GetPDGPidPos() const;
  Double_t               GetMassPos() const;

  AliFemtoThreeVector   *GetTrueMomentumNeg() const;
  AliFemtoLorentzVector *GetEmissionPointNeg() const;
  Int_t                  GetPDGPidNeg() const;
  Double_t               GetMassNeg() const;

  Int_t                  GetOrigin() const;

  void                   SetTrueMomentum(AliFemtoThreeVector *aMom);
  void                   SetTrueMomentum(const AliFemtoThreeVector& aMom);
  void                   SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz);
    
    void                   SetMotherMomentum(AliFemtoThreeVector *aMom);
    void                   SetMotherMomentum(const AliFemtoThreeVector& aMom);
    void                   SetMotherMomentum(Double_t aPx, Double_t aPy, Double_t aPz);
    
    
  void                   SetEmissionPoint(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPoint(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPid(Int_t aPid);
  void                   SetMotherPdgCode(Int_t motherPdg);
  void                   SetMass(Double_t aMass);

  void                   SetTrueMomentumPos(AliFemtoThreeVector *aMom);
  void                   SetTrueMomentumPos(const AliFemtoThreeVector& aMom);
  void                   SetTrueMomentumPos(Double_t aPx, Double_t aPy, Double_t aPz);
  void                   SetEmissionPointPos(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPointPos(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPointPos(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPidPos(Int_t aPid);
  void                   SetMassPos(Double_t aMass);

  void                   SetTrueMomentumNeg(AliFemtoThreeVector *aMom);
  void                   SetTrueMomentumNeg(const AliFemtoThreeVector& aMom);
  void                   SetTrueMomentumNeg(Double_t aPx, Double_t aPy, Double_t aPz);
  void                   SetEmissionPointNeg(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPointNeg(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPointNeg(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPidNeg(Int_t aPid);
  void                   SetMassNeg(Double_t aMass);

  void                   SetOrigin(Int_t origin);

// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* Clone() const;

protected:
  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const;

  AliFemtoThreeVector   *fTrueMomentum;  ///< True (simulated) momentum
  AliFemtoThreeVector   *fMotherMomentum;///< Momentum of mother particle
  AliFemtoLorentzVector *fEmissionPoint; ///< Emission point coordinates
  Int_t                  fPDGPid;        ///< True PID of the particle
  Int_t                  fMotherPdg;     ///< PDG code of particle's mother
  Double_t               fMass;          ///< True particle mass

  //daughter particles
  AliFemtoThreeVector   *fTrueMomentumPos;  ///< True (simulated) momentum of positive daughter
  AliFemtoLorentzVector *fEmissionPointPos; ///< Emission point coordinates of positive daughter
  Int_t                  fPDGPidPos;        ///< True PID of positive daughter
  Double_t               fMassPos;          ///< True particle mass of positive daughter

  AliFemtoThreeVector   *fTrueMomentumNeg;  ///< True (simulated) momentum of negative daughter
  AliFemtoLorentzVector *fEmissionPointNeg; ///< Emission point coordinates of negative daughter
  Int_t                  fPDGPidNeg;        ///< True PID of negative daughter
  Double_t               fMassNeg;          ///< True particle mass of negative daughter

  Int_t                  fOrigin;           ///< Origin of particles, 0 - physical primary, 1 - secondary from weak decay, 2 - secondary from material
};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::Clone() const
{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
