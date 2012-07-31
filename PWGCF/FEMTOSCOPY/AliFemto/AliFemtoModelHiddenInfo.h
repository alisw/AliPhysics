////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelHiddenInfo - the hidden info for model calculations         ///
/// Stores information needed for the weight generation - the true           ///
/// simulated momenta, freeze-out coordinates from model and particle PID    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELHIDDENINFO_H
#define ALIFEMTOMODELHIDDENINFO_H

#include <TH1D.h>
#include "AliFemtoTypes.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoHiddenInfo.h"

class AliFemtoModelHiddenInfo : public AliFemtoHiddenInfo{

public:
  AliFemtoModelHiddenInfo();
  AliFemtoModelHiddenInfo(const AliFemtoModelHiddenInfo &aInfo);
  virtual ~AliFemtoModelHiddenInfo();

  AliFemtoModelHiddenInfo& operator=(const AliFemtoModelHiddenInfo& aInfo);
  AliFemtoThreeVector   *GetTrueMomentum() const;
  AliFemtoLorentzVector *GetEmissionPoint() const;
  Int_t                  GetPDGPid() const;
  Double_t               GetMass() const;

  AliFemtoThreeVector   *GetTrueMomentumPos() const;
  AliFemtoLorentzVector *GetEmissionPointPos() const;
  Int_t                  GetPDGPidPos() const;
  Double_t               GetMassPos() const;

  AliFemtoThreeVector   *GetTrueMomentumNeg() const;
  AliFemtoLorentzVector *GetEmissionPointNeg() const;
  Int_t                  GetPDGPidNeg() const;
  Double_t               GetMassNeg() const;

  void                   SetTrueMomentum(AliFemtoThreeVector *aMom);
  void                   SetTrueMomentum(const AliFemtoThreeVector& aMom);
  void                   SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz);
  void                   SetEmissionPoint(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPoint(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPid(Int_t aPid);
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

// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* Clone() const;
  
 protected:
  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const;

  AliFemtoThreeVector   *fTrueMomentum;  // True (simulated) momentum
  AliFemtoLorentzVector *fEmissionPoint; // Emission point coordinates
  Int_t                  fPDGPid;        // True PID of the particle
  Double_t               fMass;          // True particle mass

  //daughter particles
  AliFemtoThreeVector   *fTrueMomentumPos;  // True (simulated) momentum of positive daughter
  AliFemtoLorentzVector *fEmissionPointPos; // Emission point coordinates of positive daughter
  Int_t                  fPDGPidPos;        // True PID of positive daughter
  Double_t               fMassPos;          // True particle mass of positive daughter

  AliFemtoThreeVector   *fTrueMomentumNeg;  // True (simulated) momentum of negative daughter
  AliFemtoLorentzVector *fEmissionPointNeg; // Emission point coordinates of negative daughter
  Int_t                  fPDGPidNeg;        // True PID of negative daughter
  Double_t               fMassNeg;          // True particle mass of negative daughter
};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::Clone() const{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
