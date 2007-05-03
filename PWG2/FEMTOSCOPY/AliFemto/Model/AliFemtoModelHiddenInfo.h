////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelHiddenInfo - the hidden info for model calculations         ///
/// Stores information needed for the weight generation - the true           ///
/// simulated momenta, freeze-out coordinates from model and particle PID    ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoModelHiddenInfo_hh
#define AliFemtoModelHiddenInfo_hh

#include <TH1D.h>
#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoThreeVector.h"
#include "Infrastructure/AliFemtoLorentzVector.h"
#include "Base/AliFemtoHiddenInfo.h"

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

  void                   SetTrueMomentum(AliFemtoThreeVector *aMom);
  void                   SetTrueMomentum(const AliFemtoThreeVector& aMom);
  void                   SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz);
  void                   SetEmissionPoint(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPoint(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPid(Int_t aPid);
  void                   SetMass(Double_t aMass);

// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* getParticleHiddenInfo() const;
  virtual AliFemtoHiddenInfo* clone() const;

 protected:
  AliFemtoThreeVector   *fTrueMomentum;  // True (simulated) momentum
  AliFemtoLorentzVector *fEmissionPoint; // Emission point coordinates
  Int_t                  fPDGPid;        // True PID of the particle
  Double_t               fMass;
};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoModelHiddenInfo::clone() const{
  // return exact copy of this hidden info
  return getParticleHiddenInfo();
}

#endif
