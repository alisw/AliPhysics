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
//K                                                                          M//
//K//////////////////////////////////////////////////////////////////////////M//
#ifndef ALIFEMTOMODELALLHIDDENINFO_H
#define ALIFEMTOMODELALLHIDDENINFO_H

#include <TH1D.h>
#include "AliFemtoTypes.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoHiddenInfo.h"
#include "AliFemtoModelHiddenInfo.h"

class AliFemtoModelAllHiddenInfo : public AliFemtoModelHiddenInfo{

public:
  AliFemtoModelAllHiddenInfo();
  AliFemtoModelAllHiddenInfo(const AliFemtoModelAllHiddenInfo &aInfo);
  virtual ~AliFemtoModelAllHiddenInfo();

  AliFemtoModelAllHiddenInfo& operator=(const AliFemtoModelAllHiddenInfo& aInfo);

  AliFemtoLorentzVector *GetTrueMomentumMother() const;
  AliFemtoLorentzVector *GetEmissionPointMother() const;
  Int_t                  GetPDGPidMother() const;
  AliFemtoLorentzVector *GetTrueMomentumChild1() const;
  AliFemtoLorentzVector *GetTrueMomentumChild2() const;
  Int_t                  GetPDGPidChild1() const;
  Int_t                  GetPDGPidChild2() const;

  void                   SetTrueMomentumMother(AliFemtoLorentzVector *aMomMother);
  void                   SetTrueMomentumMother(const AliFemtoLorentzVector& aMomMother);
  void                   SetTrueMomentumMother(Double_t aMotherPx, Double_t aMotherPy, Double_t aMotherPz, Double_t aMotherE);
  void                   SetEmissionPointMother(AliFemtoLorentzVector *aPos);
  void                   SetEmissionPointMother(const AliFemtoLorentzVector& aPos);
  void                   SetEmissionPointMother(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT);
  void                   SetPDGPidMother(Int_t aPidMother);
  void                   SetTrueMomentumChild1(AliFemtoLorentzVector *aMomChild1);
  void                   SetTrueMomentumChild1(const AliFemtoLorentzVector& aMomChild1);
  void                   SetTrueMomentumChild1(Double_t aChild1Px, Double_t aChild1Py, Double_t aChild1Pz, Double_t aChild1E);
  void                   SetTrueMomentumChild2(AliFemtoLorentzVector *aMomChild2);
  void                   SetTrueMomentumChild2(const AliFemtoLorentzVector& aMomChild2);
  void                   SetTrueMomentumChild2(Double_t aChild2Px, Double_t aChild2Py, Double_t aChild2Pz, Double_t aChild2E);
  void                   SetPDGPidChild1(Int_t aPidChild1);
  void                   SetPDGPidChild2(Int_t aPidChild2);


// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoModelHiddenInfo* Clone() const;
 protected:
  virtual AliFemtoModelHiddenInfo* GetParticleHiddenInfo() const;

  AliFemtoLorentzVector *fTrueMomentumMother;  // True (simulated) momentum of Mother (100.,100.,100.,100. if Mother -1)
  AliFemtoLorentzVector *fEmissionPointMother; // Emission point coordinates Mother;
  Int_t                  fPDGPidMother;        // True PID of the particle mother
  AliFemtoLorentzVector *fTrueMomentumChild1;  // True (simulated) momentum of Child1  (200.,200.,200.,200. if Child1 -1)
  AliFemtoLorentzVector *fTrueMomentumChild2;  // True (simulated) momentum of Child2  (200.,200.,200.,200. if Child2 -1)
  Int_t                  fPDGPidChild1;        // True PID of the particle child1 (-1 if is not)
  Int_t                  fPDGPidChild2;        // True PID of the particle child2 (-1 if is not)

};
//_______________________________________
inline AliFemtoModelHiddenInfo* AliFemtoModelAllHiddenInfo::Clone() const{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
