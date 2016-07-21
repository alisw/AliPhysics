////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHISIMPLEWITHCORR_H
#define ALIFEMTOCORRFCTNDETADPHISIMPLEWITHCORR_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections : public AliFemtoCorrFctn {
public:
  enum ParticleType {kNoCorrection=0, kPion=1, kKaon=2, kProton=3, kAll=4, kPionMinus=5, kKaonMinus=6, kProtonMinus=7, kLambda=8};

  AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections(const AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections();

  AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections& operator=(const AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void SetParticleTypes(ParticleType partType1, ParticleType partType2);
  void SetParticle1Type(ParticleType partType);
  void SetParticle2Type(ParticleType partType);

  void SetReadHiddenInfo(bool read);

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fDPhiDEtaNumerator;          // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function
  TH2D *fDPhiDEtaHiddenNumerator;          // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaHiddenDenominator;        // Denominator of dEta dPhi function

  double fphiL;
  double fphiT;
  
  ParticleType part1;
  ParticleType part2;

  bool fReadHiddenInfo;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections, 1)
#endif
};


#endif
