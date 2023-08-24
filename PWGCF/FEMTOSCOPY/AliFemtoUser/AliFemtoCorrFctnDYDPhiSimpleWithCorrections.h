/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                 //
// AliFemtoCorrFctnDYDPhiSimpleWithCorrections - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        		       //
// and rapidity (y) difference                                              			       //
//                                                                     			       //
//                                                                     			       //
/////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDYDPHISIMPLEWITHCORR_H
#define ALIFEMTOCORRFCTNDYDPHISIMPLEWITHCORR_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"


class AliFemtoCorrFctnDYDPhiSimpleWithCorrections : public AliFemtoCorrFctn {
public:
  enum ParticleType {kNoCorrection=0, kPion=1, kKaon=2, kProton=3, kAll=4, kPionMinus=5, kKaonMinus=6, kProtonMinus=7, kLambda=8, kLambdaMinus=9,  kXiMinus=10, kXiPlus=11};

  AliFemtoCorrFctnDYDPhiSimpleWithCorrections(const char* title, const int& aPhiBins, const int& aYBins, const double& mass1, const double& mass2);
  AliFemtoCorrFctnDYDPhiSimpleWithCorrections(const AliFemtoCorrFctnDYDPhiSimpleWithCorrections& aCorrFctn);
  virtual ~AliFemtoCorrFctnDYDPhiSimpleWithCorrections();

  AliFemtoCorrFctnDYDPhiSimpleWithCorrections& operator=(const AliFemtoCorrFctnDYDPhiSimpleWithCorrections& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void SetParticleTypes(ParticleType partType1, ParticleType partType2);
  void SetParticle1Type(ParticleType partType);
  void SetParticle2Type(ParticleType partType);
  /*void SetParticleTypes(int part, double mass);
  
  if(part==1)
     fMass1=mass;
  else 
     fMass2=mass;
     */
 // void SetReadHiddenInfo(bool read);

  void WriteHistos();
  virtual TList* GetOutputList();

  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDYDPhiSimpleWithCorrections(*this); }

private:
  TH2D *fDPhiDYNumerator;            // Numerator of dY dPhi function
  TH2D *fDPhiDYDenominator;          // Denominator of dY dPhi function

  TH1D *fPhi;
  TH1D *fY;

  double fphiL;
  double fphiT;
  
  ParticleType part1;
  ParticleType part2;
  
  
  double fMass1;
  double fMass2;
  /*
  TH2D *fDPhiDYNumerator;          // Numerator of DY dPhi function
  TH2D *fDPhiDYDenominator;        // Denominator of DY dPhi function
  TH2D *fDPhiDYHiddenNumerator;          // Numerator of DY dPhi function from MC
  TH2D *fDPhiDYHiddenDenominator;        // Denominator of DY dPhi function from MC

  TH2D *fDPhiDYHiddenPrimaryNumerator;          // Numerator of DY dPhi function from MC, physical primaries
  TH2D *fDPhiDYHiddenPrimaryDenominator;        // Denominator of DY dPhi function from MC, physical primaries

  TH2D *fDPhiDYHiddenSecWeakNumerator;          // Numerator of DY dPhi function from MC, secondaries from weak decays
  TH2D *fDPhiDYHiddenSecWeakDenominator;        // Denominator of DY dPhi function from MC, secondaries from weak decays

  TH2D *fDPhiDYHiddenSecMatNumerator;          // Numerator of DY dPhi function from MC, secondaries from material
  TH2D *fDPhiDYHiddenSecMatDenominator;        // Denominator of DY dPhi function from MC, secondaries from material

  TH2D *fDPhiDYHiddenPrimaryNumeratorData;          // Numerator of DY dPhi function from MC, physical primaries, filled with data values
  TH2D *fDPhiDYHiddenPrimaryDenominatorData;        // Denominator of DY dPhi function from MC, physical primaries, filled with data values

  TH2D *fDPhiDYHiddenSecWeakNumeratorData;          // Numerator of DY dPhi function from MC, secondaries from weak decays, filled with data values
  TH2D *fDPhiDYHiddenSecWeakDenominatorData;        // Denominator of DY dPhi function from MC, secondaries from weak decays, filled with data values

  TH2D *fDPhiDYHiddenSecMatNumeratorData;          // Numerator of DY dPhi function from MC, secondaries from material, filled with data values
  TH2D *fDPhiDYHiddenSecMatDenominatorData;        // Denominator of DY dPhi function from MC, secondaries from material, filled with data values
//
  double fphiL;
  double fphiT;

  int fYBins;
  int fPhiBins;

  TString ftitle;

  ParticleType part1;
  ParticleType part2;
  
  double fMass1;
  double fMass2;

  bool fReadHiddenInfo;
*/
#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDYDPhiSimpleWithCorrections, 1)
#endif
};
inline void AliFemtoCorrFctnDYDPhiSimpleWithCorrections::SetParticleTypes(ParticleType partType1, ParticleType partType2) {   part1 = partType1;
   part2 = partType2;}
   
inline void AliFemtoCorrFctnDYDPhiSimpleWithCorrections::SetParticle1Type(ParticleType partType) {   part1= partType;}
inline void AliFemtoCorrFctnDYDPhiSimpleWithCorrections::SetParticle2Type(ParticleType partType) {   part2= partType;}
  

#endif
