///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoESDTrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                              //
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Marek Chojnacki (WUT), mchojnacki@knf.pw.edu.pl               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOXITRACKCUT_H
#define ALIFEMTOXITRACKCUT_H

#include "AliFemtoV0TrackCut.h"

#include "TH1D.h"

class AliFemtoXiTrackCut : public AliFemtoV0TrackCut
{
  public:
  enum XiType {kXiMinus = 0, kXiPlus=1, kAll = 99, kXiMinusMC = 101, kXiPlusMC=102};
  typedef enum XiType AliFemtoXiType;


  AliFemtoXiTrackCut();
  virtual ~AliFemtoXiTrackCut();

  AliFemtoXiTrackCut(const AliFemtoXiTrackCut& aCut);
  AliFemtoXiTrackCut& operator=(const AliFemtoXiTrackCut& aCut);
  virtual AliFemtoXiTrackCut* Clone();

  virtual bool Pass(const AliFemtoXi* aXi);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtXi;}

  bool IsPionNSigmaBac(float mom, float nsigmaTPCPi, float nsigmaTOFPi);  //!!!!! There is no need for this, as we can call AliFemtoV0TrackCut::IsPionNSigma

  void SetEtaXi(double x);
  void SetPtXi(double min, double max);
  void SetChargeXi(int x);  //!!!!!!!!!!!!!!!!!!!!!! To be deleted!!!  See comment at member fCharge
  void SetMaxDecayLengthXi(double x);
  void SetEtaBac(double x);
  void SetPtBac(double min, double max);
  void SetTPCnclsBac(int x);
  void SetNdofBac(double x);
  void SetStatusBac(unsigned long x);
  void SetInvariantMassXi(double min, double max);
  void SetMaxDcaXi(double x);
  void SetMinDcaXiBac(double x);
  void SetMaxDcaXiDaughters(double x);
  void SetMinCosPointingAngleXi(double min);
  void SetMinCosPointingAngleV0toXi(double min);
  void SetParticleTypeXi(short x);

  //-----The fMinvPurityAidHistoXi is built immediately before the (final) invariant mass cut, and thus may be used to calculate the purity of the Xi collection
  void SetMinvPurityAidHistoXi(const char* name, const char* title, const int& nbins, const float& aInvMassMin, const float& aInvMassMax);  //set the Minv histogram attributes and automatically sets flag fBuildPurityAidXi=true
  TH1D* GetMinvPurityAidHistoXi();
  //--If initiated, fMinvPurityAidHistoXi will be in the output list after all other Xi cut monitors (pass and fail)
  virtual TList *GetOutputList();  //include fMinvPurityAidHistoXi and fMinvPurityAidHistoV0 in the output list if they exist

 protected:   // here are the quantities I want to cut on...

  double fMaxEtaXi;
  double fMinPtXi;
  double fMaxPtXi;
  int fChargeXi;  //!!!!!!!!!!!!!!!!!! To be deleted!!!  Currently, this member is not used anywhere in the code.
		  // It seems redundant to set both fParticleTypeXi and fChargeXi.  Also, removing this member will
		  // eliminate any error caused by a user, for instance, setting fParticleTypeXi=kXiMinus and fCharge=+1,
		  // especially since fCharge is currently initiated in the constructor to 

  double fMaxEtaBac;
  double fMinPtBac;
  double fMaxPtBac;
  int fTPCNclsBac;
  double fNdofBac;
  unsigned long fStatusBac;
  double fMaxDcaXi;
  double fMinDcaXiBac;
  double fMaxDcaXiDaughters;
  double fMinCosPointingAngleXi;
  double fMinCosPointingAngleV0toXi;
  double fMaxDecayLengthXi;
  double fInvMassXiMin;
  double fInvMassXiMax;
  short  fParticleTypeXi;

  bool fBuildPurityAidXi;
  TH1D* fMinvPurityAidHistoXi;


#ifdef __ROOT__
  ClassDef(AliFemtoXiTrackCut, 1);
#endif

};

inline TH1D* AliFemtoXiTrackCut::GetMinvPurityAidHistoXi() {return fMinvPurityAidHistoXi;}
#endif
