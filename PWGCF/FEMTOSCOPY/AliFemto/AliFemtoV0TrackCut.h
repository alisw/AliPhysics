///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoESDTrackCut: A basic track cut that used information from     //
// ALICE ESD to accept or reject the track.                              //  
// Enables the selection on charge, transverse momentum, rapidity,       //
// pid probabilities, number of ITS and TPC clusters                     //
// Author: Marek Chojnacki (WUT), mchojnacki@knf.pw.edu.pl               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOV0TRACKCUT_H
#define ALIFEMTOV0TRACKCUT_H

#include "AliFemtoTrackCut.h"

class AliFemtoV0TrackCut : public AliFemtoParticleCut 
{
  public:
  enum V0Type {kLambda = 0, kAntiLambda=1, kAll=99, kLambdaMC=101, kAntiLambdaMC=102};
  typedef enum V0Type AliFemtoV0Type;


  AliFemtoV0TrackCut();
  virtual ~AliFemtoV0TrackCut();

  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtV0;}
  
  void SetInvariantMassLambda(double,double);
  void SetMinDaughtersToPrimVertex(double);
  void SetMaxDcaV0Daughters(double);
  void SetMaxDcaV0(double);
  void SetMaxCosPointingAngle(double);
  void SetParticleType(short);
  void SetEta(double);
  void SetPt(double,double);
  void SetEtaDaughters(float);
  void SetTPCnclsDaughters(int);
  void SetNdofDaughters(int);
  void SetStatusDaughters(unsigned long);
  void SetPtDaughters(float,float);
  void SetOnFlyStatus(bool);

  //----n sigma----
  bool IsKaonTPCdEdxNSigma(float mom, float nsigmaK);
  bool IsKaonTOFNSigma(float mom, float nsigmaK);
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);
  
 private:   // here are the quantities I want to cut on...

  double            fInvMassLambdaMin;   //invariant mass lambda min
  double            fInvMassLambdaMax;   //invariant mass lambda max
  double            fMinDcaDaughtersToVert; //DCA of daughters to primary vertex
  double            fMaxDcaV0Daughters;     //Max DCA of v0 daughters at Decay vertex
  double            fMaxDcaV0;
  
  double            fMaxCosPointingAngle;
  short             fParticleType; //0-lambda
  double            fEta;
  double            fPtMin;
  double            fPtMax;
  bool              fOnFlyStatus;

  float fMaxEtaDaughters;			            // Eta of positive daughter
  int   fTPCNclsDaughters;			            // No. of cls of pos daughter
  int   fNdofDaughters;			                    // No. of degrees of freedom of the pos. daughter track
  unsigned long fStatusDaughters;			    // Status (tpc refit, its refit...)
  float fPtMinDaughters;
  float fPtMaxDaughters;

#ifdef __ROOT__ 
  ClassDef(AliFemtoV0TrackCut, 1)
#endif

};


#endif

