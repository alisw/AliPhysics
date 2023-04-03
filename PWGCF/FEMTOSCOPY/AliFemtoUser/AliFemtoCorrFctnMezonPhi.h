////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnMezonPhi - A correlation function that analyzes             //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference. Idealy particles come from mezon phi  //
// decay																	  //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
// Edits: Daniel Rodak                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNMEZONPHI_H
#define ALIFEMTOCORRFCTNMEZONPHI_H

#include "AliFemtoCorrFctn.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class AliFemtoCorrFctnMezonPhi : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnMezonPhi(const char* title, const int& aBins, const double& aMin, const double& aMax, const double& aMass1, const double& aMass2, const int& sections, const double& pMin, const double& pMax);
  AliFemtoCorrFctnMezonPhi(const AliFemtoCorrFctnMezonPhi& aCorrFctn);
  virtual ~AliFemtoCorrFctnMezonPhi();

  AliFemtoCorrFctnMezonPhi& operator=(const AliFemtoCorrFctnMezonPhi& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);
  //virtual void AddPairPt(AliFemtoPair* aPair);

  virtual void Finish();


  void SetParticleMasses(double mass1, double mass2);
  void SetParticle1Mass(double mass);
  void SetParticle2Mass(double mass);


  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnMezonPhi(*this); }

protected:

  TH1D **fNumInvMass; //invariant mass plot same events
  TH1D **fDenInvMass; //invariant mass plot mixed events
  TH1D **fNumTransvMom; //transver momentum of same events parent particle
  TH1D **fDenTransvMom; //transver momentum of mixed events parent particle
  TH2D *fNumInvMass2D;
  TH2D *fDenInvMass2D;
  
  int fBins;
  int fSections;
  double fMin;
  double fMax;
  double fpMin;
  double fpMax;
  double fMass1;
  double fMass2;
  double fIntersection;
  
  TString fTitle;



#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnMezonPhi, 1)
#endif
};


#endif
