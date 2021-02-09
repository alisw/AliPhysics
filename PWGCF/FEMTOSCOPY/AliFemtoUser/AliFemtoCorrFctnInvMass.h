////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnInvMass - A correlation function that analyzes      //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNINVMASS_H
#define ALIFEMTOCORRFCTNINVMASS_H

#include "AliFemtoCorrFctn.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

class AliFemtoCorrFctnInvMass : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnInvMass(const char* title, const int& aBins, const double& aMin, const double& aMax, const double& aMass1, const double& aMass2);
  AliFemtoCorrFctnInvMass(const AliFemtoCorrFctnInvMass& aCorrFctn);
  virtual ~AliFemtoCorrFctnInvMass();

  AliFemtoCorrFctnInvMass& operator=(const AliFemtoCorrFctnInvMass& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();


  void SetParticleMasses(double mass1, double mass2);
  void SetParticle1Mass(double mass);
  void SetParticle2Mass(double mass);


  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnInvMass(*this); }

protected:

  TH1D *fNumInvMass; //invariant mass plot same events
  TH1D *fDenInvMass; //invariant mass plot mixed events
  
  int fBins;
  double fMin;
  double fMax;
  double fMass1;
  double fMass2;
  
  TString fTitle;



#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnInvMass, 1)
#endif
};


#endif
