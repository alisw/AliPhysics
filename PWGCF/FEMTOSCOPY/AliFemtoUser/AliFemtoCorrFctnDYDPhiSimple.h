////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDYDPhiSimple - A correlation function that analyzes        //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and rapidity (y) difference                                                //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch,                                  //
//          Piotr Modzelewski Piotr.Mateusz.Modzelewski@cern.ch               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDYDPHISIMPLE_H
#define ALIFEMTOCORRFCTNDYDPHISIMPLE_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDYDPhiSimple : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnDYDPhiSimple(const char* title, const int& aPhiBins, const int& aYBins, const double& mass1, const double& mass2);
  AliFemtoCorrFctnDYDPhiSimple(const AliFemtoCorrFctnDYDPhiSimple& aCorrFctn);
  virtual ~AliFemtoCorrFctnDYDPhiSimple();

  AliFemtoCorrFctnDYDPhiSimple& operator=(const AliFemtoCorrFctnDYDPhiSimple& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();

  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDYDPhiSimple(*this); }

private:

  TH2D *fDPhiDYNumerator;            // Numerator of dY dPhi function
  TH2D *fDPhiDYDenominator;          // Denominator of dY dPhi function

  TH1D *fPhi;
  TH1D *fY;

  double fphiL;
  double fphiT;

  double fMass1;
  double fMass2;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDYDPhiSimple, 1)
#endif
};


#endif

