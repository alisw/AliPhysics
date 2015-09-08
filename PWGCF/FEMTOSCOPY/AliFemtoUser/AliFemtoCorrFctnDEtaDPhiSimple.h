////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiSimple - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHISIMPLE_H
#define ALIFEMTOCORRFCTNDETADPHISIMPLE_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiSimple : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kEta=2};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnDEtaDPhiSimple(char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoCorrFctnDEtaDPhiSimple(const AliFemtoCorrFctnDEtaDPhiSimple& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiSimple();

  AliFemtoCorrFctnDEtaDPhiSimple& operator=(const AliFemtoCorrFctnDEtaDPhiSimple& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();


  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fDPhiDEtaNumerator;          // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function

  double fphiL;
  double fphiT;
  


#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiSimple, 1)
#endif
};


#endif
