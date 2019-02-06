////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnPairFractions - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Malgorzata Janik majanik@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNPAIRFRACTIONS_H
#define ALIFEMTOCORRFCTNPAIRFRACTIONS_H

#include "TH1F.h"
#include "TH2F.h"

#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnPairFractions : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kEta=2};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnPairFractions(const char* title);
  AliFemtoCorrFctnPairFractions(const AliFemtoCorrFctnPairFractions& aCorrFctn);
  virtual ~AliFemtoCorrFctnPairFractions();

  AliFemtoCorrFctnPairFractions& operator=(const AliFemtoCorrFctnPairFractions& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnPairFractions(*this); }

  void SetDoDEtaDPhiMaps(bool dodedp=true);

private:
  TH1F *fPairFractions;
  TH1F *fPairFractionsDen;

  double fphiL;
  double fphiT;

  bool detadphi;
  TH2F *fPairFractionsDEtaDPhi[7];
  TH2F *fPairFractionsDenDEtaDPhi[7];


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnPairFractions, 1);
  /// \endcond
#endif
};


#endif
