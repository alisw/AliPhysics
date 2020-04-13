////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle mass minvariant with various mass assumptions                     //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNPAIRSFORCORRFIT_H
#define ALIFEMTOCORRFCTNPAIRSFORCORRFIT_H

#include "TH1F.h"
#include "TNtuple.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnPairsForCorrFit : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnPairsForCorrFit(const char* title);
  AliFemtoCorrFctnPairsForCorrFit(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn);
  virtual ~AliFemtoCorrFctnPairsForCorrFit();

  AliFemtoCorrFctnPairsForCorrFit& operator=(const AliFemtoCorrFctnPairsForCorrFit& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnPairsForCorrFit(*this); }
 
private:

  TNtuple* mNtuple;
  TH1F* hKstar; 

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnPairsForCorrFit, 1)
#endif
};


#endif

