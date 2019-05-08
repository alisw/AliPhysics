////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle mass minvariant with various mass assumptions                     //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNPTKSTARMONITOR_H
#define ALIFEMTOCORRFCTNPTKSTARMONITOR_H

#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnPtKstar : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnPtKstar(const char* title);
  AliFemtoCorrFctnPtKstar(const AliFemtoCorrFctnPtKstar& aCorrFctn);
  virtual ~AliFemtoCorrFctnPtKstar();

  AliFemtoCorrFctnPtKstar& operator=(const AliFemtoCorrFctnPtKstar& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnPtKstar(*this); }

private:

  TH2D *fPtKstar;   // pt vs k* 
  TH2D *fPtKstarDen;   // pt vs k* mixed pairs

  TH2D *fPtKstar2part;   // pt vs k* 
  TH2D *fPtKstarDen2part;   // pt vs k* mixed pairs 

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnPtKstar, 1)
#endif
};


#endif

