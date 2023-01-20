////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle mass minvariant with various mass assumptions                     //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNPTKSTARBETA_H
#define ALIFEMTOCORRFCTNPTKSTARBETA_H

#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnPtKstarBeta : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnPtKstarBeta(const char* title);
  AliFemtoCorrFctnPtKstarBeta(const AliFemtoCorrFctnPtKstarBeta& aCorrFctn);
  virtual ~AliFemtoCorrFctnPtKstarBeta();

  AliFemtoCorrFctnPtKstarBeta& operator=(const AliFemtoCorrFctnPtKstarBeta& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnPtKstarBeta(*this); }

  void SetBtPtKStart(bool BtPtKStart);

private:

  TH2D *fPtKstar;   // pt vs k* 
  TH2D *fPtKstarDen;   // pt vs k* mixed pairs

  TH2D *fPtKstar2part;   // pt vs k* 
  TH2D *fPtKstarDen2part;   // pt vs k* mixed pairs 

  TH2D *fPtBeta[10];   // pt vs k* part1
  TH2D *fPtBeta2part[10];  //pt vs k* part2

  TH2D *fPtkT[10];   // pt vs k* part1
  TH2D *fPtkT2part[10];  //pt vs k* part2
  
  bool fBtPtKStart;

  
  TH2D *fPairPtKstar2part;   // pair pt vs k* 
  TH2D *fPairPtKstarDen2part;   // pair pt vs k* mixed pairs 

  TH2D *fPtKstar_kT[10];   // pt vs k* for kT 
  TH2D *fPtKstarDen_kT[10];   // pt vs k* mixed pairs for kT
  TH2D *fPtKstar2part_kT[10];   // pt vs k* for kT
  TH2D *fPtKstarDen2part_kT[10];   // pt vs k* mixed pairs for kT
  TH2D *fPairPtKstar2part_kT[10];   // pair pt vs k* for kT
  TH2D *fPairPtKstarDen2part_kT[10];   // pair pt vs k* mixed pairs for kT 

  
  TH2D *fKstarBetaT;  //k* vs BetaT
  
#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnPtKstarBeta, 1)
#endif
};


#endif

