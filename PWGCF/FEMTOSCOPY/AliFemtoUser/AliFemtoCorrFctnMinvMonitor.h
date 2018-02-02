////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle mass minvariant with various mass assumptions                     //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//          Anna Zaborowska azaborow@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNMINVMONITOR_H
#define ALIFEMTOCORRFCTNMINVMONITOR_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnMinvMonitor : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnMinvMonitor(const char* title);
  AliFemtoCorrFctnMinvMonitor(const AliFemtoCorrFctnMinvMonitor& aCorrFctn);
  virtual ~AliFemtoCorrFctnMinvMonitor();

  AliFemtoCorrFctnMinvMonitor& operator=(const AliFemtoCorrFctnMinvMonitor& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnMinvMonitor(*this); }

private:

  TH1D *fMinveeFail;   // ee mass assumption - failed pairs
  TH1D *fMinvee;   // ee mass assumption - passed pairs
  TH1D *fMinv2piFail;   // 2 pi mass assumption - failed pairs
  TH1D *fMinv2pi;   // 2 pi mass assumption - passed pairs
  TH1D *fMinvppiFail;   // p pi mass assumption - failed pairs
  TH1D *fMinvppi;   // p pi mass assumption - passed pairs

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnMinvMonitor, 1)
#endif
};


#endif

