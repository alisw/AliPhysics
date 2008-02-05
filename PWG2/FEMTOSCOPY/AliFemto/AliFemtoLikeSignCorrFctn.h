////////////////////////////////////////////////////////////////////////////////
/// AliFemtoLikeSignCorrFctn - the pure virtual base class for the like sign ///
/// correlation function. All like sign correlation functions  must inherit  ///
/// from this one                                                            ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoLikeSignCorrFctn_hh
#define AliFemtoLikeSignCorrFctn_hh

class AliFemtoPair;
#include "AliFemtoCorrFctn.h"

class AliFemtoLikeSignCorrFctn : public AliFemtoCorrFctn {

  friend class AliFemtoLikeSignAnalysis;

public:
  AliFemtoLikeSignCorrFctn(){/* no-op */};
  AliFemtoLikeSignCorrFctn(const AliFemtoLikeSignCorrFctn& aCorrFctn);
  virtual ~AliFemtoLikeSignCorrFctn(){/* no-op */};

  virtual void AddLikeSignPositivePair(const AliFemtoPair* aPair) = 0;
  virtual void AddLikeSignNegativePair(const AliFemtoPair* aPair) = 0;

  virtual AliFemtoLikeSignCorrFctn* Clone() { return 0;}
  virtual TList* GetOutputList() = 0;

  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
};
//________________________________________
inline AliFemtoLikeSignCorrFctn::AliFemtoLikeSignCorrFctn(const AliFemtoLikeSignCorrFctn& c) { fyAnalysis =0; }

#endif
