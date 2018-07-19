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
  AliFemtoLikeSignCorrFctn():AliFemtoCorrFctn() {/* no-op */};
  AliFemtoLikeSignCorrFctn(const AliFemtoLikeSignCorrFctn& aCorrFctn);
  virtual ~AliFemtoLikeSignCorrFctn(){/* no-op */};
  AliFemtoLikeSignCorrFctn& operator=(const AliFemtoLikeSignCorrFctn& aCorrFctn);

  virtual void AddLikeSignPositivePair(const AliFemtoPair* aPair) = 0;
  virtual void AddLikeSignNegativePair(const AliFemtoPair* aPair) = 0;

  virtual TList* GetOutputList() = 0;

};
//________________________________________
inline AliFemtoLikeSignCorrFctn::AliFemtoLikeSignCorrFctn(const AliFemtoLikeSignCorrFctn& /* c */):AliFemtoCorrFctn() { fyAnalysis =0; }
inline AliFemtoLikeSignCorrFctn& AliFemtoLikeSignCorrFctn::operator=(const AliFemtoLikeSignCorrFctn& aCorrFctn) {   if (this != &aCorrFctn) { AliFemtoCorrFctn::operator=(aCorrFctn); } return *this; }


#endif
