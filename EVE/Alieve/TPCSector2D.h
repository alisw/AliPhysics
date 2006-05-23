// $Header$

#ifndef ALIEVE_TPCSector2D_H
#define ALIEVE_TPCSector2D_H

#include "TPCSectorViz.h"


namespace Alieve {

class TPCSector2DEditor;
class TPCSector2DGL;

class TPCSector2D : public TPCSectorViz
{
  friend class TPCSector2DGL;
  friend class TPCSector2DEditor;

protected:
  Bool_t      fShowMax;
  Bool_t      fAverage;

  Bool_t      fUseTexture;

public:
  TPCSector2D(const Text_t* n="TPCSector2D", const Text_t* t=0);
  virtual ~TPCSector2D();

  void SetShowMax(Bool_t sm)  { fShowMax  = sm;  IncRTS(); }
  void SetAverage(Bool_t avg) { fAverage  = avg; IncRTS(); }

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option="");

  ClassDef(TPCSector2D, 1); // Visualization of TPC raw-data in 2D
}; // endclass TPCSector2D

}

#endif
