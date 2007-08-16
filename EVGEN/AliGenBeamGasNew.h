#ifndef ALIGENBEAMGASNEW_H
#define ALIGENBEAMGASNEW_H

/* $Id$ */

#include "AliGenCocktail.h"

class AliGenBeamGasNew : public AliGenCocktail
{
 public:
  AliGenBeamGasNew();
  AliGenBeamGasNew(const AliGenBeamGasNew& rhs);
  //  AliGenBeamGasNew& operator=(const AliGenBeamGasNew&);
  virtual ~AliGenBeamGasNew();
  void SetTimeWindow(Float_t twindow) {fTwindow = twindow;}
  virtual void Generate();
  void VertexInternal();
  virtual void Init();

 protected:
  Float_t fItime;   // time of bg-interaction
  Float_t fTwindow; // time-window in which tpc-gate is open

  ClassDef(AliGenBeamGasNew,1);

};
#endif
