// $Header$

#ifndef ALIEVE_TPCSector3D_H
#define ALIEVE_TPCSector3D_H

#include <Alieve/TPCSectorViz.h>
#include <Alieve/TPCSectorData.h>
#include <Reve/BoxSet.h>


namespace Alieve {

class TPCSector3D: public TPCSectorViz
{
  friend class TPCSector3DEditor;

protected:
  void MakeFrameBox();
  void LoadPadrow(TPCSectorData::RowIterator& iter,
                  Float_t sx, Float_t sy, Float_t pw, Float_t ph);
  void UpdateBoxes();


  Reve::BoxSet fBoxSet;
  Reve::Box    fFrameBox;
  Float_t      fDriftVel;
  Float_t      fZStep;

public:
  TPCSector3D(const Text_t* n="TPCSector3D", const Text_t* t=0);
  virtual ~TPCSector3D();

  virtual UInt_t IncRTS();

  void SetDriftVel(Float_t v) { fDriftVel = v; IncRTS(); }
  void SetZStep(Float_t step) { fZStep = step; IncRTS(); }

  virtual void SetRnrFrame(Bool_t rf);

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option="");

  ClassDef(TPCSector3D, 1);
}; // endclass TPCSector3D

}

#endif
