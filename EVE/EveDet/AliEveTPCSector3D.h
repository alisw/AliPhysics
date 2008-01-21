// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TPCSector3D_H
#define ALIEVE_TPCSector3D_H

#include <EveDet/AliEveTPCSectorViz.h>
#include <EveDet/AliEveTPCSectorData.h>

#include <TEveBoxSet.h>
#include <TEvePointSet.h>


class AliEveTPCSector3D : public AliEveTPCSectorViz
{
  friend class AliEveTPCSector3DEditor;
  friend class AliEveTPCSector3DGL;

protected:
  void LoadPadrow(AliEveTPCSectorData::RowIterator& iter,
                  Float_t sx, Float_t sy, Float_t pw, Float_t ph);
  void UpdateBoxes();
  void SetupPointSetArray();

  TEveBoxSet        fBoxSet;
  TEvePointSetArray fPointSetArray;
  Float_t             fPointFrac;
  Float_t             fPointSize;
  Bool_t              fPointSetOn;
  Int_t               fPointSetMaxVal;

  Float_t             fDriftVel;
  Float_t             fZStep;

public:
  AliEveTPCSector3D(const Text_t* n="AliEveTPCSector3D", const Text_t* t=0);
  virtual ~AliEveTPCSector3D();

  void SetPointFrac(Float_t f) { fPointFrac = f; IncRTS(); }
  void SetPointSize(Float_t s) { fPointSize = s; }

  void SetDriftVel(Float_t v) { fDriftVel = v; IncRTS(); }
  void SetZStep(Float_t step) { fZStep = step; IncRTS(); }

  virtual void SetRnrFrame(Bool_t rf);

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option="");

  ClassDef(AliEveTPCSector3D, 1);
}; // endclass AliEveTPCSector3D

#endif
