// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSector3D_H
#define AliEveTPCSector3D_H

#include <AliEveTPCSectorViz.h>
#include <AliEveTPCSectorData.h>

#include <TEveBoxSet.h>
#include <TEvePointSet.h>

//------------------------------------------------------------------------------
// AliEveTPCSector3D
//
// Visualization of TPC raw-data in 3D.

class AliEveTPCSector3D : public AliEveTPCSectorViz
{
  friend class AliEveTPCSector3DEditor;
  friend class AliEveTPCSector3DGL;

public:
  AliEveTPCSector3D(const Text_t* n="AliEveTPCSector3D", const Text_t* t=0);
  virtual ~AliEveTPCSector3D() {}

  void SetPointFrac(Float_t f) { fPointFrac = f; IncRTS(); }
  void SetPointSize(Float_t s) { fPointSize = s; }

  void SetDriftVel(Float_t v) { fDriftVel = v; IncRTS(); }
  void SetZStep(Float_t step) { fZStep = step; IncRTS(); }

  void SetCutOnEta(Bool_t cut) { fCutOnEta = cut; IncRTS(); }
  void SetEtaMax(Float_t max) { fEtaMax = max; IncRTS(); }
  void SetEtaMin(Float_t min) { fEtaMin = min; IncRTS(); }

  Float_t GetCutOnEta() { return fCutOnEta; }
  Float_t GetEtaMax() { return fEtaMax; }
  Float_t GetEtaMin() { return fEtaMin; }

  virtual void SetRnrFrame(Bool_t rf);

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option="");

protected:
  void LoadPadrow(AliEveTPCSectorData::RowIterator& iter,
                  Float_t sx, Float_t sy, Float_t pw, Float_t ph);
  void UpdateBoxesAndPoints();
  void SetupPointSetArray();

  TEveBoxSet          fBoxSet;          // BoxSet used to display digits as boxes.
  TEvePointSetArray   fPointSetArray;   // PointSet used to display digits as points.
  Float_t             fPointFrac;       // Fraction of signal range shown as points.
  Float_t             fPointSize;       // Point size.
  Bool_t              fPointSetOn;      // PointSet initialized.
  Int_t               fPointSetMaxVal;  // Maximum signal value for data in pointset.

  Float_t             fDriftVel;        // Drift velocity for 'z' coordinate.
  Float_t             fZStep;           // Z width of a time-bin.

  Bool_t              fCutOnEta;        // Bool used to apply eta cut
  Float_t             fEtaMin;          // Minimum eta value for eta cut
  Float_t             fEtaMax;          // Maximum eta value for eta cut

  ClassDef(AliEveTPCSector3D, 0); // Visualization of TPC raw-data in 3D.
};

#endif
