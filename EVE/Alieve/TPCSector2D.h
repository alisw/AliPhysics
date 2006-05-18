// $Header$

#ifndef ALIEVE_TPCSector2D_H
#define ALIEVE_TPCSector2D_H

#include <Reve/RenderElement.h>

#include <TNamed.h> 
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TGeometry.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TAtt3D.h>
#include <TAttBBox.h>


namespace Alieve {

class TPCData;

class TPCSector2DEditor;
class TPCSector2DGL;

class TPCSector2D : public TNamed, public TAtt3D, public TAttBBox, public Reve::RenderElement
{
  friend class TPCSector2DGL;
  friend class TPCSector2DEditor;

private:
  void Init();

protected:
  TPCData*    fTPCData; 

  // These change data representation:
  Int_t       fSectorID; 
  Bool_t      fShowMax;
  Bool_t      fAverage;
  Int_t       fMinTime;     
  Int_t       fMaxTime;
  Short_t     fThreshold;
  Int_t       fMaxVal;

  Bool_t      fRnrInn;
  Bool_t      fRnrOut1;
  Bool_t      fRnrOut2;

  Bool_t      fRnrFrame;
  Bool_t      fUseTexture;
  Color_t     fFrameCol;

  Double_t    fMatrix[16];
  Bool_t      fTrans;
  UInt_t      fRTS;       //! Rendering TimeStamp

public:
  TPCSector2D(const Text_t* n="TPCSector2D", const Text_t* t=0, Color_t col=2) : 
    TNamed(n,t), Reve::RenderElement(fFrameCol), fFrameCol(col), fRTS(1)
  { Init(); }

  virtual ~TPCSector2D();

  void SetDataSource(TPCData* data);
  virtual void SetSectorID(Int_t id);

  void SetRnrInn(Bool_t r)  { fRnrInn = r; ++fRTS; }
  void SetRnrOut1(Bool_t r) { fRnrOut1 = r; ++fRTS; }
  void SetRnrOut2(Bool_t r) { fRnrOut2 = r; ++fRTS; }

  void SetShowMax(Bool_t sm)  { fShowMax  = sm; ++fRTS; }
  void SetAverage(Bool_t avg) { fAverage  = avg; ++fRTS; }
  void SetMinTime(Int_t mt)   { fMinTime  = mt; ++fRTS; }
  void SetMaxTime(Int_t mt)   { fMaxTime  = mt; ++fRTS; }
  void SetThreshold(Short_t t) { fThreshold =  t; ++fRTS; }
  void SetMaxVal(Int_t mv)    { fMaxVal   = mv; ++fRTS; }

  virtual void ComputeBBox();

  virtual void Paint(Option_t* option="");

  virtual void SetTrans(Bool_t t);

  virtual Bool_t CanEditMainColor()  { return true; }

  ClassDef(TPCSector2D, 1);
}; // endclass TPCSector2D

}

#endif
