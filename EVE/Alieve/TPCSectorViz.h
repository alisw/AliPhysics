// $Header$

#ifndef ALIEVE_TPCSectorViz_H
#define ALIEVE_TPCSectorViz_H

#include <Reve/RenderElement.h>

#include <TNamed.h> 
#include <TAtt3D.h>
#include <TAttBBox.h>


namespace Alieve {

class TPCData; class TPCSectorData;

class TPCSectorVizEditor;
class TPCSector2D;  class TPCSector2DEditor;  class TPCSector2DGL;
class TPCSector3D;  class TPCSector3DEditor;  class TPCSector3DGL;

class TPCSectorViz : public TNamed, public TAtt3D, public TAttBBox,
                     public Reve::RenderElement
{
  friend class TPCSectorVizEditor;
  friend class TPCSector2D;
  friend class TPCSector2DEditor;
  friend class TPCSector2DGL;
  friend class TPCSector3D;
  friend class TPCSector3DEditor;
  friend class TPCSector3DGL;

protected:
  TPCData*    fTPCData; 
  Int_t       fSectorID;

  Int_t       fMinTime;     
  Int_t       fMaxTime;
  Short_t     fThreshold;
  Int_t       fMaxVal;

  Bool_t      fRnrInn;
  Bool_t      fRnrOut1;
  Bool_t      fRnrOut2;

  Color_t     fFrameColor;
  Bool_t      fRnrFrame;
  Bool_t      fTrans;
  Double_t    fMatrix[16];
  UInt_t      fRTS;       //! Rendering TimeStamp

  void SetupColor(Int_t val, UChar_t* pix) const;

public:
  TPCSectorViz(const Text_t* n="TPCSectorViz", const Text_t* t=0);
  virtual ~TPCSectorViz();

  virtual UInt_t IncRTS()           { return ++fRTS; }
  virtual Bool_t CanEditMainColor() { return true; }

  void SetDataSource(TPCData* data);
  void SetSectorID(Int_t id);

  TPCData*       GetData()     const { return fTPCData; }
  Int_t          GetSectorID() const { return fSectorID; }
  TPCSectorData* GetSectorData() const;

  void SetMinTime(Int_t mt)    { fMinTime   = mt; IncRTS(); }
  void SetMaxTime(Int_t mt)    { fMaxTime   = mt; IncRTS(); }
  void SetThreshold(Short_t t) { fThreshold =  t; IncRTS(); }
  void SetMaxVal(Int_t mv)     { fMaxVal    = mv; IncRTS(); }

  void SetRnrInn(Bool_t r)     { fRnrInn  = r; IncRTS(); }
  void SetRnrOut1(Bool_t r)    { fRnrOut1 = r; IncRTS(); }
  void SetRnrOut2(Bool_t r)    { fRnrOut2 = r; IncRTS(); }

  void SetFrameColor(Color_t col)     { fFrameColor = col; IncRTS(); }
  virtual void SetRnrFrame(Bool_t rf) { fRnrFrame = rf;  IncRTS(); }
  void SetTrans(Bool_t t);

  ClassDef(TPCSectorViz, 1); // Base-class for TPC raw-data visualization
}; // endclass TPCSectorViz

}

#endif
