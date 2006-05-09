#ifndef ALIEVE_TPCSegment_H
#define ALIEVE_TPCSegment_H

#include <Reve/RenderElement.h>

#include <Alieve/TPCDigitsInfo.h>

#include <TNamed.h> 
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TGeometry.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TAtt3D.h>
#include <TAttBBox.h>


namespace Alieve {

  class TPCSegmentEditor;
  class TPCSegmentGL;

  class TPCSegment : public TNamed, public TAtt3D, public TAttBBox, public Reve::RenderElement
  {
    friend class TPCSegmentGL;
    friend class TPCSegmentEditor;
  private:
    void Init();

  protected:
    TPCDigitsInfo*      fInfo; 

    Bool_t              fRnrFrame;
    Bool_t              fUseTexture;
    Color_t             fFrameCol;

    // These change data representation:
    Int_t               fID; 
    Bool_t		fShowMax;
    Int_t               fMinTime;     
    Int_t               fMaxTime;
    Short_t             fTreshold;
    Int_t               fMaxVal;

    Double_t            fMatrix[16];
    Bool_t              fTrans;
    UInt_t              fRTS;       //! Rendering TimeStamp

  public:
    TPCSegment(const Text_t* n="TPCSegment", const Text_t* t=0, Color_t col=2) : 
      TNamed(n,t), Reve::RenderElement(fFrameCol), fFrameCol(col), fRTS(1)
    { Init(); }


    virtual ~TPCSegment();

    void SetInfo(TPCDigitsInfo* diginfo);
    virtual void SetSegmentID(Int_t id);

    void SetShowMax(Bool_t sm)  { fShowMax  = sm; ++fRTS; }
    void SetMinTime(Int_t mt)   { fMinTime  = mt; ++fRTS; }
    void SetMaxTime(Int_t mt)   { fMaxTime  = mt; ++fRTS; }
    void SetTreshold(Short_t t) { fTreshold =  t; ++fRTS; }
    void SetMaxVal(Int_t mv)    { fMaxVal   = mv; ++fRTS; }

    virtual void ComputeBBox();

    virtual void Paint(Option_t* option = "");

    virtual void SetTrans(Bool_t t);

    virtual Bool_t CanEditMainColor()  { return true; }
    ClassDef(TPCSegment, 1);
  }; // endclass TPCSegment

}

#endif
