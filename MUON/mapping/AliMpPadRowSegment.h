// $Id$
// ---------------------------------------------------------------
// Category: sector
//
// Class AliMpPadRowSegment
// ------------------------
// Class describing a pad row segment composed of the 
// the identic pads.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef M_PAD_ROW_SEGMENT_H
#define M_PAD_ROW_SEGMENT_H

#include <TObject.h>
#include <TVector2.h>

class AliMpPadRow;
class AliMpMotif;

class AliMpPadRowSegment : public TObject
{
  public:
    AliMpPadRowSegment(AliMpPadRow* padRow, AliMpMotif* motif, Int_t motifPositionId,
                   Int_t nofPads);
    AliMpPadRowSegment();
    virtual ~AliMpPadRowSegment();

    // methods
    virtual Double_t  LeftBorderX() const;
    virtual Double_t  RightBorderX() const;
    virtual Double_t  HalfSizeY() const;

    // get methods
    virtual AliMpPadRow*  GetPadRow() const;
    virtual AliMpMotif*   GetMotif() const;    
    virtual Int_t     GetMotifPositionId() const;
            Int_t     GetNofPads() const {return fNofPads;}     

    // set methods
    void  SetOffsetX(Double_t offsetX);  

  private:
    // methods
    Double_t  FirstPadCenterX() const;
    Double_t  LastPadCenterX() const;
    Double_t  FirstPadBorderX() const;
    Double_t  LastPadBorderX() const;

    // data members
    Int_t     fNofPads;  //number of pads
    Double_t  fOffsetX;  //the x position of the right border
    AliMpPadRow*  fPadRow;   //the pad row containing this segment 
    AliMpMotif*   fMotif;    //the motif 
    Int_t     fMotifPositionId;  // the motif position id
    
  ClassDef(AliMpPadRowSegment,1)  //Row segment
};

#endif //M_PAD_ROW_SEGMENT_H

