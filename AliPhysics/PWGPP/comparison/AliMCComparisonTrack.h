#ifndef AliMCComparisonTrack_h
#define AliMCComparisonTrack_h

//-------------------------------------------------------------------------
//
// A calss for keeping the MC track information used in  
// the comparison tasks by:  Andrei.Zalite@cern.ch
//
//-------------------------------------------------------------------------

#include "TObject.h"


class AliMCComparisonTrack: public TObject 
{
  public:
    AliMCComparisonTrack();
    virtual ~AliMCComparisonTrack() {}
    
    void SetLabel(Int_t Label)
      {fMCLabel = Label;}
    void SetPDG(Int_t PDG) 
      {fMCPdg = PDG;}
    void SetPz(Float_t PZ)
      {fPz = PZ;}
    void SetPt(Float_t Ptrans)
      {fPt = Ptrans;}
    void SetPhi(Float_t PhiAngle)
      {fPhi = PhiAngle;}
    void SetLocalX(Float_t localX)
      {fLocalX =localX;}
    void SetLocalY(Float_t localY)
      {fLocalY =localY;}  
    void SetZ(Float_t Zcoor)
      {fZ = Zcoor;}
      
    Int_t GetLabel() const
      {return fMCLabel;}  
    Int_t GetPDG() const
      {return fMCPdg;}
    Float_t GetPz() const
      {return fPz;}  
    Float_t GetPt() const
      {return fPt;}
    Float_t GetPhi() const
      {return fPhi;}
    Float_t GetLocalX() const
      {return fLocalX;}
    Float_t GetLocalY() const
      {return fLocalY;}      
    Float_t GetZ() const
      {return fZ;}
  
  private:
    Int_t fMCLabel; // track label
    Int_t fMCPdg;   // PDG particle code
    Float_t fPz;    // z-component of momentum
    Float_t fPt;    // component of momentum in the transverse plane
    Float_t fPhi;   // phi angle of the particle direction
    Float_t fLocalX; // local x-coordinate
    Float_t fLocalY; // local y-coordinate
    Float_t fZ;      // z-coordinate
    
    ClassDef(AliMCComparisonTrack, 2); 
};


#endif
