#ifndef AliMCComparisonTrack_cxx
#define AliMCComparisonTrack_cxx

#include "TObject.h"


class AliMCComparisonTrack: public TObject 
{
  public:
    AliMCComparisonTrack();
    virtual ~AliMCComparisonTrack() {}
    
    void SetLabel(Int_t Label)
      {MCLabel = Label;}
    void SetPDG(Int_t PDG) 
      {MCPdg = PDG;}
    void SetPz(Float_t PZ)
      {Pz = PZ;}
    void SetPt(Float_t Ptrans)
      {Pt = Ptrans;}
    void SetPhi(Float_t PhiAngle)
      {Phi = PhiAngle;}
    void SetLocalX(Float_t localX)
      {LocalX =localX;}
    void SetLocalY(Float_t localY)
      {LocalY =localY;}  
    void SetZ(Float_t Zcoor)
      {Z = Zcoor;}
      
    Int_t GetLabel()
      {return MCLabel;}  
    Int_t GetPDG()
      {return MCPdg;}
    Float_t GetPz()
      {return Pz;}  
    Float_t GetPt()
      {return Pt;}
    Float_t GetPhi()
      {return Phi;}
    Float_t GetLocalX()
      {return LocalX;}
    Float_t GetLocalY()
      {return LocalY;}      
    Float_t GetZ()
      {return Z;}
  
  private:
      // track label
    Int_t MCLabel;
      // PDG particle code
    Int_t MCPdg;
      // z-component of momentum
    Float_t Pz;
      // component of momentum in the transverse plane
    Float_t Pt; 
      // phi angle of the particle direction
    Float_t Phi;
      // local x-coordinate
    Float_t LocalX;
      // local y-coordinate
    Float_t LocalY;
      // z-coordinate
    Float_t Z;
    
    ClassDef(AliMCComparisonTrack, 1); 
};


#endif
