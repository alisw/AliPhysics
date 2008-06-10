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
    void SetTPCMomentum(Float_t Px, Float_t Py, Float_t Pz)
      {TPCPx = Px; TPCPy = Py; TPCPz = Pz;}
      
    Int_t GetLabel()
      {return MCLabel;}  
    Int_t GetPDG()
      {return MCPdg;}
    Float_t GetTPCPx()
      {return TPCPx;}
    Float_t GetTPCPy()
      {return TPCPy;}  
    Float_t GetTPCPz()
      {return TPCPz;}  
  
  private:
      // track label
    Int_t MCLabel;
      // particle PDG
    Int_t MCPdg;
      // momentum at the entrance of the TPC
    Float_t TPCPx;
    Float_t TPCPy;
    Float_t TPCPz;
    
    ClassDef(AliMCComparisonTrack, 1); 
};


#endif
