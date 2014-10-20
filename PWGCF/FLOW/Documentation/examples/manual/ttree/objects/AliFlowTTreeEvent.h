#ifndef AliFlowTTreeEVENT_H
#define AliFlowTTreeEVENT_H

#include <TObject.h>

class AliFlowTTreeEvent : public TObject{

 public:

  AliFlowTTreeEvent();
  
  void  SetRun(Int_t run)       {fRun   = run;}
  void  SetV0M(Float_t V0M)     {fV0M   = V0M;}
  void  SetTRK(Float_t TRK)     {fTRK   = TRK;}
  void  SetZvtx(Float_t Zvtx)   {fZvtx  = Zvtx;}

  Int_t         GetRun() const  {return fRun;}
  Float_t       GetV0M() const  {return fV0M;}
  Float_t       GetTRK() const  {return fTRK;}
  Float_t       GetZvtx() const {return fZvtx;}

 private:
  Int_t     fRun;        // run number
  Float_t   fV0M;        // centrality V0
  Float_t   fTRK;        // centrality TRK
  Float_t   fZvtx;       // rec vertex

  virtual ~AliFlowTTreeEvent(); // default destructor

  ClassDef(AliFlowTTreeEvent, 1);    // Help class
};

#endif
	
