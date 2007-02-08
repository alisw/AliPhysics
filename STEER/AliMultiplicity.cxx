#include "AliMultiplicity.h"

ClassImp(AliMultiplicity)

//______________________________________________________________________
AliMultiplicity::AliMultiplicity():
  TObject(),
  fNtracks(0),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fNsingle(0),
  fThsingle(0),
  fPhisingle(0)

{
  // Default Constructor
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(Int_t ntr, Float_t *t,  Float_t *ph, Float_t *df, Int_t ns, Float_t *ts, Float_t *ps):
  TObject(),
  fNtracks(ntr),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fNsingle(ns),
  fThsingle(0),
  fPhisingle(0)
{
// Standard constructor
  if(ntr>0){
    fTh = new Float_t [ntr];
    fPhi = new Float_t [ntr];
    fDeltPhi = new Float_t [ntr];
    for(Int_t i=0;i<fNtracks;i++){
      fTh[i]=t[i];
      fPhi[i]=ph[i];
      fDeltPhi[i]=df[i];
    }
  }
  if(ns>0){
    fThsingle = new Float_t [ns];
    fPhisingle = new Float_t [ns];
    for(Int_t i=0;i<fNsingle;i++){
      fThsingle[i]=ts[i];
      fPhisingle[i]=ps[i];
    }
  }
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(const AliMultiplicity& m):
  TObject(m),
  fNtracks(m.fNtracks),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fNsingle(m.fNsingle),
  fThsingle(0),
  fPhisingle(0)
{
  // copy constructor

  Duplicate(m);

}

//______________________________________________________________________
AliMultiplicity &AliMultiplicity::operator=(const AliMultiplicity& m){
  // assignment operator
  if(this == &m)return *this;
  ((TObject *)this)->operator=(m);

  if (fTh)delete [] fTh;
  if(fPhi)delete [] fPhi;
  if(fDeltPhi)delete [] fDeltPhi;
  if(fThsingle)delete [] fThsingle;
  if(fPhisingle)delete [] fPhisingle;
  Duplicate(m);

  return *this;
}

//______________________________________________________________________
void AliMultiplicity::Duplicate(const AliMultiplicity& m){
  // used by copy constructor and assignment operator
  fNtracks = m.fNtracks;
  if(fNtracks>0){
    fTh = new Float_t [fNtracks];
    fPhi = new Float_t [fNtracks];
    fDeltPhi = new Float_t [fNtracks];
  }
  else {
    fTh = 0;
    fPhi = 0;
    fDeltPhi = 0;
  }
  fNsingle = m.fNsingle;
  if(fNsingle>0){
    fThsingle = new Float_t [fNsingle];
    fPhisingle = new Float_t [fNsingle];
  }
  else {
    fThsingle = 0;
    fPhisingle = 0;
  }
  memcpy(fTh,m.fTh,fNtracks*sizeof(Float_t));
  memcpy(fPhi,m.fPhi,fNtracks*sizeof(Float_t));
  memcpy(fDeltPhi,m.fDeltPhi,fNtracks*sizeof(Float_t));
  memcpy(fThsingle,m.fThsingle,fNsingle*sizeof(Float_t));
  memcpy(fPhisingle,m.fPhisingle,fNsingle*sizeof(Float_t));
}

//______________________________________________________________________
AliMultiplicity::~AliMultiplicity(){
  // Destructor
  if (fTh)delete [] fTh;
  if(fPhi)delete [] fPhi;
  if(fDeltPhi)delete [] fDeltPhi;
  if(fThsingle)delete [] fThsingle;
  if(fPhisingle)delete [] fPhisingle;
}


