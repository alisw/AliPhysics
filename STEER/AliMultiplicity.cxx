#include "AliMultiplicity.h"

ClassImp(AliMultiplicity)

//______________________________________________________________________
AliMultiplicity::AliMultiplicity():TObject() {
  // Default Constructor
  fNtracks = 0;
  fTh = 0;
  fPhi = 0;
  fDeltPhi = 0;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(Int_t ntr, Float_t *t,  Float_t *ph, Float_t *df):TObject() {
// Standard constructor
  fNtracks = ntr;
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
  else {
    fTh = 0;
    fPhi = 0;
    fDeltPhi = 0;
  }
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(const AliMultiplicity& m):TObject(m){
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
  memcpy(fTh,m.fTh,fNtracks*sizeof(Float_t));
  memcpy(fPhi,m.fPhi,fNtracks*sizeof(Float_t));
  memcpy(fDeltPhi,m.fDeltPhi,fNtracks*sizeof(Float_t));
}

//______________________________________________________________________
AliMultiplicity::~AliMultiplicity(){
  // Destructor
  if (fTh)delete [] fTh;
  if(fPhi)delete [] fPhi;
  if(fDeltPhi)delete [] fDeltPhi;
}


