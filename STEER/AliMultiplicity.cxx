#include <string.h>
#include "AliMultiplicity.h"

ClassImp(AliMultiplicity)

//______________________________________________________________________
AliMultiplicity::AliMultiplicity():
  TObject(),
  fNtracks(0),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fLabels(0),
  fNsingle(0),
  fThsingle(0),
  fPhisingle(0)

{
  // Default Constructor
  fFiredChips[0] = -1;
  fFiredChips[1] = -1;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(Int_t ntr, Float_t *t,  Float_t *ph, Float_t *df, Int_t *labels, Int_t ns, Float_t *ts, Float_t *ps):
  TObject(),
  fNtracks(ntr),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fLabels(0),
  fNsingle(ns),
  fThsingle(0),
  fPhisingle(0)
{
// Standard constructor
  if(ntr>0){
    fTh = new Float_t [ntr];
    fPhi = new Float_t [ntr];
    fDeltPhi = new Float_t [ntr];
    fLabels = new Int_t[ntr];
    for(Int_t i=0;i<fNtracks;i++){
      fTh[i]=t[i];
      fPhi[i]=ph[i];
      fDeltPhi[i]=df[i];
      fLabels[i] = labels[i];
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
  fFiredChips[0] = -1;
  fFiredChips[1] = -1;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(const AliMultiplicity& m):
  TObject(m),
  fNtracks(m.fNtracks),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fLabels(0),
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

  if(fTh)delete [] fTh;fTh = 0;
  if(fPhi)delete [] fPhi;fPhi = 0; 
  if(fDeltPhi)delete [] fDeltPhi;fDeltPhi = 0; 
  if(fLabels)delete [] fLabels;fLabels = 0;
  if(fThsingle)delete [] fThsingle;fThsingle = 0;
  if(fPhisingle)delete [] fPhisingle;fPhisingle = 0;
  Duplicate(m);

  return *this;
}

//______________________________________________________________________
void AliMultiplicity::Duplicate(const AliMultiplicity& m){
  // used by copy constructor and assignment operator
  fNtracks = m.fNtracks;
  if(fNtracks>0){
    fTh = new Float_t[fNtracks];
    fPhi = new Float_t[fNtracks];
    fDeltPhi = new Float_t[fNtracks];
    fLabels = new Int_t[fNtracks];
  }
  else {
    fTh = 0;
    fPhi = 0;
    fDeltPhi = 0;
    fLabels = 0;
  }
  fNsingle = m.fNsingle;
  if(fNsingle>0){
    fThsingle = new Float_t[fNsingle];
    fPhisingle = new Float_t[fNsingle];
  }
  else {
    fThsingle = 0;
    fPhisingle = 0;
  }
  memcpy(fTh,m.fTh,fNtracks*sizeof(Float_t));
  memcpy(fPhi,m.fPhi,fNtracks*sizeof(Float_t));
  memcpy(fDeltPhi,m.fDeltPhi,fNtracks*sizeof(Float_t));
  if (m.fLabels)   memcpy(fLabels,m.fLabels,fNtracks*sizeof(Int_t));
  if (m.fThsingle) memcpy(fThsingle,m.fThsingle,fNsingle*sizeof(Float_t));
  if (m.fPhisingle) memcpy(fPhisingle,m.fPhisingle,fNsingle*sizeof(Float_t));

  fFiredChips[0] = m.fFiredChips[0];
  fFiredChips[1] = m.fFiredChips[1];
}

//______________________________________________________________________
AliMultiplicity::~AliMultiplicity(){
  // Destructor
  if(fTh)delete [] fTh;fTh = 0;
  if(fPhi)delete [] fPhi;fPhi = 0; 
  if(fDeltPhi)delete [] fDeltPhi;fDeltPhi = 0; 
  if(fLabels)delete [] fLabels;fLabels = 0;
  if(fThsingle)delete [] fThsingle;fThsingle = 0;
  if(fPhisingle)delete [] fPhisingle;fPhisingle = 0;

}


