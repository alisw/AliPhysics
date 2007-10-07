#include <string.h>
#include "AliMultiplicity.h"

ClassImp(AliMultiplicity)

//______________________________________________________________________
AliMultiplicity::AliMultiplicity():
  TObject(),
  fNtracks(0),
  fNsingle(0),
  fLabels(0),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fThsingle(0),
  fPhisingle(0)
{
  // Default Constructor
  fFiredChips[0] = 0;
  fFiredChips[1] = 0;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(Int_t ntr, Float_t *t,  Float_t *ph, Float_t *df, Int_t *labels, Int_t ns, Float_t *ts, Float_t *ps):
  TObject(),
  fNtracks(ntr),
  fNsingle(ns),
  fLabels(0),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
  fThsingle(0),
  fPhisingle(0)
{
// Standard constructor
  if(ntr>0){
    fLabels = new Int_t[ntr];
    fTh = new Double_t [ntr];
    fPhi = new Double_t [ntr];
    fDeltPhi = new Double_t [ntr];
    for(Int_t i=0;i<fNtracks;i++){
      fTh[i]=t[i];
      fPhi[i]=ph[i];
      fDeltPhi[i]=df[i];
      fLabels[i] = labels[i];
    }
  }
  if(ns>0){
    fThsingle = new Double_t [ns];
    fPhisingle = new Double_t [ns];
    for(Int_t i=0;i<fNsingle;i++){
      fThsingle[i]=ts[i];
      fPhisingle[i]=ps[i];
    }
  }
  fFiredChips[0] = 0;
  fFiredChips[1] = 0;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(const AliMultiplicity& m):
  TObject(m),
  fNtracks(m.fNtracks),
  fNsingle(m.fNsingle),
  fLabels(0),
  fTh(0),
  fPhi(0),
  fDeltPhi(0),
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
    fTh = new Double_t[fNtracks];
    fPhi = new Double_t[fNtracks];
    fDeltPhi = new Double_t[fNtracks];
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
    fThsingle = new Double_t[fNsingle];
    fPhisingle = new Double_t[fNsingle];
  }
  else {
    fThsingle = 0;
    fPhisingle = 0;
  }
  if(m.fTh)memcpy(fTh,m.fTh,fNtracks*sizeof(Double_t));
  if(m.fPhi)memcpy(fPhi,m.fPhi,fNtracks*sizeof(Double_t));
  if(m.fDeltPhi)memcpy(fDeltPhi,m.fDeltPhi,fNtracks*sizeof(Double_t));
  if(m.fLabels)memcpy(fLabels,m.fLabels,fNtracks*sizeof(Int_t));
  if(m.fThsingle)memcpy(fThsingle,m.fThsingle,fNsingle*sizeof(Double_t));
  if(m.fPhisingle)memcpy(fPhisingle,m.fPhisingle,fNsingle*sizeof(Double_t));

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


