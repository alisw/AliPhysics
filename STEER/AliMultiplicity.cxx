#include <string.h>
#include "AliMultiplicity.h"
#include "AliLog.h"

ClassImp(AliMultiplicity)

//______________________________________________________________________
AliMultiplicity::AliMultiplicity():
  TObject(),
  fNtracks(0),
  fNsingle(0),
  fLabels(0),
  fLabelsL2(0),
  fTh(0),
  fPhi(0),
  fDeltTh(0),
  fDeltPhi(0),
  fThsingle(0),
  fPhisingle(0),
  fFastOrFiredChips(1200),
  fClusterFiredChips(1200)
{
  // Default Constructor
  fFiredChips[0] = 0;
  fFiredChips[1] = 0;
  for(Int_t ilayer = 0; ilayer < 6; ilayer++)fITSClusters[ilayer] = 0;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(Int_t ntr, Float_t *th,  Float_t *ph, Float_t *dth, Float_t *dph, Int_t *labels, Int_t* labelsL2, Int_t ns, Float_t *ts, Float_t *ps, Short_t nfcL1, Short_t nfcL2, const TBits & fFastOr):
  TObject(),
  fNtracks(ntr),
  fNsingle(ns),
  fLabels(0),
  fLabelsL2(0),
  fTh(0),
  fPhi(0),
  fDeltTh(0),
  fDeltPhi(0),
  fThsingle(0),
  fPhisingle(0),
  fFastOrFiredChips(1200),
  fClusterFiredChips(1200)
{
// Standard constructor
  if(ntr>0){
    fLabels = new Int_t[ntr];
    fLabelsL2 = new Int_t[ntr];
    fTh = new Double_t [ntr];
    fPhi = new Double_t [ntr];
    fDeltTh = new Double_t [ntr];
    fDeltPhi = new Double_t [ntr];
    for(Int_t i=0;i<fNtracks;i++){
      fTh[i]=th[i];
      fPhi[i]=ph[i];
      fDeltTh[i]=dth[i];
      fDeltPhi[i]=dph[i];
      fLabels[i] = labels[i];
      fLabelsL2[i] = labelsL2[i];
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
  fFiredChips[0] = nfcL1;
  fFiredChips[1] = nfcL2;
  fFastOrFiredChips = fFastOr;
  for(Int_t ilayer = 0; ilayer < 6; ilayer++)fITSClusters[ilayer] = 0;
}

//______________________________________________________________________
AliMultiplicity::AliMultiplicity(const AliMultiplicity& m):
  TObject(m),
  fNtracks(m.fNtracks),
  fNsingle(m.fNsingle),
  fLabels(0),
  fLabelsL2(0),
  fTh(0),
  fPhi(0),
  fDeltTh(0),
  fDeltPhi(0),
  fThsingle(0),
  fPhisingle(0),
  fFastOrFiredChips(1200),
  fClusterFiredChips(1200)
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
  if(fDeltTh)delete [] fDeltTh;fDeltTh= 0; 
  if(fDeltPhi)delete [] fDeltPhi;fDeltPhi = 0; 
  if(fLabels)delete [] fLabels;fLabels = 0;
  if(fLabelsL2)delete [] fLabelsL2;fLabelsL2 = 0;
  if(fThsingle)delete [] fThsingle;fThsingle = 0;
  if(fPhisingle)delete [] fPhisingle;fPhisingle = 0;
  Duplicate(m);

  return *this;
}

void AliMultiplicity::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliMultiplicity *robj = dynamic_cast<AliMultiplicity*>(&obj);
  if(!robj)return; // not an AliMultiplicity
  *robj = *this;

}


//______________________________________________________________________
void AliMultiplicity::Duplicate(const AliMultiplicity& m){
  // used by copy constructor and assignment operator
  fNtracks = m.fNtracks;
  if(fNtracks>0){
    fTh = new Double_t[fNtracks];
    fPhi = new Double_t[fNtracks];
    fDeltTh = new Double_t[fNtracks];
    fDeltPhi = new Double_t[fNtracks];
    fLabels = new Int_t[fNtracks];
    fLabelsL2 = new Int_t[fNtracks];
  }
  else {
    fTh = 0;
    fPhi = 0;
    fDeltTh = 0;
    fDeltPhi = 0;
    fLabels = 0;
    fLabelsL2 = 0;
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
  if(m.fDeltTh)memcpy(fDeltTh,m.fDeltTh,fNtracks*sizeof(Double_t));
  if(m.fDeltPhi)memcpy(fDeltPhi,m.fDeltPhi,fNtracks*sizeof(Double_t));
  if(m.fLabels)memcpy(fLabels,m.fLabels,fNtracks*sizeof(Int_t));
  if(m.fLabelsL2)memcpy(fLabelsL2,m.fLabelsL2,fNtracks*sizeof(Int_t));
  if(m.fThsingle)memcpy(fThsingle,m.fThsingle,fNsingle*sizeof(Double_t));
  if(m.fPhisingle)memcpy(fPhisingle,m.fPhisingle,fNsingle*sizeof(Double_t));

  fFiredChips[0] = m.fFiredChips[0];
  fFiredChips[1] = m.fFiredChips[1];
  for(Int_t ilayer = 0; ilayer < 6; ilayer++){
    fITSClusters[ilayer] = m.fITSClusters[ilayer];
  }

  

  fFastOrFiredChips = m.fFastOrFiredChips;
  fClusterFiredChips = m.fClusterFiredChips;
}

//______________________________________________________________________
AliMultiplicity::~AliMultiplicity(){
  // Destructor
  if(fTh)delete [] fTh;fTh = 0;
  if(fPhi)delete [] fPhi;fPhi = 0; 
  if(fDeltTh)delete [] fDeltTh;fDeltTh = 0; 
  if(fDeltPhi)delete [] fDeltPhi;fDeltPhi = 0; 
  if(fLabels)delete [] fLabels;fLabels = 0;
  if(fLabelsL2)delete [] fLabelsL2;fLabelsL2 = 0;
  if(fThsingle)delete [] fThsingle;fThsingle = 0;
  if(fPhisingle)delete [] fPhisingle;fPhisingle = 0;

}

//______________________________________________________________________
void AliMultiplicity::SetLabel(Int_t i, Int_t layer, Int_t label)
{
    if(i>=0 && i<fNtracks) {
	if (layer == 0) {
	    fLabels[i] = label;
	    return;
	} else if (layer == 1) {
	    if (fLabelsL2) {
		fLabelsL2[i] = label;
		return;
	    }
	}
    }
    Error("SetLabel","Invalid track number %d or layer %d",i,layer);
}

//______________________________________________________________________
UInt_t AliMultiplicity::GetNumberOfITSClusters(Int_t layMin, Int_t layMax) const {

  if(layMax < layMin) {
    AliError("layer min > layer max");
    return 0;
  }
  if(layMin < 0) {
    AliError("layer min < 0");
    return 0;
  }
  if(layMax < 0) {
    AliError("layer max > 0");
    return 0;
  }

  Int_t sum=0; 
  for (Int_t i=layMin; i<=layMax; i++) sum+=fITSClusters[i]; 
  return sum; 

}
