#include "ITSModule.h"

#include <AliITSdigitSPD.h>
#include <AliITSdigitSDD.h>
#include <AliITSdigitSSD.h>

#include <TStyle.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

Short_t ITSModule::fgSDDThreshold  = 5;
Short_t ITSModule::fgSDDMaxVal     = 80;
Short_t ITSModule::fgSSDThreshold  = 2;
Short_t ITSModule::fgSSDMaxVal     = 100;

ClassImp(ITSModule)

/**************************************************************************/

ITSModule::ITSModule(const Text_t* n, const Text_t* t, Color_t col) :
  Reve::RenderElement(fFrameColor),
  QuadSet(n, t),
  fInfo(0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0),
  fFrameColor(col)
{}

ITSModule::ITSModule(Int_t id, ITSDigitsInfo* info, Color_t col) :
  Reve::RenderElement(fFrameColor),
  QuadSet(Form("ITS module %d", id)),
  fInfo  (0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0),
  fFrameColor(col)
{
  SetDigitsInfo(info);
  SetID(id);
}

ITSModule::~ITSModule()
{
  if(fInfo) fInfo->DecRefCount();
}

/**************************************************************************/

void ITSModule::SetMainColor(Color_t col)
{
  Reve::RenderElement::SetMainColor(col);
  if(!fQuads.empty()) {
    fQuads.front().ColorFromIdx(col);
  }
}

/**************************************************************************/

void ITSModule::SetDigitsInfo(ITSDigitsInfo* info)
{
  if(fInfo) fInfo->DecRefCount();
  fInfo = info;
  if(fInfo) fInfo->IncRefCount();
}

void ITSModule::SetID(Int_t id)
{
  static const Exc_t eH("ITSModule::SetID ");

  if(fInfo == 0)
    throw(eH + "ITSDigitsInfo not set.");

  if (id < fInfo->fGeom->GetStartSPD() || id > fInfo->fGeom->GetLastSSD())
    throw(eH + Form("%d is not valid. ID range from %d to %d", id,
		     fInfo->fGeom->GetStartSPD(), fInfo->fGeom->GetLastSSD()));

  fID = id;
  InitModule();
}

/**************************************************************************/

void ITSModule::InitModule()
{
  fInfo->fGeom->GetModuleId(fID,fLayer,fLadder,fDet);
  TString strLadder = "Ladder";
  TString strSensor = "Sensor";
  TString symname;
  Int_t id, nsector, nstave, nladder, rest;
  
  if (fID <= fInfo->fGeom->GetLastSPD()) {
    symname+=strLadder;
    if (fID<80) {
      nsector = fID/8;
      rest=fID-nsector*8;
      nstave=1;
      if (rest<4) nstave=0;
      rest-=nstave*4;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 0;
      fDx = fInfo->fSegSPD->Dx()*0.00005;
      fDz = 3.48; 
      fDy = fInfo->fSegSPD->Dy()*0.00005;
    } else {
      id=fID-80;
      nsector = id/8;
      rest=id-nsector*8;
      nstave=1;
      if (rest<4) nstave=0;
      rest-=nstave*4;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 0;
      fDx = fInfo->fSegSPD->Dx()*0.00005;
      fDz = 3.48; 
      fDy = fInfo->fSegSPD->Dy()*0.00005;
      fDetID = 0;
      fDx = fInfo->fSegSPD->Dx()*0.00005;
      fDz = 3.48; 
      fDy = fInfo->fSegSPD->Dy()*0.00005;
    }
  }
  else if (fID <= fInfo->fGeom->GetLastSDD()) {
    symname+=strSensor;
    if (fID<324) {
      id = fID-240;
      nladder = id/6;
      rest=id-nladder*6;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 1;
      fDx = fInfo->fSegSDD->Dx()*0.0001;
      fDz = fInfo->fSegSDD->Dz()*0.00005;
      fDy = fInfo->fSegSDD->Dy()*0.00005;
    } else {
      id = fID-324;
      nladder = id/8;
      rest=id-nladder*8;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 1;
      fDx = fInfo->fSegSDD->Dx()*0.0001;
      fDz = fInfo->fSegSDD->Dz()*0.00005;
      fDy = fInfo->fSegSDD->Dy()*0.00005;
    }
  }
  else {
    symname+=strSensor;
    if (fID<1248) {
      id = fID-500;
      nladder = id/22;
      rest=id-nladder*22;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 2;
      fInfo->fSegSSD->SetLayer(fLayer);  
      fDx = fInfo->fSegSSD->Dx()*0.00005;
      fDz = fInfo->fSegSSD->Dz()*0.00005;
      fDy = fInfo->fSegSSD->Dy()*0.00005;
    } else {
      id = fID-1248;
      nladder = id/25;
      rest=id-nladder*25;
      symname+=rest;
      SetName(symname.Data());
      fDetID = 2;
      fInfo->fSegSSD->SetLayer(fLayer);  
      fDx = fInfo->fSegSSD->Dx()*0.00005;
      fDz = fInfo->fSegSSD->Dz()*0.00005;
      fDy = fInfo->fSegSSD->Dy()*0.00005;
    }
  }

  LoadQuads();  
  ComputeBBox();
  SetTrans();
}

void ITSModule::LoadQuads()
{
  // printf("its module load quads \n");
  Float_t x = fDx;
  Float_t z = fDz;
  Bool_t aboveThreshold = false;

  // Module frame in xy plane
  fQuads.push_back(Reve::Quad(fFrameColor));
  Float_t dy = -0.;
  Float_t* p = fQuads.back().vertices;
  p[0] = -x;  p[1] =  dy; p[2]  = -z;
  p[3] = -x;  p[4] =  dy; p[5]  =  z;
  p[6] =  x;  p[7] =  dy; p[8]  =  z;
  p[9] =  x;  p[10] = dy; p[11] = -z;

  // Digits
  TClonesArray *digits;
  Int_t ndigits;
  Float_t dpx,dpz; 
  Int_t i,j;
  digits  = fInfo->GetDigits(fID, fDetID );
  ndigits = digits->GetEntriesFast(); 
  Int_t n_col = gStyle->GetNumberOfColors();

  switch(fDetID) {

  case 0: { // SPD
    aboveThreshold = true;
    AliITSsegmentationSPD* seg =  fInfo->fSegSPD; 
    AliITSdigitSPD *d=0;

    for (Int_t k=0; k<ndigits; k++) {
      d=(AliITSdigitSPD*)digits->UncheckedAt(k);
      j = d->GetCoord1();
      i = d->GetCoord2();
      x  = -seg->Dx()/2 + seg->Dpx(0) *i;
      x *=  0.0001;
      fInfo->GetSPDLocalZ(j,z);
      dpx = seg->Dpx(i)*0.0001;
      dpz = seg->Dpz(j)*0.0001;

      fQuads.push_back(Reve::Quad(7));
      Float_t* p = fQuads.back().vertices;
      p[0] = x;        p[1] = 0.; p[2]  = z;
      p[3] = x;        p[4] = 0.; p[5]  = z + dpz;
      p[6] = x + dpx;  p[7] = 0.; p[8]  = z + dpz;
      p[9] = x + dpx;  p[10] =0.; p[11] = z;
    }
    break;
  }

  case 1: { // SDD
    AliITSsegmentationSDD* seg =  fInfo->fSegSDD; 
    AliITSdigitSDD *d=0;
    x = 2*fDx;
    z = 2*fDz;
    for (Int_t k=0; k<ndigits; k++) {
      d=(AliITSdigitSDD*)digits->UncheckedAt(k);

      if (d->GetSignal() > fgSDDThreshold) {
	j = d->GetCoord1();
	i = d->GetCoord2();
	aboveThreshold = true;
	seg->DetToLocal(i,j,x,z);
	dpx = seg->Dpx(i)*0.0001;
	dpz = seg->Dpz(j)*0.0001;

	Int_t ci = gStyle->GetColorPalette
	  (TMath::Min(n_col - 1,
		      (n_col*(d->GetSignal() - fgSDDThreshold))/(fgSDDMaxVal - fgSDDThreshold)));
	fQuads.push_back(Reve::Quad(ci, p));
	Float_t* p = fQuads.back().vertices;
	p[0] = x;        p[1] = 0.; p[2]  = z;
	p[3] = x;        p[4] = 0.; p[5]  = z + dpz;
	p[6] = x + dpx;  p[7] = 0.; p[8]  = z + dpz;
	p[9] = x + dpx;  p[10] =0.; p[11] = z;
      }
    }
    break;
  }

  case 2: { // SSD
    AliITSsegmentationSSD* seg = fInfo->fSegSSD; 
    AliITSdigitSSD *d=0;
    Float_t ap,an,a;
    seg->Angles(ap,an);
    for (Int_t k=0; k<ndigits; k++) {
      d=(AliITSdigitSSD*)digits->UncheckedAt(k);
      if(d->GetSignal() > fgSSDThreshold){
	aboveThreshold = true;
	j = d->GetCoord1();
	i = d->GetCoord2();
	seg->DetToLocal(i,j,x,z);

	if( d->GetCoord1() == 1) {
	  a = ap;
	}
	else {
	  a = -an;
	}
     	fQuads.push_back(Reve::Quad());
	Int_t ci = gStyle->GetColorPalette
	  (TMath::Min(n_col - 1,
		      (n_col*(d->GetSignal() - fgSSDThreshold))/(fgSSDMaxVal - fgSSDThreshold)));

	fQuads.back().ColorFromIdx(ci);
	Float_t* p = fQuads.back().vertices;
        
	p[0] = x-TMath::Tan(a)*fDz;  p[1] =  0; p[2]  = -fDz;
	p[3] = x+TMath::Tan(a)*fDz;  p[4] =  0; p[5]  = fDz ;
	p[6] = x+TMath::Tan(a)*fDz;  p[7] =  0; p[8]  = fDz  ;
	p[9] = x-TMath::Tan(a)*fDz;  p[10] = 0; p[11] = -fDz;
	//	printf("%3d -> %3d -> %8x\n", d->GetSignal(), ci, fQuads.back().color);
      }
    }
    break;
  }

  }
}

/**************************************************************************/

void ITSModule::SetTrans()
{
  Double_t pos[3];
  Double_t rot[9];
  fInfo->fGeom->GetTrans(fID,pos);
  fInfo->fGeom->GetRotMatrix(fID,rot);
  Double_t *s, *d;

  // column major ii
  s = &rot[0]; d = &fMatrix[0];
  d[0] = s[0]; d[1] = s[3]; d[2] = s[6]; d[3] = 0;
  s = &rot[1]; d = &fMatrix[4];
  d[0] = s[0]; d[1] = s[3]; d[2] = s[6]; d[3] = 0;
  s = &rot[2]; d = &fMatrix[8];
  d[0] = s[0]; d[1] = s[3]; d[2] = s[6]; d[3] = 0;
  s = &pos[0]; d = &fMatrix[12];
  d[0] = s[0]; d[1] = s[1]; d[2] = s[2]; d[3] = 1;

  fTrans = true;
}

/**************************************************************************/

void ITSModule::Print(Option_t* ) const
{
  printf("ID %d, layer %d, ladder %d, det %d \n", fID, fLayer, fLadder, fDetID);
}
