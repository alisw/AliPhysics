#include "AliITSUMatLUT.h"
#include "AliLog.h"
#include "AliTrackerBase.h"
#include <TString.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom.h>


//___________________________________________________________
AliITSUMatLUT::AliITSUMatLUT() 
  :fRMin(0)
  ,fRMax(0)
  ,fDRInv(-1)
  ,fDR(0)
  ,fNBins(0) 
{
  // def c-tor
  for (int i=kNParTypes;i--;) fData[i]=0;
}

//___________________________________________________________
AliITSUMatLUT::~AliITSUMatLUT() 
{
  // d-tor
  for (int i=kNParTypes;i--;) delete fData[i];
}

//___________________________________________________________
AliITSUMatLUT::AliITSUMatLUT(Double_t rmin,Double_t rmax,Int_t nbin)
  :fRMin(rmin)
  ,fRMax(rmax)
  ,fDRInv(0)
  ,fDR(0)
  ,fNBins(nbin)
{
  //
  if (rmin<0 || rmax-rmin<1e-4 || nbin<1) AliFatal(Form("Illegal parameters Rmin:%f Rmax:%f Nbins:%d",rmin,rmax,nbin));
  fDRInv = fNBins/(fRMax-fRMin);
  fDR    = (fRMax-fRMin)/fNBins;
  for (int i=kNParTypes;i--;) {
    fData[i] = new Double_t[fNBins];
    memset(fData[i],0,fNBins*sizeof(Double_t));
  }
  //
}

//___________________________________________________________
AliITSUMatLUT::AliITSUMatLUT(const AliITSUMatLUT& src)
  :TObject(src)
  ,fRMin(src.fRMin)
  ,fRMax(src.fRMax)
  ,fDRInv(src.fDRInv)
  ,fDR(src.fDR)
  ,fNBins(src.fNBins)
{
  //
  for (int i=kNParTypes;i--;) {
    if (src.fData[i]) {
      fData[i] = new Double_t[fNBins];
      memcpy(fData[i],src.fData[i],fNBins*sizeof(Double_t));
    }
  }
  //
}

//___________________________________________________________
AliITSUMatLUT & AliITSUMatLUT::operator=(const AliITSUMatLUT& src)
{
  // copy 
  if (this == &src) return *this;
  this->~AliITSUMatLUT();
  new(this) AliITSUMatLUT(src);
  return *this;
  //  
}

//___________________________________________________________
void AliITSUMatLUT::FillData(Int_t ntest, Double_t zmin,Double_t zmax)
{
  // filla material data 
  double start[3],stop[3],parStep[7];
  if (fNBins<1) AliFatal("Limits are not set");
  if (ntest<1 || zmin>zmax) AliFatal(Form("Wrong parameters Ntest:%d Zmin:%f Zmax:%f",ntest,zmin,zmax));
  double dr = (fRMax-fRMin)/fNBins;
  AliInfo(Form("Building material table for %.3f<R<%.3f %.3f<Z<%.3f in %d bins using %d tracks",fRMin,fRMax,zmin,zmax,fNBins,ntest));
  const double kAngEps = 1e-4; // tiny slope to avoid tracks strictly normal to Z axis
  double *tmpAcc = new double[fNBins*kNParTypes];
  for (int itst=ntest;itst--;) {
    double parInt[kNParTypes]={0};
    double r   = fRMin;
    double phi = gRandom->Rndm()*TMath::Pi()*2;
    double cs  = TMath::Cos(phi);
    double sn  = TMath::Sin(phi);
    double angz = 2*(gRandom->Rndm()-0.5)*kAngEps;
    stop[0] = r*cs;
    stop[1] = r*sn;
    stop[2] = zmin + gRandom->Rndm()*(zmax-zmin);
    Bool_t fail = kFALSE;
    for (int ir=0;ir<fNBins;ir++) {
      r += dr;
      for (int i=3;i--;) start[i] = stop[i];
      stop[0] = r*cs;
      stop[1] = r*sn;
      stop[2] += dr*angz; // to avoid tracks normal to axis
      AliTrackerBase::MeanMaterialBudget(start,stop, parStep);
      if (parStep[1]>999) {fail = kTRUE; printf("fail\n"); break;}
      //
      parInt[kParX2X0] += parStep[1];
      parInt[kParRhoL] += parStep[0]*parStep[4];
      //
      for (int ip=kNParTypes;ip--;) tmpAcc[ir*kNParTypes+ip] = parInt[ip];
    }
    if (fail) {itst++; continue;} // propagation failed
    for (int ir=0;ir<fNBins;ir++) for (int ip=kNParTypes;ip--;) fData[ip][ir] += tmpAcc[ir*kNParTypes+ip];
  }
  //
  for (int ip=kNParTypes;ip--;) for (int ir=fNBins;ir--;) fData[ip][ir] /= ntest;
  delete[] tmpAcc;
  //  Print();
  //
}

//___________________________________________________________
void AliITSUMatLUT::Print(Option_t*) const
{
  // print data 
  printf("Average material budget in %d bins for %.4f<R<%.4f\n",fNBins,fRMin,fRMax);
  printf("  # :  rMin : rMax  \t  X2X0  (  incr  )\t  RhoL  (  incr  )\n");
  for (int i=0;i<fNBins;i++) {
    double r = fRMin + fDR*i;
    printf("%4d:%7.3f:%7.3f\t%.6f(%.6f)\t%.6f(%.6f)\n",i,r,r+fDR,
	   fData[kParX2X0][i],fData[kParX2X0][i]-(i==0 ? 0:fData[kParX2X0][i-1]),
	   fData[kParRhoL][i],fData[kParRhoL][i]-(i==0 ? 0:fData[kParRhoL][i-1]));
  }

}

//___________________________________________________________
TH1* AliITSUMatLUT::GetHisto(const Option_t* option, const Char_t *name) const
{
  // extract data to histo
  if (fNBins<1) return 0;
  TString nms = name;
  if (nms.IsNull()) nms = "matLUT";
  TH1F* h = new TH1F(nms.Data(),nms.Data(),fNBins,fRMin,fRMax);
  TString opts = option;
  opts.ToLower();
  Bool_t diff = opts.Contains("d");
  Int_t par = opts.Contains("rhol") ? kParRhoL : kParX2X0;
  double valPrev = 0;
  for (int i=0;i<fNBins;i++) {
    double val = fData[par][i] - valPrev;
    if (diff) valPrev = fData[par][i];
    h->SetBinContent(i+1, val);
  }
  return h;
}

//___________________________________________________________
Double_t AliITSUMatLUT::GetMatBudget(const Double_t *pnt0, const Double_t *pnt1, Double_t *ret) const
{
  // query the mat.budget between 2 points
  double r0 = TMath::Sqrt(pnt0[0]*pnt0[0] + pnt0[1]*pnt0[1]);
  double r1 = TMath::Sqrt(pnt1[0]*pnt1[0] + pnt1[1]*pnt1[1]);
  if (r1<r0) {double t=r1;r1=r0;r0=t;}
  double dr = r1-r0;
  if (dr<1e-4) {
    for (int i=kNParTypes;i--;) ret[i] = 0;
    return 0;
  }
  //
  double di = pnt1[0]-pnt0[0];
  double dst = di*di;
  di = pnt1[1]-pnt0[1];
  dst += di*di;
  di = pnt1[2]-pnt0[2];
  dst += di*di;
  dst = TMath::Sqrt(dst);
  double dstf = dst/dr;
  //
  double par0[kNParTypes],par1[kNParTypes];
  GetData(r0,par0);
  GetData(r1,par1);
  for (int i=kNParTypes;i--;) ret[i] = dstf*(par1[i]-par0[i]);
  return dst;
  //
}

//___________________________________________________________
void AliITSUMatLUT::GetData(Double_t r, Double_t* params) const
{
  // get integrated parameters from rMin to r
  double binF = (r-fRMin)*fDRInv;
  if (binF<0) binF=0;
  else if (binF>=fNBins) binF = fNBins-1;
  int   bin   = int(binF);
  double frac = 1.-(binF-bin);
  for (int i=kNParTypes;i--;) {
    double prev = bin ? fData[i][bin-1] : 0;
    params[i] = fData[i][bin] - frac*(fData[i][bin]-prev);
  }
  //
}
