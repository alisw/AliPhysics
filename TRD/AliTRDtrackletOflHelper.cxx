#include "TObjArray.h"
#include "TMath.h"

#include "AliLog.h"
#include "AliMathBase.h"

#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDpadPlane.h"
#include "AliTRDtrackletOflHelper.h"

ClassImp(AliTRDtrackletOflHelper)
//___________________________________________________________________
AliTRDtrackletOflHelper::AliTRDtrackletOflHelper() 
  :TObject()
  ,fRow(-1)
  ,fClusters(NULL)
  ,fPadPlane(NULL)
{
// Default constructor
  fCol[0]=144; fCol[1]=0;
  fTBrange[0]=1; fTBrange[1]=21;
}

//___________________________________________________________________
AliTRDtrackletOflHelper::AliTRDtrackletOflHelper(const AliTRDtrackletOflHelper &ref) 
  :TObject(ref)
  ,fRow(-1)
  ,fClusters(NULL)
  ,fPadPlane(ref.fPadPlane)
{
// Copy constructor
  fRow = ref.fRow; fCol[0]=ref.fCol[0]; fCol[1]=ref.fCol[1];
  fTBrange[0] = ref.fTBrange[0];fTBrange[1] = ref.fTBrange[1];
  Int_t n(0);
  if(ref.fClusters){
    if((n = ref.fClusters->GetEntriesFast())) {
      fClusters = new TObjArray(n);
      for(Int_t ic(n); ic--;) fClusters->AddAt(ref.fClusters->At(ic), ic);
    }
  }
}

//___________________________________________________________________
AliTRDtrackletOflHelper& AliTRDtrackletOflHelper::operator=(const AliTRDtrackletOflHelper &rhs)
{
  if(this != &rhs){
    TObject::operator=(rhs);
    fPadPlane  = rhs.fPadPlane;
    fRow       = rhs.fRow;
    fCol[0]    = rhs.fCol[0];
    fCol[1]    = rhs.fCol[1];
    fTBrange[0]= rhs.fTBrange[0];
    fTBrange[1]= rhs.fTBrange[1];
    if(rhs.fClusters){ 
      Int_t n(rhs.fClusters->GetEntriesFast());
      if(!fClusters) fClusters = new TObjArray(n);
      if(fClusters->GetEntriesFast() != n){ 
        fClusters->Clear();
        fClusters->Expand(n);
      }
      for(Int_t ic(n); ic--;)  fClusters->AddAt(rhs.fClusters->At(ic), ic);
    } else {
      if(fClusters) delete fClusters;
    }
  }
  return *this;
}

//___________________________________________________________________
AliTRDtrackletOflHelper::~AliTRDtrackletOflHelper()
{
// Clean helper  
  if(fClusters) fClusters->Clear();
  delete fClusters;
}

//___________________________________________________________________
Int_t AliTRDtrackletOflHelper::Expand(TObjArray *cls, Int_t *mark, Int_t groupId)
{
// Allocate new clusters for this helper.  If a subset of clusters is to be allocated 
// this can be specified via the identifier "groupId" which have to be compared 
// with individual elements in the array "mark".
// 
// Return total no. of clusters in helper

  if(!fPadPlane || !fClusters ){
    AliError("Helper not initialized."); 
    return 0;
  }
  Int_t ncl(fClusters->GetEntriesFast());
  Int_t mcl(cls->GetEntriesFast());
  fClusters->Expand(mcl + ncl);
  fCol[0]=144; fCol[1]=0;
  AliTRDcluster *c(NULL);
  for(Int_t icl(0), jcl(ncl); icl<mcl; icl++){
    if(!(c=(AliTRDcluster*)cls->At(icl))) continue;
    if(mark && mark[icl]!=groupId) continue;
    fClusters->AddAt(c, jcl++);
    Int_t col(c->GetPadCol());
    Short_t *sig(c->GetSignals());
    for(Int_t icol(0), jcol(col-3); icol<7; icol++, jcol++){
      if(sig[icol]==0) continue;
      if(jcol<fCol[0]) fCol[0]=jcol;
      if(jcol>fCol[1]) fCol[1]=jcol;
    }
  }

  AliDebug(1, Form("Segment[%d] Clusters[%2d] pad[%3d %3d].", groupId, fClusters->GetEntriesFast(), fCol[0], fCol[1]));
  return fClusters->GetEntriesFast();
}


//___________________________________________________________________
Int_t AliTRDtrackletOflHelper::Init(AliTRDpadPlane *p, TObjArray *cls, Int_t *mark, Int_t groupId)
{
// Allocate clusters for this helper.  If a subset of clusters is to be allocated 
// this can be specified via the identifier "groupId" which have to be compared 
// with individual elements in the array "mark".
// 
// Return no. of clusters allocated

  if(!p){
    AliError("PadPlane not initialized."); return 0;
  }
  if(!cls){
    AliError("Cluster array not initialized."); return 0;
  }
  Int_t ncl(cls->GetEntriesFast());
  if(ncl<2){ 
    AliDebug(1, Form("Segment[%d] failed n[%d].", groupId, ncl)); return 0;
  }
  Int_t mcl(ncl);
  if(mark){
    mcl = 0;
    for(Int_t icl(ncl); icl--;) 
      if(cls->At(icl) && mark[icl]==groupId) mcl++;
  }  
  if(mcl<2){ 
    AliDebug(1, Form("Segment[%d] failed n[%d] in group.", groupId, mcl)); return 0;
  }
  if(!fClusters) fClusters = new TObjArray(mcl);
  else{ 
    fClusters->Clear();
    fClusters->Expand(mcl);
  }
  fCol[0]=144; fCol[1]=0; fRow=-1;
  AliTRDcluster *c(NULL);
  for(Int_t icl(0), jcl(0); icl<ncl; icl++){
    if(!(c=(AliTRDcluster*)cls->At(icl))) continue;
    if(mark && mark[icl]!=groupId) continue;
    fClusters->AddAt(c, jcl++);
    if(fRow<0) fRow = c->GetPadRow();
    Int_t col(c->GetPadCol());
    Short_t *sig(c->GetSignals());
    for(Int_t icol(0), jcol(col-3); icol<7; icol++, jcol++){
      if(sig[icol]==0) continue;
      if(jcol<fCol[0]) fCol[0]=jcol;
      if(jcol>fCol[1]) fCol[1]=jcol;
    }
  }
  fPadPlane = p;
  
  AliDebug(1, Form("Segment[%d] Clusters[%2d] pad[%3d %3d].", groupId, fClusters->GetEntriesFast(), fCol[0], fCol[1]));
  return fClusters->GetEntriesFast();
}

//___________________________________________________________________
Int_t AliTRDtrackletOflHelper::ClassifyTopology()
{
// Classify topology and return classification code
// 0 - normal tracklet
// 1 - delta ray candidate
// 2 - secondary candidate
// 3 - "elephant" candidate
// 4 - unknown topology

  if(!fClusters){ 
    AliError("Helper not initialized. Missing clusters.");
    return kUnknown;
  }
  Int_t ncl(fClusters->GetEntries());
  if(!ncl){ 
    AliError("Helper not initialized. No cluster allocated.");
    return kUnknown;
  }
  const Int_t kRange(22); // DEFINE based on vd 
  
  // compute local occupancy
  Int_t localOccupancy[AliTRDseedV1::kNtb]; memset(localOccupancy, 0, AliTRDseedV1::kNtb*sizeof(Int_t));
  AliTRDcluster *c(NULL);
  Double_t sy[kNcls]; Int_t mcl(0);
  for(Int_t icl(ncl), time; icl--;){
    c = (AliTRDcluster*)fClusters->At(icl);
    time = c->GetPadTime();
    if(time==0 || time>=kRange){
      sy[icl] = -1; // mark clusters outsde drift volume
      continue;
    }
    // protect against wrong error param.
    if(c->GetSigmaY2() < 1.e-5) sy[icl] = 0.02;
    else sy[icl] = TMath::Min(TMath::Sqrt(c->GetSigmaY2()), 0.04);
    localOccupancy[time]++;
    mcl++;
  }
  if(!mcl){
    AliWarning("No clusters in the active area.");
    return kUnknown;
  } 
  //compute number of 0 bins and high occupancy
  Int_t goodOccupancy(0), highOccupancy(0), lowOccupancy(0);
  for(Int_t itb(1); itb<kRange; itb++){
    switch(localOccupancy[itb]){ 
      case 0: lowOccupancy++; break;
      case 1: goodOccupancy++; break;
      default: highOccupancy++; break;
    }
  }
  AliDebug(2, Form("H[%2d | %5.2f%%] L[%2d | %5.2f%%] N[%2d | %5.2f%% | %5.2f%%]",
          highOccupancy, 100.*highOccupancy/kRange, 
          lowOccupancy, 100.*lowOccupancy/kRange, 
          goodOccupancy, 100.*goodOccupancy/kRange, 100.*goodOccupancy/mcl));
  
  // filter
  Double_t dy[kNcls];
  if(goodOccupancy==mcl) return kNormal;
  else if(Double_t(goodOccupancy)/kRange > 0.9){
    if(highOccupancy == 0 ) {
      if(!FitPSR(dy)) return kUnknown;
      Int_t nin(0);
      for(Int_t idy(ncl); idy--;){
        if(sy[idy]<0.) continue;
        if(dy[idy] > 3*sy[idy]) continue;
        nin ++;
      }
      if(Double_t(nin)/mcl > 0.8) return kNormal;
      else return kUnknown;
    } else return kDeltaRay;
  } else return kUnknown;
}


//___________________________________________________________________
void AliTRDtrackletOflHelper::FindSolidCls(Bool_t *mark, Int_t *q)
{
//  Find clusters produced by large fluctuations of energy deposits
//  Largest charge and well separation from neighbors

  Int_t ntb(AliTRDseedV1::kNtb);
  Int_t idx[ntb+1];
  TMath::Sort(ntb, q, idx, kTRUE);
  Int_t qmax = Int_t(0.3*q[idx[0]]);
  mark[0] = kFALSE;
  for(Int_t icl(ntb-5); icl<ntb; icl++) mark[icl] = kFALSE;
  for(Int_t icl(0); icl<ntb; icl++){
    Int_t jcl(idx[icl]);
    if(!mark[jcl]) continue;
    if(q[jcl-1]>q[jcl] || q[jcl+1]>q[jcl]){
      mark[jcl] = kFALSE;
      continue;
    }
    if(q[jcl] < qmax){
      mark[jcl] = kFALSE;
      continue;
    }
    for(Int_t kcl=TMath::Max(0, jcl-2); kcl<jcl+3; kcl++){
      if(kcl==jcl) continue;
      mark[kcl] = kFALSE;
    }
  }
}

//___________________________________________________________________
Bool_t AliTRDtrackletOflHelper::FitPSR(Double_t dy[200], Bool_t useSolid)
{
// Fit tracklet in Pad System of Reference [PSR] to avoid uncertainty related to 
// Lorentz angle correction

  Bool_t mark[200]; memset(mark, 0, 200*sizeof(Bool_t));
  Int_t q[200]; memset(q, 0, 200*sizeof(Int_t));
  Int_t ncl(fClusters->GetEntries());
  AliTRDcluster *c(NULL);
  for(Int_t icl(ncl); icl--;){
    c = (AliTRDcluster*)fClusters->At(icl);
    if(c->GetPadRow() != fRow) continue;
    Int_t time(c->GetPadTime());
    mark[time] = kTRUE; q[time] = Int_t(c->GetQ());
  }
  if(useSolid) FindSolidCls(mark, q);
  
  Double_t 
    x[AliTRDseedV1::kNtb], y[AliTRDseedV1::kNtb], sy[AliTRDseedV1::kNtb], 
    xf[AliTRDseedV1::kNtb], yf[AliTRDseedV1::kNtb], syf[AliTRDseedV1::kNtb], 
    par[3];
  Int_t jcl(0), kcl(0);
  for(Int_t icl(0); icl<ncl; icl++){
    c = (AliTRDcluster*)fClusters->At(icl);
    if(c->GetPadRow() != fRow) continue;
    Int_t col(c->GetPadCol());
    Int_t time(c->GetPadTime());
    Double_t center(c->GetCenter());
    Double_t cw(fPadPlane->GetColSize(col));
    //Double_t corr = AliTRDcluster::GetYcorr(AliTRDgeometry::GetLayer(det), center);
    y[jcl] = fPadPlane->GetColPos(col) + (.5 + center)*cw /*+ corr*/;
    x[jcl] = c->GetX();
    sy[jcl]= TMath::Sqrt(c->GetSigmaY2());
    if(mark[time]){
      yf[kcl] = y[jcl];
      xf[kcl] = x[jcl];
      syf[kcl]= useSolid ? 0.5/time:sy[jcl];
      kcl++;
    } 
    jcl++;
  }
  Fit(kcl, xf, yf, syf, par);

  for(Int_t icl(0); icl<jcl; icl++){
    Double_t dx(x[icl] - par[2]);
    dy[icl] = y[icl] - (par[0] + par[1]*dx);
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliTRDtrackletOflHelper::Fit(Int_t n, Double_t *x, Double_t *y, Double_t *sy, Double_t *par, Double_t sCut,  Double_t *cov)
{
// Iterative robust tracklet fit 
  if(n<3) return kFALSE;
  //compute <x>
  if(Int_t(par[2])==21122012){ 
    par[2] = 0.; 
  } else {
    par[2] = 0.; for(Int_t ic(n); ic--;) par[2] += x[ic]; par[2] /= n;
  }
  AliTRDtrackerV1::AliTRDLeastSquare &f=Fitter();
  for(Int_t iter(0); iter<3; iter++){
    f.Reset();
    Int_t jp(0);
    for(Int_t ip(0); ip<n; ip++){
      Double_t dx(x[ip]-par[2]);
      Double_t dy(y[ip] - (par[0] + par[1]*dx));
      if(iter && TMath::Abs(dy)>sCut*sy[ip]) continue;
      f.AddPoint(&dx, y[ip], sy[ip]);
      jp++;
    }
    if(jp<3) continue;
    if(!f.Eval()) continue;
    par[0]=f.GetFunctionParameter(0);
    par[1]=f.GetFunctionParameter(1);
    if(cov) f.GetCovarianceMatrix(cov);
    AliDebugGeneral("AliTRDtrackletOflHelper::Fit()", 2, Form("Iter[%d] Ncl[%2d/%2d] par[%f %f %f]", iter, jp, n, par[0], par[1], par[2]));
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliTRDtrackletOflHelper::Fit(Double_t *par, Double_t sCut) const
{
// Wrapper for clusters attach to this for static Fit function  
//
  Int_t n(fClusters->GetEntriesFast()), jc(0), dr(0);
  Double_t corr = TMath::Tan(TMath::DegToRad()*fPadPlane->GetTiltingAngle())*
                  fPadPlane->GetLengthIPad();
  Double_t x[kNcls], y[kNcls], sy[kNcls];
  AliTRDcluster *c(NULL);
  for(Int_t ic(0); ic<n; ic++){
    c = (AliTRDcluster*)fClusters->At(ic);
    if(!c->IsInChamber()) continue;
    dr = c->GetPadRow() - fRow;
    x[jc] = c->GetX();
    y[jc] = c->GetY()+corr*dr;
    sy[jc]= c->GetSigmaY2()>0?(TMath::Min(TMath::Sqrt(c->GetSigmaY2()), 0.08)):0.08;
    jc++;
  }
  return Fit(jc, x, y, sy, par, sCut);
}

//___________________________________________________________________
void AliTRDtrackletOflHelper::GetColSignals(Int_t col, Int_t adc[32], Bool_t mainRow) const
{
  memset(adc, 0, 32*sizeof(Int_t));
  if(col<fCol[0] || col>fCol[1]) return;
  
  AliTRDcluster *cc(NULL);
  for(Int_t ic(fClusters->GetEntriesFast()); ic--;){
    cc = (AliTRDcluster*)fClusters->At(ic);
    if((mainRow && cc->GetPadRow()!=fRow) ||
       (!mainRow && cc->GetPadRow()==fRow)) continue;
    Short_t *sig = cc->GetSignals();
    Int_t padcol(cc->GetPadCol());
    Int_t time(cc->GetPadTime());
    for(Int_t icol(0), jcol(padcol-3); icol<7; icol++, jcol++){
      if(jcol!=col) continue;
      adc[time]+=sig[icol];
    }
  }
}

//___________________________________________________________________
Int_t AliTRDtrackletOflHelper::GetRMS(Double_t &r, Double_t &m, Double_t &s, Double_t &xm) const
{
// Calculate Rotation[r], Mean y[m] (at mean radial position [xm]) and Sigma[s] (of a gaussian distribution in the tracklet [SR])
// for clusters attach to this helper. The Rotation and Mean are calculated without tilt correction option.
// It returns the number of clusters in 1 sigma cut.
  
  Int_t n(fClusters->GetEntriesFast());
  if(n==2){
    return n;
  }
  
  Double_t par[3] = {0., 0., 0.};
  Fit(par);
  xm= par[2];
  r = par[1];

  Double_t corr = TMath::Tan(TMath::DegToRad()*fPadPlane->GetTiltingAngle())*
                  fPadPlane->GetLengthIPad();
  Double_t y[kNcls];
  AliTRDcluster *c(NULL);
  for(Int_t ic(n); ic--;){
    c = (AliTRDcluster*)fClusters->At(ic);
    Double_t x(c->GetX() - xm);
    Int_t dr(c->GetPadRow() - fRow);
    y[ic] = c->GetY()+corr*dr - (par[0] + par[1]*x);
  }
  Double_t m1(0.);
  AliMathBase::EvaluateUni(n, y, m1, s, 0);
  m = par[0] + m1;
  Int_t n0(0);
  for(Int_t ic(n); ic--;){
    c = (AliTRDcluster*)fClusters->At(ic);
    Double_t sy = c->GetSigmaY2()>0?(TMath::Min(TMath::Sqrt(c->GetSigmaY2()), 0.08)):0.08;
    if(TMath::Abs(y[ic]-m1) <= sy) n0++;
  }
  return n0;
}

//___________________________________________________________________
Double_t AliTRDtrackletOflHelper::GetSyMean() const
{
// Calculate mean uncertainty of clusters

  Double_t sym(0.);
  Int_t n(fClusters->GetEntriesFast());
  AliTRDcluster *c(NULL);
  for(Int_t ic(n); ic--;){
    c = (AliTRDcluster*)fClusters->At(ic);
    sym += c->GetSigmaY2()>0?(TMath::Min(TMath::Sqrt(c->GetSigmaY2()), 0.08)):0.08;
  }
  return sym/=n;
}

//___________________________________________________________________
Int_t AliTRDtrackletOflHelper::Segmentation(Int_t n, Double_t *x, Double_t *y, Int_t *Index)
{
// Segmentation of clusters in the tracklet roads
//
// The user supply the coordinates of the clusters (x,y) and their number "n". Also
// The array "Index" has to be allocated by the user with a size equal or larger than "n".
// On return the function returns the number of segments found and the array "Index" i-th element 
// is filled with the index of the tracklet segment to which the i-th cluster was assigned too.
// 
// Observation
// The parameter which controls the segmentation is set inside the function "kGapSize" and it is used 
// for both "x" and "y" segmentations. An improvement can be a parametrization as function of angle of 
// incidence.
//
// author:
// Alex Bercuci <a.bercuci@gsi.de>

  if(!n || !x || !y || !Index){
    AliErrorGeneral("AliTRDtrackletOflHelper::Segmentation()", "One of the input arrays non initialized.");
    return 0;
  }
  const Int_t kBuffer = 200;
  if(n>kBuffer){
    AliWarningGeneral("AliTRDtrackletOflHelper::Segmentation()", Form("Input array size %d exceed buffer %d. Truncate.", n, kBuffer));
    n = kBuffer;
  }
  const Double_t kGapSize(0.2); // cm
  Int_t ng(0),
        nc(0);
  Double_t xx[kBuffer], dy;
  Int_t idx[kBuffer+1], jdx[kBuffer], kdx[kBuffer];
  TMath::Sort(n, y, idx);
  for(Int_t iy(0); iy<n; iy++){
    dy = iy>0?(TMath::Abs(y[idx[iy-1]]-y[idx[iy]])):0.;
    if(dy>kGapSize){
      TMath::Sort(nc, xx, jdx);
      for(Int_t ic(0), jc0, jc1; ic<nc; ic++){
        jc0 = ic>0?kdx[jdx[ic-1]]:0;
        jc1 = kdx[jdx[ic]];
        dy = TMath::Abs(y[jc0] - y[jc1]);
        if(ic && dy>kGapSize) ng++;
        Index[jc1] = ng;
        AliDebugGeneral("AliTRDtrackletOflHelper::Segmentation()", 4, Form("  y[%2d]=%+f x[%+f] %2d -> %2d -> %2d ng[%2d]", jc1, y[jc1], x[jc1], ic, jdx[ic], jc1, ng));
      }
      ng++;
      nc=0;
    }
    xx[nc] = x[idx[iy]];
    kdx[nc]= idx[iy];
    AliDebugGeneral("AliTRDtrackletOflHelper::Segmentation()", 4, Form("y[%2d]=%+f -> %2d", idx[iy], y[idx[iy]], nc));
    nc++;
  }
  if(nc){
    TMath::Sort(nc, xx, jdx);
    for(Int_t ic(0), jc0, jc1; ic<nc; ic++){
      jc0 = ic>0?kdx[jdx[ic-1]]:0;
      jc1 = kdx[jdx[ic]];
      dy = TMath::Abs(y[jc0] - y[jc1]);
      if(ic && dy>kGapSize) ng++;
      Index[jc1] = ng;
      AliDebugGeneral("AliTRDtrackletOflHelper::Segmentation()", 4, Form("  y[%2d]=%+f x[%+f|%+f] %2d -> %2d -> %2d ng[%2d]\n", jc1, y[jc1], xx[jdx[ic]], x[jc1], ic, jdx[ic], jc1, ng));
    }
    ng++;
  }
  return ng;
}

//___________________________________________________________________
void AliTRDtrackletOflHelper::SetTbRange(Float_t t0, Float_t vd)
{
// Set first time bin and total number of time bins corresponding to clusters 
// in chamber based on the calibrated info "t0" and drift velocity "vd"

  // TO CHECK
  fTBrange[0] = Int_t(t0*0.1);
  fTBrange[1] = Int_t(vd*0.1);
}

#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
//___________________________________________________________________
void AliTRDtrackletOflHelper::View(TVirtualPad *vpad)
{
// Visualization support. Draw this tracklet segment  


  Int_t row(-1), det(-1);
  Int_t n(fClusters->GetEntriesFast());
  AliInfo(Form("Processing Clusters[%2d] Class[%d].", n ,ClassifyTopology()));
  
  // prepare drawing objects
  Int_t ncols(fPadPlane->GetNcols());
  Double_t cw(fPadPlane->GetColSize(1));
  TH2 *h2 = new TH2I("h2pm", ";y_{Local} [cm]; pad time;Charge", ncols, fPadPlane->GetColPos(1)-cw, fPadPlane->GetColPos(ncols-1)+cw, 31, -30.5, 0.5);
  h2->SetMarkerColor(kWhite);h2->SetMarkerSize(2.);
  h2->SetFillColor(9);
  
  TGraph *gcls = new TGraph(n);
  gcls->SetMarkerColor(kBlack);gcls->SetMarkerStyle(20);
  TGraph *gclsLC = new TGraph(n);
  gclsLC->SetMarkerColor(kBlack);gclsLC->SetMarkerStyle(28);
  TGraph *gclsSC = new TGraph(n);
  gclsSC->SetMarkerColor(kBlack);gclsSC->SetMarkerStyle(4);gclsSC->SetMarkerSize(1.5);
  TGraph *gFP = new TGraph(n);
  gFP->SetMarkerColor(kBlack);gFP->SetLineWidth(2);
  
  // fill signal data
  Bool_t map[200]; memset(map, 0, 200*sizeof(Bool_t));
  Int_t qa[200]; memset(qa, 0, 200*sizeof(Int_t));
  AliTRDcluster *cc(NULL); Int_t tm = 30; Double_t ym=0.;
  for(Int_t ic(n); ic--;){
    cc = (AliTRDcluster*)fClusters->At(ic);
    if(cc->GetPadRow() != fRow) continue; //TODO extend for 2 pad rows
    Short_t *sig = cc->GetSignals();
    Int_t col(cc->GetPadCol());
    det = cc->GetDetector();
    row = cc->GetPadRow();
    Int_t time(cc->GetPadTime());
    map[time] = kTRUE; qa[time] = (Int_t)cc->GetQ();
    for(Int_t ipad(0), jpad(col-2); ipad<7; ipad++, jpad++){
      Int_t q = (Int_t)h2->GetBinContent(jpad, 31-time);
      h2->SetBinContent(jpad, 31-time, q+sig[ipad]);
    }
    Double_t y0 = fPadPlane->GetColPos(col) + (.5 + cc->GetCenter())*cw;
    //h2->GetXaxis()->GetBinCenter(pad)+cc->GetCenter();
    gcls->SetPoint(ic, y0, -time);
    if(time <= tm) {ym = cc->GetY() - y0; tm = time;}
  }

  // draw special clusters (solid and Lorentz corrected and fit) 
  FindSolidCls(map, qa);
  Double_t ddy[200]; FitPSR(ddy, kTRUE);
  Double_t dt, dy;
  for(Int_t ic(0), jc(0); ic<n; ic++){
    cc = (AliTRDcluster*)fClusters->At(ic);
    if(cc->GetPadRow() != fRow) continue; //TODO extend for 2 pad rows
    gcls->GetPoint(ic, dy, dt);
    gclsLC->SetPoint(ic, cc->GetY()-ym, dt);
    gFP->SetPoint(ic, dy-ddy[ic], dt);
    if(map[cc->GetPadTime()]) gclsSC->SetPoint(jc++, dy, dt);
  }
  
  
  // prepare frame histo
  h2->SetName(Form("h2s%03d%02d", det, row));
  Int_t binSrt(fCol[0]+1), binSop(fCol[1]+1);
  TH1 *h1(NULL);
  
  // show everything
  if(!vpad){ 
    vpad = new TCanvas(Form("c%03d%02d", det, row), Form("D %03d [%02d_%d_%d] R[%02d]", det, AliTRDgeometry::GetSector(det), AliTRDgeometry::GetStack(det), AliTRDgeometry::GetLayer(det), row), 700, 500);
  }
  TVirtualPad *pp(NULL);
  vpad->Divide(2,1,2.e-5,2.e-5);
  pp = vpad->cd(1); pp->SetRightMargin(0.0001);pp->SetTopMargin(0.1);pp->SetBorderMode(0); pp->SetFillColor(kWhite);
  h1 = h2->ProjectionY(); h1->GetYaxis()->SetTitle("Total Charge");
  h1->SetTitle(Form("D%03d[%02d_%d_%d] R[%02d]", 
    det,AliTRDgeometry::GetSector(det), AliTRDgeometry::GetStack(det), AliTRDgeometry::GetLayer(det), row));
  TAxis *ax(h1->GetXaxis()); for(Int_t ib(0), jb(1); ib<ax->GetNbins(); ib++, jb++) ax->SetBinLabel(jb, Form("%d", -Int_t(ax->GetBinCenter(jb))));
  h1->Draw("hbar3");

  pp = vpad->cd(2);
  pp->SetRightMargin(0.15);pp->SetLeftMargin(0.0001);pp->SetTopMargin(0.1);
  pp->SetBorderMode(0); pp->SetFillColor(kWhite);
  ax=h2->GetXaxis();
  ax->SetRange(TMath::Max(1, binSrt-1), TMath::Min(ncols, binSop+1));
  h2->Draw("coltextz");
  gPad->Update();
  TGaxis *axis = new TGaxis(gPad->GetUxmin(),
                    gPad->GetUymax(),
                    gPad->GetUxmax(),
                    gPad->GetUymax(),
                    TMath::Max(0, binSrt-2)-0.5, TMath::Min(ncols-1, binSop)+0.5, 510,"-L");

  axis->SetNdivisions(103+binSop-binSrt);
  axis->SetTitle("pad col");
  axis->Draw();

  TLegend *leg = new TLegend(0.01, 0.1, 0.84, 0.21);
  leg->SetBorderSize(1); leg->SetFillColor(kYellow-9);leg->SetTextSize(0.04);
  gcls->Draw("p");  leg->AddEntry(gcls, "Cls. in Pad [SR]", "p");
  gclsLC->Draw("p"); leg->AddEntry(gclsLC, "Lorentz Corr. Cls.", "p");
  gclsSC->Draw("p"); leg->AddEntry(gclsSC, "Solid Cls.", "p");
  gFP->Draw("l"); leg->AddEntry(gFP, "Fit in Pad [SR]", "l");
  leg->Draw();
}

//___________________________________________________________________
AliTRDtrackerV1::AliTRDLeastSquare& AliTRDtrackletOflHelper::Fitter()
{
  static AliTRDtrackerV1::AliTRDLeastSquare f;
  return f;
}
