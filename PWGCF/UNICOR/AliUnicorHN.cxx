/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// multidimensional histogram 
// Physically, the data are kept in one single one-dimensional histogram. 
// Projecting on 1, 2, and n-1 dimensions is implemented.
// The histogram can be saved on root file as the one-dimensional data 
// histogram and the axes, thus eternal forward compatibility is ensured. 
//=============================================================================

#include <cmath>
#include <stdlib.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include "AliUnicorHN.h"

ClassImp(AliUnicorHN)

//=============================================================================
AliUnicorHN::AliUnicorHN(const char *nam, Int_t ndim, TAxis **ax) 
  : TH1D(nam, nam, Albins(ndim,ax), 0.5, Albins(ndim,ax)+0.5), 
    fNdim(ndim) {

  // constructor
 
  // Above, we just have managed to create a single one-dimensional histogram 
  // with number of bins equal to the product of the numbers of bins in all 
  // dimensions. For easy inspection the histogram range was set to -0.5,n-0.5.

  for (int i=0; i<fNdim; i++) ax[i]->Copy(fAxis[i]); 
  for (int i=0; i<fNdim; i++) fAxis[i].SetName(Form("axis%d",i));
  
  // for speed reasons, number of bins of each axis is stored on fNbins
  // and the running product of the latter on fMbins
  // so fMbins = {...,fNbins[fNdim-2]*fNbins[fNdim-1],fNbins[fNdim-1],1}
  
  for (int i=0; i<fgkMaxNdim; i++) fNbins[i] = fAxis[i].GetNbins();
  for (int i=0; i<fgkMaxNdim; i++) fMbins[i] = 1;
  for (int i=fNdim-1; i>0; i--) fMbins[i-1] = fMbins[i]*fNbins[i];
  printf("   %d-dimensional histogram %s with %d bins created\n",fNdim,nam,GetNbinsX());
}
//=============================================================================
AliUnicorHN* AliUnicorHN::Retrieve(const char *filnam, const char *nam) { 

  // retrieve a multidimensional histogram from file

  TFile *f = TFile::Open(filnam,"read");
  if (!f) return 0;
  if (!f->cd(nam)) {f->Close(); return 0;}
  TH1D *hist = (TH1D*) gDirectory->Get("histo");
  if (!hist) {printf("No histogram histo on file %s in directory %s\n",filnam,nam); return 0;}
  hist->SetDirectory(gROOT);
  TAxis *ax[fgkMaxNdim]={0};
  int n=0;
  while ((ax[n] = (TAxis *) gDirectory->Get(Form("axis%d",n)))) n++; 
  f->Close();
  
  if (hist->GetNbinsX()!=Albins(n,ax)) {
    printf("number of bins of histo %d differs from product of nbins of axes %d\n",
	   hist->GetNbinsX(),Albins(n,ax));
    return 0;
  }

  // derive the name from the path (ccc from aaa/bbb/ccc)
  TString path=nam;
  TObjArray *ar = path.Tokenize("/");
  TObjString *os = (TObjString*) ar->Last();
  const char *lastnam = os->GetString().Data();

  AliUnicorHN *hi = new AliUnicorHN(lastnam,n,ax);
  //  *((TH1D*) hi) = *hist;
  hist->Copy(*((TH1D*)hi));
  hi->SetName(lastnam);
  delete hist;
  return hi;
}
//=============================================================================
Int_t AliUnicorHN::Albins(Int_t n, TAxis **ax) {

  // Calculate product of nbins of ax[0]...ax[n-1]
  // (= total number of bins of the multidimensional histogram to be made). 

  Int_t nbins = 1;
  //  while (n--) nbins *= ax[n]->GetNbins();
  for (int i=0; i<n; i++) nbins *= ax[i]->GetNbins();
  return nbins;
}
//=============================================================================
Int_t AliUnicorHN::MulToOne(const Int_t * const k) const {

  // Calculate the 1-dim index n from n-dim indices k[fNdim].
  // Valid k[i] should be between 0 and fNbins[i]-1.
  // Valid result will be between 0 and GetNbinsX()-1. 
  // Return -1 if under- or overflow in any dimension.

  Int_t n = 0;
  for (int i=0; i<fNdim; i++) {
    if (k[i]<0) return -1;
    if (k[i]>=fNbins[i]) return -1;
    n += fMbins[i]*k[i];
  }
  return n;
}
//=============================================================================
Int_t AliUnicorHN::MulToOne(Double_t *x) {

  // Calculate the 1-dim index n from n-dim vector x, representing the 
  // abscissa of the n-dim histogram. The result will be between 0 and 
  // GetNbinsX()-1. 

  Int_t k[fgkMaxNdim] = {0};
  for (int i=0; i<fNdim; i++) k[i] = fAxis[i].FindFixBin(x[i])-1;
  return MulToOne(k);
}
//=============================================================================
void AliUnicorHN::OneToMul(Int_t n, Int_t *k) const {

  // Calculate the n-dim indices k[fNdim] from 1-dim index n.
  // Valid n should be between 0 and GetNbinsX()-1. 
  // Valid results will be between 0 and fNbins[i]-1.

  div_t idi; // integer division structure 
  for (int i=0; i<fNdim; i++) {
    idi  = div(n,fMbins[i]);
    k[i] = idi.quot;  // quotient
    n    = idi.rem;   // remainder
  }
}
//=============================================================================
Int_t AliUnicorHN::Fill(Double_t *xx, Double_t w) {

  // Fill the histogram. The array xx holds the abscissa information, w is the 
  // weigth. The 1-dim histogram is filled using the standard Fill method 
  // (direct access to the arrays was tried and was not faster). 

  int nbin = MulToOne(xx);
  if (nbin == -1) return 0;
  return TH1D::Fill(nbin+1,w); 
}
//=============================================================================
Int_t AliUnicorHN::Fill(Double_t x0, Double_t x1, ...) {

  // Fill the histogram. Arguments are passed as doubles rather than array. 

  va_list ap;
  Double_t xx[fgkMaxNdim] = {x0, x1};
  va_start(ap,x1);
  for (int i=2; i<fNdim; i++) xx[i] = va_arg(ap,Double_t);
  Double_t weigth = va_arg(ap,Double_t);
  va_end(ap);
  return Fill(xx,weigth);
}
//=============================================================================
Int_t AliUnicorHN::Save() const {

  // Save the 1-dim histo and the axes in a subdirectory on file. This might 
  // not be the most elegant way but it is very simple and backward and forward 
  // compatible. 

  Int_t nbytes = 0;
  TH1D histo(*this);
  histo.SetName("histo"); // hadd bug fix; actually, does not cost memory, strange

  TDirectory *dest = gDirectory->mkdir(GetName());
  if (dest) {
    dest->cd();
    nbytes += histo.Write("histo");
    for (int i=0; i<fNdim; i++) nbytes += fAxis[i].Write();
    printf("   %s stored in %s\n",GetName(),dest->GetPath());
    dest->cd("..");
  }
  return nbytes;
}
//=============================================================================
AliUnicorHN *AliUnicorHN::ProjectAlong(const char *nam, Int_t dim, Int_t first, Int_t last) {

  // Reduce dimension dim by summing up its bins between first and last. 
  // Use root convention: bin=1 is the first bin, bin=nbins is the last. 
  // last=0 means till the last bin
  // Return the resulting fNdim-1 dimensional histogram. 

  if (dim<0 || dim>fNdim-1) return 0;
  if (last<=0) last = fNbins[dim];

  // create new (reduced) histogram

  TAxis *ax[fgkMaxNdim-1] = {0};
  int n=0;
  for (int i=0; i<fNdim; i++) if (i!=dim) ax[n++] = GetAxis(i);
  const char *name = strlen(nam)? nam : Form("%s_wo%d",GetName(),dim);
  const char *eame = Form("%s_error",name); 
  AliUnicorHN *his = new AliUnicorHN(name,fNdim-1,ax); // result histogram
  AliUnicorHN *eis = new AliUnicorHN(eame,fNdim-1,ax); // temporary storage for errors

  // sum up the content and errors squared

  int k[fgkMaxNdim] = {0}; // old hist multiindex
  int m[fgkMaxNdim] = {0}; // new hist multiindex
  for (int i=0; i<GetNbinsX(); i++) {
    OneToMul(i,k);
    if (k[dim]+1<first) continue;
    if (k[dim]+1>last) continue;
    n = 0;
    for (int j=0; j<fNdim; j++) if (j!=dim) m[n++] = k[j];
    n = his->MulToOne(m);
    his->AddBinContent(n+1,GetBinContent(i+1));
    eis->AddBinContent(n+1,GetBinError(i+1)*GetBinError(i+1));
  }

  // combine content and errors in one histogram

  for (int i=0; i<his->GetNbinsX(); i++) {
    his->SetBinError(i+1,sqrt(eis->GetBinContent(i+1)));
  }

  his->SetLineColor(this->GetLineColor());
  his->SetFillColor(this->GetFillColor());
  his->SetMarkerColor(this->GetMarkerColor());
  his->SetMarkerStyle(this->GetMarkerStyle());

  // some cleanup

  delete eis;
  return his;
}
//=============================================================================
TH1D *AliUnicorHN::ProjectOn(const char *nam, Int_t dim, const Int_t * const first, 
		     const Int_t * const last) const {

  // Project on dimension dim. Use only bins between first[i] and last[i]. 
  // Use root convention: bin=1 is the first bin, bin=nbins is the last. 
  // first[i]=0 means from the first bin
  // last[i]=0 means till the last bin

  if (dim<0 || dim>fNdim-1) return 0;

  // calculate the projection; lowest index 0

  double *yy = new double[fNbins[dim]]; // value
  double *ey = new double[fNbins[dim]]; // error
  for (int i=0; i<fNbins[dim]; i++) yy[i]=0;
  for (int i=0; i<fNbins[dim]; i++) ey[i]=0;
  Int_t *k = new Int_t[fNdim];
  for (int i=0; i<GetNbinsX(); i++) {
    OneToMul(i,k);
    int isgood = 1; // bin within the range?
    for (int j=0; j<fNdim; j++) {
      if (first) if (first[j]) if (k[j]+1<first[j]) {isgood = 0; break;}
      if (last)  if (last[j])  if (k[j]+1>last[j])  {isgood = 0; break;}
    }
    if (isgood) {
      yy[k[dim]]+=GetBinContent(i+1);
      ey[k[dim]]+=GetBinError(i+1)*GetBinError(i+1);
    }
  }

  // make the projection histogram
  // use name nam if specified; otherwise generate one

  TH1D *his;
  const char *name = strlen(nam)? nam : Form("%s_proj%d",GetName(),dim);
  if (fAxis[dim].IsVariableBinSize()) 
    his = new TH1D(name,name,fNbins[dim],fAxis[dim].GetXbins()->GetArray());
  else 
    his = new TH1D(name,name,fNbins[dim],fAxis[dim].GetXmin(),fAxis[dim].GetXmax());
  his->SetXTitle(fAxis[dim].GetTitle());
  // or his->GetXaxis()->ImportAttributes(ax);
  his->Sumw2();
  his->SetLineColor(this->GetLineColor());
  his->SetFillColor(this->GetFillColor());
  his->SetMarkerColor(this->GetMarkerColor());
  his->SetMarkerStyle(this->GetMarkerStyle());
  for (int i=0; i<his->GetNbinsX(); i++) {
    his->SetBinContent(i+1,yy[i]);
    his->SetBinError(i+1,sqrt(ey[i]));
  }

  // some cleanup

  delete [] yy;
  delete [] ey;
  delete [] k;
  //  if (name!=nam) delete [] name;

  return his;
}
//=============================================================================
TH1D *AliUnicorHN::ProjectOn(const char *nam, Int_t dim, const Double_t * const first, 
		     const Double_t * const last) const {

  // Project on dimension dim. Use only bins between first[i] and last[i]. 

  if (dim<0 || dim>fNdim-1) return 0;
  Int_t kfirst[fgkMaxNdim] = {0};
  Int_t klast[fgkMaxNdim] = {0};
  for (int i=0; i<fNdim; i++) {
    kfirst[i] = fAxis[i].FindFixBin(first[i]);
    klast[i] = fAxis[i].FindFixBin(last[i]);
  }
  return ProjectOn(nam,dim,kfirst,klast);
}
//=============================================================================
TH2D *AliUnicorHN::ProjectOn(const char *nam, Int_t dim0, Int_t dim1, const Int_t * 
		     const first, const Int_t * const last) const {

  // Project on dim1 vs dim0. Use only bins between first[i] and last[i]. 
  // Use root convention: bin=1 is the first bin, bin=nbins is the last. 
  // first[i]=0 means from the first bin
  // last[i]=0 means till the last bin

  if (dim0<0 || dim0>fNdim-1) return 0;
  if (dim1<0 || dim1>fNdim-1) return 0;

  // calculate the projection

  double **yy = new double*[fNbins[dim0]]; // value
  double **ey = new double*[fNbins[dim0]]; // error
  for (int i=0; i<fNbins[dim0]; i++) yy[i] = new double[fNbins[dim1]]; 
  for (int i=0; i<fNbins[dim0]; i++) ey[i] = new double[fNbins[dim1]]; 
  for (int i=0; i<fNbins[dim0]; i++) for (int j=0; j<fNbins[dim1]; j++) yy[i][j]=0;
  for (int i=0; i<fNbins[dim0]; i++) for (int j=0; j<fNbins[dim1]; j++) ey[i][j]=0;
  Int_t *k = new Int_t[fNdim];
  for (int i=0; i<GetNbinsX(); i++) {
    OneToMul(i,k);
    int isgood = 1; // bin within the range?
    for (int j=0; j<fNdim; j++) {
      if (first) if (first[j]) if (k[j]+1<first[j]) {isgood = 0; break;}
      if (last)  if (last[j])  if (k[j]+1>last[j])  {isgood = 0; break;}
    }
    if (isgood) {
      yy[k[dim0]][k[dim1]]+=GetBinContent(i+1);
      ey[k[dim0]][k[dim1]]+=GetBinError(i+1)*GetBinError(i+1);
    }
  }

  // make the projection histogram
  // use name nam if specified; otherwise generate one

  TH2D *his=0;
  const char *name = strlen(nam)? nam : Form("%s_proj%dvs%d",GetName(),dim1,dim0);
  if (fAxis[dim0].IsVariableBinSize() && fAxis[dim1].IsVariableBinSize()) 
    his = new TH2D(name,name,
		   fNbins[dim0],fAxis[dim0].GetXbins()->GetArray(),
		   fNbins[dim1],fAxis[dim1].GetXbins()->GetArray());
  else if (!fAxis[dim0].IsVariableBinSize() && fAxis[dim1].IsVariableBinSize()) 
    his = new TH2D(name,name,
		   fNbins[dim0],fAxis[dim0].GetXmin(),fAxis[dim0].GetXmax(),
		   fNbins[dim1],fAxis[dim1].GetXbins()->GetArray());
  else if (fAxis[dim0].IsVariableBinSize() && !fAxis[dim1].IsVariableBinSize()) 
    his = new TH2D(name,name,
		   fNbins[dim0],fAxis[dim0].GetXbins()->GetArray(),
		   fNbins[dim1],fAxis[dim1].GetXmin(),fAxis[dim1].GetXmax());
  else if (!fAxis[dim0].IsVariableBinSize() && !fAxis[dim1].IsVariableBinSize()) 
    his = new TH2D(name,name,
		   fNbins[dim0],fAxis[dim0].GetXmin(),fAxis[dim0].GetXmax(),
		   fNbins[dim1],fAxis[dim1].GetXmin(),fAxis[dim1].GetXmax());
  his->SetXTitle(fAxis[dim0].GetTitle());
  his->SetYTitle(fAxis[dim1].GetTitle());
  // or his->GetXaxis()->ImportAttributes(ax0);
  // or his->GetYaxis()->ImportAttributes(ax1);
  his->Sumw2();
  his->SetLineColor(this->GetLineColor());
  his->SetFillColor(this->GetFillColor());
  his->SetMarkerColor(this->GetMarkerColor());
  his->SetMarkerStyle(this->GetMarkerStyle());
  for (int i=0; i<his->GetNbinsX(); i++) 
  for (int j=0; j<his->GetNbinsY(); j++) {
    his->SetBinContent(i+1,j+1,yy[i][j]);
    his->SetBinError(i+1,j+1,sqrt(ey[i][j]));
  }

  // some cleanup

  for (int i=0; i<fNbins[dim0]; i++) delete [] yy[i];
  for (int i=0; i<fNbins[dim0]; i++) delete [] ey[i];
  delete [] yy;
  delete [] ey;
  delete [] k;
  //  if (name!=nam) delete [] name;

  return his;
}
//=============================================================================
TH2D *AliUnicorHN::ProjectOn(const char *nam, Int_t dim0, Int_t dim1, const Double_t * 
		     const first, const Double_t * const last) const {

  // Project on dim1 vs dim0. Use only bins between first[i] and last[i]. 

  if (dim0<0 || dim0>fNdim-1) return 0;
  if (dim1<0 || dim1>fNdim-1) return 0;
  Int_t kfirst[fgkMaxNdim] = {0};
  Int_t klast[fgkMaxNdim] = {0};
  for (int i=0; i<fNdim; i++) {
    kfirst[i] = fAxis[i].FindFixBin(first[i]);
    klast[i] = fAxis[i].FindFixBin(last[i]);
  }
  return ProjectOn(nam,dim0,dim1,kfirst,klast);
}
//=============================================================================
