#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TPaveStats.h>
#endif

TH1* GetBaseHisto(TPad* pad=0);
TFrame* GetFrame(TPad* pad=0);
TPaveStats* SetStatPad(TH1* hst,float x1,float x2,float y1,float y2,Int_t stl=-1,Int_t col=-1);
TPaveStats* GetStatPad(TH1* hst);
void SetHStyle(TH1* hst,int col=kRed,int mark=20,float mrsize=0.7);
void SetGStyle(TGraph* hst,int col=kRed,int mark=20,float mrsize=0.7);
TH1* Cumulate(TH1* histo, Bool_t doErr=kTRUE, const char* copyName=0);
TLatex* AddLabel(const char*txt,float x=0.1,float y=0.9,int color=kBlack,float size=0.04);
void wAv(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);
void wSum(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);


void SaveCanvas(){
  // This dummy method allows loading save canvas from the command line as root SaveCanvas.C+
  std::cout << "Save Canvas lib loaded" << std::endl;
  
}

void SaveCanvas(TCanvas* canv,const char* path="canv",const Option_t *option="ecg")
{
  TString name;
  TString opt = option; opt.ToLower();
  //
  TString name0 = path;
  if (name0.IsNull()) name0 = "defCanv";
  //
  TFrame* fr = GetFrame();
  if (fr) fr->SetBorderMode(0);
  if (opt.Contains("c")) {
    name = name0; name+=".C";
    canv->Print(name.Data());
  }
  //
  if (opt.Contains("e")) {
    name = name0; name+=".eps";
    canv->Print(name.Data());
  }
  //
  if (opt.Contains("g")) {
    name = name0; name+=".gif";
    canv->Print(name.Data());
  }
  //
  if (opt.Contains("p")) {
    name = name0; name+=".ps";
    canv->Print(name.Data());
  }
  //
}


TLatex* AddLabel(const char*txt,float x,float y,int color,float size)
{
  TLatex* lt = new TLatex(x,y,txt); 
  lt->SetNDC(); 
  lt->SetTextColor(color);
  lt->SetTextSize(size);
  lt->Draw();
  return lt;
}

TH1* GetBaseHisto(TPad* pad)
{
  if (!pad) pad = (TPad*)gPad;
  if (!pad) return 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TH1* hst=0;
  for (int i=0;i<size;i++) {
    TObject* obj = lst->At(i);
    if (!obj) continue;
    if (obj->InheritsFrom("TH1")) {hst = (TH1*)obj; break;}
  }
  return hst;
}

TFrame* GetFrame(TPad* pad)
{
  if (!pad) pad = (TPad*)gPad;
  if (!pad) return 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TFrame* frm=0;
  for (int i=0;i<size;i++) {
    TObject* obj = lst->At(i);
    if (!obj) continue;
    if (obj->InheritsFrom("TFrame")) {frm = (TFrame*)obj; break;}
  }
  return frm;
}

TPaveStats* GetStatPad(TH1* hst)
{
  TList *lst = hst->GetListOfFunctions();
  if (!lst) return 0;
  int nf = lst->GetSize();
  for (int i=0;i<nf;i++) {
    TPaveStats *fnc = (TPaveStats*) lst->At(i);
    if (fnc->InheritsFrom("TPaveStats")) return fnc;
  }
  return 0;
  //
}


TPaveStats* SetStatPad(TH1* hst,float x1,float x2,float y1,float y2, Int_t stl, Int_t col)
{
  TPaveStats* pad = GetStatPad(hst);
  if (!pad) return 0;
  pad->SetX1NDC( x1 );
  pad->SetX2NDC( x2 );
  pad->SetY1NDC( y1 );
  pad->SetY2NDC( y2 );
  if (stl>=0) pad->SetFillStyle(stl);
  if (col>=0) pad->SetFillColor(col);
  //
  gPad->Modified();
  return pad;
}

void SetHStyle(TH1* hst,int col,int mark,float mrsize)
{
  hst->SetLineColor(col);
  hst->SetMarkerColor(col);
  //  hst->SetFillColor(col);
  hst->SetMarkerStyle(mark);
  hst->SetMarkerSize(mrsize);
  TList *lst = hst->GetListOfFunctions();
  if (lst) {
    int nf = lst->GetSize();
    for (int i=0;i<nf;i++) {
      TObject *fnc = lst->At(i);
      if (fnc->InheritsFrom("TF1")) {
	((TF1*)fnc)->SetLineColor(col);
	((TF1*)fnc)->SetLineWidth(1);
	((TF1*)fnc)->ResetBit(TF1::kNotDraw);
      }
      else if (fnc->InheritsFrom("TPaveStats")) {
	((TPaveStats*)fnc)->SetTextColor(col);
      }
    }
  }
}

void SetGStyle(TGraph* hst,int col,int mark,float mrsize)
{
  hst->SetLineColor(col);
  hst->SetMarkerColor(col);
  hst->SetFillColor(col);
  hst->SetMarkerStyle(mark);
  hst->SetMarkerSize(mrsize);
  TList *lst = hst->GetListOfFunctions();
  if (lst) {
    int nf = lst->GetSize();
    for (int i=0;i<nf;i++) {
      TObject *fnc = lst->At(i);
      if (fnc->InheritsFrom("TF1")) {
	((TF1*)fnc)->SetLineColor(col);
	((TF1*)fnc)->SetLineWidth(1);
	((TF1*)fnc)->ResetBit(TF1::kNotDraw);
      }
      else if (fnc->InheritsFrom("TPaveStats")) {
	((TPaveStats*)fnc)->SetTextColor(col);
      }
    }
  }
}

TH1* Cumulate(TH1* histo, Bool_t doErr, const char* copyName)
{
  // create cumulative histo
  TString nname = copyName;
  if (nname.IsNull()) {
    nname = histo->GetName();
    nname += "_cml";
  }
  TH1* cml = (TH1*) histo->Clone(nname.Data());
  int nb = histo->GetNbinsX();
  double sm = 0;
  double sme = 0;
  //
  for (int i=1;i<=nb;i++) {
    sm += histo->GetBinContent(i);
    cml->SetBinContent(i,sm);
    if (!doErr) continue;
    double ee = histo->GetBinError(i);
    sme += ee*ee;
    cml->SetBinError(i, sme>0 ? TMath::Sqrt(sme):0.);
  }
  return cml;
}


void wAv(double v1,double v2, double err1,double err2, double* wv,double *we)
{
  // weighted average
  double sum=0,err=0;
  if (err1<=0 || err2<=0) {
    sum = v1+v2;
  }
  else {
    sum = v1/(err1*err1) + v2/(err2*err2);
    err = 1./(err1*err1) + 1./(err2*err2);
    sum /= err;
    err = 1./TMath::Sqrt(err);
  }
  printf("wAv %+e(%+e) | %+e(%+e) -> %+e +- %e\n",v1,err1,v2,err2,sum,err);
  if (wv) *wv = sum;
  if (we) *we = err;
}


void wSum(double v1,double v2, double err1,double err2, double* wv,double *we)
{
  // sum with errors
  double sum=0,err=0;
  sum = v1+v2;
  if (err1>0 && err2>0) err = TMath::Sqrt(err1*err1 + err2*err2);
  printf("wSum %+e(%+e) + %+e(%+e) -> %+e +- %e\n",v1,err1,v2,err2,sum,err);
  if (wv) *wv = sum;
  if (we) *we = err;
}

//_________________________________________________
Double_t* defLogAx(double xMn,double xMx, int nbin)
{
  // get array for log axis
  if (xMn<=0 || xMx<=xMn || nbin<2) {
    printf("Wrong axis request: %f %f %d\n",xMn,xMx,nbin);
    return 0;
  }
  double dx = log(xMx/xMn)/nbin;
  double *xax = new Double_t[nbin+1];
  for (int i=0;i<=nbin;i++) xax[i]= xMn*exp(dx*i);
  return xax;
}

//_________________________________________________
TH2* manip2d(TH2* histo, Bool_t cumul, const char* newname, Bool_t over, Bool_t under)
{
  // normalize each XYbin to total of X bin
  // cumul=ktRUE: cumulative of each X bin, normalized to 1
  // over, under - include over/under flows
  //
  if (!histo) return 0;
  if (!histo->InheritsFrom("TH2")) {printf("Not a TH2\n"); return 0;}
  //
  TString nmstr = newname;
  if (nmstr.IsNull()) nmstr = Form("%s_nrm%s",histo->GetName(),cumul ? "_cml" : "");
  TH2* cln = (TH2*)histo->Clone(nmstr.Data());
  //
  TAxis* xax = histo->GetXaxis();
  TAxis* yax = histo->GetYaxis();  
  int yb0 = under ? 0:1;
  int yb1 = over  ? yax->GetNbins() : 1+yax->GetNbins();;
  //
  for (int ix=0;ix<=xax->GetNbins()+1;ix++) {
    double totb = 0;
    for (int iy=0;iy<=yax->GetNbins()+1;iy++) {
      double val = histo->GetBinContent(ix,iy);
      totb  += val;
    }
    if (totb<=0) continue;
    //
    // normalize X bin to 1    
    for (int iy=yb0;iy<=yb1;iy++) {
      double val = histo->GetBinContent(ix,iy);
      double err = histo->GetBinError(ix,iy);
      cln->SetBinContent(ix,iy,val/totb);
      cln->SetBinError(ix,iy,err/totb);      
    }
    //
    if (!cumul) continue;
    //
    for (int iy=yb0+1;iy<=yb1;iy++) {
      double val = cln->GetBinContent(ix,iy);
      double err = cln->GetBinError(ix,iy);
      double valP = cln->GetBinContent(ix,iy-1);
      double errP = cln->GetBinError(ix,iy-1);
      err = err*err+errP*errP;
      err = err >0 ? TMath::Sqrt(err) : 0;
      cln->SetBinContent(ix,iy,val + valP);
      cln->SetBinError(ix,iy,err);      
    }
  }
  //
  return cln;
}
