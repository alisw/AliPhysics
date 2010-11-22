/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercialf purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id: AliTRDclusterResolution.cxx */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD cluster error parameterization                                        //
//                                                                           //
// This class is designed to produce the reference plots for a detailed study//
// and parameterization of TRD cluster errors. The following effects are taken//
// into account :                                                            //
//   - dependence with the total charge of the cluster                       //
//   - dependence with the distance from the center pad. This is monitored 
// for each layer individually since the pad size varies with layer
//   - dependence with the drift length - here the influence of anisochronity 
// and diffusion are searched
//   - dependence with the distance to the anode wire - anisochronity effects
//   - dependence with track angle (for y resolution)
// The correlation between effects is taken into account. 
// 
// Since magnetic field plays a very important role in the TRD measurement 
// the ExB correction is forced by the setter function SetExB(Int_t). The 
// argument is the detector index, if none is specified all will be 
// considered.
// 
// Two cases are of big importance.
//   - comparison with MC
//   - comparison with Kalman fit. In this case the covariance matrix of the
// Kalman fit are needed.
// 
// The functionalities implemented in this class are based on the storage 
// class AliTRDclusterInfo.
// 
// The Method
// ----------
// 
// The method to disentangle s_y and s_x is based on the relation (see also fig.)
// BEGIN_LATEX
// #sigma^{2} = #sigma^{2}_{y} + tg^{2}(#alpha_{L})*#sigma^{2}_{x_{d}} + tg^{2}(#phi-#alpha_{L})*(#sigma^{2}_{x_{d}}+#sigma^{2}_{x_{c}})
// END_LATEX
// with
// BEGIN_LATEX
// #sigma^{2}_{x_{c}} #approx 0 
// END_LATEX
// we suppose the chamber is well calibrated for t_{0} and aligned in
// radial direction. 
//
// Clusters can be radially shifted due to three causes:
//   - globally shifted - due to residual misalignment/miscalibration(t0)
//   - locally shifted - due to different local drift velocity from the mean
//   - randomly shifted - due to neighboring (radial direction) clusters 
// charge induced by asymmetry of the TRF.
//
// We estimate this effects by the relations:
// BEGIN_LATEX
// #mu_{y} = tg(#alpha_{L})*#Delta x_{d}(...) + tg(#phi-#alpha_{L})*(#Delta x_{c}(...) + #Delta x_{d}(...))
// END_LATEX
// where
// BEGIN_LATEX
// #Delta x_{d}(...) = (<v_{d}> + #delta v_{d}(x_{d}, d)) * (t + t^{*}(Q))
// END_LATEX
// and we specified explicitely the variation of drift velocity parallel 
// with the track (x_{d}) and perpendicular to it due to anisochronity (d).
// 
// For estimating the contribution from asymmetry of TRF the following
// parameterization is being used
// BEGIN_LATEX
// t^{*}(Q) = #delta_{0} * #frac{Q_{t+1} - Q_{t-1}}{Q_{t-1} + Q_{t} + Q_{t+1}}
// END_LATEX
//
//
// Clusters can also be r-phi shifted due to:
//   - wrong PRF or wrong cuts at digits level
//The following correction is applied :
// BEGIN_LATEX
// <#Delta y> = a + b * sin(c*y_{pw})
// END_LATEX

// The Models
//
//   Parameterization against total charge
//
// Obtained for B=0T at phi=0. All other effects integrated out.
// BEGIN_LATEX
// #sigma^{2}_{y}(Q) = #sigma^{2}_{y}(...) + b(#frac{1}{Q} - #frac{1}{Q_{0}}) 
// END_LATEX
// For B diff 0T the error of the average ExB correction error has to be subtracted !! 
//
//   Parameterization Sx
//
// The parameterization of the error in the x direction can be written as
// BEGIN_LATEX
// #sigma_{x} = #sigma_{x}^{||} + #sigma_{x}^{#perp}
// END_LATEX
//
// where the parallel component is given mainly by the TRF width while 
// the perpendicular component by the anisochronity. The model employed for 
// the parallel is gaus(0)+expo(3) with the following parameters
// 1  C   5.49018e-01   1.23854e+00   3.84540e-04  -8.21084e-06
// 2  M   7.82999e-01   6.22531e-01   2.71272e-04  -6.88485e-05
// 3  S   2.74451e-01   1.13815e+00   2.90667e-04   1.13493e-05
// 4  E1  2.53596e-01   1.08646e+00   9.95591e-05  -2.11625e-05
// 5  E2 -2.40078e-02   4.26520e-01   4.67153e-05  -2.35392e-04
//
// and perpendicular to the track is pol2 with the parameters
//
// Par_0 = 0.190676 +/- 0.41785
// Par_1 = -3.9269  +/- 7.49862
// Par_2 = 14.7851  +/- 27.8012
//
//   Parameterization Sy
//
// The parameterization of the error in the y direction along track uses
// BEGIN_LATEX
// #sigma_{y}^{||} = #sigma_{y}^{0} -a*exp(1/(x-b))
// END_LATEX
//
// with following values for the parameters:
// 1  sy0 2.60967e-01   2.99652e-03   7.82902e-06  -1.89636e-04
// 2  a  -7.68941e+00   1.87883e+00   3.84539e-04   9.38268e-07
// 3  b  -3.41160e-01   7.72850e-02   1.63231e-05   2.51602e-05
//
//==========================================================================
// Example how to retrive reference plots from the task
// void steerClErrParam(Int_t fig=0)
// {
//   gSystem->Load("libANALYSIS.so");
//   gSystem->Load("libTRDqaRec.so");
// 
//   // initialize DB manager
//   AliCDBManager *cdb = AliCDBManager::Instance();
//   cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//   cdb->SetRun(0);
//   // initialize magnetic field.
//   AliMagFCheb *field=new AliMagFCheb("Maps","Maps", 2, 1., 10., AliMagFCheb::k5kG);
//   AliTracker::SetFieldMap(field, kTRUE);
// 
//   AliTRDclusterResolution *res = new AliTRDclusterResolution();
//   res->SetMCdata();
//   res->Load("TRD.TaskClErrParam.root");
//   res->SetExB();  
//   res->SetVisual(); 
//   //res->SetSaveAs();
//   res->SetProcessCharge(kFALSE);
//   res->SetProcessCenterPad(kFALSE);
//   //res->SetProcessMean(kFALSE);
//   res->SetProcessSigma(kFALSE);
//   if(!res->PostProcess()) return;
//   new TCanvas;
//   res->GetRefFigure(fig);
// }
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDclusterResolution.h"
#include "AliTRDresolution.h"
#include "AliTRDinfoGen.h"
#include "info/AliTRDclusterInfo.h"
#include "info/AliTRDeventInfo.h"

#include "AliTRDcalibDB.h"
#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"
#include "AliTRDCommonParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"

#include "AliESDEvent.h"
#include "AliCDBManager.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TGeoGlobalMagField.h"
#include <TGeoMatrix.h>
#include "TObjArray.h"
#include "TTree.h"
#include "TH2I.h"
#include "TH3S.h"
#include "TAxis.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TLine.h"

ClassImp(AliTRDclusterResolution)

const Float_t AliTRDclusterResolution::fgkTimeBinLength = 1./ AliTRDCommonParam::Instance()->GetSamplingFrequency();
//_______________________________________________________
AliTRDclusterResolution::AliTRDclusterResolution()
  : AliTRDrecoTask()
  ,fCanvas(NULL)
  ,fInfo(NULL)
  ,fResults(NULL)
  ,fSubTaskMap(0)
  ,fUseCalib(7)
  ,fDet(-1)
  ,fCol(-1)
  ,fRow(-1)
  ,fExB(0.)
  ,fDt(0.)
  ,fDl(0.)
  ,fVdrift(1.5)
  ,fT0(0.)
  ,fGain(1.)
  ,fXch(0.)
  ,fZch(0.)
  ,fH(0.)
  ,fDyRange(1.3)
  ,fLy(0)
  ,fT(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
{
// Constructor
  SetNameTitle("ClErrCalib", "Cluster Error Parameterization");
  memset(fR, 0, 4*sizeof(Float_t));
  memset(fP, 0, 4*sizeof(Float_t));
}

//_______________________________________________________
AliTRDclusterResolution::AliTRDclusterResolution(const char *name)
  : AliTRDrecoTask(name, "Cluster Error Parameterization")
  ,fCanvas(NULL)
  ,fInfo(NULL)
  ,fResults(NULL)
  ,fSubTaskMap(0)
  ,fUseCalib(7)
  ,fDet(-1)
  ,fCol(-1)
  ,fRow(-1)
  ,fExB(0.)
  ,fDt(0.)
  ,fDl(0.)
  ,fVdrift(1.5)
  ,fT0(0.)
  ,fGain(1.)
  ,fXch(0.)
  ,fZch(0.)
  ,fH(0.)
  ,fDyRange(1.3)
  ,fLy(0)
  ,fT(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
{
// Constructor

  memset(fR, 0, 4*sizeof(Float_t));
  memset(fP, 0, 4*sizeof(Float_t));

  // By default register all analysis
  // The user can switch them off in his steering macro
  SetProcess(kYRes);
  SetProcess(kYSys);
  SetProcess(kMean);
  SetProcess(kSigm);
}

//_______________________________________________________
AliTRDclusterResolution::~AliTRDclusterResolution()
{
// Destructor

  if(fCanvas) delete fCanvas;
  if(fResults){
    fResults->Delete();
    delete fResults;
  }
}

//_______________________________________________________
void AliTRDclusterResolution::UserCreateOutputObjects()
{
// Build and post histo container.
// Actual population of the container with histo is done in function Histos.

  if(!fContainer) fContainer = new TObjArray(kNtasks);
 //fContainer->SetOwner(kTRUE);
  PostData(1, fContainer);
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::GetRefFigure(Int_t ifig)
{
// Steering function to retrieve performance plots

  if(!fResults) return kFALSE;
  TLegend *leg = NULL;
  TList *l = NULL;
  TObjArray *arr = NULL;
  TTree *t = NULL;
  TH2 *h2 = NULL;TH1 *h1 = NULL;
  TGraphErrors *gm(NULL), *gs(NULL), *gp(NULL);
  switch(ifig){
  case kYRes:
    if(!(arr = (TObjArray*)fResults->At(kYRes))) break;
    if(!(gm = (TGraphErrors*)arr->At(0))) break;
    if(!(gs = (TGraphErrors*)arr->At(1))) break;
    if(!(gp = (TGraphErrors*)arr->At(2))) break;
    leg= new TLegend(.7, .7, .9, .95);
    leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    gs->Draw("apl"); leg->AddEntry(gs, "Sigma / Resolution", "pl");
    gs->GetHistogram()->GetYaxis()->SetRangeUser(-50., 700.);
    gs->GetHistogram()->SetXTitle("Q [a.u.]");
    gs->GetHistogram()->SetYTitle("y - x tg(#alpha_{L}) [#mum]");
    gm->Draw("pl");leg->AddEntry(gm, "Mean / Systematics", "pl");
    gp->Draw("pl");leg->AddEntry(gp, "Abundance / Probability", "pl");
    leg->Draw();
    return kTRUE;
  case kYSys:
    if(!(arr = (TObjArray*)fResults->At(kYSys))) break;
    gPad->Divide(2, 1); l = gPad->GetListOfPrimitives();
    ((TVirtualPad*)l->At(0))->cd();
    ((TTree*)arr->At(0))->Draw(Form("y:t>>h(%d, -0.5, %f, 51, -.51, .51)",AliTRDseedV1::kNtb, AliTRDseedV1::kNtb-0.5),
            "m[0]*(ly==0&&abs(m[0])<1.e-1)", "colz");
    ((TVirtualPad*)l->At(1))->cd();
    leg= new TLegend(.7, .7, .9, .95);
    leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    leg->SetHeader("TRD Plane"); 
    for(Int_t il = 1; il<=AliTRDgeometry::kNlayer; il++){
      if(!(gm = (TGraphErrors*)arr->At(il))) return kFALSE;
      gm->Draw(il>1?"pc":"apc"); leg->AddEntry(gm, Form("%d", il-1), "pl");
      if(il>1) continue;
      gm->GetHistogram()->SetXTitle("t_{drift} [tb]");
      gm->GetHistogram()->SetYTitle("#sigma_{y}(x|cen=0) [#mum]");
      gm->GetHistogram()->GetYaxis()->SetRangeUser(150., 500.);
    }
    leg->Draw();
    return kTRUE;
  case kSigm:
    if(!(t = (TTree*)fResults->At(kSigm))) break;
    t->Draw("z:t>>h2x(23, 0.1, 2.4, 25, 0., 2.5)","sx*(1)", "lego2fb");
    h2 = (TH2F*)gROOT->FindObject("h2x");
    printf("  const Double_t sx[24][25]={\n");
    for(Int_t ix=1; ix<=h2->GetNbinsX(); ix++){
      printf("    {");
      for(Int_t iy=1; iy<h2->GetNbinsY(); iy++){
        printf("%6.4f ", h2->GetBinContent(ix, iy));
      }
      printf("%6.4f},\n", h2->GetBinContent(ix, h2->GetNbinsY()));
    }
    printf("  };\n");
    gPad->Divide(2, 1, 1.e-5, 1.e-5); l = gPad->GetListOfPrimitives();
    ((TVirtualPad*)l->At(0))->cd();
    h1 = h2->ProjectionX("hsx_pxx"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#sigma_{x}> [#mum]");
    h1->SetXTitle("t_{drift} [#mus]");
    h1->GetXaxis()->SetRange(2, AliTRDseedV1::kNtb-1); h1->Draw("pc");

    t->Draw("z:t>>h2y(23, 0.1, 2.4, 25, 0., 2.5)","sy*(1)", "lego2fb");
    h2 = (TH2F*)gROOT->FindObject("h2y");
    printf("  const Double_t sy[24][25]={\n");
    for(Int_t ix=1; ix<=h2->GetNbinsX(); ix++){
      printf("    {");
      for(Int_t iy=1; iy<h2->GetNbinsY(); iy++){
        printf("%6.4f ", h2->GetBinContent(ix, iy));
      }
      printf("%6.4f},\n", h2->GetBinContent(ix, h2->GetNbinsY()));
    }
    printf("  };\n");
    ((TVirtualPad*)l->At(1))->cd();
    h1 = h2->ProjectionX("hsy_pxx"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#sigma_{y}> [#mum]");
    h1->SetXTitle("t_{drift} [#mus]");
    h1->GetXaxis()->SetRange(2, AliTRDseedV1::kNtb-1); h1->Draw("pc");
    return kTRUE;
  case kMean:
    if(!(t = (TTree*)fResults->At(kMean))) break;
    if(!t->Draw(Form("z:t>>h2x(%d, -0.5, %3.1f, %d, 0., 2.5)", 
      AliTRDseedV1::kNtb, AliTRDseedV1::kNtb-0.5, kND),
      "dx*(1)", "goff")) break;
    h2 = (TH2F*)gROOT->FindObject("h2x");
    printf("  const Double_t dx[%d][%d]={\n", AliTRDseedV1::kNtb, kND);
    for(Int_t ix=1; ix<=h2->GetNbinsX(); ix++){
      printf("    {");
      for(Int_t iy=1; iy<h2->GetNbinsY(); iy++){
        printf("%+6.4e, ", h2->GetBinContent(ix, iy));
      }
      printf("%+6.4e},\n", h2->GetBinContent(ix, h2->GetNbinsY()));
    }
    printf("  };\n");
    gPad->Divide(2, 2, 1.e-5, 1.e-5); l = gPad->GetListOfPrimitives();
    ((TVirtualPad*)l->At(0))->cd();
    h2->Draw("lego2fb");
    ((TVirtualPad*)l->At(2))->cd();
    h1 = h2->ProjectionX("hdx_pxx"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#deltax> [#mum]");
    h1->SetXTitle("t_{drift} [tb]");
    //h1->GetXaxis()->SetRange(2, AliTRDseedV1::kNtb-1); 
    h1->Draw("pc");

    if(!t->Draw(Form("z:t>>h2y(%d, -0.5, %3.1f, %d, 0., 2.5)", 
      AliTRDseedV1::kNtb, AliTRDseedV1::kNtb-0.5, kND),
      "dy*(1)", "goff")) break;
    h2 = (TH2F*)gROOT->FindObject("h2y");
    printf("  const Double_t dy[%d][%d]={\n", AliTRDseedV1::kNtb, kND);
    for(Int_t ix=1; ix<=h2->GetNbinsX(); ix++){
      printf("    {");
      for(Int_t iy=1; iy<h2->GetNbinsY(); iy++){
        printf("%+6.4e ", h2->GetBinContent(ix, iy));
      }
      printf("%+6.4e},\n", h2->GetBinContent(ix, h2->GetNbinsY()));
    }
    printf("  };\n");
    ((TVirtualPad*)l->At(1))->cd();
    h2->Draw("lego2fb");
    ((TVirtualPad*)l->At(3))->cd();
    h1 = h2->ProjectionX("hdy_pxx"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#deltay> [#mum]");
    h1->SetXTitle("t_{drift} [tb]");
    //h1->GetXaxis()->SetRange(2, AliTRDseedV1::kNtb-1); 
    h1->Draw("pc");

    return kTRUE;
  default:
    break;
  }
  AliWarning("No container/data found.");
  return kFALSE;
}

//_______________________________________________________
TObjArray* AliTRDclusterResolution::Histos()
{
// Retrieve histograms array if already build or build it

  if(!fContainer){
    fContainer = new TObjArray(kNtasks);
    //fContainer->SetOwner(kTRUE);
  }
  if(fContainer->GetEntries() == kNtasks) return fContainer;

  TH3S *h3(NULL);TH2I *h2(NULL);
  TObjArray *arr(NULL);
  if(!HasGlobalPosition() && !LoadGlobalChamberPosition()) return NULL;
  Float_t tgt(fZch/fXch), htgt(fH*tgt);
  
  // SYSTEMATIC PLOTS
  fContainer->AddAt(arr = new TObjArray(3), kYSys);
  arr->SetName("SysY");
  // systematic plot on pw and q (dydx=ExB+h*dzdx)
  if(!(h3=(TH3S*)gROOT->FindObject(Form("Sys%s%03d", (HasMCdata()?"MC":"") ,fDet)))) {
    h3 = new TH3S(
      Form("Sys%s%03d", (HasMCdata()?"MC":""),fDet),
      Form(" Det[%d] Col[%d] Row[%d];log q [a.u.];#deltay [pw];#Delta y[cm]", fDet, fCol, fRow),
      45,   2., 6.5,            // log(q) [a.u.]
      25, -.51, .51,            // y [pw]
      60, -fDyRange, fDyRange); // dy [cm]
  } h3->Reset();
  arr->AddAt(h3, 0);
  // systematic plot on tb (only for dydx = h*tgt + exb and MPV q)
  if(!(h2 = (TH2I*)gROOT->FindObject(Form("SysTb%s%03d", (HasMCdata()?"MC":""), fDet)))){
    h2 = new TH2I(Form("SysTb%s%03d", (HasMCdata()?"MC":""), fDet),
    Form(" Det[%d] Col[%d] Row[%d];t [time bin];#Delta y[cm]", fDet, fCol, fRow),
    AliTRDseedV1::kNtb, -.5, AliTRDseedV1::kNtb-0.5,   // t [tb]
    60, -fDyRange, fDyRange);                          // dy [cm]
  } h2->Reset();
  arr->AddAt(h2, 1);
  // systematic plot on tgp and tb (for MPV q)
  if(!(h3=(TH3S*)gROOT->FindObject(Form("SysTbTgp%s%03d", (HasMCdata()?"MC":""), fDet)))){
    h3 = new TH3S(
      Form("SysTbTgp%s%03d", (HasMCdata()?"MC":""), fDet),
      Form(" Det[%d];t [time bin];tg(#phi) - h*tg(#theta) %s;#Delta y[cm]", fDet, fExB>1.e-5?"- tg(#alpha_{L})":""),
      AliTRDseedV1::kNtb, -.5, AliTRDseedV1::kNtb-0.5,   // t [tb]
      36, fExB-.18, fExB+.18,                            // tgp-h tgt-tg(aL)
      60, -fDyRange, fDyRange);                          // dy
  } h3->Reset();
  arr->AddAt(h3, 2);

  //  RESOLUTION/PULLS PLOTS
  fContainer->AddAt(arr = new TObjArray(6), kYRes);
  arr->SetName("ResY");
  // resolution plot on pw and q (for dydx=0 && B=0) np = 3 and for tb in [5, 20]
  TObjArray *arrt(NULL);
  arr->AddAt(arrt = new TObjArray(16), 0);
  arrt->SetName("PwQvsX");
  for(Int_t it(0); it<=15; it++){
    if(!(h3=(TH3S*)gROOT->FindObject(Form("Res%s%03d%02d", (HasMCdata()?"MC":"") ,fDet, it)))) {
      h3 = new TH3S(
        Form("Res%s%03d%02d", (HasMCdata()?"MC":""),fDet, it),
        Form(" Det[%d] TB[%d];log q [a.u];#deltay [pw];#Delta y[cm]", fDet, it+5),
        4,   2., 6.,            // log(q) [a.u]
        5, -.51, .51,            // y [pw]
        60, -fDyRange, fDyRange); // dy
    } h3->Reset();
    arrt->AddAt(h3, it);
  }
  // Pull plot on pw and q (for dydx=0 && B=0)
  if(!(h3=(TH3S*)gROOT->FindObject(Form("Pull%s%03d", (HasMCdata()?"MC":""), fDet)))){
    h3 = new TH3S(
      Form("Pull%s%03d", (HasMCdata()?"MC":""), fDet),
      Form(" Det[%d] Col[%d] Row[%d];log q [a.u.];#deltay [pw];#Delta y[cm]/#sigma_{y}", fDet, fCol, fRow),
      4,   2., 6.,            // log(q) [a.u]
      5, -.51, .51,            // y [pw]
      60, -4., 4.);             // dy/sy
  } h3->Reset();
  arr->AddAt(h3, 1);
  // resolution/pull plot on tb (for dydx=0 && B=0 && MPV q)
  if(!(h3 = (TH3S*)gROOT->FindObject(Form("ResPullTb%s%03d", (HasMCdata()?"MC":""), fDet)))){
    h3 = new TH3S(Form("ResPullTb%s%03d", (HasMCdata()?"MC":""), fDet),
    Form(" Det[%d] Col[%d] Row[%d];t [time bin];#Delta y[cm];#Delta y/#sigma_{y}", fDet, fCol, fRow),
    AliTRDseedV1::kNtb, -.5, AliTRDseedV1::kNtb-0.5,   // t [tb]
    60, -fDyRange, fDyRange,                           // dy [cm]
    60, -4., 4.);                                      // dy/sy
  } h3->Reset();
  arr->AddAt(h3, 2);
  // resolution plot on pw and q (for dydx=0 && B=0) np = 2
  if(!(h3=(TH3S*)gROOT->FindObject(Form("Res2%s%03d", (HasMCdata()?"MC":"") ,fDet)))) {
    h3 = new TH3S(
      Form("Res2%s%03d", (HasMCdata()?"MC":""),fDet),
      Form(" Det[%d] Col[%d] Row[%d];log q [a.u];#deltay [pw];#Delta y[cm]", fDet, fCol, fRow),
      4,   2., 6.,            // log(q) [a.u]
      5, -.51, .51,            // y [pw]
      60, -fDyRange, fDyRange); // dy
  } h3->Reset();
  arr->AddAt(h3, 3);
  // resolution plot on pw and q (for dydx=0 && B=0) np = 4
  if(!(h3=(TH3S*)gROOT->FindObject(Form("Res4%s%03d", (HasMCdata()?"MC":"") ,fDet)))) {
    h3 = new TH3S(
      Form("Res4%s%03d", (HasMCdata()?"MC":""),fDet),
      Form(" Det[%d] Col[%d] Row[%d];log q [a.u];#deltay [pw];#Delta y[cm]", fDet, fCol, fRow),
      4,   2., 6.,            // log(q) [a.u]
      5, -.51, .51,            // y [pw]
      60, -fDyRange, fDyRange); // dy
  } h3->Reset();
  arr->AddAt(h3, 4);
  // systemtic plot of tb on pw and q (for dydx=0 && B=0)
  if(!(h3=(TH3S*)gROOT->FindObject(Form("SysTbPwQ%s%03d", (HasMCdata()?"MC":"") ,fDet)))) {
    h3 = new TH3S(
      Form("SysTbPwQ%s%03d", (HasMCdata()?"MC":""),fDet),
      Form(" Det[%d] Col[%d] Row[%d];log q [a.u];#deltay [pw];t [time bin]", fDet, fCol, fRow),
      4,   2., 6.,            // log(q) [a.u]
      5, -.51, .51,            // y [pw]
      AliTRDseedV1::kNtb, -.5, AliTRDseedV1::kNtb-0.5);   // t [tb]
  } h3->Reset();
  arr->AddAt(h3, 5);



  fContainer->AddAt(arr = new TObjArray(AliTRDseedV1::kNtb), kSigm);
  arr->SetName("Resolution");
  for(Int_t it=0; it<AliTRDseedV1::kNtb; it++){
    if(!(h3=(TH3S*)gROOT->FindObject(Form("hr%s%03d_t%02d", (HasMCdata()?"MC":""), fDet, it)))){
      h3 = new TH3S(
        Form("hr%s%03d_t%02d", (HasMCdata()?"MC":""), fDet, it),
        Form(" Det[%d] t_{drift}(%2d)[bin];h*tg(#theta);tg(#phi);#Delta y[cm]", fDet, it),
        35, htgt-0.0035, htgt+0.0035, // h*tgt
        36, fExB-.18, fExB+.18,       // tgp
        60, -fDyRange, fDyRange);     // dy
    } h3->Reset();
    arr->AddAt(h3, it);
  }
  return fContainer;
}

//_______________________________________________________
void AliTRDclusterResolution::UserExec(Option_t *)
{
// Fill container histograms

  if(!(fInfo = dynamic_cast<TObjArray *>(GetInputData(1)))){
    AliError("Cluster array missing.");
    return;
  }
  if(!(fEvent = dynamic_cast<AliTRDeventInfo*>(GetInputData(2)))){
    AliError("Event Info missing.");
    return;
  }
  if(!IsCalibrated()){
    LoadCalibration();
    if(!IsCalibrated()){
      AliFatal("Loading the calibration settings failed. Check OCDB access.");
      return;
    }
    fEvent->Print();
  }
  if(!fContainer->GetEntries()) Histos();

  AliDebug(2, Form("Clusters[%d]", fInfo->GetEntriesFast()));

  Int_t det, t, np;
  Float_t x, y, z, q, dy, dydx, dzdx, cov[3], covcl[3];
  TH3S *h3(NULL); TH2I *h2(NULL);

  // define limits around ExB for which x contribution is negligible
  const Float_t kAroundZero = 3.5e-2; //(+- 2 deg)

  TObjArray *arr0 = (TObjArray*)fContainer->At(kYSys);
  TObjArray *arr1 = (TObjArray*)fContainer->At(kYRes);
  TObjArray *arr10 = (TObjArray*)arr1->At(0);
  TObjArray *arr2 = (TObjArray*)fContainer->At(kSigm);

  const AliTRDclusterInfo *cli = NULL;
  TIterator *iter=fInfo->MakeIterator();
  while((cli=dynamic_cast<AliTRDclusterInfo*>((*iter)()))){
    if((np = cli->GetNpads())>4) continue;
    cli->GetCluster(det, x, y, z, q, t, covcl);

    // select cluster according to detector region if specified
    if(fDet>=0 && fDet!=det) continue;
    if(fCol>=0 && fRow>=0){
      Int_t c,r;
      cli->GetCenterPad(c, r);
      if(TMath::Abs(fCol-c) > 5) continue;
      if(TMath::Abs(fRow-r) > 2) continue;
    }
    dy = cli->GetResolution();
    AliDebug(4, Form("det[%d] tb[%2d] q[%4.0f Log[%6.4f]] np[%d] dy[%7.2f][um] ypull[%5.2f]", det, t, q, TMath::Log(q), np, 1.e4*dy, dy/TMath::Sqrt(covcl[0])));
    
    cli->GetGlobalPosition(y, z, dydx, dzdx, &cov[0]);
    Float_t pw(cli->GetYDisplacement());

    // systematics as a function of pw and log(q)
    // only for dydx = exB + h*dzdx
    if(TMath::Abs(dydx-fExB-fH*dzdx) < kAroundZero){
      h3 = (TH3S*)arr0->At(0);
      h3->Fill(TMath::Log(q), pw, dy);
    }
    // resolution/pull as a function of pw and log(q)
    // only for dydx = 0, ExB=0
    if(TMath::Abs(fExB) < kAroundZero &&
       TMath::Abs(dydx) < kAroundZero &&
       t>=5 && t<=20 ){
      switch(np){
      case 3: // MPV np
        h3 = (TH3S*)arr10->At(t-5);
        h3->Fill(TMath::Log(q), pw, dy);
        h3 = (TH3S*)arr1->At(5);
        h3->Fill(TMath::Log(q), pw, t);
        break;
      case 2: // Min np
        h3 = (TH3S*)arr1->At(3);
        h3->Fill(TMath::Log(q), pw, dy);
        break;
      case 4: // Max np
        h3 = (TH3S*)arr1->At(4);
        h3->Fill(TMath::Log(q), pw, dy);
        break;
      }
      h3 = (TH3S*)arr1->At(1);
      h3->Fill(TMath::Log(q), pw, dy/TMath::Sqrt(covcl[0]));
    }

    // do not use problematic clusters in resolution analysis
    // TODO define limits as calibration aware (gain) !!
    //if(!AcceptableGain(fGain)) continue;
    if(q<20. || q>250.) continue;

    // systematic as a function of time bin
    // only for dydx = exB + h*dzdx and MPV q
    if(TMath::Abs(dydx-fExB-fH*dzdx) < kAroundZero){
      h2 = (TH2I*)arr0->At(1);
      h2->Fill(t, dy);
    }
    // systematic as function of tb and tgp
    // only for MPV q
    h3 = (TH3S*)arr0->At(2);
    h3->Fill(t, dydx, dy);

    // resolution/pull as a function of time bin
    // only for dydx = 0, ExB=0 and MPV q
    if(TMath::Abs(fExB) < kAroundZero &&
       TMath::Abs(dydx) < kAroundZero){
      h3 = (TH3S*)arr1->At(2);
      h3->Fill(t, dy, dy/TMath::Sqrt(covcl[0]));
    }

    // resolution as function of tb, tgp and h*tgt
    // only for MPV q
    ((TH3S*)arr2->At(t))->Fill(fH*dzdx, dydx, dy);
  }
}


//_______________________________________________________
Bool_t AliTRDclusterResolution::PostProcess()
{
// Steer processing of various cluster resolution dependences :
//
//   - process resolution dependency cluster charge
//   if(HasProcess(kYRes)) ProcessCharge();
//   - process resolution dependency on y displacement
//   if(HasProcess(kYSys)) ProcessCenterPad();
//   - process resolution dependency on drift legth and drift cell width
//   if(HasProcess(kSigm)) ProcessSigma();
//   - process systematic shift on drift legth and drift cell width
//   if(HasProcess(kMean)) ProcessMean();

  if(!fContainer) return kFALSE;
  if(!IsCalibrated()){
    AliError("Not calibrated instance.");
    return kFALSE;
  }
  TObjArray *arr = NULL;
  TTree *t=NULL;
  if(!fResults){
    TGraphErrors *g = NULL;
    fResults = new TObjArray(kNtasks);
    fResults->SetOwner();
    fResults->AddAt(arr = new TObjArray(3), kYRes);
    arr->SetOwner();
    arr->AddAt(g = new TGraphErrors(), 0);
    g->SetLineColor(kBlue); g->SetMarkerColor(kBlue);
    g->SetMarkerStyle(7); 
    arr->AddAt(g = new TGraphErrors(), 1);
    g->SetLineColor(kRed); g->SetMarkerColor(kRed);
    g->SetMarkerStyle(23); 
    arr->AddAt(g = new TGraphErrors(), 2);
    g->SetLineColor(kGreen); g->SetMarkerColor(kGreen);
    g->SetMarkerStyle(7); 

    // pad center dependence
    fResults->AddAt(arr = new TObjArray(AliTRDgeometry::kNlayer+1), kYSys);
    arr->SetOwner();
    arr->AddAt(
    t = new TTree("cent", "dy=f(y,x,ly)"), 0);
    t->Branch("ly", &fLy, "ly/B");
    t->Branch("t", &fT, "t/F");
    t->Branch("y", &fY, "y/F");
    t->Branch("m", &fR[0], "m[2]/F");
    t->Branch("s", &fR[2], "s[2]/F");
    t->Branch("pm", &fP[0], "pm[2]/F");
    t->Branch("ps", &fP[2], "ps[2]/F");
    for(Int_t il=1; il<=AliTRDgeometry::kNlayer; il++){
      arr->AddAt(g = new TGraphErrors(), il);
      g->SetLineColor(il); g->SetLineStyle(il);
      g->SetMarkerColor(il);g->SetMarkerStyle(4); 
    }


    fResults->AddAt(t = new TTree("sigm", "dy=f(dw,x,dydx)"), kSigm);
    t->Branch("t", &fT, "t/F");
    t->Branch("x", &fX, "x/F");
    t->Branch("z", &fZ, "z/F");
    t->Branch("sx", &fR[0], "sx[2]/F");
    t->Branch("sy", &fR[2], "sy[2]/F");


    fResults->AddAt(t = new TTree("mean", "dy=f(dw,x,dydx - h dzdx)"), kMean);
    t->Branch("t", &fT, "t/F");
    t->Branch("x", &fX, "x/F");
    t->Branch("z", &fZ, "z/F");
    t->Branch("dx", &fR[0], "dx[2]/F");
    t->Branch("dy", &fR[2], "dy[2]/F");
  } else {
    TObject *o = NULL;
    TIterator *iter=fResults->MakeIterator();
    while((o=((*iter)()))) o->Clear(); // maybe it is wrong but we should never reach this point
  }
  
  // process resolution dependency on charge
  if(HasProcess(kYRes)) ProcessCharge();
  
  // process resolution dependency on y displacement
  if(HasProcess(kYSys)) ProcessNormalTracks();

  // process resolution dependency on drift legth and drift cell width
  if(HasProcess(kSigm)) ProcessSigma();

  // process systematic shift on drift legth and drift cell width
  if(HasProcess(kMean)) ProcessMean();

  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::LoadCalibration()
{
// Retrieve calibration parameters from OCDB, drift velocity and t0 for the detector region specified by
// a previous call to AliTRDclusterResolution::SetCalibrationRegion().

  AliCDBManager *cdb = AliCDBManager::Instance(); // check access OCDB
  if(cdb->GetRun() < 0){
    AliError("OCDB manager not properly initialized");
    return kFALSE;
  }
  // check magnetic field
  if(!TGeoGlobalMagField::Instance() || !TGeoGlobalMagField::Instance()->IsLocked()){
    AliError("Magnetic field not available.");
    return kFALSE;
  }

  AliTRDcalibDB *fCalibration  = AliTRDcalibDB::Instance();
  AliTRDCalROC  *fCalVdriftROC(fCalibration->GetVdriftROC(fDet>=0?fDet:0))
               ,*fCalT0ROC(fCalibration->GetT0ROC(fDet>=0?fDet:0));
  const AliTRDCalDet  *fCalVdriftDet = fCalibration->GetVdriftDet();
  const AliTRDCalDet  *fCalT0Det = fCalibration->GetT0Det();

  if(IsUsingCalibParam(kVdrift)){
    fVdrift = fCalVdriftDet->GetValue(fDet>=0?fDet:0);
    if(fCol>=0 && fRow>=0) fVdrift*= fCalVdriftROC->GetValue(fCol, fRow);
  }
  fExB    = AliTRDCommonParam::Instance()->GetOmegaTau(fVdrift);
  AliTRDCommonParam::Instance()->GetDiffCoeff(fDt, fDl, fVdrift);
  if(IsUsingCalibParam(kT0)){
    fT0     = fCalT0Det->GetValue(fDet>=0?fDet:0);
    if(fCol>=0 && fRow>=0) fT0 *= fCalT0ROC->GetValue(fCol, fRow);
  }
  if(IsUsingCalibParam(kGain)) fGain = (fCol>=0 && fRow>=0)?fCalibration-> GetGainFactor(fDet, fCol, fRow):fCalibration-> GetGainFactorAverage(fDet);

  SetBit(kCalibrated);

  AliInfo(Form("Calibration D[%3d] Col[%3d] Row[%2d] : \n   t0[%5.3f] vd[%5.3f] gain[%5.3f] ExB[%f] DiffT[%f] DiffL[%f]", fDet, fCol, fRow, fT0, fVdrift, fGain, fExB, fDt, fDl));

  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::LoadGlobalChamberPosition()
{
// Retrieve global chamber position from alignment
// a previous call to AliTRDclusterResolution::SetCalibrationRegion() is mandatory.

  TGeoHMatrix *matrix(NULL);
  Double_t loc[] = {0., 0., 0.}, glb[] = {0., 0., 0.};
  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  if(!(matrix= geo->GetClusterMatrix(fDet))) {
    AliFatal(Form("Could not get transformation matrix for detector %d.", fDet));
    return kFALSE;
  }
  matrix->LocalToMaster(loc, glb);
  fXch = glb[0]+ AliTRDgeometry::AnodePos()-.5*AliTRDgeometry::AmThick() - AliTRDgeometry::DrThick();
  AliTRDpadPlane *pp(geo->GetPadPlane(fDet));
  fH = TMath::Tan(pp->GetTiltingAngle()*TMath::DegToRad());
  fZch = 0.;
  if(fRow>=0){
    fZch = pp->GetRowPos(fRow)+0.5*pp->GetLengthIPad();
  }else{
    Int_t  nrows(pp->GetNrows());
    Float_t zmax(pp->GetRow0()),
            zmin(zmax - 2 * pp->GetLengthOPad()
                - (nrows-2) * pp->GetLengthIPad()
                - (nrows-1) * pp->GetRowSpacing());
    fZch = 0.5*(zmin+zmax);
  }
  
  AliInfo(Form("Global pos. D[%3d] Col[%3d] Row[%2d] : \n   x[%6.2f] z[%6.2f] h[%+6.2f].", fDet, fCol, fRow, fXch, fZch, fH));
  SetBit(kGlobal);
  return kTRUE;
}

//_______________________________________________________
void AliTRDclusterResolution::SetCalibrationRegion(Int_t det, Int_t col, Int_t row)
{
// Set calibration region in terms of detector and pad. 
// By default detector 0 mean values are considered.

  if(det>=0 && det<AliTRDgeometry::kNdet){ 
    fDet = det;
    if(col>=0 && row>=0){
      // check pad col/row for detector
      AliTRDgeometry *geo = AliTRDinfoGen::Geometry();
      AliTRDpadPlane *pp(geo->GetPadPlane(fDet));
      if(fCol>=pp->GetNcols() ||
        fRow>=pp->GetNrows()){
        AliWarning(Form("Pad coordinates col[%d] or row[%d] incorrect for det[%d].\nLimits are max col[%d] max row[%d]. Reset to default", fCol, fRow, fDet, pp->GetNcols(), pp->GetNrows()));
        fCol = -1; fRow=-1;
      } else {
        fCol = col;
        fRow = row;
      }
    }
  } else {
    AliFatal(Form("Detector index outside range [0 %d].", AliTRDgeometry::kNdet));
  }
  return;
}

//_______________________________________________________
void AliTRDclusterResolution::SetVisual()
{
  if(fCanvas) return;
  fCanvas = new TCanvas("clResCanvas", "Cluster Resolution Visualization", 10, 10, 600, 600);
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessCharge()
{
// Resolution as a function of cluster charge.
//
// As described in the function ProcessCenterPad() the error parameterization for clusters for phi = a_L can be 
// written as:
// BEGIN_LATEX
// #sigma_{y}^{2} = #sigma_{y}^{2}|_{B=0} + tg^{2}(#alpha_{L})*#sigma_{x}^{2}
// END_LATEX
// with the contribution in case of B=0 given by:
// BEGIN_LATEX
// #sigma_{y}|_{B=0} = #sigma_{diff}*Gauss(0, s_{ly}) + #delta_{#sigma}(q)
// END_LATEX
// which further can be simplified to:
// BEGIN_LATEX
// <#sigma_{y}|_{B=0}>(q) = <#sigma_{y}> + #delta_{#sigma}(q)
// <#sigma_{y}> = #int{f(q)#sigma_{y}dq}
// END_LATEX
// The results for s_y and f(q) are displayed below:
//Begin_Html
//<img src="TRD/clusterQerror.gif">
//End_Html
// The function has to extended to accomodate gain calibration scalling and errors.
//
// Author
// Alexandru Bercuci <A.Bercuci@gsi.de>



  TObjArray *arr(NULL);
  if(!(arr = (TObjArray*)fContainer->At(kYSys))) {
    AliError("Missing systematic container");
    return;
  }
  TH3S *h3s(NULL);
  if(!(h3s = (TH3S*)arr->At(0))){
    AliError("Missing systematic histo");
    return;
  }
  // PROCESS SYSTEMATIC
  Float_t tmin(6.5), tmax(20.5), tmed(0.5*(tmin+tmax));
  TGraphErrors *g[2]; TH1 *h(NULL);
  g[0] = new TGraphErrors();
  g[0]->SetMarkerStyle(24);g[0]->SetMarkerColor(kBlue);g[0]->SetLineColor(kBlue);
  g[1] = new TGraphErrors();
  g[1]->SetMarkerStyle(24);g[1]->SetMarkerColor(kRed);g[1]->SetLineColor(kRed);
  // define model for systematic shift vs pw
  TF1 fm("fm", "[0]+[1]*sin(x*[2])", -.45,.45);
  fm.SetParameter(0, 0.); fm.SetParameter(1, 1.e-2); fm.FixParameter(2, TMath::TwoPi());
  fm.SetParNames("#deltay", "#pm#delta", "2*#pi");
  h3s->GetXaxis()->SetRange(tmin, tmax);
  if(!AliTRDresolution::Process((TH2*)h3s->Project3D("zy"), g))return;
  g[0]->Fit(&fm, "QR");
  if(fCanvas){
    g[0]->Draw("apl");
    fCanvas->Modified(); fCanvas->Update();
    h = g[0]->GetHistogram();
    h->SetTitle(fm.GetTitle());
    h->GetXaxis()->SetTitle("pw");h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle("#Delta y[cm]");h->GetYaxis()->CenterTitle();
    if(IsSaveAs()) fCanvas->SaveAs(Form("D%03d_SysNormTrack_pw.gif", fDet));
    else gSystem->Sleep(100);
  }

  // define model for systematic shift vs tb
  TF1 fx("fx", "[0]+0.1*[1]*(x-[2])", tmin, tmax);
  fx.SetParNames("#deltay", "#deltay/t", "<t>");
  fx.FixParameter(2, tmed);
  h3s->GetXaxis()->UnZoom();
  if(!AliTRDresolution::Process((TH2*)h3s->Project3D("zx"), g)) return;
  g[0]->Fit(&fx, "Q", "", tmin, tmax);
  if(fCanvas){
    g[0]->Draw("apl");
    fCanvas->Modified(); fCanvas->Update();
    h = g[0]->GetHistogram();
    h->SetTitle(fx.GetTitle());
    h->GetXaxis()->SetTitle("t [tb]");h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle("#Delta y[cm]");h->GetYaxis()->CenterTitle();
    if(IsSaveAs()) fCanvas->SaveAs(Form("D%03d_SysNormTrack_tb.gif", fDet));
    else gSystem->Sleep(100);
  }

  TH3S *h3(NULL);
  if(!(h3 = (TH3S*)fContainer->At(kYRes))) {
    AliWarning("Missing dy=f(Q) histo");
    return;
  }
  TF1 f("f", "gaus", -.5, .5);
  TAxis *ax(NULL);
  TH1 *h1(NULL);

  // compute mean error on x
  Double_t s2x = 0.; 
  for(Int_t ix=5; ix<AliTRDseedV1::kNtb; ix++){
    // retrieve error on the drift length
    s2x += AliTRDcluster::GetSX(ix);
  }
  s2x /= (AliTRDseedV1::kNtb-5); s2x *= s2x;
  //Double_t exb2 = fExB*fExB;

  arr = (TObjArray*)fResults->At(kYRes);
  TGraphErrors *gqm = (TGraphErrors*)arr->At(0);
  TGraphErrors *gqs = (TGraphErrors*)arr->At(1);
  TGraphErrors *gqp = (TGraphErrors*)arr->At(2);
  Double_t q, n = 0., entries;
  ax = h3->GetXaxis();
  for(Int_t ix=1; ix<=ax->GetNbins(); ix++){
    q = TMath::Exp(ax->GetBinCenter(ix));
    ax->SetRange(ix, ix);
    h1 = h3->Project3D("y");
    entries = h1->GetEntries();
    if(entries < 150) continue;
    h1->Fit(&f, "Q");

    // Fill sy^2 = f(q)
    Int_t ip = gqm->GetN();
    gqm->SetPoint(ip, q, 1.e4*f.GetParameter(1));
    gqm->SetPointError(ip, 0., 1.e4*f.GetParError(1));

    // correct sigma for ExB effect
    gqs->SetPoint(ip, q, 1.e4*f.GetParameter(2)/**f.GetParameter(2)-exb2*s2x)*/);
    gqs->SetPointError(ip, 0., 1.e4*f.GetParError(2)/**f.GetParameter(2)*/);

    // save probability
    n += entries;
    gqp->SetPoint(ip, q, entries);
    gqp->SetPointError(ip, 0., 0./*TMath::Sqrt(entries)*/);
  } 

  // normalize probability and get mean sy
  Double_t sm = 0., sy;
  for(Int_t ip=gqp->GetN(); ip--;){
    gqp->GetPoint(ip, q, entries);
    entries/=n;
    gqp->SetPoint(ip, q, 1.e4*entries);
    gqs->GetPoint(ip, q, sy);
    sm += entries*sy;
  }

  // error parametrization s(q) = <sy> + b(1/q-1/q0)
  TF1 fq("fq", "[0] + [1]/x", 20., 250.);
  gqs->Fit(&fq/*, "W"*/);
  printf("sm=%f [0]=%f [1]=%f\n", 1.e-4*sm, fq.GetParameter(0), fq.GetParameter(1));
  printf("  const Float_t sq0inv = %f; // [1/q0]\n", (sm-fq.GetParameter(0))/fq.GetParameter(1));
  printf("  const Float_t sqb    = %f; // [cm]\n", 1.e-4*fq.GetParameter(1));
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::ProcessNormalTracks()
{
// Resolution as a function of y displacement from pad center and drift length.
//
// Since the error parameterization of cluster r-phi position can be written as (see AliTRDcluster::SetSigmaY2()):
// BEGIN_LATEX
// #sigma_{y}^{2} = (#sigma_{diff}*Gauss(0, s_{ly}) + #delta_{#sigma}(q))^{2} + tg^{2}(#alpha_{L})*#sigma_{x}^{2} + tg^{2}(#phi-#alpha_{L})*#sigma_{x}^{2}+[tg(#phi-#alpha_{L})*tg(#alpha_{L})*x]^{2}/12
// END_LATEX
// one can see that for phi = a_L one gets the following expression:
// BEGIN_LATEX
// #sigma_{y}^{2} = #sigma_{y}^{2}|_{B=0} + tg^{2}(#alpha_{L})*#sigma_{x}^{2}
// END_LATEX
// where we have explicitely marked the remaining term in case of absence of magnetic field. Thus one can use the 
// previous equation to estimate s_y for B=0 and than by comparing in magnetic field conditions one can get the s_x.
// This is a simplified method to determine the error parameterization for s_x and s_y as compared to the one 
// implemented in ProcessSigma(). For more details on cluster error parameterization please see also 
// AliTRDcluster::SetSigmaY2()
// 
// The representation of dy=f(y_cen, x_drift| layer) can be also used to estimate the systematic shift in the r-phi 
// coordinate resulting from imperfection in the cluster shape parameterization. From the expresion of the shift derived 
// in ProcessMean() with phi=exb one gets: 
// BEGIN_LATEX
// <#Delta y>= <#delta x> * (tg(#alpha_{L})-h*dz/dx) + <#delta y - #delta x * tg(#alpha_{L})>
// <#Delta y>(y_{cen})= -h*<#delta x>(x_{drift}, q_{cl}) * dz/dx + #delta y(y_{cen}, ...)
// END_LATEX
// where all dependences are made explicit. This last expression can be used in two ways:
//   - by average on the dz/dx we can determine directly dy (the method implemented here) 
//   - by plotting as a function of dzdx one can determine both dx and dy components in an independent method.
//Begin_Html
//<img src="TRD/clusterYcorr.gif">
//End_Html
// Author
// Alexandru Bercuci <A.Bercuci@gsi.de>

  TObjArray *arr(NULL);
  TH3S *h3r(NULL), *h3t(NULL);
  if(!(arr= (TObjArray*)fContainer->At(kYRes))) {
    AliError("Missing resolution container");
    return kFALSE;
  }
  if(!(h3r = (TH3S*)arr->At(0))){
    AliError("Missing resolution pw/q histo");
    return kFALSE;
  } else if(!(Int_t)h3r->GetEntries()){
    AliError("Empty resolution pw/q histo");
    return kFALSE;
  }
  if(!(h3t = (TH3S*)arr->At(2))){
    AliError("Missing resolution t histo");
    return kFALSE;
  } else if(!(Int_t)h3t->GetEntries()){
    AliError("Empty resolution t histo");
    return kFALSE;
  }

  // local variables
  Double_t x(0.), y(0.), ex(0.), ey(0.);
  Float_t tmin(6.5), tmax(20.5), tmed(0.5*(tmin+tmax));
  TGraphErrors *g[2]; TH1 *h(NULL);
  g[0] = new TGraphErrors();
  g[0]->SetMarkerStyle(24);g[0]->SetMarkerColor(kBlue);g[0]->SetLineColor(kBlue);
  g[1] = new TGraphErrors();
  g[1]->SetMarkerStyle(24);g[1]->SetMarkerColor(kRed);g[1]->SetLineColor(kRed);

  // PROCESS RESOLUTION VS TB
  TF1 fsx("fsx", "[0]*[0]+[1]*[1]*[2]*0.1*(x-[3])", tmin, tmax);
  fsx.SetParNames("#sqrt{<#sigma^{2}(prf, q)>}(t_{med})", "D_{T}", "v_{drift}", "t_{med}");
  fsx.FixParameter(1, fDt);
  fsx.SetParameter(2, fVdrift);
  fsx.FixParameter(3, tmed);
  if(!AliTRDresolution::Process((TH2*)h3r->Project3D("yx"), g)) return kFALSE;
  for(Int_t ip(0); ip<g[1]->GetN(); ip++){
    g[1]->GetPoint(ip, x, y);ex = g[1]->GetErrorX(ip); ey = g[1]->GetErrorY(ip);
    g[1]->SetPoint(ip, x, y*y);g[1]->SetPointError(ip, ex, 2*y*ey);
  }
  g[1]->Fit(&fsx, "Q", "", tmin, tmax);
  if(fCanvas){
    g[1]->Draw("apl");
    fCanvas->Modified(); fCanvas->Update();
    h = g[1]->GetHistogram();
    h->SetTitle(fsx.GetTitle());
    h->GetXaxis()->SetTitle("t [tb]");h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle("#sigma^{2} (y) [cm^{2}]");h->GetYaxis()->CenterTitle();
    if(IsSaveAs()) fCanvas->SaveAs(Form("D%03d_ResNormTrack_tb.gif", fDet));
    else gSystem->Sleep(100);
  }

  // define model for resolution vs pw
  TF1 fg("fg", "gaus", -.5, .5); fg.FixParameter(1, 0.);
  TF1 fs("fs", "[0]*[0]*exp(-1*(x/[1])**2)+[2]", -.5, .5);
  fs.SetParNames("<#sigma^{max}(q,prf)>_{q}", "#sigma(pw)", "D_{T}^{2}*<x>");
  h3r->GetXaxis()->SetRange(tmin, tmax);
  if(!AliTRDresolution::Process((TH2*)h3r->Project3D("zy"), g, 200)) return kFALSE;
  for(Int_t ip(0); ip<g[1]->GetN(); ip++){
    g[1]->GetPoint(ip, x, y); ex = g[1]->GetErrorX(ip); ey = g[1]->GetErrorY(ip);
    g[1]->SetPoint(ip, x, y*y);g[1]->SetPointError(ip, ex, 2.*y*ey);
  }
  g[1]->Fit(&fg, "QR");
  fs.SetParameter(0, TMath::Sqrt(fg.GetParameter(0)));
  fs.SetParameter(1, fg.GetParameter(2));
  Float_t sdiff(fDt*fDt*fsx.GetParameter(2)*tmed*0.1);
  fs.SetParameter(2, sdiff);
  fs.SetParLimits(2, 0.1*sdiff, 1.9*sdiff);
  g[1]->Fit(&fs, "QR");
  if(fCanvas){
    g[1]->Draw("apl");
    fCanvas->Modified(); fCanvas->Update();
    h = g[1]->GetHistogram();
    h->SetTitle(fs.GetTitle());
    h->GetXaxis()->SetTitle("pw");h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle("#sigma^{2} (y) [cm^{2}]");h->GetYaxis()->CenterTitle();
    if(IsSaveAs()) fCanvas->SaveAs(Form("D%03d_ResNormTrack_pw.gif", fDet));
    else gSystem->Sleep(100);
  }

  AliDebug(2, Form("<s(q,prf)>[mum] = %7.3f", 1.e4*TMath::Sqrt(fsx.Eval(0.))));
  AliDebug(2, Form("<s(q)>[mum] = %7.3f", 1.e4*TMath::Sqrt(fs.Eval(-0.5)-fs.GetParameter(2))));
  AliDebug(2, Form("<s(x)>[mum] = %7.3f(prf) %7.3f(diff)", 1.e4*TMath::Sqrt(fs.GetParameter(2)), 1.e4*TMath::Sqrt(sdiff)));

  // define model for resolution vs q
  TF1 fq("fq", "[0]*[0]*exp(-1*[1]*(x-[2])**2)+[2]", 2.5, 5.5);
  fq.SetParNames("<#sigma^{max}(q,prf)>_{prf}", "slope","mean", "D_{T}^{2}*<x>");
  if(!AliTRDresolution::Process((TH2*)h3t->Project3D("yx"), g)) return kFALSE;
  for(Int_t ip(0); ip<g[1]->GetN(); ip++){
    g[1]->GetPoint(ip, x, y); ex = g[1]->GetErrorX(ip); ey = g[1]->GetErrorY(ip);
    g[1]->SetPoint(ip, x, y*y);g[1]->SetPointError(ip, ex, 2.*y*ey);
  }
  fq.SetParameter(0, 8.e-2); fq.SetParLimits(0, 0., 1.);
  fq.SetParameter(1, 1.); //fq.SetParLimits(1, -1., 0.);
  fq.SetParameter(3, sdiff); fq.SetParLimits(3, 0.1*sdiff, 1.9*sdiff);
  g[1]->Fit(&fq, "QR");
//   AliDebug(2, Form("<sq>[mum] = %7.3f", 1.e4*TMath::Sqrt(fs.Eval(-0.5)-fs.GetParameter(2)));
//   AliDebug(2, Form("<sx>[mum] = %7.3f(prf) %7.3f(diff)", 1.e4*TMath::Sqrt(fs.Eval(-0.5)-fs.GetParameter(2)), 1.e4*TMath::Sqrt(sdiff)));
  if(fCanvas){
    g[1]->Draw("apl");
    fCanvas->Modified(); fCanvas->Update();
    h = g[1]->GetHistogram();
    h->SetTitle(fs.GetTitle());
    h->GetXaxis()->SetTitle("log(q) [a.u.]");h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle("#sigma^{2} (y) [cm^{2}]");h->GetYaxis()->CenterTitle();
    if(IsSaveAs()) fCanvas->SaveAs(Form("D%03d_ResNormTrack_q.gif", fDet));
    else gSystem->Sleep(100);
  }
  return kTRUE;
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessSigma()
{
// As the r-phi coordinate is the only one which is measured by the TRD detector we have to rely on it to
// estimate both the radial (x) and r-phi (y) errors. This method is based on the following assumptions. 
// The measured error in the y direction is the sum of the intrinsic contribution of the r-phi measurement
// with the contribution of the radial measurement - because x is not a parameter of Alice track model (Kalman).
// BEGIN_LATEX
// #sigma^{2}|_{y} = #sigma^{2}_{y*} + #sigma^{2}_{x*}   
// END_LATEX
// In the general case 
// BEGIN_LATEX
// #sigma^{2}_{y*} = #sigma^{2}_{y} + tg^{2}(#alpha_{L})#sigma^{2}_{x_{drift}}   
// #sigma^{2}_{x*} = tg^{2}(#phi - #alpha_{L})*(#sigma^{2}_{x_{drift}} + #sigma^{2}_{x_{0}} + tg^{2}(#alpha_{L})*x^{2}/12)
// END_LATEX
// where we have explicitely show the lorentz angle correction on y and the projection of radial component on the y
// direction through the track angle in the bending plane (phi). Also we have shown that the radial component in the
// last equation has twp terms, the drift and the misalignment (x_0). For ideal geometry or known misalignment one 
// can solve the equation
// BEGIN_LATEX
// #sigma^{2}|_{y} = tg^{2}(#phi - #alpha_{L})*(#sigma^{2}_{x} + tg^{2}(#alpha_{L})*x^{2}/12)+ [#sigma^{2}_{y} + tg^{2}(#alpha_{L})#sigma^{2}_{x}]
// END_LATEX
// by fitting a straight line:
// BEGIN_LATEX
// #sigma^{2}|_{y} = a(x_{cl}, z_{cl}) * tg^{2}(#phi - #alpha_{L}) + b(x_{cl}, z_{cl})
// END_LATEX
// the error parameterization will be given by:
// BEGIN_LATEX
// #sigma_{x} (x_{cl}, z_{cl}) = #sqrt{a(x_{cl}, z_{cl}) - tg^{2}(#alpha_{L})*x^{2}/12}
// #sigma_{y} (x_{cl}, z_{cl}) = #sqrt{b(x_{cl}, z_{cl}) - #sigma^{2}_{x} (x_{cl}, z_{cl}) * tg^{2}(#alpha_{L})}
// END_LATEX
// Below there is an example of such dependency. 
//Begin_Html
//<img src="TRD/clusterSigmaMethod.gif">
//End_Html
//
// The error parameterization obtained by this method are implemented in the functions AliTRDcluster::GetSX() and
// AliTRDcluster::GetSYdrift(). For an independent method to determine s_y as a function of drift length check the 
// function ProcessCenterPad(). One has to keep in mind that while this method return the mean s_y over the distance
// to pad center distribution the other method returns the *STANDARD* value at center=0 (maximum). To recover the 
// standard value one has to solve the obvious equation:
// BEGIN_LATEX
// #sigma_{y}^{STANDARD} = #frac{<#sigma_{y}>}{#int{s exp(s^{2}/#sigma) ds}}
// END_LATEX
// with "<s_y>" being the value calculated here and "sigma" the width of the s_y distribution calculated in 
// ProcessCenterPad().
//  
// Author
// Alexandru Bercuci <A.Bercuci@gsi.de>

  TObjArray *arr = (TObjArray*)fContainer->At(kSigm);
  if(!arr){
    AliWarning("Missing dy=f(x_d, d_w) container");
    return;
  }

  // init visualization
  TGraphErrors *ggs = NULL;
  TGraph *line = NULL;
  if(fCanvas){
    ggs = new TGraphErrors();
    line = new TGraph();
    line->SetLineColor(kRed);line->SetLineWidth(2);
  }

  // init logistic support
  TF1 f("f", "gaus", -.5, .5);
  TLinearFitter gs(1,"pol1");
  TH1 *hFrame=NULL;
  TH1D *h1 = NULL; TH3S *h3=NULL;
  TAxis *ax = NULL;
  Double_t exb2 = fExB*fExB;
  AliTRDcluster c;
  TTree *t = (TTree*)fResults->At(kSigm);
  for(Int_t ix=0; ix<AliTRDseedV1::kNtb; ix++){
    if(!(h3=(TH3S*)arr->At(ix))) continue;
    c.SetPadTime(ix);
    fX = c.GetXloc(fT0, fVdrift);
    fT = c.GetLocalTimeBin(); // ideal
    printf(" pad time[%d] local[%f]\n", ix, fT);
    for(Int_t iz=1; iz<=h3->GetXaxis()->GetNbins(); iz++){
      ax = h3->GetXaxis();
      ax->SetRange(iz, iz);
      fZ = ax->GetBinCenter(iz);

      // reset visualization
      if(fCanvas){ 
        new(ggs) TGraphErrors();
        ggs->SetMarkerStyle(7);
      }
      gs.ClearPoints();

      for(Int_t ip=1; ip<=h3->GetYaxis()->GetNbins(); ip++){
        ax = h3->GetYaxis();
        ax->SetRange(ip, ip); 
        Double_t tgl = ax->GetBinCenter(ip);
        // finish navigation in the HnSparse

        //if(TMath::Abs(dydx)>0.18) continue;
        Double_t tgg = (tgl-fExB)/(1.+tgl*fExB);
        Double_t tgg2 = tgg*tgg;

        h1 = (TH1D*)h3->Project3D("z");
        Int_t entries = (Int_t)h1->Integral();
        if(entries < 50) continue;
        //Adjust(&f, h1);
        h1->Fit(&f, "QN");

        Double_t s2  = f.GetParameter(2)*f.GetParameter(2);
        Double_t s2e = 2.*f.GetParameter(2)*f.GetParError(2);
        // Fill sy^2 = f(tg^2(phi-a_L))
        gs.AddPoint(&tgg2, s2, s2e);

        if(!ggs) continue;
        Int_t jp = ggs->GetN();
        ggs->SetPoint(jp, tgg2, s2);
        ggs->SetPointError(jp, 0., s2e);
      }
      // TODO here a more robust fit method has to be provided
      // for which lower boundaries on the parameters have to 
      // be imposed. Unfortunately the Minuit fit does not work 
      // for the TGraph in the case of B not 0.
      if(gs.Eval()) continue;

      fR[0] = gs.GetParameter(1) - fX*fX*exb2/12.;
      AliDebug(3, Form("  s2x+x2=%f ang=%f s2x=%f", gs.GetParameter(1), fX*fX*exb2/12., fR[0]));
      fR[0] = TMath::Max(fR[0], Float_t(4.e-4)); 

      // s^2_y  = s0^2_y + tg^2(a_L) * s^2_x
      // s0^2_y = f(D_L)*x + s_PRF^2 
      fR[2]= gs.GetParameter(0)-exb2*fR[0];
      AliDebug(3, Form("  s2y+s2x=%f s2y=%f", fR[0], fR[2]));
      fR[2] = TMath::Max(fR[2], Float_t(2.5e-5)); 
      fR[0] = TMath::Sqrt(fR[0]);
      fR[1] = .5*gs.GetParError(1)/fR[0];
      fR[2] = TMath::Sqrt(fR[2]);
      fR[3] = gs.GetParError(0)+exb2*exb2*gs.GetParError(1);
      t->Fill();
      AliDebug(2, Form("xd=%4.2f[cm] sx=%6.1f[um] sy=%5.1f[um]", fX, 1.e4*fR[0], 1.e4*fR[2]));

      if(!fCanvas) continue;
      fCanvas->cd(); fCanvas->SetLogx(); //fCanvas->SetLogy();
      if(!hFrame){ 
        fCanvas->SetMargin(0.15, 0.01, 0.1, 0.01);
        hFrame=new TH1I("hFrame", "", 100, 0., .3);
        hFrame->SetMinimum(0.);hFrame->SetMaximum(.005);
        hFrame->SetXTitle("tg^{2}(#phi-#alpha_{L})");
        hFrame->SetYTitle("#sigma^{2}y[cm^{2}]");
        hFrame->GetYaxis()->SetTitleOffset(2.);
        hFrame->SetLineColor(1);hFrame->SetLineWidth(1);
        hFrame->Draw();
      } else hFrame->Reset();
      Double_t xx = 0., dxx=.2/50;
      for(Int_t ip=0;ip<50;ip++){ 
        line->SetPoint(ip, xx,  gs.GetParameter(0)+xx*gs.GetParameter(1)); 
        xx+=dxx;
      }
      ggs->Draw("pl"); line->Draw("l");
      fCanvas->Modified(); fCanvas->Update();
      if(IsSaveAs()) fCanvas->SaveAs(Form("Figures/ProcessSigma_z[%5.3f]_x[%5.3f].gif", fZ, fX));
      else gSystem->Sleep(100);
    }
  }
  return;
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessMean()
{
// By this method the cluster shift in r-phi and radial directions can be estimated by comparing with the MC.
// The resolution of the cluster corrected for pad tilt with respect to MC in the r-phi (measuring) plane can be 
// expressed by:
// BEGIN_LATEX
// #Delta y=w - y_{MC}(x_{cl})
// w = y_{cl}^{'} + h*(z_{MC}(x_{cl})-z_{cl})
// y_{MC}(x_{cl}) = y_{0} - dy/dx*x_{cl}
// z_{MC}(x_{cl}) = z_{0} - dz/dx*x_{cl}
// y_{cl}^{'} = y_{cl}-x_{cl}*tg(#alpha_{L})
// END_LATEX
// where x_cl is the drift length attached to a cluster, y_cl is the r-phi coordinate of the cluster measured by
// charge sharing on adjacent pads and y_0 and z_0 are MC reference points (as example the track references at 
// entrance/exit of a chamber). If we suppose that both r-phi (y) and radial (x) coordinate of the clusters are 
// affected by errors we can write
// BEGIN_LATEX
// x_{cl} = x_{cl}^{*} + #delta x 
// y_{cl} = y_{cl}^{*} + #delta y 
// END_LATEX 
// where the starred components are the corrected values. Thus by definition the following quantity
// BEGIN_LATEX
// #Delta y^{*}= w^{*} - y_{MC}(x_{cl}^{*})
// END_LATEX
// has 0 average over all dependency. Using this decomposition we can write:
// BEGIN_LATEX
// <#Delta y>=<#Delta y^{*}> + <#delta x * (dy/dx-h*dz/dx) + #delta y - #delta x * tg(#alpha_{L})>
// END_LATEX
// which can be transformed to the following linear dependence:
// BEGIN_LATEX
// <#Delta y>= <#delta x> * (dy/dx-h*dz/dx) + <#delta y - #delta x * tg(#alpha_{L})>
// END_LATEX
// if expressed as function of dy/dx-h*dz/dx. Furtheremore this expression can be plotted for various clusters
// i.e. we can explicitely introduce the diffusion (x_cl) and drift cell - anisochronity (z_cl) dependences. From 
// plotting this dependence and linear fitting it with:
// BEGIN_LATEX
// <#Delta y>= a(x_{cl}, z_{cl}) * (dy/dx-h*dz/dx) + b(x_{cl}, z_{cl})
// END_LATEX
// the systematic shifts will be given by:
// BEGIN_LATEX
// #delta x (x_{cl}, z_{cl}) = a(x_{cl}, z_{cl})
// #delta y (x_{cl}, z_{cl}) = b(x_{cl}, z_{cl}) + a(x_{cl}, z_{cl}) * tg(#alpha_{L})
// END_LATEX
// Below there is an example of such dependency. 
//Begin_Html
//<img src="TRD/clusterShiftMethod.gif">
//End_Html
//
// The occurance of the radial shift is due to the following conditions 
//   - the approximation of a constant drift velocity over the drift length (larger drift velocities close to 
//     cathode wire plane)
//   - the superposition of charge tails in the amplification region (first clusters appear to be located at the 
//     anode wire)
//   - the superposition of charge tails in the drift region (shift towards anode wire)
//   - diffusion effects which convolute with the TRF thus enlarging it
//   - approximate knowledge of the TRF (approximate measuring in test beam conditions) 
// 
// The occurance of the r-phi shift is due to the following conditions 
//   - approximate model for cluster shape (LUT)
//   - rounding-up problems
//
// The numerical results for ideal simulations for the radial and r-phi shifts are displayed below and used 
// for the cluster reconstruction (see the functions AliTRDcluster::GetXcorr() and AliTRDcluster::GetYcorr()). 
//Begin_Html
//<img src="TRD/clusterShiftX.gif">
//<img src="TRD/clusterShiftY.gif">
//End_Html
// More details can be found in the presentation given during the TRD
// software meeting at the end of 2008 and beginning of year 2009, published on indico.cern.ch.
// 
// Author 
// Alexandru Bercuci <A.Bercuci@gsi.de>


 
  TObjArray *arr = (TObjArray*)fContainer->At(kMean);
  if(!arr){
    AliWarning("Missing dy=f(x_d, d_w) container");
    return;
  }

  // init logistic support
  TF1 f("f", "gaus", -.5, .5);
  TF1 line("l", "[0]+[1]*x", -.15, .15);
  TGraphErrors *gm = new TGraphErrors();
  TH1 *hFrame=NULL;
  TH1D *h1 = NULL; TH3S *h3 =NULL;
  TAxis *ax = NULL;

  AliDebug(1, Form("Calibrate for Det[%3d] t0[%5.3f] vd[%5.3f]", fDet, fT0, fVdrift));

  AliTRDcluster c;
  TTree *t = (TTree*)fResults->At(kMean);
  for(Int_t ix=0; ix<AliTRDseedV1::kNtb; ix++){
    if(!(h3=(TH3S*)arr->At(ix))) continue;
    c.SetPadTime(ix);
    fX = c.GetXloc(fT0, fVdrift);
    fT = c.GetLocalTimeBin();
    for(Int_t iz=1; iz<=h3->GetXaxis()->GetNbins(); iz++){
      ax = h3->GetXaxis();
      ax->SetRange(iz, iz);
      fZ = ax->GetBinCenter(iz);

      // reset fitter
      new(gm) TGraphErrors();
      gm->SetMarkerStyle(7);

      for(Int_t ip=1; ip<=h3->GetYaxis()->GetNbins(); ip++){
        ax = h3->GetYaxis();
        ax->SetRange(ip, ip); 
        Double_t tgl = ax->GetBinCenter(ip);
        // finish navigation in the HnSparse

        h1 = (TH1D*)h3->Project3D("z");
        Int_t entries = (Int_t)h1->Integral();
        if(entries < 50) continue;
        //Adjust(&f, h1);
        h1->Fit(&f, "QN");

        // Fill <Dy> = f(dydx - h*dzdx)
        Int_t jp = gm->GetN();
        gm->SetPoint(jp, tgl, f.GetParameter(1));
        gm->SetPointError(jp, 0., f.GetParError(1));
      }
      if(gm->GetN()<10) continue;

      gm->Fit(&line, "QN");
      fR[0] = line.GetParameter(1); // dx
      fR[1] = line.GetParError(1);
      fR[2] = line.GetParameter(0) + fExB*fR[0]; // xs = dy - tg(a_L)*dx
      t->Fill();
      AliDebug(2, Form("tb[%02d] xd=%4.2f[cm] dx=%6.2f[um] dy=%6.2f[um]", ix, fX, 1.e4*fR[0], 1.e4*fR[2]));
      if(!fCanvas) continue;

      fCanvas->cd();
      if(!hFrame){ 
        fCanvas->SetMargin(0.1, 0.02, 0.1, 0.01);
        hFrame=new TH1I("hFrame", "", 100, -.3, .3);
        hFrame->SetMinimum(-.1);hFrame->SetMaximum(.1);
        hFrame->SetXTitle("tg#phi-htg#theta");
        hFrame->SetYTitle("#Delta y[cm]");
        hFrame->GetYaxis()->SetTitleOffset(1.5);
        hFrame->SetLineColor(1);hFrame->SetLineWidth(1);
        hFrame->Draw();
      } else hFrame->Reset();
      gm->Draw("pl"); line.Draw("same");
      fCanvas->Modified(); fCanvas->Update();
      if(IsSaveAs()) fCanvas->SaveAs(Form("Figures/ProcessMean_Z[%5.3f]_TB[%02d].gif", fZ, ix));
      else gSystem->Sleep(100);
    }
  }
}
