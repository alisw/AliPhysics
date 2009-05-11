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
#include "info/AliTRDclusterInfo.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"

#include "AliLog.h"
#include "AliTracker.h"
#include "AliCDBManager.h"

#include "TROOT.h"
#include "TObjArray.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TH2I.h"
#include "TH3S.h"
#include "TTree.h"
#include "TMath.h"
#include "TLinearFitter.h"

#include "TCanvas.h"
#include "TSystem.h"

ClassImp(AliTRDclusterResolution)

const Float_t AliTRDclusterResolution::fgkTimeBinLength = 1./ AliTRDCommonParam::Instance()->GetSamplingFrequency();
//_______________________________________________________
AliTRDclusterResolution::AliTRDclusterResolution(const char *name, const char *title)
  : AliTRDrecoTask(name, title)
  ,fCanvas(0x0)
  ,fInfo(0x0)
  ,fResults(0x0)
  ,fAt(0x0)
  ,fStatus(0)
  ,fDet(-1)
  ,fExB(0.)
  ,fVdrift(0.)
  ////////////
  ,fLy(0)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
{
  memset(fR, 0, 4*sizeof(Float_t));
  // time drift axis
  fAt = new TAxis(kNTB, 0., kNTB*fgkTimeBinLength);

  // By default register all analysis
  // The user can switch them off in his steering macro
  SetProcess(kQRes);
  SetProcess(kCenter);
  SetProcess(kMean);
  SetProcess(kSigm);
}

//_______________________________________________________
AliTRDclusterResolution::~AliTRDclusterResolution()
{
  if(fCanvas) delete fCanvas;
  if(fAt) delete fAt;
  if(fResults){
    fResults->Delete();
    delete fResults;
  }
}

//_______________________________________________________
void AliTRDclusterResolution::ConnectInputData(Option_t *)
{
  fInfo = dynamic_cast<TObjArray *>(GetInputData(0));
}

//_______________________________________________________
void AliTRDclusterResolution::CreateOutputObjects()
{
  OpenFile(0, "RECREATE");
  fContainer = Histos();
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::GetRefFigure(Int_t ifig)
{
  if(!fResults) return kFALSE;
  
  TList *l = 0x0;
  TObjArray *arr = 0x0;
  TH2 *h2 = 0x0;TH1 *h1 = 0x0;
  TGraphErrors *gm(0x0), *gs(0x0), *gp(0x0);
  switch(ifig){
  case kQRes:
    if(!(arr = (TObjArray*)fResults->At(kQRes))) break;
    if(!(gm = (TGraphErrors*)arr->At(0))) break;
    if(!(gs = (TGraphErrors*)arr->At(1))) break;
    if(!(gp = (TGraphErrors*)arr->At(2))) break;
    gs->Draw("apl");
    gs->GetHistogram()->GetYaxis()->SetRangeUser(-.01, .6);
    gs->GetHistogram()->SetXTitle("Q [a.u.]");
    gs->GetHistogram()->SetYTitle("#sigma_{y} / #mu_{y} [mm] / freq");
    gm->Draw("pl");
    gp->Draw("pl");
    return kTRUE;
  case kCenter:
    if(!(arr = (TObjArray*)fResults->At(kCenter))) break;
    gPad->Divide(3, 1); l = gPad->GetListOfPrimitives();
    for(Int_t ipad = 3; ipad--;){
      if(!(h2 = (TH2F*)arr->At(ipad))) return kFALSE;
      ((TVirtualPad*)l->At(ipad))->cd();
      h2->Draw("lego2fb");
    }
    return kTRUE;
  case kSigm:
    if(!(arr = (TObjArray*)fResults->At(kSigm))) break;
    gPad->Divide(2, 1); l = gPad->GetListOfPrimitives();
    if(!(h2 = (TH2F*)arr->At(0))) return kFALSE;
    ((TVirtualPad*)l->At(0))->cd();
    h1 = h2->ProjectionY("hsx_pyy"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#sigma_{x}> [#mum]");
    h1->GetXaxis()->SetRange(2, kNTB-1); h1->Draw("pc");

    if(!(h2 = (TH2F*)arr->At(1))) return kFALSE;
    ((TVirtualPad*)l->At(1))->cd();
    h1 = h2->ProjectionY("hsy_pyy"); h1->Scale(1.e4/kND); h1->SetMarkerStyle(24);
    h1->SetYTitle("<#sigma_{y}> [#mum]");
    h1->GetXaxis()->SetRange(2, kNTB-1); h1->Draw("pc");
    return kTRUE;
  case kMean:
    if(!(arr = (TObjArray*)fResults->At(kMean))) break;
    gPad->Divide(2, 1);  l = gPad->GetListOfPrimitives();
    ((TVirtualPad*)l->At(0))->cd();
    if(!(gm = (TGraphErrors*)arr->At(0))) return kFALSE;
    gm->Draw("apl");
    gm->GetHistogram()->SetXTitle("t_{drift} [#mus]");
    gm->GetHistogram()->SetYTitle("dx [#mum]");

    ((TVirtualPad*)l->At(1))->cd();
    if(!(gm = (TGraphErrors*)arr->At(1))) return kFALSE;
    gm->Draw("apl");
    gm->GetHistogram()->SetXTitle("t_{drift} [#mus]");
    gm->GetHistogram()->SetYTitle("dy [#mum]");

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
  if(fContainer) return fContainer;
  fContainer = new TObjArray(kNtasks);
  //fContainer->SetOwner(kTRUE);

  TH2I *h2 = 0x0;
  TH3S *h3 = 0x0;
  TObjArray *arr = 0x0;

  fContainer->AddAt(h2 = new TH2I("Charge", "dy=f(q)", 50, 2.2, 7.5, 60, -.3, .3), kQRes);
  h2->SetXTitle("log(q) [a.u.]");
  h2->SetYTitle("#Delta y[cm]");
  h2->SetZTitle("entries");

  fContainer->AddAt(arr = new TObjArray(AliTRDgeometry::kNlayer), kCenter);
  arr->SetName("Center");
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(h3=(TH3S*)gROOT->FindObject(Form("hc_l%1d", il)))){ 
      h3 = new TH3S(
        Form("hc_l%1d", il), 
        Form(" ly [%d]", il), 
        kNTB, fAt->GetBinLowEdge(1), fAt->GetBinUpEdge(kNTB),   // x
        51, -.51, .51, // y 
        60, -.3, .3); // dy
      h3->SetXTitle("x [#mus]");
      h3->SetYTitle("y [pw]");
      h3->SetZTitle("#Delta y[cm]");
    } h3->Reset();
    arr->AddAt(h3, il);
  }

  fContainer->AddAt(arr = new TObjArray(kNTB), kSigm);
  arr->SetName("Resolution");
  for(Int_t ix=0; ix<kNTB; ix++){
    if(!(h3=(TH3S*)gROOT->FindObject(Form("hr_x%02d", ix)))){
      h3 = new TH3S(
        Form("hr_x%02d", ix), 
        Form(" t_{drift}(%3.1f-%3.1f)[#mus]", fAt->GetBinLowEdge(ix+1), fAt->GetBinUpEdge(ix+1)), 
        kND, 0., 2.5,   // z 
        35, -.35, .35, // tgp 
        60, -.3, .3); // dy
      h3->SetXTitle("z [mm]");
      h3->SetYTitle("tg#phi");
      h3->SetZTitle("#Delta y[cm]");
    }
    arr->AddAt(h3, ix);
  }

  fContainer->AddAt(arr = new TObjArray(kNTB), kMean);
  arr->SetName("Systematics");
  for(Int_t ix=0; ix<kNTB; ix++){
    if(!(h3=(TH3S*)gROOT->FindObject(Form("hs_x%02d", ix)))){
      h3 = new TH3S(
        Form("hs_x%02d", ix), 
        Form(" t_{drift}(%3.1f-%3.1f)[#mus]", fAt->GetBinLowEdge(ix+1), fAt->GetBinUpEdge(ix+1)), 
        kND, 0., 2.5,   // z 
        35, -.35, .35, // tgp-h tgt 
        60, -.3, .3); // dy
      h3->SetXTitle("z [mm]");
      h3->SetYTitle("tg(#phi) - h*tg(#theta)");
      h3->SetZTitle("#Delta y[cm]");
    }
    arr->AddAt(h3, ix);
  }

  return fContainer;
}

//_______________________________________________________
void AliTRDclusterResolution::Exec(Option_t *)
{
  if(!HasExB()) AliWarning("ExB was not set. Call SetExB() before running the task.");

  Int_t det, t;
  Float_t x, y, z, q, dy, dydx, dzdx, cov[3], covcl[3];
  TH2I *h2 = 0x0;
  TH3S *h3 = 0x0;

  // define limits around ExB for which x contribution is negligible
  const Float_t kDtgPhi = 3.5e-2; //(+- 2 deg)

  TObjArray *arr0 = (TObjArray*)fContainer->At(kCenter);
  TObjArray *arr1 = (TObjArray*)fContainer->At(kSigm);
  TObjArray *arr2 = (TObjArray*)fContainer->At(kMean);

  const AliTRDclusterInfo *cli = 0x0;
  TIterator *iter=fInfo->MakeIterator();
  while((cli=dynamic_cast<AliTRDclusterInfo*>((*iter)()))){
    cli->GetCluster(det, x, y, z, q, t, covcl);
    if(fDet>=0 && fDet!=det) continue;
    
    dy = cli->GetResolution();
    cli->GetGlobalPosition(y, z, dydx, dzdx, &cov[0]);

    // resolution as a function of cluster charge
    // only for phi equal exB 
    if(TMath::Abs(dydx-fExB) < kDtgPhi){
      h2 = (TH2I*)fContainer->At(kQRes);
      h2->Fill(TMath::Log(q), dy);
    }

    // do not use problematic clusters in resolution analysis
    // TODO define limits as calibration aware (gain) !!
    if(q<20. || q>250.) continue;

    x = (t+.5)*fgkTimeBinLength; // conservative approach !!

    // resolution as a function of y displacement from pad center
    // only for phi equal exB
    if(TMath::Abs(dydx-fExB) < kDtgPhi){
      h3 = (TH3S*)arr0->At(AliTRDgeometry::GetLayer(det));
      h3->Fill(x, cli->GetYDisplacement(), dy);
    }

    Int_t ix = fAt->FindBin(x);
    if(ix==0 || ix == fAt->GetNbins()+1){
      AliWarning(Form("Drift time %3.1f outside allowed range", x));
      continue;
    }

    // fill histo for resolution (sigma)
    ((TH3S*)arr1->At(ix-1))->Fill(10.*cli->GetAnisochronity(), dydx, dy);

    // fill histo for systematic (mean)
    ((TH3S*)arr2->At(ix-1))->Fill(10.*cli->GetAnisochronity(), dydx-cli->GetTilt()*dzdx, dy);  
  }
  PostData(0, fContainer);
}


//_______________________________________________________
Bool_t AliTRDclusterResolution::PostProcess()
{
  if(!fContainer) return kFALSE;
  if(!HasExB()) AliWarning("ExB was not set. Call SetExB() before running the post processing.");
  
  TObjArray *arr = 0x0;
  TTree *t=0x0;
  if(!fResults){
    TGraphErrors *g = 0x0;
    fResults = new TObjArray(kNtasks);
    fResults->SetOwner();
    fResults->AddAt(arr = new TObjArray(3), kQRes);
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
    fResults->AddAt(t = new TTree("cent", "dy=f(y,x,ly)"), kCenter);
    t->Branch("ly", &fLy, "ly/B");
    t->Branch("x", &fX, "x/F");
    t->Branch("y", &fY, "y/F");
    t->Branch("m", &fR[0], "m[2]/F");
    t->Branch("s", &fR[2], "s[2]/F");


    fResults->AddAt(t = new TTree("sigm", "dy=f(dw,x,dydx)"), kSigm);
    t->Branch("x", &fX, "x/F");
    t->Branch("z", &fZ, "z/F");
    t->Branch("sx", &fR[0], "sx[2]/F");
    t->Branch("sy", &fR[2], "sy[2]/F");


    fResults->AddAt(t = new TTree("mean", "dy=f(dw,x,dydx - h dzdx)"), kMean);
    t->Branch("x", &fX, "x/F");
    t->Branch("z", &fZ, "z/F");
    t->Branch("dx", &fR[0], "dx[2]/F");
    t->Branch("dy", &fR[2], "dy[2]/F");
  } else {
    TObject *o = 0x0;
    TIterator *iter=fResults->MakeIterator();
    while((o=((*iter)()))) o->Clear(); // maybe it is wrong but we should never reach this point
  }
  
  // precalculated value of tg^2(alpha_L)
  Double_t exb2 = fExB*fExB;
  // square of the mean value of sigma drift length.
  // has to come from previous calibration 
  //Double_t sxd2 = 1.;// [mm^2]

  printf("ExB[%e] ExB2[%e]\n", fExB, exb2);

  // process resolution dependency on charge
  if(HasProcess(kQRes)) ProcessCharge();
  
  // process resolution dependency on y displacement
  if(HasProcess(kCenter)) ProcessCenterPad();

  // process resolution dependency on drift legth and drift cell width
  if(HasProcess(kSigm)) ProcessSigma();

  // process systematic shift on drift legth and drift cell width
  if(HasProcess(kMean)) ProcessMean();

  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDclusterResolution::SetExB(Int_t det, Int_t col, Int_t row)
{
  // check OCDB
  AliCDBManager *cdb = AliCDBManager::Instance();
  if(cdb->GetRun() < 0){
    AliError("OCDB manager not properly initialized");
    return kFALSE;
  }

  // check magnetic field
  if(TMath::Abs(AliTracker::GetBz()) < 1.e-10){
    AliWarning("B=0. Magnetic field may not be initialized. Continue if you know what you are doing !");
  }

  // set reference detector if any
  if(det>=0 && det<AliTRDgeometry::kNdet) fDet = det;
  else det = 0;

  AliTRDcalibDB *fCalibration  = AliTRDcalibDB::Instance();
  AliTRDCalROC  *fCalVdriftROC = fCalibration->GetVdriftROC(det);
  const AliTRDCalDet  *fCalVdriftDet = fCalibration->GetVdriftDet();

  fVdrift = fCalVdriftDet->GetValue(det) * fCalVdriftROC->GetValue(col, row);
  fExB   = AliTRDCommonParam::Instance()->GetOmegaTau(fVdrift);
  SetBit(kExB);
  return kTRUE;
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
  TH2I *h2 = 0x0;
  if((h2 = (TH2I*)fContainer->At(kQRes))) {
    AliWarning("Missing dy=f(Q) histo");
    return;
  }
  TF1 f("f", "gaus", -.5, .5);
  TAxis *ax = 0x0;
  TH1D *h1 = 0x0;

  TObjArray *arr = (TObjArray*)fResults->At(kQRes);
  TGraphErrors *gqm = (TGraphErrors*)arr->At(0);
  TGraphErrors *gqs = (TGraphErrors*)arr->At(1);
  TGraphErrors *gqp = (TGraphErrors*)arr->At(2);
  Double_t q, n = 0., entries;
  ax = h2->GetXaxis();
  for(Int_t ix=1; ix<=ax->GetNbins(); ix++){
    q = TMath::Exp(ax->GetBinCenter(ix));
    if(q<20. || q>250.) continue; // ?!

    h1 = h2->ProjectionY("py", ix, ix);
    entries = h1->GetEntries();
    if(entries < 50) continue;
    Adjust(&f, h1);
    h1->Fit(&f, "Q");

    // Fill sy^2 = f(q)
    Int_t ip = gqm->GetN();
    gqm->SetPoint(ip, q, 10.*f.GetParameter(1));
    gqm->SetPointError(ip, 0., 10.*f.GetParError(1));

    // correct sigma for ExB effect
    gqs->SetPoint(ip, q, 1.e1*f.GetParameter(2)/**f.GetParameter(2)-exb2*sxd2*/);
    gqs->SetPointError(ip, 0., 1.e1*f.GetParError(2)/**f.GetParameter(2)*/);

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
    gqp->SetPoint(ip, q, entries);
    gqs->GetPoint(ip, q, sy);
    sm += entries*sy;
  }

  // error parametrization s(q) = <sy> + b(1/q-1/q0)
  TF1 fq("fq", "[0] + [1]/x", 20., 250.);
  gqs->Fit(&fq);
  printf("sy(Q) :: sm[%f] b[%f] 1/q0[%f]\n", sm, fq.GetParameter(1), (sm-fq.GetParameter(0))/fq.GetParameter(1));
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessCenterPad()
{
  TObjArray *arr = (TObjArray*)fContainer->At(kCenter);
  if(!arr) {
    AliWarning("Missing dy=f(y | x, ly) container");
    return;
  }
  TF1 f("f", "gaus", -.5, .5);
  TH1D *h1 = 0x0; TH3S *h3=0x0;
  TTree *t = (TTree*)fResults->At(kCenter);
  TAxis *ax = 0x0;
  for(Int_t il=0; il<arr->GetEntriesFast(); il++){
    if(!(h3 = (TH3S*)arr->At(il))) continue;

    fLy = il;
    for(Int_t ix=1; ix<=h3->GetXaxis()->GetNbins(); ix++){
      ax = h3->GetXaxis();
      ax->SetRange(ix, ix);     
      fX  = ax->GetBinCenter(ix);
      //printf("  x[%2d]=%4.2f\n", ix, fX);
      for(Int_t iy=1; iy<=h3->GetYaxis()->GetNbins(); iy++){
        ax = h3->GetYaxis();
        ax->SetRange(iy, iy);
        fY  = ax->GetBinCenter(iy);
        //printf("    y[%2d]=%5.2f\n", iy, fY);
        // finish navigation in the HnSparse

        h1 = (TH1D*)h3->Project3D("z");
        Int_t entries = (Int_t)h1->Integral();
        if(entries < 50) continue;
        //Adjust(&f, h1);
        h1->Fit(&f, "QN");

    
        // Fill sy,my=f(y_w,x,ly)
        fR[0] = f.GetParameter(1); fR[1] = f.GetParError(1);
        fR[2] = f.GetParameter(2); fR[3] = f.GetParError(2);

        //printf("ly[%d] x[%3.1f] y[%+5.2f] m[%5.3f] s[%5.3f] \n", fLy, fX, fY, fR[0], fR[2]);
        t->Fill();
      }
    }
    if(!fCanvas) continue;

    t->Draw("y:x>>h(23, 0.1, 2.4, 51, -.51, .51)",
            Form("m[0]*(ly==%d&&abs(m[0])<1.e-1)", fLy),
            "lego2fb");
    fCanvas->Modified(); fCanvas->Update();
    if(IsSaveAs()) fCanvas->SaveAs(Form("Figures/ProcessCenter_ly[%d].gif", fLy));
    else gSystem->Sleep(100);
  }
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessSigma()
{
  TObjArray *arr = (TObjArray*)fContainer->At(kSigm);
  if(!arr){
    AliWarning("Missing dy=f(x_d, d_w) container");
    return;
  }

  // init visualization
  TGraphErrors *ggs = 0x0;
  TGraph *line = 0x0;
  if(fCanvas){
    ggs = new TGraphErrors();
    line = new TGraph();
    line->SetLineColor(kRed);line->SetLineWidth(2);
  }

  // init logistic support
  TF1 f("f", "gaus", -.5, .5);
  TLinearFitter gs(1,"pol1");
  TH1 *hFrame=0x0;
  TH1D *h1 = 0x0; TH3S *h3=0x0;
  TAxis *ax = 0x0;
  Double_t exb2 = fExB*fExB;

  TTree *t = (TTree*)fResults->At(kSigm);
  for(Int_t ix=0; ix<kNTB; ix++){
    if(!(h3=(TH3S*)arr->At(ix))) continue;
    fX = fAt->GetBinCenter(ix+1);
    
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
        Int_t ip = ggs->GetN();
        ggs->SetPoint(ip, tgg2, s2);
        ggs->SetPointError(ip, 0., s2e);
      }
      if(gs.Eval()) continue;

      // s^2_x = s0^2_x - x^2*tg^2(a_L)/12
      fR[0] = gs.GetParameter(1)/* - x*x*exb2/12.*/;
      if(fR[0]<0.) continue; 
      fR[0] = TMath::Sqrt(fR[0]);
      fR[1] = .5*gs.GetParError(1)/fR[0];

      // s^2_y  = s0^2_y + tg^2(a_L) * s^2_x
      // s0^2_y = f(D_L)*x + s_PRF^2 
      fR[2]= gs.GetParameter(0)/*-exb2*sx*/;
      if(fR[1] <0.) continue;
      fR[2] = TMath::Sqrt(fR[2]);
      fR[3] = gs.GetParError(0)+exb2*exb2*gs.GetParError(1);
      t->Fill();

      if(!fCanvas) continue;
      fCanvas->cd(); fCanvas->SetLogx(); //fCanvas->SetLogy();
      if(!hFrame){ 
        hFrame=new TH1I("hFrame", "", 100, 0., .3);
        hFrame->SetMinimum(0.);hFrame->SetMaximum(.005);
        hFrame->SetXTitle("tg^{2}(#phi-#alpha_{L})");
        hFrame->SetYTitle("#sigma^{2}y[cm^{2}]");
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

      printf("    xd=%4.1f[cm] sx=%5.3e[cm] sy=%5.3e[cm]\n", fX, TMath::Sqrt(fR[0]), TMath::Sqrt(fR[1]));
    }
  }
  return;
}

//_______________________________________________________
void AliTRDclusterResolution::ProcessMean()
{
  TObjArray *arr = (TObjArray*)fContainer->At(kMean);
  if(!arr){
    AliWarning("Missing dy=f(x_d, d_w) container");
    return;
  }

  // init logistic support
  TF1 f("f", "gaus", -.5, .5);
  TF1 line("l", "[0]+[1]*x", -.15, .15);
  TGraphErrors *gm = new TGraphErrors();
  TH1 *hFrame=0x0;
  TH1D *h1 = 0x0; TH3S *h3 =0x0;
  TAxis *ax = 0x0;

  TTree *t = (TTree*)fResults->At(kMean);
  for(Int_t ix=0; ix<kNTB; ix++){
    if(!(h3=(TH3S*)arr->At(ix))) continue;
    fX = fAt->GetBinCenter(ix+1);
  
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
        Int_t ip = gm->GetN();
        gm->SetPoint(ip, tgl, f.GetParameter(1));
        gm->SetPointError(ip, 0., f.GetParError(1));
        ip++;
      }
      if(gm->GetN()<4) continue;

      gm->Fit(&line, "QN");
      fR[0] = line.GetParameter(1); // dx
      fR[1] = line.GetParError(1);
      fR[2] = line.GetParameter(0) + fExB*fR[0]; // xs = dy - tg(a_L)*dx
      t->Fill();

      if(!fCanvas) continue;
      fCanvas->cd();
      if(!hFrame){ 
        hFrame=new TH1I("hFrame", "", 100, -.3, .3);
        hFrame->SetMinimum(-.1);hFrame->SetMaximum(.1);
        hFrame->SetXTitle("tg#phi-htg#theta");
        hFrame->SetYTitle("#Deltay[cm]");
        hFrame->SetLineColor(1);hFrame->SetLineWidth(1);
        hFrame->Draw();
      } else hFrame->Reset();
      gm->Draw("pl"); line.Draw("same");
      fCanvas->Modified(); fCanvas->Update();
      if(IsSaveAs()) fCanvas->SaveAs(Form("Figures/ProcessMean_Z[%5.3f]_X[%5.3f].gif", fZ, fX));
      else gSystem->Sleep(100);
      printf("    xd=%4.2f[cm] dx=%5.3e[cm] dy=%5.3e[cm]\n", fX, fR[0], fR[2]);
    }
  }
  
  // draw shift results
  //t->Draw("z:x>>h(24, 0, 2.4, 25, 0, 2.5)", "dx*(abs(dx)<1.e-2)", "lego2fb");
  //t->Draw("z:x>>h(24, 0, 2.4, 25, 0, 2.5)", "dy*(abs(dx)<1.e-2)", "lego2fb");

}
