#include "AliTRDclusterResolution.h"
#include "AliTRDtrackInfo/AliTRDclusterInfo.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

#include "AliLog.h"

#include "TObjArray.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2I.h"
#include "TMath.h"

ClassImp(AliTRDclusterResolution)


//_______________________________________________________
AliTRDclusterResolution::AliTRDclusterResolution()
  : AliTRDrecoTask("ClErrParam", "Cluster Error Parametrization")
  ,fInfo(0x0)
  ,fResults(0x0)
  ,fAt(0x0)
  ,fAd(0x0)
  ,fExB(0.)
{
  fAt = new TAxis(kNTB, -0.075, (kNTB-.5)*.15);
  fAd = new TAxis(kND, 0., .25);
}

//_______________________________________________________
AliTRDclusterResolution::~AliTRDclusterResolution()
{
  if(fAd) delete fAd;
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
  
  TObjArray *arr = 0x0;
  TH2 *h2 = 0x0;
  TGraphErrors *gm(0x0), *gs(0x0);
  switch(ifig){
  case kQRes:
    if(!(arr = (TObjArray*)fResults->At(kQRes))) break;
    if(!(gm = (TGraphErrors*)arr->At(0))) break;
    if(!(gs = (TGraphErrors*)arr->At(1))) break;
    gs->Draw("apl");
    gs->GetHistogram()->SetXTitle("Log(Q) [a.u.]");
    gs->GetHistogram()->SetYTitle("#sigma_{y}^{2} [mm^{2}] [a.u.]");
    return kTRUE;
  case kYRes:
    if(!(arr = (TObjArray*)fResults->At(kYRes))) break;
    if(!(gm = (TGraphErrors*)arr->At(0))) break;
    if(!(gs = (TGraphErrors*)arr->At(1))) break;
    gs->Draw("apl");
    gs->GetHistogram()->SetXTitle("y [w]");
    gs->GetHistogram()->SetYTitle("#sigma_{y}^{2} [mm^{2}] [a.u.]");
    return kTRUE;
  case kSXRes:
    if(!(h2 = (TH2F*)fResults->At(kSXRes))) break;
    h2->Draw("lego2");
    return kTRUE;
  case kSYRes:
    if(!(h2 = (TH2F*)fResults->At(kSYRes))) break;
    h2->Draw("lego2");
    return kTRUE;
  default:
    break;
  }
  
  return kFALSE;
}

//_______________________________________________________
TObjArray* AliTRDclusterResolution::Histos()
{
  if(fContainer) return fContainer;
  fContainer = new TObjArray(3);
  //fContainer->SetOwner(kTRUE);

  TH2I *h2 = 0x0;
  TObjArray *arr = 0x0;

  fContainer->AddAt(h2 = new TH2I("h_q", "", 50, 2.2, 7.5, 100, -.5, .5), kQRes);
  h2->SetXTitle("log(q) [a.u.]");
  h2->SetYTitle("#Delta y[cm]");
  h2->SetZTitle("entries");

  Double_t w = 0.;
  AliTRDgeometry geo;
  fContainer->AddAt(arr = new TObjArray(AliTRDgeometry::kNlayer), kYRes);
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    w = .5*geo.GetPadPlane(ily, 2)->GetWidthIPad();
    arr->AddAt(h2 = new TH2I(Form("h_y%d", ily), Form("Ly[%d]", ily), 50, -w, w, 100, -.5, .5), ily);
    h2->SetXTitle("y_{w} [w]");
    h2->SetYTitle("#Delta y[cm]");
    h2->SetZTitle("entries");
  }

  fContainer->AddAt(arr = new TObjArray(kN), kSXRes);
  Int_t ih = 0;
  for(Int_t id=1; id<=fAd->GetNbins(); id++){
    for(Int_t it=1; it<=fAt->GetNbins(); it++){
      arr->AddAt(h2 = new TH2I(Form("h_d%02dt%02d", id, it), Form("d_{wire}(%3.1f-%3.1f)[mm] x_{drift}(%5.2f-%5.2f)[mm]", 10.*fAd->GetBinCenter(id)- 5.*fAd->GetBinWidth(id), 10.*fAd->GetBinCenter(id)+ 5.*fAd->GetBinWidth(id), 10.*fAt->GetBinCenter(it)- 5.*fAt->GetBinWidth(it), 10.*fAt->GetBinCenter(it)+ 5.*fAt->GetBinWidth(it)), 30, -.15, .15, 100, -.5, .5), ih++);
      h2->SetXTitle("tg#phi");
      h2->SetYTitle("#Delta y[cm]");
      h2->SetZTitle("entries");
    }
  }
  return fContainer;
}

//_______________________________________________________
void AliTRDclusterResolution::Exec(Option_t *)
{
  Int_t det, t;
  Float_t x, y, z, q, dy, dydx, dzdx, cov[3], covcl[3];
  TH2I *h2 = 0x0;

  TObjArray *arr0 = (TObjArray*)fContainer->At(kYRes);
  TObjArray *arr1 = (TObjArray*)fContainer->At(kSXRes);

  const AliTRDclusterInfo *cli = 0x0;
  TIterator *iter=fInfo->MakeIterator();
  while((cli=dynamic_cast<AliTRDclusterInfo*>((*iter)()))){
    dy = cli->GetResolution();
    Int_t it = fAt->FindBin(cli->GetDriftLength());
    if(it==0 || it == fAt->GetNbins()+1){
      AliWarning(Form("Drift length %f outside allowed range", cli->GetDriftLength()));
      continue;
    }
    Int_t id = fAd->FindBin(cli->GetAnisochronity());
    if(id==0 || id == fAd->GetNbins()+1){
      AliWarning(Form("Distance to anode %f outside allowed range", cli->GetAnisochronity()));
      continue;
    }
    if(!(h2 = (TH2I*)arr1->At((id-1)*kNTB+it-1))){
      AliWarning(Form("Missing histo at index idx[%3d] [id[%2d] it[%2d]] xd[%f] d[%f]\n", (id-1)*kNTB+it-1, id, it, cli->GetDriftLength(), cli->GetAnisochronity()));
      continue;
    }

    cli->GetGlobalPosition(y, z, dydx, dzdx, &cov[0]);
    h2->Fill(dydx, cov[0]!=0. ? dy/TMath::Sqrt(cov[0]) : dy);

    // resolution as a function of:
    //   - cluster charge and
    //   - y displacement
    // only for phi equal exB 
    if(TMath::Abs(dydx-fExB)<.01){
      cli->GetCluster(det, x, y, z, q, t, covcl);
      h2 = (TH2I*)fContainer->At(kQRes);
      h2->Fill(TMath::Log(q), dy);
      
      h2 = (TH2I*)arr0->At(AliTRDgeometry::GetLayer(det));
      h2->Fill(cli->GetYDisplacement(), dy);
    }
  }
  PostData(0, fContainer);
}


//_______________________________________________________
Bool_t AliTRDclusterResolution::PostProcess()
{
  if(!fContainer) return kFALSE;
  
  TObjArray *arr = 0x0;
  TH2 *h2 = 0x0; 
  if(!fResults){
    TGraphErrors *g = 0x0;
    fResults = new TObjArray(4);
    fResults->SetOwner();
    fResults->AddAt(arr = new TObjArray(2), kQRes);
    arr->SetOwner();
    arr->AddAt(g = new TGraphErrors(), 0);
    g->SetLineColor(kBlue); g->SetMarkerColor(kBlue);
    g->SetMarkerStyle(7); 
    arr->AddAt(g = new TGraphErrors(), 1);
    g->SetLineColor(kRed); g->SetMarkerColor(kRed);
    g->SetMarkerStyle(23); 

    fResults->AddAt(arr = new TObjArray(2), kYRes);
    arr->SetOwner();
    arr->AddAt(g = new TGraphErrors(), 0);
    g->SetLineColor(kBlue); g->SetMarkerColor(kBlue);
    g->SetMarkerStyle(7); 
    arr->AddAt(g = new TGraphErrors(), 1);
    g->SetLineColor(kRed); g->SetMarkerColor(kRed);
    g->SetMarkerStyle(23); 

    fResults->AddAt(h2 = new TH2F("hSX", "", 
      fAd->GetNbins(), fAd->GetXmin(), fAd->GetXmax(), 
      fAt->GetNbins(), fAt->GetXmin(), fAt->GetXmax()), kSXRes);
    h2->SetXTitle("d [mm]");
    h2->SetYTitle("x [mm]");
    h2->SetZTitle("#sigma_{x}^{2} [mm^{2}]");
    fResults->AddAt(h2 = (TH2F*)h2->Clone("hSY"), kSYRes);
    h2->SetZTitle("#sigma_{y}^{2} [mm^{2}]");
  } else {
    TObject *o = 0x0;
    TIterator *iter=fResults->MakeIterator();
    while((o=((*iter)()))) o->Clear(); // maybe it is wrong but we should never reach this point
  }

  TF1 f("f", "gaus", -.5, .5);
  TF1 sig("sig", "pol1", 0., .0225);

  TAxis *ax = 0x0;
  TH1D *h1 = 0x0;

  // process resolution dependency on charge
  if((h2 = (TH2I*)fContainer->At(kQRes))) {
    TObjArray *arr = (TObjArray*)fResults->At(kQRes);
    TGraphErrors *gqm = (TGraphErrors*)arr->At(0);
    TGraphErrors *gqs = (TGraphErrors*)arr->At(1);
    ax = h2->GetXaxis();
    for(Int_t ix=1; ix<=ax->GetNbins(); ix++){
      Float_t logq = ax->GetBinCenter(ix);
      h1 = h2->ProjectionY("py", ix, ix);
      if(h1->GetEntries() < 50) continue;
      Adjust(&f, h1);
      h1->Fit(&f, "Q");
  
      // Fill sy^2 = f(log(q))
      Int_t ip = gqm->GetN();
      gqm->SetPoint(ip, logq, 10.*f.GetParameter(1));
      gqm->SetPointError(ip, 0., 10.*f.GetParError(1));
      gqs->SetPoint(ip, logq, 1.e2*f.GetParameter(2)*f.GetParameter(2));
      gqs->SetPointError(ip, 0., 2.e2*f.GetParameter(2)*f.GetParError(2));
    } 
  } else AliWarning("Missing dy=f(Q) histo");

  // process resolution dependency on y displacement
  if((arr = (TObjArray*)fContainer->At(kYRes))) {
    for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
      if(!(h2 = (TH2I*)arr->At(ily))) continue;
      TObjArray *arrg = (TObjArray*)fResults->At(kYRes);
      TGraphErrors *gym = (TGraphErrors*)arrg->At(0);
      TGraphErrors *gys = (TGraphErrors*)arrg->At(1);
      ax = h2->GetXaxis();
      for(Int_t ix=1; ix<=ax->GetNbins(); ix++){
        Float_t yd = ax->GetBinCenter(ix);
        h1 = h2->ProjectionY("py", ix, ix);
        if(h1->GetEntries() < 50) continue;
        Adjust(&f, h1);
        h1->Fit(&f, "Q");
    
        // Fill sy^2 = f(log(q))
        Int_t ip = gym->GetN();
        gym->SetPoint(ip, yd, 10.*f.GetParameter(1));
        gym->SetPointError(ip, 0., 10.*f.GetParError(1));
        gys->SetPoint(ip, yd, 1.e2*f.GetParameter(2)*f.GetParameter(2));
        gys->SetPointError(ip, 0., 2.e2*f.GetParameter(2)*f.GetParError(2));
      } 
    }
  } else AliWarning("Missing dy=f(y_d) container");

  // process resolution dependency on drift legth and drift cell width
  TGraphErrors *gm = new TGraphErrors(), *gs = new TGraphErrors();
  Float_t d(0.), x(0.);
  TH2F *hsx = (TH2F*)fResults->At(kSXRes);
  TH2F *hsy = (TH2F*)fResults->At(kSYRes);
  if((arr = (TObjArray*)fContainer->At(kSXRes))) {
    for(Int_t id=1; id<=fAd->GetNbins(); id++){
      d = fAd->GetBinCenter(id); //[mm]
      printf(" Doing d = %f[mm]\n", d);
      for(Int_t it=1; it<=fAt->GetNbins(); it++){
        x = fAt->GetBinCenter(it); //[mm]
        printf("    Doing xd = %f[mm]\n", x);
        Int_t idx = (id-1)*kNTB+it-1;
        if(!(h2 = (TH2I*)arr->At(idx))) {
          AliWarning(Form("Missing histo at index idx[%3d] [id[%2d] it[%2d]] xd[%f] d[%f]\n", idx, id, it, x, d));
          continue;
        }
  
        Int_t ip = 0;
        ax = h2->GetXaxis();
        for(Int_t ix=1; ix<=ax->GetNbins(); ix++){
          Float_t dydx = ax->GetBinCenter(ix) - fExB;
          if(dydx<0.) continue;
          h1 = h2->ProjectionY("py", ix, ix);
          if(h1->GetEntries() < 50) continue;
          Adjust(&f, h1);
          h1->Fit(&f, "Q", "", -.2, .2);
  
          // Fill sy^2 = f(tg_phi^2)
          gm->SetPoint(ip, dydx, 10.*f.GetParameter(1));
          gm->SetPointError(ip, 0., 10.*f.GetParError(1));
          gs->SetPoint(ip, dydx*dydx, 1.e2*f.GetParameter(2)*f.GetParameter(2));
          gs->SetPointError(ip, 0., 2.e2*f.GetParameter(2)*f.GetParError(2));
          ip++;
        }
        if(ip<4) continue;
        for(;ip<gm->GetN(); ip++){
          gm->RemovePoint(ip);gs->RemovePoint(ip);
        }
        gs->Fit(&sig, "Q");
        hsx->SetBinContent(id, it, sig.GetParameter(1));
        hsx->SetBinError(id, it, sig.GetParError(1));
        hsy->SetBinContent(id, it, sig.GetParameter(0));
        hsy->SetBinError(id, it, sig.GetParError(0));
      }
    }
  } else AliWarning("Missing dy=f(x_d, d_wire) container");
  delete gm; delete gs;

  return kTRUE;
}


