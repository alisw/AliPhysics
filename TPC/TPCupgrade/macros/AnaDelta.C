#include <TStyle.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <THn.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TString.h>
#include <TVectorT.h>
#include <TCanvas.h>
#include <TProfile2D.h>
#include <TGraphErrors.h>
#include <TTreeStream.h>

#include <AliExternalTrackParam.h>
#include <AliTPCComposedCorrection.h>
#include <AliTPCCorrectionLookupTable.h>

#include <AliToyMCEventGenerator.h>

/*

.L $ALICE_ROOT/TPC/Upgrade/macros/AnaDelta.C+g


*/
TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
TVectorD* MakeArbitraryBinning(const char* bins);

void DumpHn(THn *hn, TTreeSRedirector &stream);
void AnaDeltaBase(TString file, TString outDir=".");
void AnaDeltaTree(TString file, TString outFile="deltas_tree.root");

void AnaDelta(Int_t type, TString file, TString output="")
{
  switch (type) {
    case 0:
      AnaDeltaBase(file,output);
      break;
    case 1:
      AnaDeltaTree(file,output);
      break;
  }
}


void AnaDeltaBase(TString file, TString outDir)
{
  //
  //
  //

  gStyle->SetOptFit();
  
  TTreeSRedirector stream(Form("%s/deltas.root",outDir.Data()));
  gROOT->cd();
  
  TFile f(file);
  THn *hn=(THn*)f.Get("hn");

  DumpHn(hn, stream);

  delete hn;
}


void AnaDeltaTree(TString file, TString outFile)
{
  if (outFile.IsNull()) outFile="deltas_tree.root";
  TFile f(file);
  gROOT->cd();
  TTree *t = (TTree*)f.Get("delta");
  Float_t soneOverPt=0.;
  Float_t radius=0.;
  Float_t trackPhi=0.;
  Float_t trackY=0.;
  Float_t trackZ=0.;
  Float_t resRphi=0.;
  Double_t trackRes=0.;
  Float_t pointY=0.;
  Float_t pointZ=0.;
  Short_t npTRD=0.;
  Short_t event=0.;
  
  t->SetBranchAddress("soneOverPt" , &soneOverPt);
  t->SetBranchAddress("r"          , &radius);
  t->SetBranchAddress("trackPhi"   , &trackPhi);
  t->SetBranchAddress("trackY"     , &trackY);
  t->SetBranchAddress("trackZ"     , &trackZ);
  t->SetBranchAddress("resRphi"    , &resRphi);
  t->SetBranchAddress("trackRes"   , &trackRes);
  t->SetBranchAddress("pointY"     , &pointY);
  t->SetBranchAddress("pointZ"     , &pointZ);
  t->SetBranchAddress("npTRD"      , &npTRD);
  t->SetBranchAddress("event"      , &event);
  
  // make binning
  TVectorD *vR   = MakeLinBinning(10,86.,250.);
  TVectorD *vPhi = MakeLinBinning(18*8,0.,2*TMath::Pi());
  TVectorD *vZ   = MakeLinBinning(50,-250.,250.);

  const Int_t nbins=4;
  Int_t bins[nbins]    = {vR->GetNrows()-1, vPhi->GetNrows()-1, vZ->GetNrows()-1, 80};
//   Int_t bins[nbins]    = {16, 18*5, 50, 80};
  Double_t xmin[nbins] = {86. , 0.,           -250., -2.};
  Double_t xmax[nbins] = {250., 2*TMath::Pi(), 250.,  2.};
  THnF *hn = new THnF("hn", "hn", nbins, bins, xmin, xmax);

  hn->GetAxis(0)->Set(vR  ->GetNrows()-1, vR  ->GetMatrixArray());
  hn->GetAxis(1)->Set(vPhi->GetNrows()-1, vPhi->GetMatrixArray());
  hn->GetAxis(2)->Set(vZ  ->GetNrows()-1, vZ  ->GetMatrixArray());

  hn->GetAxis(0)->SetNameTitle("r","r (cm)");
  hn->GetAxis(1)->SetNameTitle("phi","#varphi");
  hn->GetAxis(2)->SetNameTitle("z","z (cm)");
  hn->GetAxis(3)->SetNameTitle("drphi","#Delta(r#varphi)");
  
  for (Int_t iev=0; iev<t->GetEntries(); ++iev) {
    t->GetEntry(iev);
    
    // cuts
    // -- on trd
    if (npTRD<2) continue;
    Double_t pt=1./TMath::Abs(soneOverPt);
    if (pt<0.8) continue;
    
    Float_t resRphiRandom = resRphi*trackRes;
    Float_t deviation     = pointY-(trackY+resRphiRandom);
    
    Double_t xx[4]={radius, trackPhi, trackZ ,deviation};
    hn->Fill(xx);
  }
  
  // do fits and fill tree
  TTreeSRedirector stream(outFile.Data());
  gROOT->cd();

  DumpHn(hn, stream);

  stream.GetFile()->cd();
  hn->Write();
  
  delete hn;
  delete vR;
  delete vPhi;
  delete vZ;
}


void AnaDeltaTree2(TString file/*, TString outDir="."*/)
{
  //
  // NOTE: not finished
  //
  TFile f(file);
  gROOT->cd();
  TTree *t = (TTree*)f.Get("delta");
  Float_t soneOverPt=0.;
  Float_t radius=0.;
  Float_t trackPhi=0.;
  Float_t trackY=0.;
  Float_t trackZ=0.;
  Float_t resRphi=0.;
  Float_t trackRes=0.;
  Float_t pointY=0.;
  Float_t pointZ=0.;
  Float_t npTRD=0.;
  Float_t event=0.;
  
  t->SetBranchAddress("soneOverPt" , &soneOverPt);
  t->SetBranchAddress("r"          , &radius);
  t->SetBranchAddress("trackPhi"   , &trackPhi);
  t->SetBranchAddress("trackY"     , &trackY);
  t->SetBranchAddress("trackZ"     , &trackZ);
  t->SetBranchAddress("resRphi"    , &resRphi);
  t->SetBranchAddress("trackRes"   , &trackRes);
  t->SetBranchAddress("pointY"     , &pointY);
  t->SetBranchAddress("pointZ"     , &pointZ);
  t->SetBranchAddress("npTRD"      , &npTRD);
  t->SetBranchAddress("event"      , &event);

  // make binning
  TVectorD *vZ   = MakeLinBinning(50,-250.,250.);
  TVectorD *vPhi = MakeLinBinning(18*8,0.,TMath::Pi());
  TVectorD *vR   = MakeLinBinning(16,86.,250.);

  TObjArray arrZ(vZ->GetNrows()-1);
  arrZ.SetOwner();
  
  for (Int_t iev=0; iev<t->GetEntries(); ++iev) {
    t->GetEntry(iev);
    
    // cuts
    // -- on trd
    if (npTRD<2) continue;

    Float_t resRphiRandom=resRphi*trackRes;

    Int_t binZ   = TMath::BinarySearch(vZ->GetNrows(),vZ->GetMatrixArray(),(Double_t)trackZ);
    Int_t binPhi = TMath::BinarySearch(vPhi->GetNrows(),vPhi->GetMatrixArray(),(Double_t)trackPhi);
    Int_t binR   = TMath::BinarySearch(vR->GetNrows(),vR->GetMatrixArray(),(Double_t)radius);

    if (binZ<0)   binZ=0;
    if (binPhi<0) binPhi=0;
    if (binR<0)   binR=0;

    TObjArray *arrPhi=(TObjArray*)arrZ.UncheckedAt(binZ);
    if (!arrPhi) {
      arrPhi=new TObjArray(vPhi->GetNrows()-1);
      arrZ.AddAt(arrPhi,binZ);
    }

    TObjArray *arrR=(TObjArray*)arrPhi->UncheckedAt(binPhi);
    if (!arrR) {
      arrR=new TObjArray(vR->GetNrows()-1);
      arrPhi->AddAt(arrR,binPhi);
    }

    TH1S *h = (TH1S*)arrR->UncheckedAt(binR);
    if (!h) {
      h = new TH1S(Form("h_%02d_%02d_%d02",binZ, binPhi, binR),
                   Form("z,phi,r: %02d,%02d,%d02; #Delta r#phi (cm)",binZ, binPhi, binR),
                   80, -2., 2.);
      arrR->AddAt(h, binR);
    }

    h->Fill(trackY+resRphiRandom-pointY);
  }

  // do fits and fill tree
}

void AnaDeltaResiduals(TString fluctuationMap, TString averageMap, TString outFile="deltas_residuals.root")
{
  //
  //
  //

  TFile fFluct(fluctuationMap);
  AliTPCCorrectionLookupTable *corrFluct = (AliTPCCorrectionLookupTable*)fFluct.Get("map");
  fFluct.Close();
  
  TFile fAverage(averageMap);
  AliTPCCorrectionLookupTable *corrAverage = (AliTPCCorrectionLookupTable*)fAverage.Get("map");
  fAverage.Close();
  
//   TObjArray *arrMaps = new TObjArray(2);
//   arrMaps->Add(corrAverage); // correction with the average Map
//   arrMaps->Add(corrFluct);   // distortion with the fluctuation Map
  
  // create the composed correction
  // if the weight are set to +1 and -1, the first map will be responsible for the distortions
  // The second map for the corrections
  // !!!!! In AliTPCComposedCorrection::GetDistortion MakeInverseIterator is called !!!!
  // for this reason we have to add the maps in the wrong order
  
//   AliTPCComposedCorrection *residualDistortion = new AliTPCComposedCorrection(arrMaps, AliTPCComposedCorrection::kQueueResidual);
  Float_t dummy=0;
//   TVectorD weights(2);
//   weights(0)=+1.;
//   weights(1)=-AliToyMCEventGenerator::GetSCScalingFactor(corrFluct, corrAverage,dummy);
//   residualDistortion->SetWeights(&weights);

  corrAverage->SetCorrScaleFactor(AliToyMCEventGenerator::GetSCScalingFactor(corrFluct, corrAverage,dummy));

  TVectorD *vR   = MakeLinBinning(10,86.,250.);
  TVectorD *vPhi = MakeLinBinning(18*8,0.,2*TMath::Pi());
  TVectorD *vZ   = MakeLinBinning(50,-250.,250.);
  
  const Int_t nbins=4;
  Int_t bins[nbins]    = {vR->GetNrows()-1, vPhi->GetNrows()-1, vZ->GetNrows()-1, 80};
  //   Int_t bins[nbins]    = {16, 18*5, 50, 80};
  Double_t xmin[nbins] = {86. , 0.,           -250., -2.};
  Double_t xmax[nbins] = {250., 2*TMath::Pi(), 250.,  2.};
  THnF *hn = new THnF("hn", "hn", nbins, bins, xmin, xmax);
  
  hn->GetAxis(0)->Set(vR  ->GetNrows()-1, vR  ->GetMatrixArray());
  hn->GetAxis(1)->Set(vPhi->GetNrows()-1, vPhi->GetMatrixArray());
  hn->GetAxis(2)->Set(vZ  ->GetNrows()-1, vZ  ->GetMatrixArray());
  
  AliExternalTrackParam vv;
  
  for (Float_t iz=-245; iz<=245; iz+=2) {
    Short_t roc=(iz>=0)?0:18;
    for (Float_t ir=86; ir<250; ir+=1) {
      for (Float_t iphi=0; iphi<TMath::TwoPi(); iphi+=0.5*TMath::DegToRad()){
        Float_t x=ir*(Float_t)TMath::Cos(iphi);
        Float_t y=ir*(Float_t)TMath::Sin(iphi);
        Float_t x3[3]    = {x,y,iz};
        Float_t x3dc[3]    = {x,y,iz};
        Float_t dx3[3]   = {0.,0.,0.};
//         residualDistortion->GetDistortion(x3,roc,dx3);
        corrFluct->DistortPoint(x3dc,roc);
        corrAverage->CorrectPoint(x3dc,roc);
        dx3[0]=x3dc[0]-x3[0];
        dx3[1]=x3dc[1]-x3[1];
        dx3[2]=x3dc[2]-x3[2];
        
        Double_t ddx3[3]={dx3[0], dx3[1], dx3[2]};
        vv.Global2LocalPosition(ddx3,iphi);

        Double_t xx[4]={ir, iphi, iz ,ddx3[1]};
        hn->Fill(xx);
        
      }
    }
  }

  TTreeSRedirector stream(outFile.Data());
  gROOT->cd();
  
  DumpHn(hn, stream);
  
  stream.GetFile()->cd();
  hn->Write();
  
  delete hn;
  delete vR;
  delete vPhi;
  delete vZ;
  
//   delete residualDistortion;
}

void DumpHn(THn *hn, TTreeSRedirector &stream)
{
  TAxis *ar   = hn->GetAxis(0);
  TAxis *aphi = hn->GetAxis(1);
  TAxis *az   = hn->GetAxis(2);

  // output Hn
  const Int_t nbins=3;
  Int_t bins[nbins]    = {1,1,1};
  Double_t xmin[nbins] = {0.,0.,0.};
  Double_t xmax[nbins] = {1.,1.,1.};
  THnF hnRes("hnRes", "hnRes", nbins, bins, xmin, xmax);

  ar  ->Copy(*hnRes.GetAxis(0));
  aphi->Copy(*hnRes.GetAxis(1));
  az  ->Copy(*hnRes.GetAxis(2));
  
  
  for (Int_t iz=0; iz<az->GetNbins(); ++iz) {
    az->SetRange(iz+1,iz+1);
    TObjArray arrFits;
    arrFits.SetName(Form("z_%02d",iz));
    arrFits.SetOwner();
    
    for (Int_t ir=0; ir<ar->GetNbins(); ++ir) {
      ar->SetRange(ir+1,ir+1);
      for (Int_t iphi=0; iphi<aphi->GetNbins(); ++iphi) {
        aphi->SetRange(iphi+1,iphi+1);

        Float_t cr       = 0.;
        Float_t cphi     = 0.;
        Float_t cz       = 0.;
        Float_t mean     = 0.;
        Float_t meanErr  = 0.;
        Float_t sigma    = 0.;
        Float_t sigmaErr = 0.;
        Int_t   entries  = 0.;
        Float_t chi2ndf  = 0.;
        Float_t mean2    = 0.;
        Float_t meanErr2 = 0.;
        Float_t rms2     = 0.;
        Float_t rmsErr2  = 0.;
        
        TH1 *hProj = hn->Projection(3);
        if (hProj->GetEntries()>1) {
          TF1 fg("fg","gaus",-2,2);
          cr   = ar->GetBinCenter(ir+1);
          cphi = aphi->GetBinCenter(iphi+1);
          cz   = az->GetBinCenter(iz+1);
          hProj->SetNameTitle(Form("h_%02d_%02d_%02d",iz, iphi, ir),
          Form("z,phi,r: %02d,%02d,%02d (%.2f, %.2f, %.2f)",iz,iphi,ir, cz, cphi, cr )
          );
          hProj->Fit(&fg,"LMQR");
          arrFits.Add(hProj);

          mean     = fg.GetParameter(1);
          meanErr  = fg.GetParError(1);
          sigma    = fg.GetParameter(2);
          sigmaErr = fg.GetParError(2);
          entries  = hProj->GetEntries();
          chi2ndf  = fg.GetChisquare()/fg.GetNDF();
          mean2    = hProj->GetMean();
          meanErr2 = hProj->GetMeanError();
          rms2     = hProj->GetRMS();
          rmsErr2  = hProj->GetRMSError();
        } else {
          delete hProj;
        }
        
        stream << "d" <<
        "ir="          << ir       <<
        "iphi="        << iphi     <<
        "iz="          << iz       <<
        "cr="          << cr       <<
        "cphi="        << cphi     <<
        "cz="          << cz       <<
        "mean="        << mean     <<
        "meanErr="     << meanErr  <<
        "sigma="       << sigma    <<
        "sigmaErr="    << sigmaErr <<
        "histMean="    << mean2    <<
        "histMeanErr=" << meanErr2 <<
        "histRMS="     << rms2    <<
        "histRMSErr="  << rmsErr2 <<
        "entries="     << entries  <<
        "chi2ndf="     << chi2ndf  <<
        "\n";

//        Double_t x[nbins]={cr, cphi, cz};
//        if (meanErr<0.3) hnRes.Fill(x,mean);
      }
    }
    stream.GetFile()->cd();
    arrFits.Write(0x0,TObject::kSingleKey);
    gROOT->cd();
  }

  stream.GetFile()->cd();
  hnRes.Write();
  gROOT->cd();
}

//______________________________________________________________________________
TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    printf("For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//______________________________________________________________________________
TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make linear binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* MakeArbitraryBinning(const char* bins)
{
  //
  // Make arbitrary binning, bins separated by a ','
  //
  TString limits(bins);
  if (limits.IsNull()){
    printf("Bin Limit string is empty, cannot add the variable");
    return 0x0;
  }
  
  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    printf("Need at leas 2 bin limits, cannot add the variable");
    delete arr;
    return 0x0;
  }
  
  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }
  
  delete arr;
  return binLimits;
}

void PlotFromTree(TTree *d, TString outDir=".")
{
  TCanvas *c=new TCanvas;
  gStyle->SetOptStat(0);
  d->SetMarkerStyle(20);
  d->SetMarkerSize(1);
  
  TProfile2D pRZ("pRZ",";z (cm); r(cm)",50,-250,250,10,85,250);
  d->Draw("entries:cr:cz>>pRZ","","profcolz");
  pRZ.GetZaxis()->UnZoom();
  c->SaveAs(Form("%s/entries_average.png",outDir.Data()));
  d->Draw("entries:cr:cz>>pRZ","iphi==2","profcolz");
  c->SaveAs(Form("%s/entries_onePhi.png",outDir.Data()));
  
  pRZ.SetMaximum(0.04);
  d->Draw("meanErr:cr:cz>>pRZ","","profcolz");
  c->SaveAs(Form("%s/meanErr_average.png",outDir.Data()));
  d->Draw("meanErr:cr:cz>>pRZ","iphi==2","profcolz");
  c->SaveAs(Form("%s/meanErr_onePhi.png",outDir.Data()));
  
  
  d->Draw("mean:cphi:cr","iz==25","colz");
  c->SaveAs(Form("%s/mean_oneZ_phi_allR.png",outDir.Data()));
  d->Draw("mean:meanErr:cphi","iz==25&&ir==2","goff");
  TGraphErrors *grmean_phi=new TGraphErrors(d->GetSelectedRows(),d->GetV3(),d->GetV1(),0,d->GetV2());
  grmean_phi->SetTitle(";#varphi;#LT#Delta r#varphi#GT");
  grmean_phi->SetMarkerStyle(20);
  grmean_phi->SetMarkerSize(1);
  grmean_phi->Draw("ap");
  c->SaveAs(Form("%s/mean_oneZ_phi_oneR.png",outDir.Data()));
  
  d->Draw("mean:cr:cphi","iz==25","colz");
  c->SaveAs(Form("%s/mean_oneZ_r_allPhi.png",outDir.Data()));
  
  d->Draw("mean:meanErr:cr","iz==25&&iphi==2","goff");
  TGraphErrors *grmean_r=new TGraphErrors(d->GetSelectedRows(),d->GetV3(),d->GetV1(),0,d->GetV2());
  grmean_r->SetTitle(";r (cm);#LT#Delta r#varphi#GT");
  grmean_r->SetMarkerStyle(20);
  grmean_r->SetMarkerSize(1);
  grmean_r->Draw("ap");
  c->SaveAs(Form("%s/mean_oneZ_r_onePhi.png",outDir.Data()));
  
  
  d->Draw("meanErr:cphi:cr","iz==25","colz");
  c->SaveAs(Form("%s/meanErr_oneZ_phi_allR.png",outDir.Data()));
  d->Draw("meanErr:cphi","iz==25&&ir==2");
  c->SaveAs(Form("%s/meanErr_oneZ_phi_oneR.png",outDir.Data()));
  
  d->Draw("meanErr:cr:cphi","iz==25","colz");
  c->SaveAs(Form("%s/meanErr_oneZ_r_allPhi.png",outDir.Data()));
  
  d->Draw("meanErr:cr","iz==25&&iphi==2");
  c->SaveAs(Form("%s/meanErr_oneZ_r_onePhi.png",outDir.Data()));
  
}


