/*

 // example usage:
 TH2D h2("h2",";p (GeV/#it{c});d#it{E}/d#it{x} (arb. units)",100,0,20,100,0,200);
 h2.Draw();
 .L $ALICE_ROOT/ANALYSIS/macros/DrawTPCPIDSplines.C
 DrawTPCPIDSplines(122374,"pass2");



*/

void DrawTPCPIDSplines(Int_t run, TString recoPass, Float_t dEdxMin=0., Float_t dEdxMax=200, Float_t pMin=0.1, Float_t pMax=20.)
{
  //
  // simple macro to draw splines using the AliPIDResponse class
  // only splines for electrons, pions, kaons and protons 
  // NOTE: This is only for simple overlaying purposes not for QA.
  //       The TPC pid depends on eta, the splines are eta averaged
  //         so an overlay could look wrong, depending on the eta
  //         of the track selection.
  //       For QA purposes the nSigma values with full eta correction
  //         should be used as done e.g. in AliAnalysisTaskPIDqa
  //

  Int_t recoPassNumber = 0;
  if (recoPass.Contains("pass1") ) {
    recoPassNumber=1;
  } else if (recoPass.Contains("pass2") ) {
    recoPassNumber=2;
  } else if (recoPass.Contains("pass3") ) {
    recoPassNumber=3;
  } else if (recoPass.Contains("pass4") ) {
    recoPassNumber=4;
  } else if (recoPass.Contains("pass5") ) {
    recoPassNumber=5;
  }

  AliPIDResponse *pidres = new AliPIDResponse(kFALSE);
  //
  pidres->SetOADBPath("$ALICE_PHYSICS/OADB");
  pidres->SetCachePID(kFALSE);
  pidres->SetUseTPCEtaCorrection(kFALSE);
  pidres->SetUseTPCMultiplicityCorrection(kFALSE);
  pidres->SetTunedOnData(kFALSE,2);
  pidres->SetCurrentAliRootRev(62720);

  AliESDEvent *ev=new AliESDEvent;
  pidres->InitialiseEvent(ev,recoPassNumber,recoPass,run);

  const  Double_t first=pMin;
  const  Double_t last=pMax;
  const  Double_t expMax=TMath::Log(last/first);
  const  Int_t nbinsX=100;

  AliTPCPIDResponse &tpcpid=pidres->GetTPCResponse();

  for (Int_t i=0; i<AliPID::kSPECIES; ++i){
    if (i==AliPID::kMuon) continue;
    TGraph *gr=new TGraph;
    gr->SetNameTitle(Form("Spline_Graph_%s",AliPID::ParticleName(i)),Form("%s splines;p (GeV/#it{c});d#it{E}/d#it{x} (arb. units)",AliPID::ParticleName(i)));
    for (Int_t ip=0; ip<nbinsX; ++ip) {
      Double_t p=first*TMath::Exp(expMax/nbinsX*(Double_t)ip);
      Double_t dEdx=tpcpid.GetExpectedSignal(p,i);
      if (p<pMin || p>pMax || dEdx<dEdxMin || dEdx>dEdxMax ) continue;
      gr->SetPoint(gr->GetN(),p,dEdx);
    }
    gr->Draw("c");
    gr->SetLineWidth(2);
  }
}


