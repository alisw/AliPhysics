TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);

Int_t lastRun=0;
AliPIDResponse *p=0x0;
AliPID dummy;

/*
  // Examples how to use this macro
  // ===| Draw graphs from the OADB -> used in Analysis |=======================
  //
  // ---| Draw pion graph and proton graph |------------------------------------
  .L $ALICE_PHYSICS/PWGPP/TPC/macros/PIDCalib/getdEdxGraph.C

  TCanvas cdEdx("cdEdx","dE/dx");
  cdEdx.SetLogx();

  TH1F hDummy("hDummy",";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)",100,0,20);
  hDummy.SetMaximum(1000);
  hDummy.Draw();

  TGraph *grPio=getdEdxGraph(165772,"pass2",AliPID::kPion);
  TGraph *grPro=getdEdxGraph(165772,"pass2",AliPID::kProton);

  grPio->SetLineColor(kRed);
  grPro->SetLineColor(kBlue);
  grPio->Draw("l");
  grPro->Draw("l");

  // ---| Example: Draw Deuteron line and +- 3sigma |---------------------------
  .L $ALICE_PHYSICS/PWGPP/TPC/macros/PIDCalib/getdEdxGraph.C

  TCanvas cdEdx("cdEdx","dE/dx");
  cdEdx.SetLogx();
  TH1F hDummy("hDummy",";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)",100,0,20);
  hDummy.SetMaximum(1000);
  hDummy.Draw();

  TGraph *grDeuteron     = getdEdxGraph(165772,"pass2",AliPID::kDeuteron);
  TGraph *grDeuteron_p3s = getdEdxGraph(165772,"pass2",AliPID::kDeuteron, kFALSE, +3);
  TGraph *grDeuteron_m3s = getdEdxGraph(165772,"pass2",AliPID::kDeuteron, kFALSE, -3);

  grDeuteron    ->SetLineColor(kRed);
  grDeuteron_m3s->SetLineColor(kRed);
  grDeuteron_p3s->SetLineColor(kRed);
  grDeuteron_m3s->SetLineStyle(kDotted);
  grDeuteron_p3s->SetLineStyle(kDotted);

  grDeuteron    ->Draw("l");
  grDeuteron_m3s->Draw("l");
  grDeuteron_p3s->Draw("l");

  // ===| Draw graphs for PID in tracking |=====================================
  //
  // ---| Draw proton graphs and +- 5sigma lines |------------------------------
  .L $ALICE_PHYSICS/PWGPP/TPC/macros/PIDCalib/getdEdxGraph.C

  TCanvas cdEdx("cdEdx","dE/dx");
  cdEdx.SetLogx();

  TH1F hDummy("hDummy",";#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} (arb. units)",100,0,20);
  hDummy.SetMaximum(1000);
  hDummy.Draw();

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  TGraph *grProtonTracking    = getdEdxGraphPIDReco(225000, AliPID::kProton)
  TGraph *grProtonTracking_p5 = getdEdxGraphPIDReco(225000, AliPID::kProton,15)
  TGraph *grProtonTracking_m5 = getdEdxGraphPIDReco(225000, AliPID::kProton,-15)

  grProtonTracking    -> SetLineColor(kBlack);
  grProtonTracking_p5 -> SetLineColor(kBlack);
  grProtonTracking_m5 -> SetLineColor(kBlack);

  grProtonTracking    -> SetLineStyle(kSolid);
  grProtonTracking_p5 -> SetLineStyle(kDashed);
  grProtonTracking_m5 -> SetLineStyle(kDashed);

  grProtonTracking    -> SetLineWidth(3);
  grProtonTracking_p5 -> SetLineWidth(3);
  grProtonTracking_m5 -> SetLineWidth(3);

  grProtonTracking   ->Draw("l")
  grProtonTracking_p5->Draw("l")
  grProtonTracking_m5->Draw("l")

 */

TGraph* getdEdxGraph(Int_t run, TString recoPass, AliPID::EParticleType particle,
                     Bool_t isMC=kFALSE, Double_t nSigma=0,
                     Double_t xmin=0.1, Double_t xmax=20., Double_t dEdxMax=1000)
{
  AliESDEvent ev;

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

  if (run!=lastRun) {
    delete p;
    p=new AliPIDResponse(isMC);
    p->SetUseTPCMultiplicityCorrection();
    p->SetOADBPath("$ALICE_PHYSICS/OADB");
    p->InitialiseEvent(&ev,recoPassNumber, recoPass, run);
    lastRun=run;
  }
  AliTPCPIDResponse &tpcpid=p->GetTPCResponse();
  tpcpid.SetCurrentEventMultiplicity(5000);
  AliESDtrack tr;
  tr.SetTPCsignal(0,0,120);

  TGraph *gr=new TGraph;
  TVectorD *bins=MakeLogBinning(200,xmin,xmax);
  Double_t xyz[3]={0.,0.,0.};
  Double_t cv[21]={0.};
  for (Int_t ibin=0; ibin<bins->GetNrows(); ++ibin) {
    const Double_t p     = (*bins)[ibin];
    Double_t pxyz[3]={p,0.,0.};
    tr.Set(xyz, pxyz, cv, 1);
    //const Double_t dEdx  = tpcpid.GetExpectedSignal(p,particle);
    //const Double_t sigma = tpcpid.GetExpectedSigma(p,100,particle);
    Double_t dEdx  = tpcpid.GetExpectedSignal(&tr,particle,AliTPCPIDResponse::kdEdxDefault,kFALSE,kTRUE);
    const Double_t sigma = tpcpid.GetExpectedSigma(&tr,particle,AliTPCPIDResponse::kdEdxDefault,kFALSE,kTRUE);
    printf("%.2f +- %.2f (%.2f) %.2f\n", dEdx, sigma, sigma/dEdx, dEdx+nSigma*sigma);
    dEdx += nSigma*sigma;
    if (dEdx>dEdxMax) continue;
    gr->SetPoint(gr->GetN(), p, dEdx);
  }

  delete bins;
  return gr;
}

//______________________________________________________________________________
TGraph* getdEdxGraphPIDReco(Int_t run, AliPID::EParticleType particle, Double_t nSigma=0,
                            Double_t xmin=0.1, Double_t xmax=20., Double_t dEdxMax=1000)
{
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man->GetDefaultStorage() && !man->GetRaw()) {
    std::cout << "ERROR: Set up CDB manager properly, with the required specific storages, e.g." << std::endl;
    std::cout << "AliCDBManager *man = AliCDBManager::Instance();" << std::endl;
    std::cout << "man->SetDefaultStorage(\"raw://\");" << std::endl;
    return 0x0;
  }

  if (run!=lastRun) {
    man->SetRun(run);

    AliGRPManager grpMan;
    grpMan.ReadGRPEntry();
    grpMan.SetMagField();

    AliMagF * field = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField());
    AliTPCcalibDB::Instance()->SetExBField(field);

    AliTPCcalibDB::Instance()->SetRun(run);
    delete p;
    p=new AliESDpid();

    AliTPCReconstructor r;
    r.GetPidSettings((AliESDpid*)p);
    lastRun=run;
  }

  return getdEdxGraph(run, "0", particle, kFALSE, nSigma, xmin, xmax, dEdxMax);
}

//______________________________________________________________________________
TF1* getPIDforRecoFunction(Int_t run, Double_t xmin=0.1, Double_t xmax=20., AliPID::EParticleType particle=AliPID::EParticleType(AliPID::kSPECIESC))
{
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man->GetDefaultStorage() && !man->GetRaw()) {
    std::cout << "ERROR: Set up CDB manager properly, e.g." << std::endl;
    std::cout << "AliCDBManager *man = AliCDBManager::Instance();" << std::endl;
    std::cout << "man->SetDefaultStorage(\"raw://\");" << std::endl;
    return 0x0;
  }

  man->SetRun(run);

  AliCDBEntry *parEntry=(AliCDBEntry*)man->Get("TPC/Calib/Parameters");
  AliTPCParamSR *par=(AliTPCParamSR*)parEntry->GetObject();

  TVectorD *vBBPID = par->GetBetheBlochParameters();
  Double_t mass=1;
  if (particle!=AliPID::kSPECIESC) {
    mass=AliPID::ParticleMass(particle);
  }
  TVectorD vParams(*vBBPID);
  vParams.ResizeTo(6);
  vParams(5)=mass;

  TF1* funcBB = new TF1("fBBPidForReco", "50*AliExternalTrackParam::BetheBlochAleph(x/[5],[0],[1],[2],[3],[4])", xmin, xmax);
  funcBB->SetParameters(vBBPID->GetMatrixArray());

  return funcBB;
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
    printf("ERROR: For Log binning xmin and xmax must be > 1e-20.");
    return 0x0;
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
void DrawNsigmaLinesPIDtracking(const int run, const float nSigma=15.f, const int nparticles=int(AliPID::kSPECIES))
{
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");

  const int colors[AliPID::kSPECIES] = {kBlack, kRed-2, kBlue-2, kMagenta, kGray+2};
  // draw n-sigma lines for all default particle species
  for (int ipart=0; ipart<nparticles; ++ipart) {
    TGraph *grTracking    = getdEdxGraphPIDReco(run, AliPID::EParticleType(ipart));
    TGraph *grTracking_p5 = getdEdxGraphPIDReco(run, AliPID::EParticleType(ipart),nSigma);
    TGraph *grTracking_m5 = getdEdxGraphPIDReco(run, AliPID::EParticleType(ipart),-1*nSigma);

    //grTracking    -> SetLineColor(colors[ipart]);
    //grTracking_p5 -> SetLineColor(colors[ipart]);
    //grTracking_m5 -> SetLineColor(colors[ipart]);

    grTracking    -> SetLineColor(kBlack);
    grTracking_p5 -> SetLineColor(kBlack);
    grTracking_m5 -> SetLineColor(kBlack);

    grTracking    -> SetLineStyle(kSolid);
    grTracking_p5 -> SetLineStyle(kDashed);
    grTracking_m5 -> SetLineStyle(kDashed);

    grTracking    -> SetLineWidth(3);
    grTracking_p5 -> SetLineWidth(3);
    grTracking_m5 -> SetLineWidth(3);

    grTracking   ->Draw("l");
    grTracking_p5->Draw("l");
    grTracking_m5->Draw("l");
  }
}

//______________________________________________________________________________
TObject* GetObjectFromPath(TObject *iterable, const TString path)
{
  if (!iterable) return 0x0;
  TObjArray *arr = path.Tokenize("/");
  for (int i=0; i<arr->GetEntriesFast(); ++i) {
    TString& name = ((static_cast<TObjString*>(arr->At(i)))->String());
    printf("%s: %s\n", iterable->GetName(), name.Data());
    if (iterable->InheritsFrom(TDirectory::Class())) {
      iterable = static_cast<TDirectory*>(iterable)->Get(name);
    }
    else {
      iterable = iterable->FindObject(name);
    }
    if (!iterable) return 0x0;
  }
  return iterable;
}
