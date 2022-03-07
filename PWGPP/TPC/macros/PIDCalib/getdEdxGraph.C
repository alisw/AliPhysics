TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);

Int_t lastRun = 0;
AliPIDResponse* p = 0x0;
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

  initPID(165772,"pass2");
  TGraph *grPio=getdEdxGraph(AliPID::kPion);
  TGraph *grPro=getdEdxGraph(AliPID::kProton);

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

  initPID(165772,"pass2");
  TGraph *grDeuteron     = getdEdxGraph(AliPID::kDeuteron);
  TGraph *grDeuteron_p3s = getdEdxGraph(AliPID::kDeuteron, 0, 0, kTRUE, kTRUE, kTRUE, +3);
  TGraph *grDeuteron_m3s = getdEdxGraph(AliPID::kDeuteron, 0, 0, kTRUE, kTRUE, kTRUE, -3);

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
  TGraph *grProtonTracking    = getdEdxGraphPIDReco(run, AliPID::kProton);
  TGraph *grProtonTracking_p5 = getdEdxGraphPIDReco(run, AliPID::kProton, 15);
  TGraph *grProtonTracking_m5 = getdEdxGraphPIDReco(run, AliPID::kProton, -15);

  grProtonTracking    -> SetLineColor(kBlack);
  grProtonTracking_p5 -> SetLineColor(kBlack);
  grProtonTracking_m5 -> SetLineColor(kBlack);

  grProtonTracking    -> SetLineStyle(kSolid);
  grProtonTracking_p5 -> SetLineStyle(kDashed);
  grProtonTracking_m5 -> SetLineStyle(kDashed);

  grProtonTracking    -> SetLineWidth(3);
  grProtonTracking_p5 -> SetLineWidth(3);
  grProtonTracking_m5 -> SetLineWidth(3);

  grProtonTracking   ->Draw("l");
  grProtonTracking_p5->Draw("l");
  grProtonTracking_m5->Draw("l");

  // if you would like to use the cvmfs OCDB, use the next line, otherwise omit it
  gSystem->Setenv("OCDB_PATH","/cvmfs/alice-ocdb.cern.ch");
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  Int_t run = 270824
  for (Int_t i=0; i<AliPID::kSPECIESC; ++i) {
    TGraph *grTracking    = getdEdxGraphPIDReco(run, AliPID::EParticleType(i));
    TGraph *grTracking_p5 = getdEdxGraphPIDReco(run, AliPID::EParticleType(i), 15);
    TGraph *grTracking_m5 = getdEdxGraphPIDReco(run, AliPID::EParticleType(i), -15);

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

  // Draw BB used for MC
  gSystem->Setenv("OCDB_PATH","/cvmfs/alice-ocdb.cern.ch");
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  Int_t run = 270824
  for (Int_t i=0; i<AliPID::kSPECIESC; ++i) {
    TF1* f = GetPIDforMCFunction(run);
    f->Draw("same");
  }

 */

void initPID(Int_t run, const TString recoPass, const TString customPIDResponse = "", const TString customEtaMaps = "", Bool_t isMC = kFALSE)
{
  Int_t recoPassNumber = 0;
  if (recoPass.Contains("pass1")) {
    recoPassNumber = 1;
  } else if (recoPass.Contains("pass2")) {
    recoPassNumber = 2;
  } else if (recoPass.Contains("pass3")) {
    recoPassNumber = 3;
  } else if (recoPass.Contains("pass4")) {
    recoPassNumber = 4;
  } else if (recoPass.Contains("pass5")) {
    recoPassNumber = 5;
  }

  if (run != lastRun) {
    delete p;
    p = new AliPIDResponse(isMC);
    p->SetOADBPath("$ALICE_PHYSICS/OADB");

    // force loading of corrections
    // if they are applied can be steered from the getdEdx function below
    p->SetUseTPCEtaCorrection(kTRUE);
    p->SetUseTPCMultiplicityCorrection(kTRUE);
    p->SetUseTPCPileupCorrection(kTRUE);

    if (!customPIDResponse.IsNull()) {
      p->SetCustomTPCpidResponseOADBFile(customPIDResponse);
    }

    if (!customEtaMaps.IsNull()) {
      p->SetCustomTPCetaMaps(customEtaMaps);
    }

    AliESDEvent ev;
    p->InitialiseEvent(&ev, recoPassNumber, recoPass, run);
    lastRun = run;
  }
}

TGraph* getdEdxGraph(AliPID::EParticleType particle, Float_t eta = 0, Float_t multiplicity = 0,
                     Bool_t useEtaCorrection = kTRUE, Bool_t useMultiplicityCorrection = kTRUE, Bool_t usePileupCorrection = kTRUE,
                     Double_t nSigma = 0, Double_t xmin = 0.1, Double_t xmax = 20., Double_t dEdxMax = 1000)
{

  AliTPCPIDResponse& tpcpid = p->GetTPCResponse();
  tpcpid.SetCurrentEventMultiplicity(multiplicity);

  AliESDtrack tr;
  tr.SetTPCsignal(0, 0, 120);

  const Double_t phi = 0.;

  TGraph* gr = new TGraph;
  TVectorD* bins = MakeLogBinning(200, xmin, xmax);
  Double_t xyz[3] = {0., 0., 0.};
  Double_t cv[21] = {0.};
  for (Int_t ibin = 0; ibin < bins->GetNrows(); ++ibin) {
    const Double_t pabs = (*bins)[ibin];
    const Double_t theta = 2 * TMath::ATan(TMath::Exp(-eta));
    const Double_t pz = pabs * std::cos(theta);
    const Double_t pt = pabs * std::sin(theta);
    const Double_t px = pt * std::sin(phi);
    const Double_t py = pt * std::cos(phi);
    Double_t pxyz[3] = {px, py, pz};
    tr.Set(xyz, pxyz, cv, 1);
    Double_t dEdx = tpcpid.GetExpectedSignal(&tr, particle, AliTPCPIDResponse::kdEdxDefault, useEtaCorrection, useMultiplicityCorrection, usePileupCorrection);
    const Double_t sigma = tpcpid.GetExpectedSigma(&tr, particle, AliTPCPIDResponse::kdEdxDefault, useEtaCorrection, useMultiplicityCorrection, usePileupCorrection);
    printf("%3d, %.2f : %.2f +- %.2f (%.2f) %.2f\n", ibin, pabs, dEdx, sigma, sigma / dEdx, dEdx + nSigma * sigma);
    dEdx += nSigma * sigma;
    if (dEdx > dEdxMax) {
      continue;
    }
    gr->SetPoint(gr->GetN(), pabs, dEdx);
  }

  delete bins;
  return gr;
}

TGraph* getdEdxGraph(Int_t run, TString recoPass, AliPID::EParticleType particle,
                     Bool_t isMC = kFALSE, Double_t nSigma = 0,
                     Double_t xmin = 0.1, Double_t xmax = 20., Double_t dEdxMax = 1000)
{
  initPID(run, recoPass, "", "", isMC);
  return getdEdxGraph(particle, 0, 0, kTRUE, kTRUE, kTRUE, nSigma, xmin, xmax, dEdxMax);
}

//______________________________________________________________________________
TGraph* getdEdxGraphPIDReco(Int_t run, AliPID::EParticleType particle, Double_t nSigma = 0,
                            Double_t xmin = 0.1, Double_t xmax = 20., Double_t dEdxMax = 1000)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->GetDefaultStorage() && !man->GetRaw()) {
    std::cout << "ERROR: Set up CDB manager properly, with the required specific storages, e.g." << std::endl;
    std::cout << "AliCDBManager *man = AliCDBManager::Instance();" << std::endl;
    std::cout << "man->SetDefaultStorage(\"raw://\");" << std::endl;
    return 0x0;
  }

  if (run != lastRun) {
    man->SetRun(run);

    AliGRPManager grpMan;
    grpMan.ReadGRPEntry();
    grpMan.SetMagField();

    AliMagF* field = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField());
    AliTPCcalibDB::Instance()->SetExBField(field);

    AliTPCcalibDB::Instance()->SetRun(run);
    delete p;
    p = new AliESDpid();

    AliTPCReconstructor r;
    r.GetPidSettings((AliESDpid*)p);
    lastRun = run;
  }

  // return getdEdxGraph(run, "0", particle, kFALSE, nSigma, xmin, xmax, dEdxMax);
  return getdEdxGraph(particle, 0, 0, kTRUE, kTRUE, kTRUE, nSigma, xmin, xmax, dEdxMax);
}

//______________________________________________________________________________
/// return Bethe Bloch function for used in MC simulation or during reconstruction
/// \param run run number
/// \param xmin minimum momentum
/// \param xmax maximum momentum
/// \param particle particle type
/// \param type 0: MC BB parameters, 1: PID for reconstruction BB parameters
TF1* GetPIDFunction(Int_t run, Double_t xmin = 0.1, Double_t xmax = 20., AliPID::EParticleType particle = AliPID::EParticleType(AliPID::kSPECIESC), Int_t type = 0)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->GetDefaultStorage() && !man->GetRaw()) {
    std::cout << "ERROR: Set up CDB manager properly, e.g." << std::endl;
    std::cout << "AliCDBManager *man = AliCDBManager::Instance();" << std::endl;
    std::cout << "man->SetDefaultStorage(\"raw://\");" << std::endl;
    return 0x0;
  }

  man->SetRun(run);

  AliCDBEntry* parEntry = (AliCDBEntry*)man->Get("TPC/Calib/Parameters");
  AliTPCParamSR* par = (AliTPCParamSR*)parEntry->GetObject();

  TVectorD* vBBPID = par->GetBetheBlochParameters();
  if (type == 1)
    vBBPID = par->GetBetheBlochParametersMC();
  Double_t mass = 1;
  if (particle != AliPID::kSPECIESC) {
    mass = AliPID::ParticleMass(particle);
  }
  TVectorD vParams(*vBBPID);
  vParams.ResizeTo(6);
  vParams(5) = mass;

  TF1* funcBB = new TF1(Form("fBBPidForReco_%d_%s", run, AliPID::ParticleShortName(particle)), "50*AliExternalTrackParam::BetheBlochAleph(x/[5],[0],[1],[2],[3],[4])", xmin, xmax);
  funcBB->SetParameters(vParams.GetMatrixArray());
  funcBB->SetNpx(500);

  return funcBB;
}

//______________________________________________________________________________
/// return Bethe Bloch function for used in during reconstruction
/// \see GetPIDFunction
TF1* GetPIDforRecoFunction(Int_t run, Double_t xmin = 0.1, Double_t xmax = 20., AliPID::EParticleType particle = AliPID::EParticleType(AliPID::kSPECIESC))
{
  return GetPIDFunction(run, xmin, xmax, particle, 0);
}

//______________________________________________________________________________
/// return Bethe Bloch function for used in MC simulation
/// \see GetPIDFunction
TF1* GetPIDforMCFunction(Int_t run, Double_t xmin = 0.1, Double_t xmax = 20., AliPID::EParticleType particle = AliPID::EParticleType(AliPID::kSPECIESC))
{
  return GetPIDFunction(run, xmin, xmax, particle, 1);
}

//______________________________________________________________________________
TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //

  // check limits
  if (xmin < 1e-20 || xmax < 1e-20) {
    printf("ERROR: For Log binning xmin and xmax must be > 1e-20.");
    return 0x0;
  }
  if (xmax < xmin) {
    Double_t tmp = xmin;
    xmin = xmax;
    xmax = tmp;
  }
  TVectorD* binLim = new TVectorD(nbinsX + 1);
  Double_t first = xmin;
  Double_t last = xmax;
  Double_t expMax = TMath::Log(last / first);
  for (Int_t i = 0; i < nbinsX + 1; ++i) {
    (*binLim)[i] = first * TMath::Exp(expMax / nbinsX * (Double_t)i);
    printf("%d, %.2f\n", i, (*binLim)[i]);
  }
  return binLim;
}

//______________________________________________________________________________
void DrawNsigmaLinesPIDtracking(const int run, const float nSigma = -1., const float xMin = 0.1, const float xMax = 20., const float dEdxMax = 5000., const int nparticles = int(AliPID::kSPECIES))
{
  AliCDBManager* man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");

  // static int lastRun = 0;
  //  get configures sigma from OCDB
  const int colors[AliPID::kSPECIESC] = {kBlack, kGray + 2, kRed + 3, kOrange + 2, kYellow + 2, kGreen + 2, kBlue - 2, kMagenta, kMagenta + 2};
  // draw n-sigma lines for all default particle species
  for (int ipart = 0; ipart < nparticles; ++ipart) {
    TGraph* grTracking = getdEdxGraphPIDReco(run, (AliPID::EParticleType)ipart, 0, xMin, xMax, dEdxMax);
    float sigma = AliTPCcalibDB::Instance()->GetParameters()->GetSigmaRangePIDinTracking();
    if (run != lastRun) {
      printf("Found sigma setting: %.1f\n", sigma);
      lastRun = run;
    }
    if (nSigma >= 0) {
      printf("Using simga setting %.2f instead of %.2f\n", nSigma, sigma);
      sigma = nSigma;
    }

    TGraph* grTracking_p5 = getdEdxGraphPIDReco(run, (AliPID::EParticleType)ipart, sigma, xMin, xMin, dEdxMax);
    TGraph* grTracking_m5 = getdEdxGraphPIDReco(run, (AliPID::EParticleType)ipart, -1 * sigma, xMin, xMin, dEdxMax);

    grTracking->SetLineColor(colors[ipart]);
    grTracking_p5->SetLineColor(colors[ipart]);
    grTracking_m5->SetLineColor(colors[ipart]);

    // grTracking    -> SetLineColor(kBlack);
    // grTracking_p5 -> SetLineColor(kBlack);
    // grTracking_m5 -> SetLineColor(kBlack);

    grTracking->SetLineStyle(kSolid);
    grTracking_p5->SetLineStyle(kDashed);
    grTracking_m5->SetLineStyle(kDashed);

    grTracking->SetLineWidth(3);
    grTracking_p5->SetLineWidth(3);
    grTracking_m5->SetLineWidth(3);

    grTracking->Draw("l");
    grTracking_p5->Draw("l");
    grTracking_m5->Draw("l");
  }
}

//______________________________________________________________________________
TObject* GetObjectFromPath(TObject* iterable, const TString path)
{
  if (!iterable)
    return 0x0;
  TObjArray* arr = path.Tokenize("/");
  for (int i = 0; i < arr->GetEntriesFast(); ++i) {
    TString& name = ((static_cast<TObjString*>(arr->At(i)))->String());
    printf("%s: %s\n", iterable->GetName(), name.Data());
    if (iterable->InheritsFrom(TDirectory::Class())) {
      iterable = static_cast<TDirectory*>(iterable)->Get(name);
    } else {
      iterable = iterable->FindObject(name);
    }
    if (!iterable)
      return 0x0;
  }
  return iterable;
}
