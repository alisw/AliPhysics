/*!
   \file AliMCSpectraWeights.cxx
   \brief Class for re-weighting particle abundances in MC simulations
   \author Patrick Huhn
   \date 25/10/2019
*/
#include "AliMCSpectraWeights.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"

ClassImp(AliMCSpectraWeights);

AliMCSpectraWeights::AliMCSpectraWeights()
  : TNamed(), // default constructor
  fHistMCGenPrimTrackParticle(0), fHistDataFractions(0), fHistMCFractions(0), fHistMCWeights(0),
  fBinsPt(0), fBinsMultCent(0), fstCollisionSystem("pp"),
  fNPartTypes(6), fstPartTypes(0), fNCentralities(0), fstCentralities(0),
  fstFileMCSpectra(""), fstFilePublished(""),
  fstSavedObjName("fMCSpectraWeights"),fstSavedListName("dNdPt_ParCor"),
  fUseMultiplicity(kTRUE), fbTaskStatus(0), fFlag(AliMCSpectraWeights::SysFlag::kNominal) {

    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
  }

// non-default contructor; way to go
AliMCSpectraWeights::AliMCSpectraWeights(const char *collisionSystem, const char* name, AliMCSpectraWeights::SysFlag flag)
  : TNamed(name, name), fHistMCGenPrimTrackParticle(0), fHistDataFractions(0), fHistMCFractions(0),
  fHistMCWeights(0), fBinsPt(0), fBinsMultCent(0),
  fstCollisionSystem("pp"), fNPartTypes(6), fstPartTypes(0), fNCentralities(0), fstCentralities(0),
  fstFileMCSpectra(""), fstFilePublished(""),
  fstSavedObjName("fMCSpectraWeights"),fstSavedListName("dNdPt_ParCor"),
  fUseMultiplicity(kTRUE), fbTaskStatus(0), fFlag(flag) {

    fstCollisionSystem = collisionSystem;
    fstCollisionSystem.ToLower();
    //setting uniform name
    if(fstCollisionSystem=="pp" || fstCollisionSystem=="p-p") fstCollisionSystem = "pp";
    else if(fstCollisionSystem=="ppb" || fstCollisionSystem=="p-pb") fstCollisionSystem="ppb";
    else if(fstCollisionSystem=="pbpb" || fstCollisionSystem=="pb-pb") fstCollisionSystem="pbpb";
    else if(fstCollisionSystem=="xexe" || fstCollisionSystem=="xe-xe") fstCollisionSystem="xexe";
    else fstCollisionSystem="pp";

    // set default Binning
    // pT binning
    Double_t bins[44] = {0.0,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,   0.45, 0.5,
      0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85,  0.9,  0.95,
      1.0,  1.1,  1.2,  1.4,  1.6,  1.8,  2.0,   2.2,  2.4,
      2.6,  2.8,  3.0,  3.2,  3.6,  4.0,  5.0,   6.0,  8.0,
      10.0, 13.0, 20.0, 30.0, 50.0, 80.0, 100.0, 200.0};
    fBinsPt = new TArrayD(44, bins);
    // if(bins) delete bins;

    // multiplicity binning
    if (fstCollisionSystem.Contains("pp")) {
      fNCentralities = 10;
      fstCentralities = new TString[fNCentralities];
      for (int icent = 0; icent < fNCentralities; icent++) {
        fstCentralities[icent] = Form("%d", icent);
      }
      Double_t mulBins[51] = {0};
      for (int i = 0; i < 51; ++i) {
        mulBins[i] = i;
      }
      fBinsMultCent = new TArrayD(51, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("ppb")){
      Double_t mulBins[301] = {0};
      for (int i = 0; i < 301; ++i) {
        mulBins[i] = i;
      }
      fBinsMultCent = new TArrayD(301, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("pbpb")) {
      Double_t mulBins[201] = {0};
      for (int i = 0; i < 201; ++i) {
        mulBins[i] = i * 25;
      }
      fBinsMultCent = new TArrayD(201, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("xexe")) {
      Double_t mulBins[201] = {0};
      for (int i = 0; i < 201; ++i) {
        mulBins[i] = i * 25;
      }
      fBinsMultCent = new TArrayD(201, mulBins);
      // if(mulBins) delete mulBins;
    } else {
      Double_t mulBins[101] = {0};
      for (int i = 0; i < 101; ++i) {
        mulBins[i] = i;
      }
      fBinsMultCent = new TArrayD(101, mulBins);
      // if(mulBins) delete mulBins;
    }
    if(!fstCollisionSystem.Contains("pp"))
    {
      TString cent[] = {"MB", "c0005", "c0510", "c0010", "c1020", "c2040", "c4060", "c6080"};
      fNCentralities = 8;
      fstCentralities = new TString[fNCentralities];
      for (int icent = 0; icent < fNCentralities; icent++) {
        fstCentralities[icent] = cent[icent];
      }
    }

    fstPartTypes = new TString[6];
    fstPartTypes[AliMCSpectraWeights::ParticleType::kPion] = "Pion";
    fstPartTypes[AliMCSpectraWeights::ParticleType::kProtons] = "Proton";
    fstPartTypes[AliMCSpectraWeights::ParticleType::kKaon] = "Kaon";
    fstPartTypes[AliMCSpectraWeights::ParticleType::kSigmaMinus] = "SigmaMinus";
    fstPartTypes[AliMCSpectraWeights::ParticleType::kSigmaPlus] = "SigmaPlus";
    fstPartTypes[AliMCSpectraWeights::ParticleType::kRest] = "Rest";

    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
    fstFilePublished = "alien:///alice/cern.ch/user/p/phuhn/AllPublishedFractions.root";
    // fstFilePublished = "/home/alidock/particle-composition-correction/data/published/ppMult/ppDataMultFractions.root";
  }

AliMCSpectraWeights::~AliMCSpectraWeights() {
  // if (fHistMCGenPrimTrackParticle)
  // delete fHistMCGenPrimTrackParticle;
  if (fHistDataFractions)
    delete fHistDataFractions;
  if (fHistMCFractions)
    delete fHistMCFractions;
  if (fHistMCWeights)
    delete fHistMCWeights;
  if (fBinsPt)
    delete fBinsPt;
  if (fBinsMultCent)
    delete fBinsMultCent;
  if (fstPartTypes)
    delete fstPartTypes;
}

void AliMCSpectraWeights::Init() {
  // Histograms
  AliMCSpectraWeights::InitHistos();

  if (fstFileMCSpectra.Length() > 5) // *.root
  {
    // printf("AliMCSpectraWeights:: Opening %s \n", fstFileMCSpectra.Data());
    TFile* fInput = TFile::Open(fstFileMCSpectra.Data());
    if(fInput)
    {
      if (fInput->GetNkeys() != 1)
      {
        if(!fstSavedListName.Contains("dNdPt_ParCor"))
          printf("AliMCSpectraWeights::WARNING: more than 1 list in the streamed file; please specify; using 1st list;\n\n");
      }
      else fstSavedListName  = fInput->GetListOfKeys()->At(0)->GetName();
      // printf("AliMCSpectraWeights:: Loading %s from list %s\n", fstSavedObjName.Data(), fstSavedListName.Data());
      TList *listMC = (TList*)fInput->Get(fstSavedListName);
      if(!listMC) {printf("AliMCSpectraWeights::ERROR: could not load list in streamed file\n");}
      else{
        fHistMCGenPrimTrackParticle = (THnF*)listMC->FindObject("fHistMCGenPrimTrackParticle");
        if(!fHistMCGenPrimTrackParticle)
        {
          printf("AliMCSpectraWeights::WARNING: Couln't get fHistMCGenPrimTrackParticle\n" );
          AliMCSpectraWeights* inWeights = (AliMCSpectraWeights*)listMC->FindObject(fstSavedObjName.Data());
          if(AliMCSpectraWeights::LoadFromAliMCSpectraWeight(inWeights))
            fbTaskStatus = AliMCSpectraWeights::TaskState::kMCSpectraObtained;
        }
        if(fHistMCGenPrimTrackParticle->GetEntries()>0) fbTaskStatus = AliMCSpectraWeights::TaskState::kMCSpectraObtained;
        // if(inWeights) delete inWeights;
      }
    }
    // else printf("AliMCSpectraWeights::WARNING: %s can not be loaded\n ", fstFileMCSpectra.Data());
    if(fInput) {fInput->Close(); delete fInput;}
  }

  // Loading measured fractions
  if(fbTaskStatus == AliMCSpectraWeights::TaskState::kMCSpectraObtained)
  {
    AliMCSpectraWeights::LoadMeasuredFractions();
    fbTaskStatus = AliMCSpectraWeights::TaskState::kDataFractionLoaded;
  }
  //Calculating weight factors
  if(fbTaskStatus == AliMCSpectraWeights::TaskState::kDataFractionLoaded)
  {
    if(AliMCSpectraWeights::CalculateMCWeights())
    {
      // printf("AliMCSpectraWeights::Calculating MC Weight factors succsessfull\n");
      fbTaskStatus = AliMCSpectraWeights::TaskState::kMCWeightCalculated;
    }
    // else printf("AliMCSpectraWeights::WARNING Calculating weight factors not succsessfull\n");
  }

  // printf("AliMCSpectraWeights: Status after init: %d\n",
      // AliMCSpectraWeights::GetTaskStatus());
}

void AliMCSpectraWeights::InitHistos() {
  // Initalizing histograms
  // histogram charged patricles pt:multcent:type
  const Int_t iNumberOfParticles = 6;
  Int_t nBinsTrackParticle[3] = {
    fBinsPt->GetSize() - 1,
    fBinsMultCent->GetSize() - 1, iNumberOfParticles};
  Double_t minTrackParticle[3] = {fBinsPt->GetAt(0),
    fBinsMultCent->GetAt(0), 0};
  Double_t maxTrackParticle[3] = {
    fBinsPt->GetAt(fBinsPt->GetSize() - 1),
    fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1),
    iNumberOfParticles};

  const Int_t iNumberOfParticlesDATA = 5;
  Int_t nBinsTrackParticleDATA[3] = {fBinsPt->GetSize() - 1, fBinsMultCent->GetSize() - 1,
    iNumberOfParticlesDATA};
  Double_t minTrackParticleDATA[3] = {fBinsPt->GetAt(0), fBinsMultCent->GetAt(0), 0};
  Double_t maxTrackParticleDATA[3] = {fBinsPt->GetAt(fBinsPt->GetSize() - 1), fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1), iNumberOfParticlesDATA};

  Int_t nBinsTrackParticleMC[3] = {fBinsPt->GetSize() - 1, fBinsMultCent->GetSize() - 1,
    iNumberOfParticles};
  Double_t minTrackParticleMC[3] = {fBinsPt->GetAt(0), fBinsMultCent->GetAt(0), 0};
  Double_t maxTrackParticleMC[3] = {fBinsPt->GetAt(fBinsPt->GetSize() - 1), fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1), iNumberOfParticles};

  const Int_t iNumberOfParticlesFRACTION = 5;
  Int_t nBinsTrackParticleFRACTION[3] = {fBinsPt->GetSize() - 1, fBinsMultCent->GetSize() - 1,
    iNumberOfParticlesFRACTION};
  Double_t maxTrackParticleFRACTION[3] = {
    fBinsPt->GetAt(fBinsPt->GetSize() - 1), static_cast<Double_t>(fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1)), static_cast<Double_t>(iNumberOfParticlesFRACTION)};

  fHistMCGenPrimTrackParticle =
    new THnF("fHistMCGenPrimTrackParticle",
        "histogram for charged particle composition", 3,
        nBinsTrackParticle, minTrackParticle, maxTrackParticle);
  fHistMCGenPrimTrackParticle->SetBinEdges(0, fBinsPt->GetArray());
  fHistMCGenPrimTrackParticle->SetBinEdges(1, fBinsMultCent->GetArray());
  fHistMCGenPrimTrackParticle->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistMCGenPrimTrackParticle->GetAxis(1)->SetTitle(
      "multiplicity or centrality");
  fHistMCGenPrimTrackParticle->GetAxis(2)->SetTitle("Particle type");
  fHistMCGenPrimTrackParticle->Sumw2();

  fHistDataFractions = new THnF(
      "fHistDataFractions", "DATA fractions histogram", 3,
      nBinsTrackParticleDATA, minTrackParticleDATA, maxTrackParticleDATA);
  fHistDataFractions->SetBinEdges(0, fBinsPt->GetArray());
  fHistDataFractions->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistMCGenPrimTrackParticle->GetAxis(1)->SetTitle(
      "multiplicity or centrality");
  fHistDataFractions->GetAxis(2)->SetTitle("Particle type");
  fHistDataFractions->Sumw2();

  fHistMCFractions = new THnF(
      "fHistMCFractions", "MC fractions histogram", 3,
      nBinsTrackParticleMC, minTrackParticleMC, maxTrackParticleMC);
  fHistMCFractions->SetBinEdges(0, fBinsPt->GetArray());
  fHistMCFractions->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistMCGenPrimTrackParticle->GetAxis(1)->SetTitle(
      "multiplicity or centrality");
  fHistMCFractions->GetAxis(2)->SetTitle("Particle type");
  fHistMCFractions->Sumw2();

  fHistMCWeights = new THnF(
      "fHistMCWeights", "MC weight histogram for charged particle composition",
      3, nBinsTrackParticleFRACTION, minTrackParticleDATA,
      maxTrackParticleFRACTION);
  fHistMCWeights->SetBinEdges(0, fBinsPt->GetArray());
  fHistMCWeights->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistMCGenPrimTrackParticle->GetAxis(1)->SetTitle(
      "multiplicity or centrality");
  fHistMCWeights->GetAxis(2)->SetTitle("Particle type");
  fHistMCWeights->Sumw2();

  // printf("AliMCSpectraWeights: init histos successful\n"); // works
}

void AliMCSpectraWeights::LoadMeasuredFractions() {
  // printf("AliMCSpectraWeights:: Loading %s\n", fstFilePublished.Data());
  //TFile *fMeasuredFile = AliDataFile::OpenOADB(fstFilePublished.Data());
  TFile *fMeasuredFile = TFile::Open(fstFilePublished.Data(), "OPEN");
  if (!fMeasuredFile) {
    printf(
        "AliMCSpectraWeights::Error: Could not load measured fractions in %s\n",
        fstFilePublished.Data());
    return;
  }

  TString stHistName("");

  for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
    int ibin = ipart;
    if(fstPartTypes[ipart].Contains("Pion")) ibin = static_cast<int>(AliMCSpectraWeights::ParticleType::kPion);
    else if(fstPartTypes[ipart].Contains("Proton")) ibin = static_cast<int>(AliMCSpectraWeights::ParticleType::kProtons);
    else if(fstPartTypes[ipart].Contains("Kaon")) ibin = static_cast<int>(AliMCSpectraWeights::ParticleType::kKaon);
    else if(fstPartTypes[ipart].Contains("SigmaPlus")) ibin = static_cast<int>(AliMCSpectraWeights::ParticleType::kSigmaPlus);
    else if(fstPartTypes[ipart].Contains("SigmaMinus")) ibin = static_cast<int>(AliMCSpectraWeights::ParticleType::kSigmaMinus);
    else continue;
    for (int icent = 0; icent < fNCentralities; icent++) {
      TString stFunction("Bylinkin");//TODO define functions with SysFlag
      TString stSystematic("");
      // switch (fFlag) {
      //   case /* value */:
      // }
      stHistName = Form("%s%s%s%s%s", fstCollisionSystem.Data(), fstPartTypes[ipart].Data(), fstCentralities[icent].Data(), stFunction.Data(), stSystematic.Data());
      // printf("AliMCSpectraWeights:: Loading histogram %s\n", stHistName.Data());
      TH1D *hist = (TH1D *)fMeasuredFile->Get(stHistName);
      if (!hist) {
        printf("AliMCSpectraWeights::Error: could not find %s \n",
            stHistName.Data());
        continue;
      }
      // else printf("AliMCSpectraWeights:: loading successful\n");

      // hist-> pt:multcent:type
      Double_t binEntry[3] = {0.};
      binEntry[2] = static_cast<Double_t>(ibin);// particle type
      if(fUseMultiplicity) binEntry[1] = AliMCSpectraWeights::GetMultFromCent(icent);
      else binEntry[1] =  static_cast<Double_t>(icent);
      for (int ipt = 0; ipt < hist->GetNbinsX(); ++ipt) {
        binEntry[0] = hist->GetBinCenter(ipt);
        if(binEntry[0] < 0) continue;
        // printf("AliMCSpectraWeights::DEBUG: Writing for particle %lf in %lf cent the momentum %lf the content %lf\n", binEntry[2], binEntry[1], binEntry[0], hist->GetBinContent(ipt));
        int ibinHist = fHistDataFractions->GetBin(binEntry);
        // printf("AliMCSpectraWeights::DEBUG: Writing in bin %d\n", ibinHist);
        fHistDataFractions->SetBinContent(ibinHist,
            hist->GetBinContent(ipt));
      }
      delete hist;
    }
  }
  fMeasuredFile->Close();
  delete fMeasuredFile;
  // printf("AliMCSpectraWeights: Load measured fractions finished\n");
}

Bool_t AliMCSpectraWeights::LoadFromAliMCSpectraWeight(
    AliMCSpectraWeights *obj)
{
  // printf("AliMCSpectraWeights::DEBUG: Loading MC histogram from input object\n");
  if(!obj) return kFALSE;

  if(fHistMCGenPrimTrackParticle) delete fHistMCGenPrimTrackParticle;
  fHistMCGenPrimTrackParticle = (THnF*)((THnF*)obj->GetHistMCGenPrimTrackParticles())->Clone("fHistMCGenPrimTrackParticleLoaded");
  if(!fHistMCGenPrimTrackParticle) //printf("AliMCSpectraWeights:: loading successful\n");
    {printf("AliMCSpectraWeights::ERROR: problem with loading from object\n"); return kFALSE;}
  // fHistMCGenPrimTrackParticle->SetName("fHistMCGenPrimTrackParticleLoaded");
  // if(!fHistMCGenPrimTrackParticle->GetEntries()>0) {printf("AliMCSpectraWeights::ERROR: loaded hist from object has zero entries\n"); return kFALSE;}
  // TH1D* pTProjection = (TH1D*)fHistMCGenPrimTrackParticle->Projection(0);
  // pTProjection->SetName("pTProjection");
  // printf("AliMCSpectraWeights:: mean pT: %lf\n", fHistMCGenPrimTrackParticle->Projection(0)->GetMean());

  // if(obj) delete obj;
  return kTRUE;
}

Bool_t AliMCSpectraWeights::CalcMCFractions(){
  if (!fHistMCGenPrimTrackParticle)
    return kFALSE;

  for (int icent = 0; icent < fNCentralities; icent++) {
    float currentMult = AliMCSpectraWeights::GetMultFromCent(icent);
    float nextMult = 0;
    float previousMult = 0;
    float dMultHigh = 0;
    float dMultLow = 0;
    if(icent<fBinsMultCent->GetSize()-1)
      nextMult = AliMCSpectraWeights::GetMultFromCent(icent+1);
    if(icent>0)
      previousMult = AliMCSpectraWeights::GetMultFromCent(icent-1);

    if(icent==0) dMultHigh = fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1);
    else dMultHigh = currentMult + (previousMult-currentMult)/2.;
    if(icent==fBinsMultCent->GetSize()-1) dMultLow = 0.;
    else dMultLow = currentMult + (currentMult-nextMult)/2.;

    fHistMCGenPrimTrackParticle->GetAxis(1)->SetRangeUser(dMultLow, dMultHigh);

    TH1D *h1pTMCAll = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
    if(!h1pTMCAll) {printf("AliMCSpectraWeights::ERROR could not create h1pTMCAll\n"); return kFALSE;}
    h1pTMCAll->SetName("h1pTMCAll");

    for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
      fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(ipart+1, ipart+1);
      TH1D *h1MCFraction = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
      h1MCFraction->SetName("h1MCFraction_tmp");
      h1MCFraction->Divide(h1pTMCAll);

      for (int ipt = 0; ipt < fHistMCFractions->GetAxis(0)->GetNbins(); ++ipt) {
        Double_t pt = fHistMCFractions->GetAxis(0)->GetBinCenter(ipt);
        // fHistMCFractions : pt-mult-ipart
        Double_t binEntry[3] = {pt,  currentMult, static_cast<Double_t>(ipart)};
        //Write to fHistMCFractions
        fHistMCFractions->SetBinContent(fHistMCFractions->GetBin(binEntry), h1MCFraction->GetBinContent(h1MCFraction->FindBin(pt)));
      }
      delete h1MCFraction;
    }
    fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(0, 0);
    delete h1pTMCAll;
  }
    return kTRUE;
}

Bool_t AliMCSpectraWeights::CorrectFractionsforRest()
{
  if(!fHistMCGenPrimTrackParticle || !fHistDataFractions) return kFALSE;

  for (int icent = 0; icent < fNCentralities; icent++) {
    float currentMult = AliMCSpectraWeights::GetMultFromCent(icent);
    float nextMult = 0;
    float previousMult = 0;
    float dMultHigh = 0;
    float dMultLow = 0;
    if(icent<fBinsMultCent->GetSize()-1)
      nextMult = AliMCSpectraWeights::GetMultFromCent(icent+1);
    if(icent>0)
      previousMult = AliMCSpectraWeights::GetMultFromCent(icent-1);
    if(icent==0) dMultHigh = fBinsMultCent->GetAt(fBinsMultCent->GetSize() - 1);
    else dMultHigh = currentMult + (previousMult-currentMult)/2.;
    if(icent==fBinsMultCent->GetSize()-1) dMultLow = 0.;
    else dMultLow = currentMult + (currentMult-nextMult)/2.;

    fHistMCGenPrimTrackParticle->GetAxis(1)->SetRangeUser(dMultLow, dMultHigh);
    TH1D *h1pTMCAll = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
    if(!h1pTMCAll) {printf("AliMCSpectraWeights::ERROR could not create h1pTMCAll\n"); return kFALSE;}
    h1pTMCAll->SetName("h1pTMCAll");
    fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(1, AliMCSpectraWeights::ParticleType::kRest);
    TH1D* h1RestCorrFactor = (TH1D*)fHistMCGenPrimTrackParticle->Projection(0);
    h1RestCorrFactor->SetName("h1RestCorrFactor");
    h1RestCorrFactor->Divide(h1pTMCAll);
    fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(0, 0);

    for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
      if(ipart == AliMCSpectraWeights::ParticleType::kRest) continue;
      for (int ipt = 0; ipt < fHistDataFractions->GetAxis(0)->GetNbins(); ++ipt) {
        Double_t pt = fHistDataFractions->GetAxis(0)->GetBinCenter(ipt);
        // fHistMCFractions : pt-mult-ipart
        Double_t binEntry[3] = {pt,  currentMult, static_cast<Double_t>(ipart)};
        double value = fHistDataFractions->GetBinContent(fHistDataFractions->GetBin(binEntry));
        fHistDataFractions->SetBinContent(fHistDataFractions->GetBin(binEntry), value*h1RestCorrFactor->GetBinContent(h1RestCorrFactor->FindBin(pt)));
      }
    }

    delete h1pTMCAll;
    delete h1RestCorrFactor;
  }

  return kTRUE;
}

Bool_t AliMCSpectraWeights::CalculateMCWeights(){
  if (!AliMCSpectraWeights::CalcMCFractions()) return kFALSE;
  if (!AliMCSpectraWeights::CorrectFractionsforRest()) return kFALSE;
  // correction of rest particles not measured in data fractions (see
  // AnalysisNote)

  for (int icent = 0; icent < fNCentralities; icent++) {
    for (int ipart = 0; ipart < fNPartTypes; ipart++) {
      if(ipart == AliMCSpectraWeights::ParticleType::kRest) continue;
      for (int ipt = 0; ipt < fBinsPt->GetSize(); ipt++) {
        Double_t pt = fHistMCWeights->GetAxis(0)->GetBinCenter(ipt);
        Double_t binEntry[3] = {pt,  AliMCSpectraWeights::GetMultFromCent(icent), static_cast<Double_t>(ipart)};
        // Double_t binEntryData[3] = {pt,  static_cast<Double_t>(icent), static_cast<Double_t>(ipart)};
        Double_t dFractionMC = fHistMCFractions->GetBinContent(fHistMCFractions->GetBin(binEntry));
        Double_t dFractionData = fHistDataFractions->GetBinContent(fHistDataFractions->GetBin(binEntry));
        if(dFractionMC!=0)
          fHistMCWeights->SetBinContent(fHistMCWeights->GetBin(binEntry), dFractionData/dFractionMC);
      }
    }
  }

  return kTRUE;
}

void AliMCSpectraWeights::CountEventMult()
{
  fMultOrCent=0;
  double eta = 0.8;
  if(fstCollisionSystem.Contains("pp")) eta = 0.5;
  AliStack* fMCStack = fMCEvent->Stack();
  for (size_t ipart = 0; ipart < fMCStack->GetNtrack(); ipart++) {
    TParticle *mcGenParticle = fMCStack->Particle(ipart);
    if(!fMCStack->IsPhysicalPrimary(ipart)) continue;// secondary rejection
    if(TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01) continue; // neutral rejection
    double pEta = mcGenParticle->Eta();
    if(TMath::Abs(pEta) > eta) continue; //acceptance cut
    if(mcGenParticle->Pt() < 0.05) continue; //TODO hard coded low pT cut
    ++fMultOrCent;
  }
}

Double_t
AliMCSpectraWeights::GetMCSpectraWeight(TParticle *mcGenParticle,
    Float_t eventMultiplicityOrCentrality) {
  Double_t weight = 1;
  if (!mcGenParticle->GetPDG()) {return 1;}
  if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01) {return 1;} // charge rejection
  // if (!mcGenParticle->IsPrimary()) {printf("particle not primary\n"); return 1;}
  Int_t particleType = AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
  if (particleType == AliMCSpectraWeights::ParticleType::kRest) {return 1;}
  if (fbTaskStatus == AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
    // rest particles can not be tuned
    Double_t icent = AliMCSpectraWeights::GetCentFromMult(eventMultiplicityOrCentrality);
    Double_t pt = mcGenParticle->Pt();
    if(pt<0.15) return 1;
    if(pt>=20) pt=19.9;
    Double_t binEntry[3] = {pt, AliMCSpectraWeights::GetMultFromCent(icent),  static_cast<Double_t>(particleType)};
    weight = fHistMCWeights->GetBinContent(fHistMCWeights->GetBin(binEntry));
    if(weight==0) weight=1;// printf("AliMCSpectraWeights:: got weight 0; return 1;\n");}
  // printf("AliMCSpectraWeights:: got weight %lf for pid %d at pt %lf\n", weight, mcGenParticle->GetPdgCode(), binEntry[0]);
  } //else printf("Status not right\n" );
  return weight;
}

Double_t
AliMCSpectraWeights::GetMCSpectraWeight(TParticle *mcGenParticle,
    AliMCEvent* mcEvent) {
    if(mcEvent!=fMCEvent)
    {
      fMCEvent=mcEvent;
      AliMCSpectraWeights::CountEventMult();
    }
    return AliMCSpectraWeights::GetMCSpectraWeight(mcGenParticle, fMultOrCent);
}

void AliMCSpectraWeights::FillMCSpectra(AliMCEvent* mcEvent) {
  if (fbTaskStatus >= AliMCSpectraWeights::TaskState::kMCSpectraObtained) return;

  if(mcEvent!=fMCEvent)
  {
    fMCEvent=mcEvent;
    AliMCSpectraWeights::CountEventMult();
  }

  AliStack* MCStack = fMCEvent->Stack();
  if (!MCStack) {printf("AliMCSpectraWeights::ERROR: fMCStack not available\n"); return;}

  for (Int_t iParticle = 0; iParticle < MCStack->GetNtrack(); iParticle++){
    TParticle *mcGenParticle = MCStack->Particle(iParticle);
    if(!mcGenParticle) {printf("AliMCSpectraWeights::ERROR: mcGenParticle  not available\n"); continue;}
    if (!mcGenParticle->GetPDG()) continue;
    if(!MCStack->IsPhysicalPrimary(iParticle)) continue;
    float partEta = mcGenParticle->Eta();
    if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01) continue;
    if(partEta > 0.8 || partEta < -0.8) continue; // apply same acceptance as in published spectra
    Int_t particleType = AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
    if(particleType<0) continue;
    Double_t binEntry[3] = {static_cast<Double_t>(mcGenParticle->Pt()),
      static_cast<Double_t>(fMultOrCent),
      static_cast<Double_t>(particleType)};
    fHistMCGenPrimTrackParticle->Fill(binEntry);
  }
}

/// Function to return Particle ID for Histograms
Int_t AliMCSpectraWeights::IdentifyMCParticle(TParticle *mcParticle) {
  // if(!mcParticle->GetPDG()) return -1;
  Int_t ipdg = TMath::Abs(
    mcParticle->GetPdgCode()); // Abs() because antiparticles are negaitve...
  if (ipdg == 211)
    return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kPion);
  if (ipdg == 321)
    return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kKaon);
  if (ipdg == 2212)
    return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kProtons);
  if (ipdg == 3222)
    return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kSigmaPlus);
  if (ipdg == 3112)
    return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kSigmaMinus);
  // if(ipdg==3334) return AliMCSpectraWeights::ParticleType::kOmegaMinus;
  // if(ipdg==3312) return AliMCSpectraWeights::ParticleType::kXiMinus;
  // if(ipdg==11) return AliMCSpectraWeights::ParticleType::kElectron;
  // if(ipdg==13) return AliMCSpectraWeights::ParticleType::kMuon;
  //printf("AliMCSpectraWeights:: pdf code of rest particle %d\n", ipdg);
  return static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kRest);
}

Double_t AliMCSpectraWeights::GetMultFromCent(int CentBin){
  if(fstCollisionSystem.Contains("pp"))
  {
    switch (CentBin) {
      case 0: return 21.3;
      case 1: return 16.5;
      case 2: return 13.5;
      case 3: return 11.5;
      case 4: return 10.1;
      case 5: return 8.45;
      case 6: return 6.72;
      case 7: return 5.4;
      case 8: return 3.9;
      case 9: return 2.26;
      default: return 0;
    }
  }
  else if(fstCollisionSystem.Contains("ppb"))//TODO include other systems
  {

  }
  else if(fstCollisionSystem.Contains("pbpb"))
  {

  }
  else if(fstCollisionSystem.Contains("xexe"))
  {

  }

  return -1;
}

Double_t AliMCSpectraWeights::GetCentFromMult(double dMult)
{
  if(fstCollisionSystem.Contains("pp"))
  {
    if(dMult > 18) return 0;
    if(dMult > 14.5) return 1;
    if(dMult > 12) return 2;
    if(dMult > 10.9) return 3;
    if(dMult > 9) return 4;
    if(dMult > 7) return 5;
    if(dMult > 5.7) return 6;
    if(dMult > 4.3) return 7;
    if(dMult > 3) return 8;
    if(dMult > 0) return 9;
  }
  else if(fstCollisionSystem.Contains("ppb"))//TODO include other systems
  {

  }
  else if(fstCollisionSystem.Contains("pbpb"))
  {

  }
  else if(fstCollisionSystem.Contains("xexe"))
  {

  }

  return -1;
}
