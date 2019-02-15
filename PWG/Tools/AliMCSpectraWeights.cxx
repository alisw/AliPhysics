/*!
  \file AliMCSpectraWeights.cxx
  \brief "Description"
  \author Patrick Huhn
  \date 17/10/2018
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
  fHistMCGenPrimTrackParticle(0), fHistDataFractions(0), fHistMCWeights(0),
  fBinsPt(0), fBinsMultCent(0), fstCollisionSystem("pp"),
  fNPartTypes(6), fPartTypes(0), fstFileMCSpectra(""), fstFilePublished(""),
  fstSavedObjName("fMCSpectraWeights"),fstSavedListName("dNdPt_test"),
  fUseMultiplicity(kFALSE), fbTaskStatus(0) {

    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
  }

// not-default contructor; way to go
AliMCSpectraWeights::AliMCSpectraWeights(const char *collisionSystem, const char* name)
  : TNamed(name, name), fHistMCGenPrimTrackParticle(0), fHistDataFractions(0),
  fHistMCWeights(0), fBinsPt(0), fBinsMultCent(0),
  fstCollisionSystem("pp"), fNPartTypes(6), fPartTypes(0),
  fstFileMCSpectra(""), fstFilePublished(""),
  fstSavedObjName("fMCSpectraWeights"),fstSavedListName("dNdPt_test"),
  fUseMultiplicity(kFALSE), fbTaskStatus(0) {

    fstCollisionSystem = collisionSystem;
    fstCollisionSystem.ToLower();

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
    if (fstCollisionSystem.Contains("pp") || fstCollisionSystem.Contains("p-p")) {
      Double_t mulBins[101] = {0};
      for (int i = 0; i < 101; ++i) {
        mulBins[i] = i;
      }
      fBinsMultCent = new TArrayD(101, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("ppb") ||
        fstCollisionSystem.Contains("p-pb")) {
      Double_t mulBins[301] = {0};
      for (int i = 0; i < 301; ++i) {
        mulBins[i] = i;
      }
      fBinsMultCent = new TArrayD(301, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("pbpb") ||
        fstCollisionSystem.Contains("pb-pb")) {
      Double_t mulBins[201] = {0};
      for (int i = 0; i < 201; ++i) {
        mulBins[i] = i * 25;
      }
      fBinsMultCent = new TArrayD(201, mulBins);
      // if(mulBins) delete mulBins;
    } else if (fstCollisionSystem.Contains("xexe") ||
        fstCollisionSystem.Contains("xe-xe")) {
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

    fPartTypes = new TString[6];
    fPartTypes[AliMCSpectraWeights::ParticleType::kPion] = "Pion";
    fPartTypes[AliMCSpectraWeights::ParticleType::kProtons] = "Proton";
    fPartTypes[AliMCSpectraWeights::ParticleType::kKaon] = "Kaon";
    fPartTypes[AliMCSpectraWeights::ParticleType::kSigmaMinus] = "SigmaMinus";
    fPartTypes[AliMCSpectraWeights::ParticleType::kSigmaPlus] = "SigmaPlus";
    fPartTypes[AliMCSpectraWeights::ParticleType::kRest] = "Rest";

    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
    fstFilePublished = "$ALICE_PHYSICS/OADB/PWGPP/data/AllPublishedFractions.root";
  }

AliMCSpectraWeights::~AliMCSpectraWeights() {
  // if (fHistMCGenPrimTrackParticle)
    // delete fHistMCGenPrimTrackParticle;
  if (fHistDataFractions)
    delete fHistDataFractions;
  if (fHistMCWeights)
    delete fHistMCWeights;
  if (fBinsPt)
    delete fBinsPt;
  if (fBinsMultCent)
    delete fBinsMultCent;
  if (fPartTypes)
    delete fPartTypes;
}

/*!
  \brief "Description"
  \param "Param description"
  \pre "Pre-conditions"
  \post "Post-conditions"
  \return "Return of the function"
  */
void AliMCSpectraWeights::Init() {
  // Histograms
  AliMCSpectraWeights::InitHistos();

  if (fstFileMCSpectra.Length() > 5) // *.root
  {
    TFile* fInput = TFile::Open(fstFileMCSpectra.Data());
    if(fInput)
    {
      if (fInput->GetNkeys() != 1)
      {if(!fstSavedListName.Contains("dNdPt_test"))
        printf("AliMCSpectraWeights::WARNING: more than 1 list in the streamed file; please specify; using 1st list;\n\n");}
      else fstSavedListName  = fInput->GetListOfKeys()->At(0)->GetName();
      printf("AliMCSpectraWeights:: Loading %s from list %s\n", fstSavedObjName.Data(), fstSavedListName.Data());
      TList *listMC = (TList*)fInput->Get(fstSavedListName);
      if(!listMC) {printf("AliMCSpectraWeights::ERROR: could not load list in streamed file\n");}
      else{
        AliMCSpectraWeights* inWeights = (AliMCSpectraWeights*)listMC->FindObject(fstSavedObjName.Data());
        if(AliMCSpectraWeights::LoadFromAliMCSpectraWeight(inWeights))
          fbTaskStatus = AliMCSpectraWeights::TaskState::kMCSpectraObtained;
        else fHistMCGenPrimTrackParticle = (THnF*)listMC->FindObject("fHistMCGenPrimTrackParticle");
        if(!fHistMCGenPrimTrackParticle) {printf("AliMCSpectraWeights::WARNING: Couln't get fHistMCGenPrimTrackParticle\n" ); return;}
        // fHistMCGenPrimTrackParticle->SetDirectory(0);
        if(fHistMCGenPrimTrackParticle->GetEntries()>0) fbTaskStatus = AliMCSpectraWeights::TaskState::kMCSpectraObtained;
        if(inWeights) delete inWeights;
      }
    }
    else printf("AliMCSpectraWeights::WARNING: %s can not be loaded\n ", fstFileMCSpectra.Data());
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
      printf("AliMCSpectraWeights::Calculating MC Weight factors succsessfull\n");
      fbTaskStatus = AliMCSpectraWeights::TaskState::kMCWeightCalculated;
    }
    else printf("AliMCSpectraWeights::WARNING Calculating weight factors not succsessfull\n");
  }

  // calculate weights
  //if (fbTaskStatus == AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
  //  if (AliMCSpectraWeights::CalculateMCWeights())
  //    fbTaskStatus = AliMCSpectraWeights::TaskState::kMCWeightCalculated;
  //}

  printf("AliMCSpectraWeights: Status after init: %d\n",
      AliMCSpectraWeights::GetTaskStatus());
}

void AliMCSpectraWeights::InitHistos() {
  // Initalizing histograms
  // TODO implement multiplicity or centrality dependence
  // histogram charged patricles pt:eta:multcent:type
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
  // TString binNameTrackParticle[iNumberOfParticles] = {
    // "Pion", "Kaon", "Proton", "SigmaPlus", "SigmaMinus", "Rest"};

  const Int_t iNumberOfParticlesDATA = 5;
  Int_t nBinsTrackParticleDATA[2] = {fBinsPt->GetSize() - 1,
    iNumberOfParticlesDATA};
  Double_t minTrackParticleDATA[2] = {fBinsPt->GetAt(0), 0};
  Double_t maxTrackParticleDATA[2] = {fBinsPt->GetAt(fBinsPt->GetSize() - 1),
    iNumberOfParticlesDATA};
  // TString binNameTrackParticleDATA[iNumberOfParticlesDATA] = {
    // "Pion", "Kaon", "Proton", "SigmaPlus", "SigmaMinus"};

  const Int_t iNumberOfParticlesFRACTION = 5;
  Int_t nBinsTrackParticleFRACTION[2] = {fBinsPt->GetSize() - 1,
    iNumberOfParticlesFRACTION};
  Double_t maxTrackParticleFRACTION[2] = {
    fBinsPt->GetAt(fBinsPt->GetSize() - 1), iNumberOfParticlesFRACTION};
  // TString binNameTrackParticleFRACTION[iNumberOfParticlesFRACTION] = {
    // "Pion", "Kaon", "Proton", "SigmaPlus", "SigmaMinus"};

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

  // for (Int_t ii = 1; ii <= fHistMCGenPrimTrackParticle->GetAxis(2)->GetNbins();
  //     ii++) {
  //   fHistMCGenPrimTrackParticle->GetAxis(2)->SetBinLabel(
  //       ii, binNameTrackParticle[ii - 1].Data());
  // }

  fHistDataFractions = new THnF(
      "fHistDataFractions", "DATA fractions histogram", 2,
      nBinsTrackParticleDATA, minTrackParticleDATA, maxTrackParticleDATA);
  fHistDataFractions->SetBinEdges(0, fBinsPt->GetArray());
  fHistDataFractions->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistDataFractions->GetAxis(1)->SetTitle("Particle type");
  fHistDataFractions->Sumw2();

  // for (Int_t ii = 1; ii <= fHistDataFractions->GetAxis(1)->GetNbins(); ii++) {
  //   fHistDataFractions->GetAxis(1)->SetBinLabel(
  //       ii, binNameTrackParticleDATA[ii - 1].Data());
  // }

  fHistMCWeights = new THnF(
      "fHistMCWeights", "MC weight histogram for charged particle composition",
      2, nBinsTrackParticleFRACTION, minTrackParticleDATA,
      maxTrackParticleFRACTION);
  fHistMCWeights->SetBinEdges(0, fBinsPt->GetArray());
  fHistMCWeights->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistMCWeights->GetAxis(1)->SetTitle("Particle type");
  fHistMCWeights->Sumw2();

  // for (Int_t ii = 1; ii <= fHistMCWeights->GetAxis(1)->GetNbins(); ii++) {
  //   fHistMCWeights->GetAxis(1)->SetBinLabel(
  //       ii, binNameTrackParticleFRACTION[ii - 1].Data());
  // }

  printf("AliMCSpectraWeights: init histos successful\n"); // works
}

void AliMCSpectraWeights::LoadMeasuredFractions() {//TODO for some reasons only pTBin=1 filled
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
    if(fPartTypes[ipart].Contains("Pion")) ibin = AliMCSpectraWeights::ParticleType::kPion;
    else if(fPartTypes[ipart].Contains("Proton")) ibin = AliMCSpectraWeights::ParticleType::kProtons;
    else if(fPartTypes[ipart].Contains("Kaon")) ibin = AliMCSpectraWeights::ParticleType::kKaon;
    else if(fPartTypes[ipart].Contains("SigmaPlus")) ibin = AliMCSpectraWeights::ParticleType::kSigmaPlus;
    else if(fPartTypes[ipart].Contains("SigmaMinus")) ibin = AliMCSpectraWeights::ParticleType::kSigmaMinus;
    else continue;
    stHistName = Form("%s%sStatBylinkin", fstCollisionSystem.Data(),
        fPartTypes[ipart].Data()); // TODO name schema
    // printf("AliMCSpectraWeights:: Loading histogram %s\n", stHistName.Data());
    TH1D *hist = (TH1D *)fMeasuredFile->Get(stHistName);
    if (!hist) {
      printf("AliMCSpectraWeights::Error: could not find %s \n",
          stHistName.Data());
      continue;
    }
    else printf("AliMCSpectraWeights:: loading successful\n");
    Double_t binEntry[2] = {0.};
    binEntry[1] = static_cast<Double_t>(ibin);
    for (int ipt = 0; ipt < hist->GetNbinsX(); ++ipt) {
      //fHistMCWeights: pT:PartType
      binEntry[0] = hist->GetBinCenter(ipt);
      if(binEntry[0] < 0) continue;
      // printf("AliMCSpectraWeights::DEBUG: Writing for particle %lf the momentum %lf the content %lf\n", binEntry[1], binEntry[0], hist->GetBinContent(ipt));

      fHistDataFractions->SetBinContent(fHistDataFractions->GetBin(binEntry),
          hist->GetBinContent(ipt));
    }
    delete hist;
  }
  fMeasuredFile->Close();
  delete fMeasuredFile;
  // printf("AliMCSpectraWeights: Load measured fractions finished\n");
}

Bool_t AliMCSpectraWeights::LoadFromAliMCSpectraWeight(
    AliMCSpectraWeights *obj)
{
  printf("AliMCSpectraWeights::DEBUG: Loading MC histogram from input object\n");
  if(!obj) return kFALSE;

  if(fHistMCGenPrimTrackParticle) delete fHistMCGenPrimTrackParticle;
  fHistMCGenPrimTrackParticle = (THnF*)((THnF*)obj->GetHistMCGenPrimTrackParticles())->Clone("fHistMCGenPrimTrackParticleLoaded");
  if(fHistMCGenPrimTrackParticle) printf("AliMCSpectraWeights:: loading successful\n");
  else {printf("AliMCSpectraWeights::ERROR: problem with loading from object\n"); return kFALSE;}
  // fHistMCGenPrimTrackParticle->SetName("fHistMCGenPrimTrackParticleLoaded");
  // if(!fHistMCGenPrimTrackParticle->GetEntries()>0) {printf("AliMCSpectraWeights::ERROR: loaded hist from object has zero entries\n"); return kFALSE;}
  // TH1D* pTProjection = (TH1D*)fHistMCGenPrimTrackParticle->Projection(0);
  // pTProjection->SetName("pTProjection");
  // printf("AliMCSpectraWeights:: mean pT: %lf\n", fHistMCGenPrimTrackParticle->Projection(0)->GetMean());

  // if(obj) delete obj;
  return kTRUE;
}

Bool_t AliMCSpectraWeights::CalculateMCWeights(){//TODO make constant ratio @large pT
  if (!fHistMCGenPrimTrackParticle || !fHistDataFractions)
    return kFALSE;
  else printf("AliMCSpectraWeights:: start calculating weight factors\n");
  // printf("AliMCSpectraWeights:: mean pT: %lf\n", fHistMCGenPrimTrackParticle->Projection(0)->GetMean());
  // correction of rest particles not measured in data fractions (see
  // AnalysisNote)
  Int_t kRestPosition = static_cast<Int_t>(AliMCSpectraWeights::ParticleType::kRest);
  // for (int i = 0; i < fNPartTypes; i++) {
    // if(fPartTypes[i].Contains("rest") || fPartTypes[i].Contains("Rest")) {kRestPosition=i; break;}
  // }
  TH1D *h1pTMCAll = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
  if(!h1pTMCAll) {printf("AliMCSpectraWeights::ERROR could not create h1pTMCAll\n"); return kFALSE;}
  h1pTMCAll->SetName("h1pTMCAll");
  fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(0, kRestPosition - 1);
  // printf("AliMCSpectraWeights:: h1pTMCAll created; now h1pTMCRestFraction\n");
  TH1D *h1pTMCRestFraction = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
  if(!h1pTMCRestFraction) {printf("AliMCSpectraWeights::ERROR could not create h1pTMCRestFraction\n"); return kFALSE;}
  h1pTMCRestFraction->SetName("h1pTMCRestFraction");
  h1pTMCRestFraction->Divide(h1pTMCAll);
  // printf("AliMCSpectraWeights:: Calculated histo for rest correction\n");
  int inPart = kRestPosition-1;//fHistMCWeights->GetAxis(1)->GetNbins();
  for (int ipart = 0; ipart < inPart; ++ipart) {
    // setting range
    // printf("AliMCSpectraWeights:: Now at %s\n", fPartTypes[ipart].Data());
    fHistDataFractions->GetAxis(1)->SetRange(ipart, ipart);
    fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(ipart, ipart);

    // printf("AliMCSpectraWeights:: Before projection\n");
    TH1D *h1DataFraction = (TH1D *)fHistDataFractions->Projection(0);
    h1DataFraction->SetName("h1DataFraction_tmp");
    TH1D *h1MCFraction = (TH1D *)fHistMCGenPrimTrackParticle->Projection(0);
    h1MCFraction->SetName("h1MCFraction_tmp");
    h1MCFraction->Divide(h1pTMCAll);
    // printf("AliMCSpectraWeights:: Before ratio of fractions\n");
    h1DataFraction->Multiply(h1pTMCRestFraction);// TODO here might be an issue due to binning
    h1DataFraction->Divide(h1MCFraction);
      // printf("AliMCSpectraWeights:: ratio succsess; filling fHistMCWeights\n");
    for (int ipt = 0; ipt < fHistMCWeights->GetAxis(0)->GetNbins(); ++ipt) {
      Double_t pt = fHistMCWeights->GetAxis(0)->GetBinCenter(ipt);
      Double_t binEntry[2] = {pt, static_cast<Double_t>(ipart)};
      // printf("AliMCSpectraWeights:: Filling for %s and pT %lf\n", fPartTypes[ipart].Data(), binEntry[0]);
      fHistMCWeights->SetBinContent(
          fHistMCWeights->GetBin(binEntry),
          h1DataFraction->GetBinContent(h1DataFraction->FindBin(pt)));
    }
    // printf("AliMCSpectraWeights:: deleting tmp hists\n");
    delete h1DataFraction;
    delete h1MCFraction;
  }
  fHistDataFractions->GetAxis(1)->SetRange(0, kRestPosition-1);
  fHistMCGenPrimTrackParticle->GetAxis(2)->SetRange(0, kRestPosition);
  // printf("AliMCSpectraWeights:: start calculating weight factors...done\n");
  return kTRUE;
}

Double_t
AliMCSpectraWeights::GetMCSpectraWeight(TParticle *mcGenParticle,
    Float_t eventMultiplicityOrCentrality) {
  Double_t weight = 1;
  if (!mcGenParticle->GetPDG())
    return 1;
  // if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01)
  //   return 1; // charge rejection
  Int_t particleType = AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
  if (fbTaskStatus >= AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
    // rest particles can not be tuned
    if (particleType > AliMCSpectraWeights::ParticleType::kSigmaMinus)
      return 1;

    Double_t binEntry[2] = {0.};
    binEntry[0] = mcGenParticle->Pt();
    if(binEntry[0] < 0.15) return 1;
    binEntry[1] = static_cast<Double_t>(particleType);
    weight = fHistMCWeights->GetBinContent(fHistMCWeights->GetBin(binEntry));
    if(weight==0) weight=1;// printf("AliMCSpectraWeights:: got weight 0; return 1;\n");}
    // printf("AliMCSpectraWeights:: got weight %lf for pid %d at pt %lf\n", weight, mcGenParticle->GetPdgCode(), binEntry[0]);
  } else {
    printf("AliMCSpectraWeights: MC spectra not obtained, yet.\n");
  }
  return weight;
}

void AliMCSpectraWeights::FillMCSpectra(AliMCEvent* mcEvent,
    Float_t eventMultiplicityOrCentrality) {
  if (fbTaskStatus == AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
    printf("AliMCSpectraWeights:: MC spectra already obtained; step skipped\n");
    return;
  }

  AliStack* fMCStack = mcEvent->Stack();
  if (!fMCStack) {printf("AliMCSpectraWeights::ERROR: fMCStack not available\n"); return;}

  for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack(); iParticle++){
    TParticle *mcGenParticle = fMCStack->Particle(iParticle);
    if(!mcGenParticle) {printf("AliMCSpectraWeights::ERROR: mcGenParticle  not available\n"); continue;}
    if (!mcGenParticle->GetPDG()) continue;
    if(!fMCStack->IsPhysicalPrimary(iParticle)) continue;
    float partEta = mcGenParticle->Eta();
    if(partEta > 0.8 || partEta < -0.8) continue; // apply same acceptance as in published spectra
    Int_t particleType = AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
    if(particleType<0) continue;
    Double_t binEntry[3] = {static_cast<Double_t>(mcGenParticle->Pt()),
      static_cast<Double_t>(eventMultiplicityOrCentrality),
      static_cast<Double_t>(particleType)};
    fHistMCGenPrimTrackParticle->Fill(binEntry);
  }
}

/// Function to return Particle ID for Histograms
Int_t AliMCSpectraWeights::IdentifyMCParticle(TParticle *mcParticle) {
  // if(!mcParticle->GetPDG()) return -1;
  if (TMath::Abs(mcParticle->GetPDG()->Charge()) < 0.01)
    return -1; // charge rejection;
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

void AliMCSpectraWeights::GetMCTrackHist(THnF *hist) {//TODO fix this function
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
  TString binNameTrackParticle[iNumberOfParticles] = {
    "Pion", "Kaon", "Proton", "SigmaPlus", "SigmaMinus", "Rest"};

  hist = new THnF("fHistMCGenPrimTrackParticleOutput",
      "histogram for charged particle composition", 3,
      nBinsTrackParticle, minTrackParticle, maxTrackParticle);
  hist->SetBinEdges(0, fBinsPt->GetArray());
  hist->SetBinEdges(1, fBinsMultCent->GetArray());
  hist->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  hist->GetAxis(1)->SetTitle("multiplicity or centrality");
  hist->GetAxis(2)->SetTitle("Particle type");
  hist->Sumw2();
  // for (Int_t ii = 1; ii <= hist->GetAxis(3)->GetNbins(); ii++) {
  //   hist->GetAxis(3)->SetBinLabel(ii, binNameTrackParticle[ii - 1].Data());
  // }

  // hist =
  // (THnF*)fHistMCGenPrimTrackParticle->Clone("fHistMCGenPrimTrackParticleOutput");
  Int_t binEntry[3] = {0};
  for (int ibin = 0; ibin < fHistMCGenPrimTrackParticle->GetNbins(); ++ibin) {
    Double_t content =
      fHistMCGenPrimTrackParticle->GetBinContent(ibin, binEntry);
    // Double_t error = fHistMCGenPrimTrackParticle->GetBinError(ibin);

    // Double_t pt =
    // fHistMCGenPrimTrackParticle->GetAxis(0)->GetBinCenter(binEntry[0]);
    // Double_t eta =
    // fHistMCGenPrimTrackParticle->GetAxis(1)->GetBinCenter(binEntry[1]);
    // Double_t cent =
    // fHistMCGenPrimTrackParticle->GetAxis(2)->GetBinCenter(binEntry[2]);
    // Double_t partType =
    // fHistMCGenPrimTrackParticle->GetAxis(3)->GetBinCenter(binEntry[3]);

    hist->SetBinContent(binEntry, content);
    // hist->SetBinError(binEntry, error);
  }
}
