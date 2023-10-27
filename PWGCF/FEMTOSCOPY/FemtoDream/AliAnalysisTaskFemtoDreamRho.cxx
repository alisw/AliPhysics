/*
 * AliAnalysisTaskFemtoDreamRho.cxx
 *
 *  Created on: 19 Jul 2023
 *      Author: M. Korwieser
 */

#include "AliAnalysisTaskFemtoDreamRho.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"

ClassImp(AliAnalysisTaskFemtoDreamRho)
    AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho()
    : AliAnalysisTaskSE(),
      fTrigger(AliVEvent::kINT7),
      fIsMC(false),
      fDoCleaning(false),
      fResonanceHistograms(),
      fFlagHistogram(),
      fpTerrorHistogram(),
      fpTCorrerrorHistogram(),
      ptHist_RhoMCTrue(),
      fHistogramPDG(),
      fOutput(),
      fEvent(),
      fTrack(),
      fRhoParticle(),
      fEventCuts(),
      fPosPionCuts(),
      fNegPionCuts(),
      fRhoCuts(),
      fPosProtonCuts(),
      fNegProtonCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize()
{
}

AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho(const char *name,
                                                           bool isMC, bool doCleaning)
    : AliAnalysisTaskSE(name),
      fTrigger(AliVEvent::kINT7),
      fIsMC(isMC),
      fDoCleaning(doCleaning),
      fResonanceHistograms(),
      fFlagHistogram(),
      fpTerrorHistogram(),
      fpTCorrerrorHistogram(),
      ptHist_RhoMCTrue(),
      fHistogramPDG(),
      fOutput(),
      fEvent(),
      fTrack(),
      fRhoParticle(),
      fEventCuts(),
      fPosPionCuts(),
      fNegPionCuts(),
      fRhoCuts(),
      fPosProtonCuts(),
      fNegProtonCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize(2000)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskFemtoDreamRho::~AliAnalysisTaskFemtoDreamRho() {}

std::map<TString, std::pair<TH1F *, TH2F *>> AliAnalysisTaskFemtoDreamRho::CreateResonanceHistograms(const std::vector<TString> &resonanceList)
{
  std::map<TString, std::pair<TH1F *, TH2F *>> histograms;

  for (const auto &resonanceName : resonanceList)
  {
    // Create TH1F for pT spectrum
    TString ptHistName = TString::Format("histPtSpectrum_%s", resonanceName.Data());
    TString ptHistTitle = TString::Format("%s pT Spectrum", resonanceName.Data());
    TH1F *ptHist = new TH1F(ptHistName, ptHistTitle, 100, 0.0, 10.0);
    ptHist->GetXaxis()->SetTitle("pT (GeV/c)");

    // Create TH2F for invariant mass vs pT
    TString massPtHistName = TString::Format("histInvariantMassPt_%s", resonanceName.Data());
    TString massPtHistTitle = TString::Format("%s Invariant Mass vs. pT", resonanceName.Data());
    TH2F *massPtHist = new TH2F(massPtHistName, massPtHistTitle, 100, 0.0, 10.0, 100, 0.0, 5.0);
    massPtHist->GetXaxis()->SetTitle("pT (GeV/c)");
    massPtHist->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    histograms[resonanceName] = std::make_pair(ptHist, massPtHist);
  }

  return histograms;
}

void AliAnalysisTaskFemtoDreamRho::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, fTrigger);
  fOutput->Add(fEvent->GetEvtCutList());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fRhoParticle = new AliFemtoDreamv0();
  fRhoParticle->SetPDGCode(fRhoCuts->GetPDGv0());
  fRhoParticle->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterPos(
      fRhoCuts->GetPDGPosDaug()); // order +sign doesnt play a role
  fRhoParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterNeg(
      fRhoCuts->GetPDGNegDaug()); // only used for MC Matching
  fRhoParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack *[fTrackBufferSize];

  if (!fEventCuts)
  {
    AliFatal("Event Cuts not set!");
  }
  fEventCuts->InitQA();
  fOutput->Add(fEventCuts->GetHistList());

  if (!fPosPionCuts)
  {
    AliFatal("Track Cuts for positive pion not set!");
  }
  fPosPionCuts->Init();
  fPosPionCuts->SetName("PosPion");
  fOutput->Add(fPosPionCuts->GetQAHists());

  if (fPosPionCuts->GetIsMonteCarlo())
  {
    fPosPionCuts->SetMCName("MCPosPion");
    fOutput->Add(fPosPionCuts->GetMCQAHists());
  }

  if (!fNegPionCuts)
  {
    AliFatal("Track Cuts for negative pion not set!");
  }
  fNegPionCuts->Init();
  fNegPionCuts->SetName("NegPion");
  fOutput->Add(fNegPionCuts->GetQAHists());
  if (fNegPionCuts->GetIsMonteCarlo())
  {
    fNegPionCuts->SetMCName("MCNegPion");
    fOutput->Add(fNegPionCuts->GetMCQAHists());
  }

  if (!fRhoCuts)
  {
    AliFatal("Cuts for the Rho not set!");
  }
  fRhoCuts->Init();
  fRhoCuts->SetName("Rho");
  fOutput->Add(fRhoCuts->GetQAHists());
  if (fRhoCuts->GetIsMonteCarlo())
  {
    fRhoCuts->SetMCName("MCRho");
    fOutput->Add(fRhoCuts->GetMCQAHists());
  }

  /*if (!fRhoParticleMCTrueCuts)
  {
    AliFatal("Cuts for the RhoMCTrue not set!");
  }
  fRhoParticleMCTrueCuts->Init();
  fRhoParticleMCTrueCuts->SetName("RhoMCTrue");
  fOutput->Add(fRhoParticleMCTrueCuts->GetQAHists());
  if (fRhoParticleMCTrueCuts->GetIsMonteCarlo())
  {
    fRhoParticleMCTrueCuts->SetMCName("MCRhoMCTrue");
    fOutput->Add(fRhoParticleMCTrueCuts->GetMCQAHists());
  }*/

  if (!fPosProtonCuts)
  {
    AliFatal("Track Cuts for Proton not set!");
  }
  fPosProtonCuts->Init();
  fPosProtonCuts->SetName("Proton");
  fOutput->Add(fPosProtonCuts->GetQAHists());
  if (fPosProtonCuts->GetIsMonteCarlo())
  {
    fPosProtonCuts->SetMCName("MCProton");
    fOutput->Add(fPosProtonCuts->GetMCQAHists());
  }

  if (!fNegProtonCuts)
  {
    AliFatal("Track Cuts for AntiProton not set!");
  }
  fNegProtonCuts->Init();
  fNegProtonCuts->SetName("AntiProton");
  fOutput->Add(fNegProtonCuts->GetQAHists());
  if (fNegProtonCuts->GetIsMonteCarlo())
  {
    fNegProtonCuts->SetMCName("MCAntiProton");
    fOutput->Add(fNegProtonCuts->GetMCQAHists());
  }

  fPairCleaner =
      new AliFemtoDreamPairCleaner(3, 0, fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());

  // Define the resonance list (modify this list as per your requirements)
  std::vector<TString> resonanceList = {"rho", "omega", "kshort", "f0", "f2"};

  // Create histograms for each resonance and store them in the map
  fResonanceHistograms = CreateResonanceHistograms(resonanceList);
  auto *fHistListResonance = new TList();
  fHistListResonance->SetName("V0MCinvestigation");
  fHistListResonance->SetOwner();

  // Add histograms to the output list
  for (auto &pair : fResonanceHistograms)
  {
    fHistListResonance->Add(pair.second.first);  // TH1F for pT spectrum
    fHistListResonance->Add(pair.second.second); // TH2F for invariant mass vs pT
  }
  fOutput->Add(fHistListResonance);

  // Create histograms for RhoMCTrue particles
  auto *fHistListRhoMCTrue = new TList();
  fHistListRhoMCTrue->SetName("RhoMCTrue");
  fHistListRhoMCTrue->SetOwner();
  // String
  TString NameIng = "RhoMCTrue";

  // Create TH1F for pT spectrum
  TString ptHistName = TString::Format("histPtSpectrum_%s", NameIng.Data());
  TString ptHistTitle = TString::Format("%s pT Spectrum", NameIng.Data());
  TH1F *ptHist_RhoMCTrue = new TH1F(ptHistName, ptHistTitle, 1000, 0.0, 3.0);
  ptHist_RhoMCTrue->GetXaxis()->SetTitle("pT (GeV/c)");

  // Create TH2F for invariant mass vs pT
  TString massPtHistName = TString::Format("histInvariantMassPt_%s", NameIng.Data());
  TString massPtHistTitle = TString::Format("%s Invariant Mass vs. pT", NameIng.Data());
  TH2F *massPtHist_RhoMCTrue = new TH2F(massPtHistName, massPtHistTitle, 100, 0.0, 10.0, 100, 0.0, 5.0);
  massPtHist_RhoMCTrue->GetXaxis()->SetTitle("pT (GeV/c)");
  massPtHist_RhoMCTrue->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

  fHistListRhoMCTrue->Add(ptHist_RhoMCTrue); // TH1F for pT spectrum
  fHistListRhoMCTrue->Add(massPtHist_RhoMCTrue);

  fOutput->Add(fHistListRhoMCTrue); // Add the histogram to your output list

  // Create the histogram to track flags
  fFlagHistogram = new TH1F("FlagHistogram", "Flag Changes in if(fIsMC && CheckV0s)", 7, 0.5, 7.5);
  fFlagHistogram->GetXaxis()->SetBinLabel(1, "useThisV0");
  fFlagHistogram->GetXaxis()->SetBinLabel(2, "sameMother");
  fFlagHistogram->GetXaxis()->SetBinLabel(3, "diffMother");
  fFlagHistogram->GetXaxis()->SetBinLabel(4, "diffMotherID");
  fFlagHistogram->GetXaxis()->SetBinLabel(5, "dauArePhysPrim");
  fFlagHistogram->GetXaxis()->SetBinLabel(6, "dauAreNotPions");
  fFlagHistogram->GetXaxis()->SetBinLabel(7, "No investigation");

  fOutput->Add(fFlagHistogram); // Add the histogram to your output list

  // Create the histogram to track rho daughter momentum
  fpTerrorHistogram = new TH1F("fpTerrorHistogram", "fpTerrorHistogram", 4000, 0., 4.);
  fpTCorrerrorHistogram = new TH2F("fpTCorrerrorHistogram", "fpTCorrerrorHistogram", 4000, 0., 4., 4000, 0., 4.);

  fOutput->Add(fpTerrorHistogram);     // Add the histogram to your output list
  fOutput->Add(fpTCorrerrorHistogram); // Add the histogram to your output list

  // Create the histogram to track flags
  fHistogramPDG = new TH1F("fHistogramPDG", "fHistogramPDG", 100000, 0, 100000);
  fHistogramPDG->GetXaxis()->SetBinLabel(211, "pions");
  fHistogramPDG->GetXaxis()->SetBinLabel(2212, "protons");

  fOutput->Add(fHistogramPDG); // Add the histogram to your output list

  PostData(1, fOutput);
}

// Helper functions

void AliAnalysisTaskFemtoDreamRho::UserExec(Option_t *)
{

  AliAODEvent *Event = static_cast<AliAODEvent *>(fInputEvent);
  if (!Event)
  {
    AliWarning("No Input Event");
  }

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent))
    return;

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track)
      continue;
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // First we want to combine all charged pions with each other in the SE in order to find Rhos
  static std::vector<AliFemtoDreamBasePart> Particles; // pi+ candidates
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles; // pi- candidates
  AntiParticles.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles; // Rho candidates
  V0Particles.clear();
  static std::vector<AliFemtoDreamv0> V0Particles_full; // Rho candidates quick fix in order to retain full information about the daughters
  V0Particles_full.clear();
  static std::vector<AliFemtoDreamBasePart> Protons; // proton candidates
  Protons.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtons; // anti-proton candidates
  AntiProtons.clear();
  // static std::vector<AliFemtoDreamBasePart> RhoMcTrue; // Rhos verified via MC information
  // RhoMcTrue.clear();

  static float massChargedPion =
      TDatabasePDG::Instance()->GetParticle(fPosPionCuts->GetPDGCode())->Mass(); // as usual to minimize uncert.
  fRhoParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // Loop to identify all charged pions & protons
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track)
      continue;
    fTrack->SetTrack(track);
    if (fPosPionCuts->isSelected(fTrack))
    {
      fTrack->SetInvMass(massChargedPion); // Since we combine these later we set the inv. mass to the PDG value
      Particles.push_back(*fTrack);
    }
    if (fNegPionCuts->isSelected(fTrack))
    {
      fTrack->SetInvMass(massChargedPion); // Since we combine these later we set the inv. mass to the PDG value
      AntiParticles.push_back(*fTrack);
    }
    if (fPosProtonCuts->isSelected(fTrack))
    {
      Protons.push_back(*fTrack);
    }
    if (fNegProtonCuts->isSelected(fTrack))
    {
      AntiProtons.push_back(*fTrack);
    }
  }

  // Construct the V0 for the Rho decay, just simple combinatorix for now
  for (const auto &posPion : Particles)
  { // Build charged pion pairs!
    for (const auto &negPion : AntiParticles)
    {
      fRhoParticle->Setv0(posPion, negPion, Event);
      if (fRhoCuts->isSelected(fRhoParticle))
      { // Check for proper Rho candidates, just Minv cut and kaon reject.
        V0Particles.push_back(*fRhoParticle);
        // V0Particles_full.push_back(*fRhoParticle);
      }
    }
  }

  // Investigate the V0s in the MC sample, assign the motherID and check M_inv distributions as a function of pT.
  bool CheckV0s = true;
  bool useThisV0 = true;
  bool sameMother = true;
  bool diffMother = true;
  bool diffMotherID = true;
  bool dauArePhysPrim = false;
  bool dauAreNotPions = true;
  const int pdgIdealDaughters = 211;

  // fFlagHistogram->Fill(1, useThisV0);
  // fFlagHistogram->Fill(2, sameMother);
  // fFlagHistogram->Fill(3, diffMother);
  // fFlagHistogram->Fill(4, diffMotherID);
  // fFlagHistogram->Fill(5, dauArePhysPrim);
  // fFlagHistogram->Fill(6, dauAreNotPions);

  // What might also be interesting is to took at the Armenteros Podolanski plot of the different contributions?
  float pTposDau = -1;
  float pTnegDau = -1;
  float pTv0 = -1;
  float mINVv0 = -1;

  // std::cout << "Check the fIsMC flag: " << fIsMC << std::endl;

  if (fIsMC && CheckV0s)
  {
    for (const auto &V0 : V0Particles) // these are the reconstructed ones, later we also take a look at the full MC sample
    {
      std::cout << "This is one of the V0s (selfreconstructed)" << std::endl;
      // AliFatal("This is one of the V0s (selfreconstructed)");

      // match IDs to mcarray
      TClonesArray *mcarray = dynamic_cast<TClonesArray *>(Event->FindListObject(
          AliAODMCParticle::StdBranchName()));
      if (!mcarray)
      {
        AliFatal("No MC Array found");
      }
      const int negID = V0.GetIDTracks().at(0);
      const int posID = V0.GetIDTracks().at(1);
      std::cout << "Check IDs: negID:" << negID << " posID:" << posID << std::endl;

      AliAODMCParticle *mcPartPos = (AliAODMCParticle *)mcarray->At(posID);
      AliAODMCParticle *mcPartNeg = (AliAODMCParticle *)mcarray->At(negID);
      if (!mcPartPos || !mcPartNeg)
      {
        useThisV0 = false;
        continue;
      }
      const int negDauPDG = mcPartNeg->GetPdgCode();
      const int posDauPDG = mcPartPos->GetPdgCode();
      std::cout << "Check PDGs: negPDG:" << negDauPDG << " posPDG:" << posDauPDG << std::endl;
      fHistogramPDG->Fill(negDauPDG);
      fHistogramPDG->Fill(posDauPDG);

      const int motherIDnegDau = mcPartNeg->GetMother();
      const int motherIDposDau = mcPartPos->GetMother();
      std::cout << "Check motherLabels: negLabelmom:" << motherIDnegDau << " posLabelmom:" << motherIDposDau << std::endl;

      if (motherIDnegDau == motherIDposDau) // same mother ID MOVE THIS BEFORE CHECKING THE PDG_MOTHER CHECK ELSE SOME CODE WILL DOUBLE
      {
        diffMotherID = false;

        if (std::abs(negDauPDG) == pdgIdealDaughters && std::abs(posDauPDG) == pdgIdealDaughters) // mc true pion pair
        {
          dauAreNotPions = false;

          // here we only look at those V0s for which we have verified with MC that both daughters are charged pions!
          // for now we only search with a depth of 1 i.e. we only ask for the mother
          if (mcPartPos->IsPhysicalPrimary() && mcPartNeg->IsPhysicalPrimary()) // pions stem from primoridal production or strong reso
          {
            dauArePhysPrim = true;

            AliAODMCParticle *mcPartPosMother = (AliAODMCParticle *)mcarray->At(motherIDposDau);
            AliAODMCParticle *mcPartNegMother = (AliAODMCParticle *)mcarray->At(motherIDnegDau);
            if (!mcPartPosMother || !mcPartNegMother)
            {
              useThisV0 = false;
              continue;
            }
            const int negDauPDGmother = mcPartPosMother->GetPdgCode();
            const int posDauPDGmother = mcPartNegMother->GetPdgCode();
            // AliFemtoDreamTrack *posDaughter = V0.GetPosDaughter(); //This will not work as the downcast from the V0 to the AliFemtoDreamBasePart removes the information about the daughters
            // AliFemtoDreamTrack *negDaughter = V0.GetNegDaughter(); //For now just take the MC momentum as these are this accessible via their MC Label.
            if (posDauPDGmother != negDauPDGmother)
            { // Should be trivially fullfilled but better be safe than sry
              sameMother = true;
              mINVv0 = V0.GetInvMass();
              pTv0 = V0.GetPt();
              // pTposDau = posDaughter->GetPt();
              // pTnegDau = negDaughter->GetPt();
              pTposDau = mcPartPos->Pt();
              pTnegDau = mcPartNeg->Pt();
            }
            else if (posDauPDGmother == negDauPDGmother)
            {
              diffMother = false;
            }
            else
            {
              AliFatal("posDauPDGmother is neither different nor equal to negDauPDGmother");
              useThisV0 = false;
            }
          }
          else // these here are V0 for which the daughers do not pass IsPhysicalPrimary()
          {
            // fill here some QC output
            dauArePhysPrim = false;
          }
        }
        else // these here are V0 which are impure
        {
          // fill here some QC output
          dauAreNotPions = true;
        }
      }
      else // these here don't have the same mother
      {
        // fill here some QC output
        diffMotherID = true;
      }
      printf("Check Flags: \n");
      printf("useThisV0: %i\n", useThisV0);
      printf("sameMother: %i\n", sameMother);
      printf("diffMother: %i\n", diffMother);
      printf("diffMotherID: %i\n", diffMotherID);
      printf("dauArePhysPrim: %i\n", dauArePhysPrim);
      printf("dauAreNotPions: %i\n", dauAreNotPions);
      // Update the histogram based on flag values
      fFlagHistogram->Fill(1, useThisV0);
      fFlagHistogram->Fill(2, sameMother);
      fFlagHistogram->Fill(3, diffMother);
      fFlagHistogram->Fill(4, diffMotherID);
      fFlagHistogram->Fill(5, dauArePhysPrim);
      fFlagHistogram->Fill(6, dauAreNotPions);
    } // end of for (const auto &V0 : V0Particles)
  }   // end of if (isMC && CheckV0s)
  else
  {
    fFlagHistogram->Fill(7, 1);
  }

  // Alternative look for the Rhos as they are in the MC
  bool fIsMCTruth = true;
  const int pdgIdealRho = 113;

  if (fIsMC && fIsMCTruth)
  {
    printf("Getting all the Rhos from MC\n");

    TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(
        Event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fArrayMCAOD)
    {
      AliFatal("No MC Array found");
    }
    int amountOfPart = fArrayMCAOD->GetEntriesFast();
    int mcpdg;
    AliFemtoDreamBasePart part;
    AliFemtoDreamBasePart part2;
    for (int iPart = 1; iPart < amountOfPart; iPart++)
    {
      AliAODMCParticle *mcPart = (AliAODMCParticle *)fArrayMCAOD->At(iPart);
      if (!(mcPart))
      {
        std::cout << "NO MC particle" << std::endl;
        continue;
      }
      if (mcPart->GetLabel() < 0)
      {
        continue;
      }
      mcpdg = mcPart->GetPdgCode();
      if (std::abs(mcpdg) == pdgIdealRho) // combine particle and antiparticle
      {
        printf("Found a Rho in MC\n");

        int firstdaughter = mcPart->GetDaughterFirst(); // The B.R. for the two pion decay is 100% hence it is okay to only check one daughter
        int lastdaughter = mcPart->GetDaughterLast();   // The B.R. for the two pion decay is 100% hence it is okay to only check one daughter

        if (firstdaughter <= amountOfPart && lastdaughter <= amountOfPart)
        {
          AliAODMCParticle *mcDaughter =
              (AliAODMCParticle *)fArrayMCAOD->At(firstdaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?
          AliAODMCParticle *mcDaughter2 =
              (AliAODMCParticle *)fArrayMCAOD->At(lastdaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?

          if (mcDaughter && mcDaughter2)
          {
            int dpdg = mcDaughter->GetPdgCode();
            double dpt = mcDaughter->Pt();
            double deta = mcDaughter->Eta();

            int dpdg2 = mcDaughter2->GetPdgCode();
            double dpt2 = mcDaughter2->Pt();
            double deta2 = mcDaughter2->Eta();
            if (std::abs(dpdg) == pdgIdealDaughters)
            {
              printf("Found a Rho with pion daughters MC\n");

              if ((dpt < 999. && dpt > 0.15) && (deta > -0.8 && deta < 0.8)) // slight acceptance cuts
              {
                printf("Daughters are in acceptance\n");

                part.SetMCParticleRePart(mcPart); // Here we set the MC values for all the kinematic variables (i.e. no detector effects)
                                                  // RhoMcTrue.push_back(part);        // Handle the MC particle like a real one, see also function above
                                                  //  The important part is to use the existing infrastructure as much as possible
                                                  //  if (fRhoParticleMCTrueCuts->isSelected(part)) //but here we need a AliFemtoDreamTrack
                                                  //{
                                                  //  here is the proper construction but for now this also only contains the acceptance cuts, basically this is needed in order to have a structured output
                                                  // }
                printf("Acceptance; pt: 0.15<%.2f<999.\n", dpt);
                fpTerrorHistogram->Fill(dpt);
                fpTCorrerrorHistogram->Fill(dpt, dpt2);
                // ptHist_RhoMCTrue->Fill(dpt);
                //   massPtHist_RhoMCTrue->Fill(, dpt);
              }
            }
          }
        }
      }
    }
  }

  // For now no default pair cleaning but may change to CleanTrackandTrack of same charge pion and proton
  if (fDoCleaning)
  {
    fPairCleaner->CleanTrackAndDecay(&Particles, &Protons, 0);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiProtons, 1);
    fPairCleaner->CleanDecay(&V0Particles, 0);
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(V0Particles);
  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  // fPairCleaner->StoreParticle(RhoMcTrue);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamRho::ResetGlobalTrackReference()
{
  for (UShort_t i = 0; i < fTrackBufferSize; i++)
  {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamRho::StoreGlobalTrackReference(
    AliAODTrack *track)
{
  // for documentation see AliFemtoDreamAnalysis

  const int trackID = track->GetID();
  if (trackID < 0)
  {
    return;
  }
  if (trackID >= fTrackBufferSize)
  {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }
  if (fGTI[trackID])
  {
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls()))
    {
      return;
    }
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
    {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}