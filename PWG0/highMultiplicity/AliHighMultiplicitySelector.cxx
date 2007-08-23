/* $Id$ */

#include "AliHighMultiplicitySelector.h"

#include <TVector3.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include <AliITSgeom.h>
#include <AliITSLoader.h>
#include <AliITSdigitSPD.h>
#include <AliITSRecPoint.h>

#include "AliPWG0Helper.h"

//
//

ClassImp(AliHighMultiplicitySelector)

AliHighMultiplicitySelector::AliHighMultiplicitySelector() :
  AliSelectorRL(),
  fChipsLayer1(0),
  fChipsLayer2(0),
  fL1vsL2(0),
  fMvsL1(0),
  fMvsL2(0),
  fChipsFired(0),
  fPrimaryL1(0),
  fPrimaryL2(0),
  fClusterZL1(0),
  fClusterZL2(0),
  fClvsL1(0),
  fClvsL2(0),
  centralRegion(kFALSE)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliHighMultiplicitySelector::~AliHighMultiplicitySelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliHighMultiplicitySelector::SlaveBegin(TTree *tree)
{
  AliSelectorRL::SlaveBegin(tree);

  fChipsFired = new TH2F("fChipsFired", ";Module;Chip;Count", 240, -0.5, 239.5, 5, -0.5, 4.5);

  fPrimaryL1 = new TNtuple("fPrimaryL1", "", "multiplicity:firedChips:chipsByPrimaries:clusters");
  fPrimaryL2 = new TNtuple("fPrimaryL2", "", "multiplicity:firedChips:chipsByPrimaries:clusters");

  fClusterZL1 = new TH1F("fClusterZL1", ";z", 400, -20, 20);
  fClusterZL2 = new TH1F("fClusterZL2", ";z", 400, -20, 20);
}

void AliHighMultiplicitySelector::Init(TTree* tree)
{
  // read the user objects

  AliSelectorRL::Init(tree);

  // enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);

    /*if (fTree->GetCurrentFile())
    {
      TString fileName(fTree->GetCurrentFile()->GetName());
      fileName.ReplaceAll("AliESDs", "geometry");

      // load geometry
      TGeoManager::Import(fileName);
      }*/
  }
}

Bool_t AliHighMultiplicitySelector::Process(Long64_t entry)
{
  //
  // processing
  //

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD, AliPWG0Helper::kMB1);

  if (!eventTriggered)
  {
    AliDebug(AliLog::kDebug, "Event not triggered. Skipping.");
    return kTRUE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  Int_t nPrim  = stack->GetNprimary();
  Int_t multiplicity21 = 0;
  Int_t multiplicity16 = 0;

  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    if (centralRegion)
    {
      if (TMath::Abs(particle->Eta()) < 1.05)
        multiplicity21++;
      if (TMath::Abs(particle->Eta()) < 0.8)
        multiplicity16++;
    }
    else
    {
      if (TMath::Abs(particle->Eta()) < 2.1)
        multiplicity21++;
      if (TMath::Abs(particle->Eta()) < 1.6)
        multiplicity16++;
    }
  }// end of mc particle

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
  {
    AliDebug(AliLog::kError, "runloader not available");
    return kFALSE;
  }

  // TDirectory::TContext restores the current directory is restored when the scope ends.
  // This helps around ROOT bug #26025 and is good behaviour anyway
  TDirectory::TContext context(0);
  AliITSLoader* loader = (AliITSLoader*) runLoader->GetLoader( "ITSLoader" );
  loader->LoadDigits("READ");
  TTree* treeD = loader->TreeD();
  if (!treeD)
  {
    AliDebug(AliLog::kError, "Could not retrieve TreeD of ITS");
    return kFALSE;
  }

  treeD->SetBranchStatus("*", 0);
  treeD->SetBranchStatus("ITSDigitsSPD.fTracks*", 1);
  treeD->SetBranchStatus("ITSDigitsSPD.fCoord1", 1);

  TClonesArray* digits = 0;
  treeD->SetBranchAddress("ITSDigitsSPD", &digits);
  if (digits);
    digits->Clear();

  // each value for both layers
  Int_t totalNumberOfFO[2];
  Int_t chipsHitByPrimaries[2];
  //Int_t chipsHitBySecondaries[2];

  for (Int_t i=0; i<2; ++i)
  {
    totalNumberOfFO[i] = 0;
    chipsHitByPrimaries[i] = 0;
    //chipsHitBySecondaries[i] = 0;
  }

  Int_t startSPD = 0; //geom->GetStartSPD();
  Int_t lastSPD  = 239; //geom->GetLastSPD();

  //printf("%d %d\n", startSPD, lastSPD);
// for (Int_t l=0; l<240; ++l) { AliITSgeomTGeo::GetModuleId(l, i, j, k); printf("%d --> %d\n", l, i); }

  // loop over modules (ladders)
  for (Int_t moduleIndex=startSPD; moduleIndex<lastSPD+1; moduleIndex++)
  {
    if (centralRegion)
    {
      if ((moduleIndex % 4) == 0 || (moduleIndex % 4) == 3)
        continue;
    }

    Int_t currentLayer = 0;
    if (moduleIndex >= 80)
      currentLayer = 1;

    treeD->GetEvent(moduleIndex);

    // get number of digits in this module
    Int_t ndigitsInModule = digits->GetEntriesFast();

    // get number of digits in each chip of the module
    Int_t ndigitsInChip[5];
    Bool_t hitByPrimary[5];
    for( Int_t iChip=0; iChip<5; iChip++)
    {
      ndigitsInChip[iChip]=0;
      hitByPrimary[iChip] = kFALSE;
    }

    // loop over digits in this module
    for (Int_t iDig=0; iDig<ndigitsInModule; iDig++)
    {
      AliITSdigitSPD* dp = (AliITSdigitSPD*) digits->At(iDig);
      Int_t column = dp->GetCoord1();
      Int_t isChip = column / 32;

      //printf("Digit %d has column %d which translates to chip %d\n", iDig, column, isChip);

      fChipsFired->Fill(moduleIndex, isChip);

      ndigitsInChip[isChip]++;

      Bool_t debug = kFALSE;

      // find out which particle caused this chip to fire
      // if we find at least one primary we consider this chip being fired by a primary
      for (Int_t track=0; track<10; ++track)
      {
        Int_t label = dp->GetTrack(track);

        if (label < 0)
          continue;

        if (debug)
          printf("track %d contains label %d\n", track, label);

        TParticle* particle = stack->Particle(label);

        if (!particle)
        {
          AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (digit loop).", label));
          continue;
        }

        if (debug)
        {
          particle->Print();
          printf("statuscode = %d, p = %f, m = %d\n", particle->GetStatusCode(), particle->P(), particle->GetFirstMother());
        }

        // TODO delta electrons should be traced back to their mother. this is e.g. solved in AliITSClusterFinderV2::CheckLabels2
        while (particle->P() < 0.02 && particle->GetStatusCode() == 0 && particle->GetFirstMother() >= 0)
        {
          label = particle->GetFirstMother();
          particle = stack->Particle(label);

          if (!particle)
            break;

          if (debug)
          {
            printf("-->\n");
            printf("statuscode = %d, p = %f, m = %d\n", particle->GetStatusCode(), particle->P(), particle->GetFirstMother());
            particle->Print();
          }
        }

        if (!particle)
        {
          AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (digit loop, finding delta electrons).", label));
          continue;
        }

        if (label > nPrim)
          continue;

        if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
          continue;

        if (debug)
          printf("This was a primary (or delta-electron of a primary)!\n");

        hitByPrimary[isChip] = kTRUE;
      }
    }

    // get number of FOs in the module
    for (Int_t ifChip=0; ifChip<5; ifChip++)
      if( ndigitsInChip[ifChip] >= 1 )
      {
        totalNumberOfFO[currentLayer]++;
        if (hitByPrimary[ifChip])
        {
          chipsHitByPrimaries[currentLayer]++;
        }
        //else
        //  chipsHitBySecondaries[currentLayer]++;
      }
  }

  //printf("Fired chips: %d %d\n", totalNumberOfFOLayer1, totalNumberOfFOLayer2);

  // now find clusters
  Int_t clustersLayer[2];
  clustersLayer[0] = 0;
  clustersLayer[1] = 0;

  loader->LoadRecPoints("READ");
  TTree* treeR = loader->TreeR();
  if (!treeR)
  {
    AliDebug(AliLog::kError, "Could not retrieve TreeR of ITS");
    return kFALSE;
  }

  // TODO branches!
  //treeR->SetBranchStatus("*", 0);

  TClonesArray* itsClusters = 0;
  treeR->SetBranchAddress("ITSRecPoints", &itsClusters);

  Int_t nTreeEntries = treeR->GetEntries();
  for (Int_t iEntry = 0; iEntry < nTreeEntries; ++iEntry)
  {
    if (!treeR->GetEvent(iEntry))
      continue;

    Int_t nClusters = itsClusters->GetEntriesFast();

    while(nClusters--)
    {
      AliITSRecPoint* cluster = (AliITSRecPoint*) itsClusters->UncheckedAt(nClusters);

      if (cluster->GetLayer() == 0)
      {
        clustersLayer[0]++;
        fClusterZL1->Fill(cluster->GetZ());
      }
      else if (cluster->GetLayer() == 1)
      {
        clustersLayer[1]++;
        fClusterZL2->Fill(cluster->GetZ());
      }
    }
  }

  fPrimaryL1->Fill(multiplicity21, totalNumberOfFO[0], chipsHitByPrimaries[0], clustersLayer[0]);
  fPrimaryL2->Fill(multiplicity16, totalNumberOfFO[1], chipsHitByPrimaries[1], clustersLayer[1]);

  return kTRUE;
}

Bool_t AliHighMultiplicitySelector::Notify()
{
  AliRunLoader* runLoader = GetRunLoader();

  if (runLoader)
  {
    AliITSLoader* loader = (AliITSLoader* )runLoader->GetLoader( "ITSLoader" );
    if (loader)
    {
      loader->UnloadDigits();
      loader->UnloadRecPoints();
    }
  }

  return AliSelectorRL::Notify();
}

void AliHighMultiplicitySelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized.");
    return;
  }

  fOutput->Add(fChipsFired);
  fOutput->Add(fPrimaryL1);
  fOutput->Add(fPrimaryL2);
  fOutput->Add(fClusterZL1);
  fOutput->Add(fClusterZL2);
}

void AliHighMultiplicitySelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fChipsFired = dynamic_cast<TH2F*> (fOutput->FindObject("fChipsFired"));
  fPrimaryL1 = dynamic_cast<TNtuple*> (fOutput->FindObject("fPrimaryL1"));
  fPrimaryL2 = dynamic_cast<TNtuple*> (fOutput->FindObject("fPrimaryL2"));
  fClusterZL1 = dynamic_cast<TH1F*> (fOutput->FindObject("fClusterZL1"));
  fClusterZL2 = dynamic_cast<TH1F*> (fOutput->FindObject("fClusterZL2"));

  if (!fClusterZL1 || !fClusterZL2 || !fChipsFired || !fPrimaryL1 || !fPrimaryL2)
  {
    AliError("Histograms not available");
    return;
  }

  WriteHistograms();
}

void AliHighMultiplicitySelector::WriteHistograms(const char* filename)
{
  TFile* file = TFile::Open(filename, "RECREATE");

  fChipsFired->Write();
  fPrimaryL1->Write();
  fPrimaryL2->Write();
  fClusterZL1->Write();
  fClusterZL2->Write();

  file->Close();
}

void AliHighMultiplicitySelector::ReadHistograms(const char* filename)
{
  TFile* file = TFile::Open(filename);

  if (!file)
    return;

  fPrimaryL1  = dynamic_cast<TNtuple*> (file->Get("fPrimaryL1"));
  fPrimaryL2  = dynamic_cast<TNtuple*> (file->Get("fPrimaryL2"));
  fChipsFired  = dynamic_cast<TH2F*> (file->Get("fChipsFired"));
  fClusterZL1  = dynamic_cast<TH1F*> (file->Get("fClusterZL1"));
  fClusterZL2  = dynamic_cast<TH1F*> (file->Get("fClusterZL2"));

  #define MULT   1001, -0.5, 1000.5
  #define BINNING_LAYER1 401, -0.5, 400.5
  #define BINNING_LAYER2 801, -0.5, 800.5

  fChipsLayer1 = new TH1F("fChipsLayer1", "Layer 1;Fired Chips;Count", BINNING_LAYER1);
  fChipsLayer2 = new TH1F("fChipsLayer2", "Layer 2;Fired Chips;Count", BINNING_LAYER2);

  fL1vsL2 = new TH2F("fL1vsL2", ";Fired Chips Layer 1;Fired Chips Layer 2", BINNING_LAYER1, BINNING_LAYER2);
  fMvsL1 = new TH2F("fMvsL1", ";true multiplicity;fired chips layer1", MULT, BINNING_LAYER1);
  fMvsL2 = new TH2F("fMvsL2", ";true multiplicity;fired chips layer2", MULT, BINNING_LAYER2);

  fClvsL1 = new TH2F("fClvsL1", ";clusters layer1;fired chips layer1", MULT, BINNING_LAYER1);
  fClvsL2 = new TH2F("fClvsL2", ";clusters layer2;fired chips layer2", MULT, BINNING_LAYER2);

  for (Int_t i = 0; i < fPrimaryL1->GetEntries(); i++)
  {
    fPrimaryL1->GetEvent(i);
    fPrimaryL2->GetEvent(i);

    Int_t multiplicity21 = (Int_t) fPrimaryL1->GetArgs()[0];
    Int_t multiplicity16 = (Int_t) fPrimaryL2->GetArgs()[0];

    Int_t totalNumberOfFO[2];
    totalNumberOfFO[0] = (Int_t) fPrimaryL1->GetArgs()[1];
    totalNumberOfFO[1] = (Int_t) fPrimaryL2->GetArgs()[1];

    Int_t chipsHitByPrimaries[2];
    chipsHitByPrimaries[0] = (Int_t) fPrimaryL1->GetArgs()[2];
    chipsHitByPrimaries[1] = (Int_t) fPrimaryL2->GetArgs()[2];

    Int_t clustersLayer[2];
    clustersLayer[0] = (Int_t) fPrimaryL1->GetArgs()[3];
    clustersLayer[1] = (Int_t) fPrimaryL2->GetArgs()[3];

    fChipsLayer1->Fill(totalNumberOfFO[0]);
    fChipsLayer2->Fill(totalNumberOfFO[1]);

    fL1vsL2->Fill(totalNumberOfFO[0], totalNumberOfFO[1]);

    fMvsL1->Fill(multiplicity21, totalNumberOfFO[0]);
    fMvsL2->Fill(multiplicity16, totalNumberOfFO[1]);

    fClvsL1->Fill(clustersLayer[0], totalNumberOfFO[0]);
    fClvsL2->Fill(clustersLayer[1], totalNumberOfFO[1]);
  }
}

TH1* AliHighMultiplicitySelector::GetTriggerEfficiency(TH2* multVsLayer, Int_t cut)
{
  //
  // returns the trigger efficiency as function of multiplicity with a given cut
  //

  //cut and multiply with x-section
  TH1* allEvents = multVsLayer->ProjectionX("fMvsL_x_total", 1, 1001);
  //allEvents->Sumw2();

  //cut and multiply with x-section
  TH1* proj = multVsLayer->ProjectionX(Form("%s_x", multVsLayer->GetName()), cut, 1001);
  //proj->Sumw2();

  //new TCanvas; allEvents->DrawCopy(); gPad->SetLogy();
  //new TCanvas; proj->DrawCopy(); gPad->SetLogy();

  // make probability distribution out of it
  // TODO binomial errors do not work??? weird...
  proj->Divide(proj, allEvents, 1, 1, "B");

  return proj;
}

TH1* AliHighMultiplicitySelector::GetXSectionCut(TH1* xSection, TH2* multVsLayer, Int_t cut)
{
  // returns the rel. cross section of the true spectrum that is measured when a cut at <cut> is performed

  TH1* proj = GetTriggerEfficiency(multVsLayer, cut);

  //new TCanvas; proj->DrawCopy(); gPad->SetLogy();

  for (Int_t i=1; i<=proj->GetNbinsX(); i++)
  {
    if (i <= xSection->GetNbinsX())
    {
      Double_t value = proj->GetBinContent(i) * xSection->GetBinContent(i);
      Double_t error = 0;

      if (value != 0)
        error = value * (proj->GetBinError(i) / proj->GetBinContent(i) + xSection->GetBinError(i) / xSection->GetBinContent(i));

      proj->SetBinContent(i, value);
      proj->SetBinError(i, error);
    }
    else
    {
      proj->SetBinContent(i, 0);
      proj->SetBinError(i, 0);
    }
  }

  //new TCanvas; proj->DrawCopy(); gPad->SetLogy();

  return proj;
}

void AliHighMultiplicitySelector::MakeGraphs2(const char* title, TH1* xSection, TH2* fMvsL)
{
  TGraph* effGraph = new TGraph;
  effGraph->SetTitle(Form("%s;Cut on fired chips;mult. where eff. >50%%", title));
  TGraph* ratioGraph = new TGraph;
  ratioGraph->SetTitle(Form("%s;Cut on fired chips;x-section_(>=eff. limit) / x-section_(total)", title));
  TGraph* totalGraph = new TGraph;
  totalGraph->SetTitle(Form("%s;Cut on fired chips;rel x-section_(>=eff. limit)", title));

  for (Int_t cut = 0; cut <= 300; cut+=50)
  {
    TH1* proj = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("clone");

    //proj->Rebin(3);
    //proj->Scale(1.0 / 3);

    new TCanvas; proj->DrawCopy();

    Int_t limitBin = proj->GetNbinsX()+1;
    while (limitBin > 1 && proj->GetBinContent(limitBin-1) > 0.5)
      limitBin--;

    Float_t limit = proj->GetXaxis()->GetBinCenter(limitBin);

    effGraph->SetPoint(effGraph->GetN(), cut, limit);

    proj = GetXSectionCut(xSection, fMvsL, cut);

    Double_t ratio = 0;
    Double_t total = 0;
    if (proj->Integral(1, 1001) > 0)
    {
      ratio = proj->Integral(proj->FindBin(limit), 1001) / proj->Integral(1, 1001);
      total = proj->Integral(proj->FindBin(limit), 1001);
    }

    ratioGraph->SetPoint(ratioGraph->GetN(), cut, ratio);
    totalGraph->SetPoint(totalGraph->GetN(), cut, total);

    Printf("Cut at %d --> trigger eff. is > 0.5 for mult. >= %.2f. That is the case for %f of the triggered, %e of all events", cut, limit, ratio, total);
  }

  TCanvas* canvas = new TCanvas(Form("%s_Efficiency", title), Form("%s_Efficiency", title), 1200, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  effGraph->Draw("A*");

  for (Int_t i=8; i<=10; ++i)
  {
    TLine* line = new TLine(0, xSection->GetMean() * i, 300, xSection->GetMean() * i);
    line->Draw();
  }

  canvas->cd(2);  ratioGraph->Draw("A*");
  canvas->cd(3);  gPad->SetLogy(); totalGraph->Draw("A*");

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
}

void AliHighMultiplicitySelector::MakeGraphs(const char* title, TH1* xSection, TH2* fMvsL, Int_t limit)
{
  // relative x-section, once we have a collision
  xSection->Scale(1.0 / xSection->Integral());

  TGraph* ratioGraph = new TGraph;
  ratioGraph->SetTitle(Form("%s;Cut on fired chips;x-section_(>=%d) / x-section_(total)", title, limit));
  TGraph* totalGraph = new TGraph;
  totalGraph->SetTitle(Form("%s;Cut on fired chips;rel x-section_(>=%d)", title, limit));

  Double_t max = 0;
  Int_t bestCut = -1;
  Double_t bestRatio = -1;
  Double_t bestTotal = -1;
  Int_t fullCut = -1;
  Double_t fullRatio = -1;
  Double_t fullTotal = -1;

  fMvsL->Sumw2();

  for (Int_t cut = 50; cut <= 300; cut+=2)
  {
    TH1* proj = GetXSectionCut(xSection, fMvsL, cut);

    Double_t ratio = 0;
    Double_t total = 0;
    if (proj->Integral(1, 1000) > 0)
    {
      ratio = proj->Integral(limit, 1000) / proj->Integral(1, 1000);
      total = proj->Integral(limit, 1000);
    }

    max = TMath::Max(max, total);

    //printf("Cut at %d: rel. x-section_(>=%d) = %e; x-section_(>=%d) / x-section_(total) = %f\n", cut, limit, total, limit, ratio);

    if (total < max * 0.9 && bestCut == -1)
    {
      bestCut = cut;
      bestRatio = ratio;
      bestTotal = total;
    }

    if (ratio == 1 && fullCut == -1)
    {
      fullCut = cut;
      fullRatio = ratio;
      fullTotal = total;
    }

    ratioGraph->SetPoint(ratioGraph->GetN(), cut, ratio);
    totalGraph->SetPoint(totalGraph->GetN(), cut, total);
  }

  if (bestCut != -1)
    printf("Best cut at %d: rel. x-section_(>=%d) = %e %%; x-section_(>=%d) / x-section_(total) = %f %%\n", bestCut, limit, bestTotal, limit, bestRatio);
  if (fullCut != -1)
    printf("100%% cut at %d: rel. x-section_(>=%d) = %e %%; x-section_(>=%d) / x-section_(total) = %f %%\n", fullCut, limit, fullTotal, limit, fullRatio);

  TCanvas* canvas = new TCanvas(Form("%s_RatioXSection_%d", title, limit), Form("%s_RatioXSection_%d", title, limit), 600, 400);
  ratioGraph->Draw("A*");
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  canvas = new TCanvas(Form("%s_TotalXSection_%d", title, limit), Form("%s_TotalXSection_%d", title, limit), 600, 400);
  totalGraph->Draw("A*");
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
}

void AliHighMultiplicitySelector::JPRPlots()
{
  /*

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->JPRPlots();

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  for (Int_t i=0; i<2; ++i)
  {
    if (!xSections[i])
      continue;

    TH1* xSection = xSections[i];
    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;
    //Int_t cut = (i == 0) ? 164 : 150; // 8 times the mean
    //Int_t cut = (i == 0) ? 178 : 166; // 9 times the mean
    Int_t cut = (i == 0) ? 190 : 182; // 10 times the mean

    // limit is N times the mean
    Int_t limit = (Int_t) (xSection->GetMean() * 10);

    // 10^28 lum --> 1.2 kHz
    // 10^31 lum --> 1200 kHz
    Float_t rate = 1200e3;

    // time in s
    Float_t lengthRun = 1e6;

    xSection->SetStats(kFALSE);
    xSection->SetTitle(""); //(i == 0) ? "SPD Layer 1" : "SPD Layer 2");
    xSection->GetXaxis()->SetTitle(Form("true multiplicity in |#eta| < %.1f", (i == 0) ? 2.0 : 1.5));
    xSection->GetXaxis()->SetRangeUser(0, (i == 0) ? 400 : 350);
    //xSection->GetYaxis()->SetTitle("relative cross section");
    xSection->GetYaxis()->SetTitleOffset(1.2);

    // relative x-section, once we have a collision
    xSection->Scale(1.0 / xSection->Integral());

    TH1* proj = GetXSectionCut(xSection, fMvsL, cut);

    Double_t ratio = 0;
    Double_t total = 0;
    if (proj->Integral(1, 1000) > 0)
    {
      ratio = proj->Integral(limit, 1000) / proj->Integral(1, 1000);
      total = proj->Integral(limit, 1000);
    }

    printf("Cut at %d: rel. x-section_(>=%d) = %e; x-section_(>=%d) / x-section_(total) = %f\n", cut, limit, total, limit, ratio);

    TCanvas* canvas = new TCanvas(Form("HMPlots_%d", i), Form("HMPlots_%d", i), 800, 600);
    canvas->SetLogy();
    xSection->DrawCopy();
    proj->SetLineColor(2);
    proj->SetStats(kFALSE);
    proj->DrawCopy("SAME");

    TLegend* legend = new TLegend(0.15, 0.15, 0.45, 0.3);
    legend->SetFillColor(0);
    legend->AddEntry(xSection, "no trigger");
    legend->AddEntry(proj, Form("FO trigger > %d chips", cut));
    legend->Draw();

    TLine* line = new TLine(limit, xSection->GetMinimum() * 0.5, limit, xSection->GetMaximum() * 2);
    line->SetLineWidth(2);
    line->Draw();

    canvas->SaveAs(Form("%s.gif", canvas->GetName()));

    TCanvas* canvas2 = new TCanvas(Form("HMPlots_%d_Random", i), Form("HMPlots_%d_Random", i), 800, 600);
    //canvas2->SetTopMargin(0.05);
    //canvas2->SetRightMargin(0.05);
    canvas2->SetLogy();
    xSection->DrawCopy("HIST");

    TLegend* legend2 = new TLegend(0.15, 0.15, 0.6, 0.3);
    legend2->SetFillColor(0);
    legend2->AddEntry(xSection, "no trigger");

    TH1* proj2 = (TH1*) proj->Clone("random");
    proj2->Reset();
    // MB lengthRun s 100 Hz
    Int_t nTrigger = (Int_t) (100 * lengthRun * proj->Integral(1, 1000));
    proj2->FillRandom(proj, nTrigger);

    TH1* proj3 = (TH1*) proj2->Clone("random_clone");
    proj3->Divide(proj);
    proj3->Fit("pol0", "0", "");
    proj2->Scale(1.0 / proj3->GetFunction("pol0")->GetParameter(0));

    /*
    proj2->DrawCopy("SAME");
    legend2->AddEntry(proj2, Form("%d evts, FO > %d chips (%d evts)", nTrigger, cut, (Int_t) (nTrigger * ratio)));
    */

    proj2 = (TH1*) proj->Clone("random2");
    proj2->Reset();
    // 10^31 lum --> 1200 kHz; lengthRun s
    nTrigger = (Int_t) (rate * proj->Integral(1, 1000) * lengthRun);
    proj2->FillRandom(proj, nTrigger);

    proj3 = (TH1*) proj2->Clone("random_clone2");
    proj3->Divide(proj);
    proj3->Fit("pol0", "0", "");
    proj2->Scale(1.0 / proj3->GetFunction("pol0")->GetParameter(0));

    proj2->SetLineColor(4);
    proj2->SetMarkerStyle(7);
    proj2->SetMarkerColor(4);
    proj2->DrawCopy("SAME P");
    //legend2->AddEntry(proj2, Form("%d evts, FO > %d chips (%d evts)", nTrigger, cut, (Int_t) (nTrigger * ratio)));
    legend2->AddEntry(proj2, Form("FO trigger > %d chips", cut));

    legend2->Draw();
    line->Draw();

    canvas2->SaveAs(Form("%s.gif", canvas2->GetName()));
    canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  }
}

void AliHighMultiplicitySelector::Ntrigger()
{
  //
  // produces a spectrum created with N triggers
  // number of triggers and thresholds for the moment fixed
  //

  /*

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Ntrigger();

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  // 10^28 lum --> 1.2 kHz
  // 10^31 lum --> 1200 kHz
  //Float_t rate = 1200e3;
  Float_t rate = 1200e3;

  // time in s
  Float_t lengthRun = 1e6;

  Int_t colors[] = { 2, 3, 4, 6, 7, 8 };
  Int_t markers[] = { 7, 2, 4, 5, 6, 27 };

  // put to 2 for second layer
  for (Int_t i=0; i<1; ++i)
  {
    if (!xSections[i])
      continue;

    TH1* xSection = xSections[i];
    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;

    Int_t nCuts = 6;
    Int_t cuts[] = { 0, 164, 178, 190, 204, 216 };
    //Int_t nCuts = 4;
    //Int_t cuts[] = { 0, 164, 190, 216 };
    // desired trigger rate in Hz
    Float_t ratePerTrigger[] = { 100, 1, 1, 1, 1, 1 };

    xSection->SetStats(kFALSE);
    xSection->SetTitle(""); //(i == 0) ? "SPD Layer 1" : "SPD Layer 2");
    xSection->GetXaxis()->SetTitle(Form("true multiplicity in |#eta| < %.1f", (i == 0) ? 2.0 : 1.5));
    xSection->GetXaxis()->SetRangeUser(0, (i == 0) ? 450 : 350);
    //xSection->GetYaxis()->SetTitle("relative cross section");
    xSection->GetYaxis()->SetTitleOffset(1.2);

    // relative x-section, once we have a collision
    xSection->Scale(1.0 / xSection->Integral());

    TCanvas* canvas2 = new TCanvas(Form("HMPlots2_%d_Random", i), Form("HMPlots2_%d_Random", i), 800, 600);
    canvas2->SetTopMargin(0.05);
    canvas2->SetRightMargin(0.05);
    canvas2->SetLogy();
    xSection->DrawCopy("HIST");

    TLegend* legend2 = new TLegend(0.15, 0.15, 0.6, 0.4);
    legend2->SetFillColor(0);
    legend2->AddEntry(xSection, "cross-section");

    for (Int_t currentCut = 0; currentCut<nCuts; ++currentCut)
    {
      Int_t cut = cuts[currentCut];

      TH1* proj = GetXSectionCut(xSection, fMvsL, cut);

      Double_t total = 0;
      if (proj->Integral(1, 1000) > 0)
        total = proj->Integral(1, 1000);

      printf("Cut at %d: rel. x-section = %e\n", cut, total);

      TH1* proj2 = (TH1*) proj->Clone("random2");
      proj2->Reset();

      // calculate downscale factor
      Float_t normalRate = rate * proj->Integral(1, 1000);
      Float_t downScale = normalRate / ratePerTrigger[currentCut];
      if (downScale < 1)
        downScale = 1;
      Long64_t nTrigger = (Long64_t) (normalRate / downScale * lengthRun);

      Printf("Normal rate is %f, downscale: %f, Simulating %lld triggers", normalRate, downScale, nTrigger);
      proj2->FillRandom(proj, nTrigger);

      Int_t removed = 0;
      for (Int_t bin=1; bin<proj2->GetNbinsX(); ++bin)
        if (proj2->GetBinContent(bin) < 5)
        {
          removed += (Int_t) proj2->GetBinContent(bin);
          proj2->SetBinContent(bin, 0);
        }

      Printf("Removed %d events", removed);

      proj2->Scale(1.0 / nTrigger * proj->Integral(1, 1000));

      proj2->SetLineColor(colors[currentCut]);
      proj2->SetMarkerStyle(markers[currentCut]);
      proj2->SetMarkerColor(colors[currentCut]);
      proj2->DrawCopy("SAME P");
      legend2->AddEntry(proj2, Form("%lld evts, FO > %d chips", nTrigger, cut));
    }

    legend2->Draw();

    canvas2->SaveAs(Form("%s.gif", canvas2->GetName()));
    canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  }
}

void AliHighMultiplicitySelector::DrawHistograms()
{
  /*

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_pythia.root");
  x->DrawHistograms();


  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing.root");
  x->DrawHistograms();

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_central.root");
  x->DrawHistograms();

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->DrawHistograms();

  */

  /*TCanvas* canvas = new TCanvas("chips", "chips", 600, 400);

  fChipsLayer2->SetLineColor(2);
  fChipsLayer2->SetStats(kFALSE);
  fChipsLayer1->SetStats(kFALSE);
  fChipsLayer2->SetTitle("");
  fChipsLayer2->DrawCopy();
  fChipsLayer1->DrawCopy("SAME");
  canvas->SaveAs("chips.gif");

  canvas = new TCanvas("L1vsL2", "L1vsL2", 600, 400);
  fL1vsL2->SetStats(kFALSE);
  fL1vsL2->DrawCopy("COLZ");
  gPad->SetLogz();
  canvas->SaveAs("L1vsL2.gif");*/

  TCanvas *canvas = new TCanvas("L1", "L1", 800, 600);
  canvas->SetTopMargin(0.05);
  canvas->SetRightMargin(0.12);
  fMvsL1->SetStats(kFALSE);
  fMvsL1->DrawCopy("COLZ");
  gPad->SetLogz();

  canvas->SaveAs("L1NoCurve.gif");
  canvas->SaveAs("L1NoCurve.eps");

  // draw corresponding theoretical curve
  TF1* func = new TF1("func", "[0]*(1-(1-1/[0])**x)", 1, 1000);
  func->SetParameter(0, 400-5*2);
  func->DrawCopy("SAME");

  canvas->SaveAs("L1.gif");

  canvas = new TCanvas("L2", "L2", 600, 400);
  //fMvsL2->GetYaxis()->SetRangeUser(0, 150);
  fMvsL2->SetStats(kFALSE);
  fMvsL2->DrawCopy("COLZ");
  gPad->SetLogz();
  func->SetParameter(0, 800-5*4);
  func->DrawCopy("SAME");
  canvas->SaveAs("L2.gif");

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (file)
  {
    TH1* xSection2 = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
    TH1* xSection15 = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

    MakeGraphs2("Layer1", xSection2, fMvsL1);
    return;

    // 5 times the mean
    //MakeGraphs("Layer1", xSection2, fMvsL1, (Int_t) (xSection2->GetMean() * 5)); //75 * 2 * 2);
    //MakeGraphs("Layer2", xSection15, fMvsL2, (Int_t) (xSection15->GetMean() * 5)); //(Int_t) (75 * 1.5 * 2));

    MakeGraphs("Layer1", xSection2, fMvsL1, (Int_t) (xSection2->GetMean() * 8));
    MakeGraphs("Layer2", xSection15, fMvsL2, (Int_t) (xSection15->GetMean() * 8));
    MakeGraphs("Layer1", xSection2, fMvsL1, (Int_t) (xSection2->GetMean() * 9));
    MakeGraphs("Layer2", xSection15, fMvsL2, (Int_t) (xSection15->GetMean() * 9));
    MakeGraphs("Layer1", xSection2, fMvsL1, (Int_t) (xSection2->GetMean() * 10));
    MakeGraphs("Layer2", xSection15, fMvsL2, (Int_t) (xSection15->GetMean() * 10));

    file->Close();
  }

  // make spread hists
  TGraph* spread = new TGraph;
  spread->SetTitle("Spread L1;true multiplicity;RMS");

  for (Int_t i=1; i<=fMvsL1->GetNbinsX(); ++i)
  {
    TH1* proj = fMvsL1->ProjectionY("proj", i, i);
    spread->SetPoint(spread->GetN(), i, proj->GetRMS());
  }

  canvas = new TCanvas("SpreadL1", "SpreadL1", 600, 400);
  spread->Draw("A*");
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  TF1* log = new TF1("log", "[0]*log([1]*x)", 1, 150);
  log->SetParLimits(0, 0, 10);
  log->SetParLimits(1, 1e-5, 10);

  spread->Fit(log, "", "", 1, 150);
  log->DrawCopy("SAME");

  TGraph* spread2 = new TGraph;
  spread2->SetTitle("Spread L2;true multiplicity;RMS");

  for (Int_t i=1; i<=fMvsL1->GetNbinsX(); ++i)
  {
    TH1* proj = fMvsL2->ProjectionY("proj", i, i);
    spread2->SetPoint(spread2->GetN(), i, proj->GetRMS());
  }

  canvas = new TCanvas("SpreadL2", "SpreadL2", 600, 400);
  spread2->Draw("A*");
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  spread2->Fit(log, "", "", 1, 150);
  log->DrawCopy("SAME");

  canvas = new TCanvas("Clusters_L1", "Clusters_L1", 600, 400);
  fClvsL1->SetStats(kFALSE);
  fClvsL1->DrawCopy("COLZ");
  gPad->SetLogz();

  func->SetParameter(0, 400-5*2);
  func->DrawCopy("SAME");

  canvas->SaveAs("Clusters_L1.gif");

  canvas = new TCanvas("Clusters_L2", "Clusters_L2", 600, 400);
  fClvsL2->SetStats(kFALSE);
  fClvsL2->DrawCopy("COLZ");
  gPad->SetLogz();
  func->SetParameter(0, 800-5*4);
  func->DrawCopy("SAME");
  canvas->SaveAs("Clusters_L2.gif");

  canvas = new TCanvas("ChipsFired", "ChipsFired", 600, 400);
  //fChipsFired->GetYaxis()->SetRangeUser(0, 150);
  fChipsFired->SetStats(kFALSE);
  fChipsFired->DrawCopy("COLZ");
  canvas->SaveAs("ChipsFired.gif");

  /*TH1F* tresholdHistL1 = new TH1F("tresholdHistL1", ";chip treshold;<n>", BINNING_LAYER1);
  TH1F* tresholdHistL2 = new TH1F("tresholdHistL2", ";chip treshold;<n>", BINNING_LAYER2);

  for (Int_t treshold = 0; treshold < 800; treshold++)
  {
    if (fPrimaryL1->Draw("multiplicity>>mult", Form("firedChips>%d", treshold), "goff") > 0)
    {
      TH1F* mult = dynamic_cast<TH1F*> (fPrimaryL1->GetHistogram());
      if (mult)
        tresholdHistL1->Fill(treshold, mult->GetMean());
    }
    if (fPrimaryL2->Draw("multiplicity>>mult", Form("firedChips>%d", treshold), "goff") > 0)
    {
      TH1F* mult = dynamic_cast<TH1F*> (fPrimaryL2->GetHistogram());
      if (mult)
        tresholdHistL2->Fill(treshold, mult->GetMean());
    }
  }

  canvas = new TCanvas("TresholdL1", "TresholdL1", 600, 400);
  tresholdHistL1->Draw();
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  canvas = new TCanvas("TresholdL2", "TresholdL2", 600, 400);
  tresholdHistL2->Draw();
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));*/

  fPrimaryL1->Draw("(chipsByPrimaries/firedChips):multiplicity>>eff1", "", "prof goff");
  fPrimaryL2->Draw("(chipsByPrimaries/firedChips):multiplicity>>eff2", "", "prof goff");

  canvas = new TCanvas("Efficiency", "Efficiency", 600, 400);
  fPrimaryL1->GetHistogram()->SetStats(kFALSE);
  fPrimaryL1->GetHistogram()->Draw();
  fPrimaryL2->GetHistogram()->SetLineColor(2);
  fPrimaryL2->GetHistogram()->SetStats(kFALSE);
  fPrimaryL2->GetHistogram()->Draw("SAME");
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  canvas = new TCanvas("ClustersZL1", "ClustersZL1", 600, 400);
  fClusterZL1->Rebin(2);
  fClusterZL1->Draw();
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  canvas = new TCanvas("ClustersZL2", "ClustersZL2", 600, 400);
  fClusterZL2->Draw();
  fClusterZL2->Rebin(2);
  canvas->SaveAs(Form("%s.gif", canvas->GetName()));
}
