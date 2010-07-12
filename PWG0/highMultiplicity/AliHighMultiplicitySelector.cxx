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
#include <TMath.h>
#include <TLegend.h>
#include <TText.h>

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
// Selector that produces plots needed for the high multiplicity analysis with the
// pixel detector
// Needs cluster file
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
  fCentralRegion(kFALSE)
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
  // create histograms

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

  static AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis;
  Bool_t eventTriggered = triggerAnalysis->IsTriggerBitFired(fESD->GetTriggerMask(), AliTriggerAnalysis::kMB1);

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

    if (fCentralRegion)
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
  if (digits)
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
    if (fCentralRegion)
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
  // get next ITS runloader

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
  // write the histograms to a file

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
  // read the data from a file and fill histograms

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

TH1* AliHighMultiplicitySelector::GetTriggerEfficiency(TH2* multVsLayer, Int_t cut, Int_t upperCut)
{
  //
  // returns the trigger efficiency as function of multiplicity with a given cut
  //

  //cut and multiply with x-section
  TH1* allEvents = multVsLayer->ProjectionX("fMvsL_x_total", 1, 1001);
  //allEvents->Sumw2();

  //cut and multiply with x-section
  TH1* proj = multVsLayer->ProjectionX(Form("%s_x", multVsLayer->GetName()), cut, upperCut);
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
  // creates graphs

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
  // more plots

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

void AliHighMultiplicitySelector::Ntrigger(Bool_t relative)
{
  //
  // produces a spectrum created with N triggers
  // number of triggers and thresholds for the moment fixed
  //

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Ntrigger();

  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Ntrigger(kFALSE);

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  // 10^28 lum --> 1.4 kHz
  // 10^31 lum --> 1400 kHz
  //Float_t rate = 1400e3;
  Float_t rate = 1.4e3;

  // time in s
  Float_t lengthRun = 1e5;

  Int_t colors[] = { 2, 3, 4, 6, 7, 8 };
  Int_t markers[] = { 7, 2, 4, 5, 6, 27 };

  // put to 2 for second layer
  for (Int_t i=0; i<1; ++i)
  {
    if (!xSections[i])
      continue;

    TH1* xSection = xSections[i];
    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;

    //Int_t nCuts = 6;
    //Int_t cuts[] = { 0, 164, 178, 190, 204, 216 };
    //Int_t nCuts = 4;
    //Int_t cuts[] = { 0, 164, 190, 216 };

    //Int_t nCuts = 4;
    //Int_t cuts[] = { 0, 114, 145, 165 };
    //Float_t ratePerTrigger[] = { 60, 13.3, 13.3, 13.3 };

    Int_t nCuts = 3;
    Int_t cuts[] = { 0, 114, 148 };

    //Int_t nCuts = 3;
    //Int_t cuts[] = { 0, 126, 162 };
    //Float_t ratePerTrigger[] = { 60, 20.0, 20.0 };

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

    if (relative)
      xSection->DrawCopy("HIST");

    TLegend* legend2 = new TLegend(0.15, 0.15, 0.6, 0.4);
    legend2->SetFillColor(0);
    legend2->AddEntry(xSection, "cross-section");

    for (Int_t currentCut = 0; currentCut<nCuts; ++currentCut)
    {
      Int_t cut = cuts[currentCut];

      TH1* triggerEff = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");

      TH1* proj = GetXSectionCut(xSection, fMvsL, cut);

      Float_t triggerLimit = 0;
      for (Int_t bin = 1; bin <= triggerEff->GetNbinsX(); bin++)
        if (triggerEff->GetBinContent(bin) < 0.5)
          triggerLimit = triggerEff->GetXaxis()->GetBinCenter(bin);

      Printf("Efficiency limit (50%%) is at multiplicity %f", triggerLimit);

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
      nTrigger = TMath::Nint(((Float_t) nTrigger) / 1000) * 1000;

      Printf("Normal rate is %f, downscale: %f, Simulating %lld triggers", normalRate, downScale, nTrigger);
      if (nTrigger == 0)
        continue;

      proj2->FillRandom(proj, nTrigger);

      Int_t removed = 0;
      for (Int_t bin=1; bin<proj2->GetNbinsX(); ++bin)
        if (proj2->GetBinContent(bin) < 5)
        {
          removed += (Int_t) proj2->GetBinContent(bin);
          proj2->SetBinContent(bin, 0);
        }

      Printf("Removed %d events", removed);

      if (relative)
        proj2->Scale(1.0 / nTrigger * proj->Integral(1, 1000));

      proj2->SetLineColor(colors[currentCut]);
      proj2->SetMarkerStyle(markers[currentCut]);
      proj2->SetMarkerColor(colors[currentCut]);

      if (relative || currentCut > 0) {
        proj2->DrawCopy("SAME P");
      } else
        proj2->DrawCopy(" P");

      TString eventStr;
      if (nTrigger > 1e6)
      {
        eventStr.Form("%lld M", nTrigger / 1000 / 1000);
      }
      else if (nTrigger > 1e3)
      {
        eventStr.Form("%lld K", nTrigger / 1000);
      }
      else
      	eventStr.Form("%lld", nTrigger);

      TString triggerStr;
      if (cut == 0)
      {
      	triggerStr = "minimum bias";
      }
      else
      	triggerStr.Form("FO > %d chips", cut);

      legend2->AddEntry(proj2, Form("%s evts, %s", eventStr.Data(), triggerStr.Data()));

      if (triggerLimit > 1)
      {
        TLine* line = new TLine(triggerLimit, proj2->GetMinimum(), triggerLimit, proj2->GetMaximum());
        line->SetLineColor(colors[currentCut]);
        line->Draw();
      }
    }

    legend2->Draw();

    canvas2->SaveAs(Form("%s.gif", canvas2->GetName()));
    canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  }
}

void AliHighMultiplicitySelector::Contamination()
{
  //
  //

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Contamination();

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  // rate = L * sigma
  // sigma ~ 80 mb (Pythia 14 TeV)
  // 10^28 lum --> 8e2 Hz
  // 10^31 lum --> 8e5 Hz
  Double_t rates[] = { 8e2, 8e3, 8e4, 8e5 };

  Int_t nCuts = 4;
  Int_t cuts[] = { 104, 134, 154, 170 };

  // put to 2 for second layer
  for (Int_t i=0; i<1; ++i)
  {
    if (!xSections[i])
      continue;

    // relative x-section, once we have a collision
    xSections[i]->Scale(1.0 / xSections[i]->Integral());

    Int_t max = xSections[i]->GetNbinsX();
    max = 500;

    Float_t* xSection = new Float_t[max];
    for (Int_t mult = 0; mult < max; mult++)
      xSection[mult] = xSections[i]->GetBinContent(mult+1);

    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;

    TGraph* graph = new TGraph;

    for (Int_t currentCut = 0; currentCut<nCuts; ++currentCut)
    {
      Int_t cut = cuts[currentCut];
      Double_t rate = rates[currentCut];
      //Double_t rate = rates[3];

      // coll. in 100 ns window
      Double_t windowSize = 100e-9;
      //Double_t windowSize = 25e-9;
      Double_t collPerWindow = windowSize * rate;
      Printf("coll/window = %f", collPerWindow);
      Double_t windowsPerSecond = 1.0 / windowSize;

      TH1* triggerEffHist = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");
      Float_t* triggerEff = new Float_t[max];
      for (Int_t mult = 0; mult < max; mult++)
        triggerEff[mult] = triggerEffHist->GetBinContent(mult+1);

      Double_t triggerRate = 0;
      for (Int_t mult = 0; mult < max; mult++)
        triggerRate += xSection[mult] * triggerEff[mult];

      triggerRate *= TMath::Poisson(1, collPerWindow) * windowsPerSecond;

      Printf("Rate for 1 collision is %f Hz", triggerRate);

      Double_t triggerRate2 = 0;
      for (Int_t mult = 0; mult < max; mult++)
        for (Int_t mult2 = mult; mult2 < max; mult2++)
          if (mult+mult2 < max)
            triggerRate2 += ((mult2 > mult) ? 2. : 1.) * xSection[mult] * xSection[mult2] * triggerEff[mult+mult2];

      triggerRate2 *= TMath::Poisson(2, collPerWindow) * windowsPerSecond;

      Printf("Rate for 2 collisions is %f Hz --> %.1f%%", triggerRate2, triggerRate2 / triggerRate * 100);

      Double_t triggerRate3 = 0;

      for (Int_t mult = 0; mult < max; mult++)
        for (Int_t mult2 = mult; mult2 < max-mult; mult2++)
          for (Int_t mult3 = 0; mult3 < max-mult-mult2; mult3++)
            //if (mult+mult2+mult3 < max)
              triggerRate3 += ((mult2 > mult) ? 2. : 1.) * xSection[mult] * xSection[mult2] * xSection[mult3] * triggerEff[mult+mult2+mult3];

      triggerRate3 *= TMath::Poisson(3, collPerWindow) * windowsPerSecond;
      //triggerRate3 *= collPerWindow * collPerWindow * rate;

      Printf("Rate for 3 collisions is %f Hz --> %.1f%%", triggerRate3, triggerRate3 / triggerRate * 100);

      Float_t totalContamination = (triggerRate2 + triggerRate3) / (triggerRate + triggerRate2 + triggerRate3);

      Printf("Total contamination is %.1f%%", totalContamination * 100);

      graph->SetPoint(graph->GetN(), cut, totalContamination);

      continue;

      Double_t triggerRate4 = 0;
      for (Int_t mult = 0; mult < max; mult++)
        for (Int_t mult2 = mult; mult2 < max-mult; mult2++)
          for (Int_t mult3 = 0; mult3 < max-mult-mult2; mult3++)
            for (Int_t mult4 = 0; mult4 < max-mult-mult2-mult3; mult4++)
              //if (mult+mult2+mult3+mult4 < max)
                triggerRate4 += ((mult2 > mult) ? 2. : 1.) * xSection[mult] * xSection[mult2] * xSection[mult3] * xSection[mult4] * triggerEff[mult+mult2+mult3+mult4];

      //triggerRate4 *= collPerWindow * collPerWindow * collPerWindow * rate;
      triggerRate4 *= TMath::Poisson(4, collPerWindow) * windowsPerSecond;

      Printf("Rate for 4 collisions is %f Hz --> %.1f%%", triggerRate4, triggerRate4 / triggerRate * 100);

      // general code for n collisions follows, however much slower...
      /*
      const Int_t maxdepth = 4;
      for (Int_t depth = 1; depth <= maxdepth; depth++) {
        Double_t triggerRate = 0;

        Int_t m[maxdepth];
        for (Int_t d=0; d<maxdepth; d++)
          m[d] = 0;

        while (m[0] < max) {
          Double_t value = 1;
          Int_t sum = 0;
          for (Int_t d=0; d<depth; d++) {
            value *= xSection[m[d]];
            sum += m[d];
          }

          if (sum < max) {
            value *= triggerEff[sum];
            triggerRate += value;
          }

          Int_t increase = depth-1;
          ++m[increase];
          while (m[increase] == max && increase > 0) {
            m[increase] = 0;
            --increase;
            ++m[increase];
          }
        }

        triggerRate *= rate * TMath::Power(collPerWindow, depth - 1);

        Printf("Rate for %d collisions is %f Hz", depth, triggerRate);
      }*/
    }

    new TCanvas; graph->Draw("AP*");
  }
}

void AliHighMultiplicitySelector::Contamination2()
{
  //
  // produces a spectrum created with N triggers
  // number of triggers and thresholds for the moment fixed
  //

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Contamination2();

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  Int_t nCuts = 4;
  Int_t cuts[] = { 104, 134, 154, 170 };

  new TCanvas;

  Int_t colors[] = { 2, 3, 4, 6, 7, 8 };
  Int_t markers[] = { 7, 2, 4, 5, 6, 27 };

  // put to 2 for second layer
  for (Int_t i=0; i<1; ++i)
  {
    if (!xSections[i])
      continue;

    // relative x-section, once we have a collision
    xSections[i]->Scale(1.0 / xSections[i]->Integral());

    Int_t max = xSections[i]->GetNbinsX();
    max = 500;

    Float_t* xSection = new Float_t[max];
    for (Int_t mult = 0; mult < max; mult++)
      xSection[mult] = xSections[i]->GetBinContent(mult+1);

    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;

    for (Int_t currentCut = 0; currentCut<nCuts; ++currentCut)
    {
      TGraph* graph = new TGraph;
      graph->SetMarkerColor(colors[currentCut]);
      graph->SetMarkerStyle(markers[currentCut]);

      Int_t cut = cuts[currentCut];

      TH1* triggerEffHist = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");
      Float_t* triggerEff = new Float_t[max];
      for (Int_t mult = 0; mult < max; mult++)
        triggerEff[mult] = triggerEffHist->GetBinContent(mult+1);

      Double_t triggerRate = 0;
      for (Int_t mult = 0; mult < max; mult++)
        triggerRate += xSection[mult] * triggerEff[mult];

      Printf("Raw value for 1 collision is %e", triggerRate);

      Double_t triggerRate2 = 0;
      for (Int_t mult = 0; mult < max; mult++)
        for (Int_t mult2 = mult; mult2 < max; mult2++)
          if (mult+mult2 < max)
            triggerRate2 += ((mult2 > mult) ? 2. : 1.) * xSection[mult] * xSection[mult2] * triggerEff[mult+mult2];

      Printf("Raw value for 2 collisions is %e", triggerRate2);

      for (Double_t doubleRate = 0; doubleRate <= 0.3; doubleRate += 0.005)
      {
        Float_t totalContamination = (triggerRate2 * doubleRate) / (triggerRate + triggerRate2 * doubleRate);

        //Printf("Total contamination is %.1f%%", totalContamination * 100);

        graph->SetPoint(graph->GetN(), doubleRate, totalContamination);
      }

      graph->Draw((currentCut == 0) ? "A*" : "* SAME");
      graph->GetXaxis()->SetRangeUser(0, 1);
    }
  }
}

void AliHighMultiplicitySelector::Contamination3()
{
  //
  // draws the contamination as function of treshold depending on a number a set of input MC and rate parameters
  //

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Contamination3();

  */

  // output file
  TFile* output = TFile::Open("contamination3.root", "RECREATE");

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx_10TeV.root");
  if (!file)
    return;
    
  TCanvas* c = new TCanvas;
  c->SetGridx();
  c->SetGridy();
  
  TLegend* legend = new TLegend(0.7, 0.2, 1, 0.5);
  legend->SetNColumns(2);
  
  TH2* dummy = new TH2F("dummy", ";Layer 1 Threshold;Contamination", 100, 95, 255, 100, 0, 1);
  dummy->SetStats(kFALSE);
  dummy->Draw();
    
  for (Int_t mc = 0; mc < 6; mc++)
  {
    TH1* xSections[2];
    TString str;
    str.Form("xSection2Ex_%d_%d", mc/3, mc%3);
    Printf("%s", str.Data());
    file->cd();
    xSections[0] = dynamic_cast<TH1*> (gFile->Get(str));
    //xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));
  
    // prob for a collision in a bunch crossing
    Int_t nRates = 1;
    //Float_t rates[] = {0.02, 0.05, 0.1, 0.15, 0.2};
    Float_t rates[] = {0.0636};
    
    // bunch crossing rate in Hz
    Float_t bunchCrossingRate = 24. * 11245.5;
  
    Int_t colors[] = { 2, 3, 4, 6, 7, 8 };
    Int_t markers[] = { 7, 2, 4, 5, 6, 27 };
  
    // put to 2 for second layer
    for (Int_t i=0; i<1; ++i)
    {
      if (!xSections[i])
        continue;
  
      // relative x-section, once we have a collision
      xSections[i]->Scale(1.0 / xSections[i]->Integral());
  
      Int_t max = xSections[i]->GetNbinsX();
      max = 500;
  
      Float_t* xSection = new Float_t[max];
      for (Int_t mult = 0; mult < max; mult++)
        xSection[mult] = xSections[i]->GetBinContent(mult+1);
  
      TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;
  
      for (Int_t currentRate = 0; currentRate<nRates; ++currentRate)
      {
        TGraph* graph = new TGraph;
        graph->SetMarkerColor(colors[currentRate]);
        graph->SetMarkerStyle(markers[currentRate]);
 
	TGraph* graph2 = new TGraph;
 
        Float_t rate = rates[currentRate];
  
        Double_t singleRate = TMath::Poisson(1, rate);
        Double_t doubleRate = TMath::Poisson(2, rate);
        Double_t tripleRate = TMath::Poisson(3, rate);
  
        Printf("single = %f, double = %f, triple = %f", singleRate, doubleRate, tripleRate);
    
        for (Int_t cut = 100; cut <= 251; cut += 10)
        {
          Printf("Cut at %d", cut);
  
          TH1* triggerEffHist = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");
          Float_t* triggerEff = new Float_t[max];
          for (Int_t mult = 0; mult < max; mult++)
            triggerEff[mult] = triggerEffHist->GetBinContent(mult+1);
    
          Double_t triggerRate = 0;
          for (Int_t mult = 0; mult < max; mult++)
            triggerRate += xSection[mult] * triggerEff[mult];
    
          //Printf("  Raw value for 1  collision is %e; Rate: %.1f Hz", triggerRate, triggerRate * singleRate * bunchCrossingRate);
    
          Double_t triggerRate2 = 0;
          for (Int_t mult = 0; mult < max; mult++)
            for (Int_t mult2 = mult; mult2 < max; mult2++)
              if (mult+mult2 < max)
                triggerRate2 += ((mult2 > mult) ? 2. : 1.) * xSection[mult] * xSection[mult2] * triggerEff[mult+mult2];
    
          //Printf("  Raw value for 2 collisions is %e; Rate: %.1f Hz", triggerRate2, triggerRate2 * doubleRate * bunchCrossingRate);
    
          Double_t triggerRate3 = 0;
          for (Int_t mult = 0; mult < max; mult++)
            for (Int_t mult2 = 0; mult2 < max; mult2++)
              for (Int_t mult3 = 0; mult3 < max; mult3++)
                if (mult+mult2+mult3 < max)
                  triggerRate3 += xSection[mult] * xSection[mult2] * xSection[mult3] * triggerEff[mult+mult2+mult3];
          
          //Printf("  Raw value for 3 collisions is %e; Rate: %.1f Hz", triggerRate3, triggerRate3 * tripleRate * bunchCrossingRate);
          
          Printf("  Rates: %.1f Hz; %.1f Hz; %.1f Hz", triggerRate * singleRate * bunchCrossingRate, triggerRate2 * doubleRate * bunchCrossingRate, triggerRate3 * tripleRate * bunchCrossingRate);
          
          Float_t totalTrigger = (triggerRate * singleRate + triggerRate2 * doubleRate + triggerRate3 * tripleRate);
          
          Printf("  Total trigger rate: %.1f Hz", totalTrigger * bunchCrossingRate);
          
          //if (totalTrigger * bunchCrossingRate > 200)
          //  continue;
    
          Float_t totalContamination = (triggerRate2 * doubleRate + triggerRate3 * tripleRate) / totalTrigger;
          //if (totalContamination > 0.99)
	  //  break;
    
          Printf("  Total contamination is %.1f%%", totalContamination * 100);
    
          graph->SetPoint(graph->GetN(), cut, totalContamination);
	  graph2->SetPoint(graph->GetN(), cut, totalTrigger * bunchCrossingRate);
        }
  
        graph->SetMarkerStyle(mc+20);
        graph->SetMarkerColor(currentRate+1);
        graph->Draw("P SAME");
        graph->GetXaxis()->SetTitle("Layer 1 threshold");
        graph->GetYaxis()->SetTitle("Contamination");
        graph->GetYaxis()->SetRangeUser(0, 1);
        
        if (currentRate == 0)
        {
          const char* legendLabel[] = { "Pythia Slope 1", "Pythia Slope 2", "Pythia Slope 3", "Phojet Slope 1", "Phojet Slope 2", "Phojet Slope 3" };
          legend->AddEntry(graph, legendLabel[mc], "P");
        }

	output->cd();
	graph->Write(Form("%s_%d_cont", str.Data(), currentRate));
	graph2->Write(Form("%s_%d_rate", str.Data(), currentRate));
      }
    }
  }
  
  output->Close();
  
  legend->Draw();
}

void AliHighMultiplicitySelector::Contamination_Reach()
{
  // plot the multiplicity reach based on the output from Contamination3()
  // for each rate case and each MC, a certain number of events is required starting counting from the highest multiplicity
  // note that the reach of different MC cannot be compared with each other

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->Contamination_Reach();

  */

  TCanvas* c = new TCanvas("c", "c", 800, 600);
  c->Divide(2, 3);
  
  // prob for a collision in a bunch crossing
  Int_t nRates = 1;
  //Float_t rates[] = {0.02, 0.05, 0.1, 0.15, 0.2};
  Float_t rates[] = {0.0636};

  // bunch crossing rate in Hz
  Float_t bunchCrossingRate = 24. * 11245.5;

  TH2* dummy = new TH2F("dummy", ";Coll/bunch crossing;multiplicity reach", 100, 0, 0.3, 100, 50, 350);
  //TH2* dummy = new TH2F("dummy", ";Coll/bunch crossing;fractional cross-section", 100, 0, 0.3, 1000, 1e-6, 0.1);
  dummy->SetStats(kFALSE);

  const char* legendLabel[] = { "Pythia Slope 1", "Pythia Slope 2", "Pythia Slope 3", "Phojet Slope 1", "Phojet Slope 2", "Phojet Slope 3" };

  TFile* mcFile = TFile::Open("crosssectionEx_10TeV.root");
  TFile* contFile = TFile::Open("contamination3.root");

  // for comparison: how many MB events can one take at the same time
  Int_t mbEvents = 2e6 * 500;

  for (Int_t mc=0; mc<6; mc++)
  {
    mcFile->cd();
    TH1* mcHist = (TH1*) gFile->Get(Form("xSection2Ex_%d_%d", mc/3, mc%3));
    mcHist->Scale(1.0 / mcHist->Integral());
    
    c->cd(mc+1);//->SetLogy();
    c->SetGridx();
    c->SetGridy();
    dummy->Draw();
    
    Int_t color = 0;
    for (Int_t requiredEvents = 300; requiredEvents <= 3000000; requiredEvents *= 10)
    {
      TGraph* reach = new TGraph;

      color++;
      if (color == 5)
        color++;
      
      Float_t requiredRate = (Float_t) requiredEvents / 1e6;
      Printf("Required rate is %f", requiredRate);
    
      // find reach without trigger
      Int_t mbReach = 1000;
      while (mcHist->Integral(mcHist->FindBin(mbReach), mcHist->GetNbinsX()) < (Float_t) requiredEvents / mbEvents && mbReach > 1)
        mbReach--;
      Printf("MB reach is %d with %f events", mbReach, mcHist->Integral(mcHist->FindBin(mbReach), mcHist->GetNbinsX()) * mbEvents);
      
      for (Int_t rate=0; rate<nRates; rate++)
      {
        contFile->cd();
        TGraph* cont = (TGraph*) gFile->Get(Form("xSection2Ex_%d_%d_%d_cont", mc/3, mc%3, rate));
        TGraph* rateh = (TGraph*) gFile->Get(Form("xSection2Ex_%d_%d_%d_rate", mc/3, mc%3, rate));
  
        Double_t singleRate = TMath::Poisson(1, rates[rate]);
        Double_t totalCollRate = singleRate * bunchCrossingRate;
        Printf("collisions/bc: %f; coll. rate: %f", singleRate, totalCollRate);
  
        // find 200 Hz limit
        Int_t low = 100;
        while (rateh->Eval(low) > 200)
          low++;
  
        // find contamination limit
        Int_t high = 100;
        while (cont->Eval(high) < 0.9 && high < 250)
          high++;
  
        Printf("MC %d Rate %d: Acceptable threshold region is %d to %d", mc, rate, low, high);
        // find reachable multiplicity; include contamination in rate calculation
        // TH1* triggerEffHist = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");
  
        // trigger efficiency as function of multiplicity in range mult <= ... <= high
        //new TCanvas; fMvsL1->Draw("COLZ");
        TH1* triggerEffHist = (TH1*) GetTriggerEfficiency(fMvsL1, low, high);
        
        //new TCanvas; triggerEffHist->DrawCopy();
        triggerEffHist->Multiply(mcHist);
        
        Float_t fractionXSection = triggerEffHist->Integral();
        Printf("The fraction of the cross-section is %f", fractionXSection);
        
        //new TCanvas; triggerEffHist->DrawCopy();
        triggerEffHist->Scale(totalCollRate);
        //new TCanvas; triggerEffHist->DrawCopy(); gPad->SetLogy();
        
        Float_t achievedRate = 0;
        Int_t mult = 1000;
        while (1)
        {
          achievedRate = triggerEffHist->Integral(triggerEffHist->FindBin(mult), triggerEffHist->GetNbinsX());
  
          if (achievedRate >= requiredRate)
            break;
      
          if (mult == 1)
            break;
      
          mult--;
        } 
  
        Printf("Achieved rate %f above multiplicity %d", achievedRate, mult);
        
        if (achievedRate < requiredRate)
        {
          Printf("Achieved rate too low");
          continue;
        }
  
        reach->SetPoint(reach->GetN(), rates[rate], mult);
        //reach->SetPoint(reach->GetN(), rates[rate], fractionXSection);
        
      
        //return;
      }
  
      reach->SetMarkerColor(color);
      reach->Draw("SAME*");
      
      TLine* line = new TLine(0, mbReach, 0.3, mbReach);
      line->SetLineColor(color);
      //line->Draw();
    }    
  
    TText* text = new TText;
    text->DrawText(0.2, 325, legendLabel[mc]);
    //return;
  }
}

void AliHighMultiplicitySelector::DrawHistograms()
{
  // draws the histograms

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

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx++
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

  TLine* line = new TLine(fMvsL1->GetXaxis()->GetXmin(), 150, fMvsL1->GetXaxis()->GetXmax(), 150);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw();

  canvas->SaveAs("L1NoCurveCut.gif");
  canvas->SaveAs("L1NoCurveCut.eps");

  return;

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

TGraph* AliHighMultiplicitySelector::IntFractRate()
{
  // A plot which shows the fractional rate above any threshold
  // as function of threshold (i.e. the integral of dSigma/dN as function of
  // N, normalised to 1 for N=0)

  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->IntFractRate();
  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return 0;

  TH1* xSection;
  xSection = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));

  TGraph* result = new TGraph;

  for (Int_t threshold = 0; threshold < 300; threshold += 2)
  {
    TH1* proj = GetXSectionCut(xSection, fMvsL1, threshold);

    //new TCanvas; proj->DrawCopy();

    Double_t integral = proj->Integral();

    Printf("Cut at %d, integral is %e", threshold, integral);

    result->SetPoint(result->GetN(), threshold, integral);
  }

  TCanvas* canvas = new TCanvas("IntFractRate", "IntFractRate", 600, 400);
  gPad->SetLogy();

  result->Draw("A*");
  result->GetXaxis()->SetTitle("threshold");
  result->GetYaxis()->SetTitle("integrated fractional rate above threshold");

  canvas->SaveAs(Form("%s.gif", canvas->GetName()));

  return result;
}

void AliHighMultiplicitySelector::MBComparison()
{
  //
  // finds the threshold from which onwards the number of found events above N times the mean
  // is higher using a high mult. trigger than just triggering with MB
  //

  /*

  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
  .L AliHighMultiplicitySelector.cxx+g
  x = new AliHighMultiplicitySelector();
  x->ReadHistograms("highmult_hijing100k.root");
  x->MBComparison();

  */

  // get x-sections
  TFile* file = TFile::Open("crosssectionEx.root");
  if (!file)
    return;

  TH1* xSections[2];
  xSections[0] = dynamic_cast<TH1*> (gFile->Get("xSection2Ex"));
  xSections[1] = dynamic_cast<TH1*> (gFile->Get("xSection15Ex"));

  // rate = L * sigma
  // sigma ~ 80 mb (Pythia 14 TeV)
  // 10^28 lum --> 8e2 Hz
  // 10^31 lum --> 8e5 Hz
  Int_t nRates = 4;
  Double_t rates[] = { 8e2, 8e3, 8e4, 8e5 };

  // threshold in number of fired chips for corresponding rate
  //Int_t cuts[] = { 104, 134, 154, 170 }; // values for 20 Hz
  Int_t cuts[] = { 82, 124, 147, 164 };    // values for 50 Hz

  // bandwidth, fractions (for MB, high mult.)
  Float_t bandwidth = 1e3;
  Float_t fractionMB = 0.5;
  Float_t fractionHM = 0.05;

  // different limits to define "interesting events"
  Int_t nLimits = 9;
  Int_t limits[] = { 0, 1, 2, 4, 6, 7, 8, 9, 10 };

  // put to 2 for second layer
  for (Int_t i=0; i<1; ++i)
  {
    if (!xSections[i])
      continue;

    TH1* xSection = xSections[i];
    TH2* fMvsL = (i == 0) ? fMvsL1: fMvsL2;

    // relative x-section, once we have a collision
    xSection->Scale(1.0 / xSection->Integral());

    xSection->SetStats(kFALSE);
    xSection->SetTitle(""); //(i == 0) ? "SPD Layer 1" : "SPD Layer 2");
    xSection->GetXaxis()->SetTitle(Form("true multiplicity in |#eta| < %.1f", (i == 0) ? 2.0 : 1.5));
    xSection->GetXaxis()->SetRangeUser(0, (i == 0) ? 450 : 350);
    xSection->GetYaxis()->SetTitleOffset(1.2);

    TCanvas* canvas = new TCanvas("MBComparison", "MBComparison", 1000, 800);
    canvas->Divide(3, 3);

    for (Int_t currentLimit = 0; currentLimit<nLimits; currentLimit++)
    {
      // limit is N times the mean
      Int_t limit = (Int_t) (xSection->GetMean() * limits[currentLimit]);
      if (limit < 1)
        limit = 1;

      TGraph* graphMB = new TGraph;
      graphMB->SetTitle(Form("Events with %d times above <n> (i.e. n >= %d)", limits[currentLimit], limit));
      graphMB->SetMarkerStyle(20);

      TGraph* graphBoth = new TGraph;
      graphBoth->SetMarkerStyle(21);

      Float_t min = bandwidth;
      Float_t max = 0;

      for (Int_t current = 0; current<nRates; ++current)
      {
        Float_t rate = rates[current];
        Int_t cut = cuts[current];

        TH1* triggerEff = (TH1*) GetTriggerEfficiency(fMvsL, cut)->Clone("triggerEff");
        TH1* proj = GetXSectionCut(xSection, fMvsL, cut);

        Float_t downScaleMB1 = rate / bandwidth;
        if (downScaleMB1 < 1)
          downScaleMB1 = 1;

        Float_t downScaleMB2 = rate / (bandwidth * fractionMB);
        if (downScaleMB2 < 1)
          downScaleMB2 = 1;

        Float_t downScaleHM = rate * proj->Integral(1, xSection->GetNbinsX()) / (bandwidth * fractionHM);
        if (downScaleHM < 1)
          downScaleHM = 1;

        Float_t rateMB1 = rate / downScaleMB1 * xSection->Integral(limit, xSection->GetNbinsX());
        Float_t rateMB2 = rate / downScaleMB2 * xSection->Integral(limit, xSection->GetNbinsX());
        Float_t rateHM = rate / downScaleHM * proj->Integral(limit, xSection->GetNbinsX());
        Float_t combinedRate = rateMB2 + rateHM;

        graphMB->SetPoint(graphMB->GetN(), rate, rateMB1);
        graphBoth->SetPoint(graphBoth->GetN(), rate, combinedRate);

        min = TMath::Min(min, TMath::Min(rateMB1, combinedRate));
        max = TMath::Max(min, TMath::Max(rateMB1, combinedRate));

        Printf("The rates for events with %d times above <n> (i.e. n >= %d) at a rate of %.2e Hz is:", limits[currentLimit], limit, rate);
        Printf("   %.2e Hz in MB-only mode", rateMB1);
        Printf("   %.2e Hz = %.2e Hz + %.2e Hz in MB + high mult. mode", combinedRate, rateMB2, rateHM);

        Printf("   The downscale factors are: %.2f %.2f %.2f", downScaleMB1, downScaleMB2, downScaleHM);

        Int_t triggerLimit = 0;
        for (Int_t bin = 1; bin <= triggerEff->GetNbinsX(); bin++)
          if (triggerEff->GetBinContent(bin) < 0.5)
            triggerLimit = (Int_t) triggerEff->GetXaxis()->GetBinCenter(bin);

        Printf("   Efficiency limit (50%%) is at multiplicity %d", triggerLimit);
        Float_t fractionGood = proj->Integral(triggerLimit, proj->GetNbinsX()) / proj->Integral();
        Printf("   %.2f %% of the events are above the trigger limit", fractionGood * 100);

        if (triggerLimit > limit)
          Printf("   WARNING: interesting events also counted inside the trigger limit");

        Printf(" ");
      }

      canvas->cd(currentLimit+1)->SetLogx();
      canvas->cd(currentLimit+1)->SetLogy();

      graphMB->Draw("AP");
      graphBoth->Draw("P SAME");

      graphMB->GetYaxis()->SetRangeUser(0.5 * min, 2 * max);
      graphMB->GetXaxis()->SetTitle("Raw rate in Hz");
      graphMB->GetYaxis()->SetTitle("Event rate in Hz");
    }
  }
}
