// $Id$
//
// Jet extractor task
// Stores jets (e.g. as AliEmcalJet objects) in trees

// Author: R. Haake

#include "AliParticleContainer.h"
#include "AliEmcalJet.h"
#include "AliAODTrack.h"
#include "TRandom3.h"
#include "AliAnalysisTaskEmcalJetExtractor.h"


ClassImp(AliAnalysisTaskEmcalJetExtractor)
ClassImp(AliBasicJet)
ClassImp(AliBasicJetConstituent)

//________________________________________________________________________
AliBasicJet::~AliBasicJet() 
{
// dummy destructor
}

//________________________________________________________________________
AliBasicJetConstituent::~AliBasicJetConstituent() 
{
// dummy destructor
}


//________________________________________________________________________
AliAnalysisTaskEmcalJetExtractor::AliAnalysisTaskEmcalJetExtractor() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetExtractor", kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetBuffer(0),
  fJetsOutput(0),
  fCounter(0),
  fRandom(0),
  fIsExtractionDefined(0),
  fExtractionType(0),
  fExtractionCriterium(0),
  fExtractionMinPt(0),
  fExtractionMaxPt(0),
  fExtractionPercentage(1.0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetExtractor::AliAnalysisTaskEmcalJetExtractor(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetBuffer(0),
  fJetsOutput(0),
  fCounter(0),
  fRandom(0),
  fIsExtractionDefined(0),
  fExtractionType(0),
  fExtractionCriterium(0),
  fExtractionMinPt(0),
  fExtractionMaxPt(0),
  fExtractionPercentage(1.0)
{
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetExtractor::~AliAnalysisTaskEmcalJetExtractor()
{
  if(fRandom) delete fRandom;
  if(fJetsOutput) delete fJetsOutput;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetExtractor::DefineExtraction(Int_t type, Int_t criterium, Double_t minPt, Double_t maxPt, Double_t percentage)
{
  if(!fIsExtractionDefined)
  {
    fIsExtractionDefined = kTRUE;
    fExtractionType = type;
    fExtractionCriterium = criterium;
    fExtractionMinPt = minPt;
    fExtractionMaxPt = maxPt;
    fExtractionPercentage = percentage;
  }
  else
    AliFatal("Tried to define extraction twice -- aborting.");
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetExtractor::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  fRandom = new TRandom3(0);
  // ### Basic container settings
  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles connected to jets
    fTracksCont       = dynamic_cast<AliTrackContainer*>(fJetsCont->GetParticleContainer());
  } else {        //no jets, just analysis tracks
    fTracksCont       = dynamic_cast<AliTrackContainer*>(GetTrackContainer(0));
  }
  if(fTracksCont) fTracksCont->SetClassName("AliAODTrack");

  // Histograms
  AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");

  PrintSettings();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetExtractor::PrintSettings()
{
  std::cout << std::endl << Form("########### Printing settings for task %s", GetName()) << std::endl;
  std::cout << "Jet container cuts:\n";
  fJetsCont->PrintCuts();
  std::cout << "---------------------------------------\n";
  std::cout << "Particle container: " << fTracksCont->GetName() << std::endl;
  std::cout << "Track filter: " << fTracksCont->GetTrackFilterType() << " (0 - no filter, 1 - custom, 2 - hybrid, 3 - TPC only)\n";
  std::cout << "---------------------------------------\n";
  std::cout << "Extraction type: " << fExtractionType << " (0 - AliEmcalJet, 1 - AliBasicJet, 2 - AliBasicJet w/constituents)\n";
  std::cout << "Extraction criterium: " << fExtractionCriterium << " (0 - MinBias, 1 - Signal jets, 2 - Background jets)\n";
  std::cout << "Extraction min pT (after background correction if rho given): " << fExtractionMinPt << "\n";
  std::cout << "Extraction max pT (after background correction if rho given): " << fExtractionMaxPt << "\n";
  std::cout << "Extraction percentage (rest is thrown away): " << fExtractionPercentage << "\n";
  std::cout << "###########" << std::endl << std::endl;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetExtractor::FillHistograms()
{
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetExtractor::ExecOnce() {
  AliAnalysisTaskEmcalJet::ExecOnce();
  if (fJetsCont && fJetsCont->GetArray() == 0) 
  {
    fJetsCont = 0;
    AliFatal("No jet container found.");
  }
  if (fTracksCont && fTracksCont->GetArray() == 0)
  {
    fTracksCont = 0;
    AliFatal("No track container found.");
  }
  if (!fIsExtractionDefined)
  {
    AliFatal("Extraction conditions not defined. Use DefineExtraction().");
  }


  fJetsOutput = new TTree("ExtractedJets", "ExtractedJets");
  if(fExtractionType==0) // AliEmcalJet
    fJetsOutput->Branch("Jets", "AliEmcalJet", &fJetBuffer, 1000);
  else if(fExtractionType==1) // AliBasicJet w/o constituents
    fJetsOutput->Branch("Jets", "AliBasicJet", &fJetBuffer, 1000);
  else if(fExtractionType==2) // AliBasicJet w/ constituents
    fJetsOutput->Branch("Jets", "AliBasicJet", &fJetBuffer, 1000);

  TList* list = dynamic_cast<TList*> (GetOutputData(1));
  list->Add(fJetsOutput);

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetExtractor::Run()
{
  Long64_t eventID = InputEvent()->GetHeader()->GetEventIdAsLong();

  fTracksCont->ResetCurrentID();
  AliAODTrack *track = static_cast<AliAODTrack*>(fTracksCont->GetNextAcceptParticle());
  // All tracks plots
  while(track) {
    FillHistogram("hTrackPhiEta", track->Phi(), track->Eta()); 
    track = static_cast<AliAODTrack*>(fTracksCont->GetNextAcceptParticle());
  }

  fJetsCont->ResetCurrentID();
  AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(); 
  while(jet) {
    // Select jets according to extraction criterium
    if( ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) < fExtractionMinPt) || ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) >= fExtractionMaxPt) )
    {  
      jet = fJetsCont->GetNextAcceptJet(); 
      continue;
    }

    // "Minimum bias" criterium
    if(fExtractionCriterium==0)
    {
      // do nothing
    }

    // Discard jets statistically
    if(fRandom->Rndm() >= fExtractionPercentage)
    {  
      jet = fJetsCont->GetNextAcceptJet(); 
      continue;
    }

    FillHistogram("hJetPhiEta", jet->Phi(), jet->Eta()); 

    // Create the jet object that will be saved to the tree
    if(fExtractionType==0) // AliEmcalJet
    {
      fJetBuffer = jet;
      fJetsOutput->Fill();
    }
    else if( (fExtractionType==1) || (fExtractionType==2) )// AliBasicJet
    {
      AliBasicJet basicJet(jet->Eta(), jet->Phi(), jet->Pt(), jet->Charge(), fJetsCont->GetJetRadius(), jet->Area(), fJetsCont->GetRhoVal(), eventID, fCent);
      if(fExtractionType==1)
        fJetBuffer = &basicJet;
      else if(fExtractionType==2) // AliBasicJet w/ constituents
      {
        for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
        {
          AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
          basicJet.AddJetConstituent(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
          fJetBuffer = &basicJet;
        }
      }
      fJetsOutput->Fill();
    }
    jet = fJetsCont->GetNextAcceptJet(); 
  }

  return kTRUE;

}

//########################################################################
// HISTOGRAM HELPER FUNCTIONS
//########################################################################

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetExtractor::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetExtractor::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetExtractor::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEmcalJetExtractor::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEmcalJetExtractor::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}
