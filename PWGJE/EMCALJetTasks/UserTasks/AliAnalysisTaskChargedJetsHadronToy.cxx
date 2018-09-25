/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <iostream>
#include <TRandom3.h>
#include <AliLog.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <AliEmcalJet.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include <AliAODTrack.h>
#include <TClonesArray.h>
#include <AliAODMCParticle.h>

#include <TFile.h>
#include <TGrid.h>
#include <TSystem.h>


#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskChargedJetsHadronToy.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskChargedJetsHadronToy)
/// \endcond
//_____________________________________________________________________________________________________
AliAnalysisTaskChargedJetsHadronToy::AliAnalysisTaskChargedJetsHadronToy() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskChargedJetsHadronToy", kTRUE), fDistributionMultiplicity(0), fDistributionPt(0), fDistributionEtaPhi(0), fMinCentrality(0), fMaxCentrality(10), fPhysicalPrimariesOnly(kTRUE), fDistributionV2(0), fDistributionV3(0), fDistributionV4(0), fDistributionV5(0), fUseMixedEvent(0), fMixedEvent_Tree(0), fMixedEvent_CurrentFile(0), fMixedEvent_CurrentFileID(-1), fMixedEvent_BaseFolder(""), fMixedEvent_TreeName("ME_tree"), fMixedEvent_CurrentEventID(0), fMixedEvent_NumTotalFiles(30), fBuffer_NumTracks(0), fBuffer_TrackPt(0), fBuffer_TrackPhi(0), fBuffer_TrackEta(0), fBuffer_TrackCharge(0), fInputArrayName(""), fOutputArrayName(""), fInputArray(0), fOutputArray(0), fRandom(), fToyCent(0), fRandomPsi3(0), fRandomPsi4(0), fRandomPsi5(0)
{
  // constructor
  fBuffer_TrackPt = new Float_t[10000];
  fBuffer_TrackPhi = new Float_t[10000];
  fBuffer_TrackEta = new Float_t[10000];
  fBuffer_TrackCharge = new Short_t[10000];
  
}


//_____________________________________________________________________________________________________
AliAnalysisTaskChargedJetsHadronToy::AliAnalysisTaskChargedJetsHadronToy(const char* name) :
  AliAnalysisTaskEmcalJet(name, kTRUE), fDistributionMultiplicity(0), fDistributionPt(0), fDistributionEtaPhi(0), fMinCentrality(0), fMaxCentrality(10), fPhysicalPrimariesOnly(kTRUE), fDistributionV2(0), fDistributionV3(0), fDistributionV4(0), fDistributionV5(0), fUseMixedEvent(0), fMixedEvent_Tree(0), fMixedEvent_CurrentFile(0), fMixedEvent_CurrentFileID(-1), fMixedEvent_BaseFolder(""), fMixedEvent_TreeName("ME_tree"), fMixedEvent_CurrentEventID(0), fMixedEvent_NumTotalFiles(30), fBuffer_NumTracks(0), fBuffer_TrackPt(0), fBuffer_TrackPhi(0), fBuffer_TrackEta(0), fBuffer_TrackCharge(0), fInputArrayName(""), fOutputArrayName(""), fInputArray(0), fOutputArray(0), fRandom(), fToyCent(0), fRandomPsi3(0), fRandomPsi4(0), fRandomPsi5(0)
{
  // constructor
  fBuffer_TrackPt = new Float_t[10000];
  fBuffer_TrackPhi = new Float_t[10000];
  fBuffer_TrackEta = new Float_t[10000];
  fBuffer_TrackCharge = new Short_t[10000];
}

//_____________________________________________________________________________________________________
AliAnalysisTaskChargedJetsHadronToy::~AliAnalysisTaskChargedJetsHadronToy()
{
  // destructor
  if(fDistributionMultiplicity)
    delete fDistributionMultiplicity;
  if(fDistributionPt)
    delete fDistributionPt;
  if(fDistributionEtaPhi)
    delete fDistributionEtaPhi;

  delete fBuffer_TrackPt;
  delete fBuffer_TrackPhi;
  delete fBuffer_TrackEta;
  delete fBuffer_TrackCharge;

}


//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fRandom = new TRandom3(0);
  AddHistogram1D<TH1D>("hTrackPt", "Tracks p_{T} distribution", "", 300, 0., 300., "p_{T} (GeV/c)", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");
  AddHistogram1D<TH1D>("hMultiplicity", "Number of tracks in acceptance vs. centrality", "", 500, 0., 5000., "N tracks","dN^{Events}/dN^{Tracks}");

  PostData(1, fOutput);
}


//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();

  // Check if input array can be loaded
  if(!fInputArrayName.IsNull())
  {
    fInputArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fInputArrayName.Data())));
    if(!fInputArray)
    {
      AliFatal(Form("Input array '%s' not found!", fInputArrayName.Data()));
    }
  }

  // Check if output arrays can be created
  if((InputEvent()->FindListObject(Form("%s", fOutputArrayName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fOutputArrayName.Data()));

  fOutputArray = new TClonesArray("AliAODTrack");
  fOutputArray->SetName(fOutputArrayName.Data());
  fInputEvent->AddObject(fOutputArray);
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronToy::Run()
{
  AssembleEvent();
  CreateQAPlots();
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::AssembleEvent()
{
  fRandomPsi3 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi3
  fRandomPsi4 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi4
  fRandomPsi5 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi5
  fToyCent    = fMinCentrality + (fMaxCentrality-fMinCentrality)*fRandom->Rndm(); // centrality value (flat from selected range)

  // Create the event from the several inputs and run the jet finder

  // ################# 1. Add input tracks (if available)
  Int_t particleCount = 0;
  if(fInputArray)
  {
    for(Int_t iPart=0; iPart<fInputArray->GetEntries(); iPart++)
    {
      // Load AOD tracks or AOD MC particles (for MC gen train)
      AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fInputArray->At(iPart));
      AliAODMCParticle* aodMCParticle = dynamic_cast<AliAODMCParticle*>(fInputArray->At(iPart));

      if(aodTrack)
      {
        new ((*fOutputArray)[particleCount]) AliAODTrack(*aodTrack);
        particleCount++;
      }
      else if(aodMCParticle)
      {
        if (fPhysicalPrimariesOnly && !aodMCParticle->IsPhysicalPrimary())
          continue;
        Float_t trackTheta = 2.*atan(exp(-aodMCParticle->Eta()));
        new ((*fOutputArray)[particleCount]) AliAODTrack();
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPt(aodMCParticle->Pt());
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPhi(aodMCParticle->Phi());
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetTheta(trackTheta); // AliAODTrack cannot set eta directly
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetCharge(aodMCParticle->Charge());
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetLabel(aodMCParticle->GetLabel());
        static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetIsHybridGlobalConstrainedGlobal();
        particleCount++;
      }
    }
  }

  // ################# 2. Create a vertex if there is none (needed by some tasks)
  if(dynamic_cast<AliESDEvent*>(InputEvent()))
  {
    if((dynamic_cast<AliESDEvent*>(InputEvent()))->GetPrimaryVertexTracks())
      if(!(dynamic_cast<AliESDEvent*>(InputEvent()))->GetPrimaryVertexTracks()->GetNContributors())
        static_cast<AliESDEvent*>(fInputEvent)->SetPrimaryVertexTracks(new AliESDVertex(0.,0., 100));
  }
  else if(dynamic_cast<AliAODEvent*>(InputEvent()))
  {
    if( (!(dynamic_cast<AliAODEvent*>(InputEvent()))->GetPrimaryVertex()) || (!(dynamic_cast<AliAODEvent*>(InputEvent()))->GetPrimaryVertex()->GetNContributors()) )
    {
      Double_t* p = new Double_t[3];
      p[0] = 0.; p[1] = 0.; p[2] = 0.; // for backwards compatibility
      AliAODVertex* vertex = new AliAODVertex(p,1.);
      vertex->SetNContributors(100);
      vertex->SetName("PrimaryVertex");
      static_cast<AliAODEvent*>(fInputEvent)->AddVertex(vertex);
    }
  }

  // ################# 3. Create toy event
  if(fUseMixedEvent) // get underlying event from mixed event files
  {
    // if input tree not loaded or index at the end, get next tree file
    if( !fMixedEvent_Tree || (fMixedEvent_CurrentEventID >= fMixedEvent_Tree->GetEntries()) )
    {
      fMixedEvent_Tree = GetNextMixedEventTree();
      if (!fMixedEvent_Tree)
      {
        AliError(Form("Could not get tree %s in file!", fMixedEvent_TreeName.Data()));
        return;
      }
    }

    fBuffer_NumTracks = 0; // Failsafe: Set num tracks to 0 if it is not read by GetEntry
    fMixedEvent_Tree->GetEntry(fMixedEvent_CurrentEventID);

    // Loop over tracks from event
    for(Int_t j=0; j<fBuffer_NumTracks; j++)
    {
      Float_t trackTheta = 2.*atan(exp(-fBuffer_TrackEta[j]));

      // Add basic particle to event
      new ((*fOutputArray)[particleCount]) AliAODTrack();
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPt(fBuffer_TrackPt[j]);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPhi(fBuffer_TrackPhi[j]);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetTheta(trackTheta); // AliAODTrack cannot set eta directly
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetCharge(fBuffer_TrackCharge[j]);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetLabel(100000 + j);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetIsHybridGlobalConstrainedGlobal();
      particleCount++;
    }

    fMixedEvent_CurrentEventID++;
  }
  else // Simple toy
  {
    Int_t multiplicity = (Int_t)fDistributionMultiplicity->GetRandom();
    for(Int_t i=0;i<multiplicity; i++)
    {
      Double_t trackPt = fDistributionPt->GetRandom();
      Double_t trackEta = 0;
      Double_t trackPhi = 0;
      static_cast<TH2*>(fDistributionEtaPhi)->GetRandom2(trackPhi, trackEta);
      Double_t trackTheta = 2.*atan(exp(-trackEta));
      Double_t trackCharge = fRandom->Rndm() - 0.5;

      if(trackCharge>0) trackCharge = 1; else trackCharge = -1;

      // Add flow to particle
      if(fDistributionV2 || fDistributionV3 || fDistributionV4 || fDistributionV5)
        trackPhi = AddFlow(trackPhi, trackPt);


      // Add basic particle to event
      new ((*fOutputArray)[particleCount]) AliAODTrack();
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPt(trackPt);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetPhi(trackPhi);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetTheta(trackTheta); // AliAODTrack cannot set eta directly
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetCharge(trackCharge);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetLabel(100000 + i);
      static_cast<AliAODTrack*>(fOutputArray->At(particleCount))->SetIsHybridGlobalConstrainedGlobal();
      particleCount++;
    }
  }

  
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::CreateQAPlots()
{
  for(Int_t iTrack=0; iTrack<fOutputArray->GetEntries(); iTrack++)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(fOutputArray->At(iTrack));
    FillHistogram("hTrackPt", track->Pt());
    FillHistogram("hTrackPhiEta", track->Phi(), track->Eta());
  }
  FillHistogram("hMultiplicity", fOutputArray->GetEntries());
}

//_____________________________________________________________________________________________________
TTree* AliAnalysisTaskChargedJetsHadronToy::GetNextMixedEventTree() 
{
  // ## Form file name
  fMixedEvent_CurrentFileID = TMath::Nint(fRandom->Rndm()*(fMixedEvent_NumTotalFiles-1));
  TString fileName;
  fileName = Form("%s%d.root", fMixedEvent_BaseFolder.Data(), fMixedEvent_CurrentFileID); 

  AliInfo(Form("Opening mixed event file: %s", fileName.Data()));
  // ## Check if file exists
  if (fileName.BeginsWith("alien://") && !gGrid)
  {
    AliInfo("Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }
  if (gSystem->AccessPathName(fileName)) {
    AliError(Form("File %s does not exist!", fileName.Data()));
    return 0;
  }

  // ## Open mixed event file and get tree object
  if (fMixedEvent_CurrentFile) fMixedEvent_CurrentFile->Close(); 

  fMixedEvent_CurrentFile = TFile::Open(fileName);
  if (!fMixedEvent_CurrentFile || fMixedEvent_CurrentFile->IsZombie())
  {
    AliError(Form("Unable to open file: %s!", fileName.Data()));
    return 0;
  }
  TTree* tree = static_cast<TTree*>(fMixedEvent_CurrentFile->Get(fMixedEvent_TreeName.Data()));
  if(tree)
  {
    tree->SetBranchAddress("NumTracks", &fBuffer_NumTracks);
    tree->SetBranchAddress("Track_Pt", fBuffer_TrackPt);
    tree->SetBranchAddress("Track_Phi", fBuffer_TrackPhi);
    tree->SetBranchAddress("Track_Eta", fBuffer_TrackEta);
    tree->SetBranchAddress("Track_Charge", fBuffer_TrackCharge);
  }

  fMixedEvent_CurrentEventID = 0;

  return tree;
}

//_____________________________________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsHadronToy::AddFlow(Double_t phi, Double_t pt)
{
  // adapted from AliFlowTrackSimple
  Double_t precisionPhi = 1e-10;
  Int_t maxNumberOfIterations  = 200;

  Double_t phi0=phi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;
  Int_t ptBin = 0;

  // Evaluate V2 for track pt/centrality
  Double_t v2 = 0;
  if(fDistributionV2)
  {
    ptBin = fDistributionV2->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV2->GetNbinsX())
      v2 = fDistributionV2->GetBinContent(fDistributionV2->GetNbinsX(), fDistributionV2->GetYaxis()->FindBin(fToyCent));
    else if(ptBin>0)
      v2 = fDistributionV2->GetBinContent(ptBin, fDistributionV2->GetYaxis()->FindBin(fToyCent));
  }

  // Evaluate V3 for track pt/centrality
  Double_t v3 = 0;
  if(fDistributionV3)
  {
    ptBin = fDistributionV3->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV3->GetNbinsX())
      v3 = fDistributionV3->GetBinContent(fDistributionV3->GetNbinsX(), fDistributionV3->GetYaxis()->FindBin(fToyCent));
    else if(ptBin>0)
      v3 = fDistributionV3->GetBinContent(ptBin, fDistributionV3->GetYaxis()->FindBin(fToyCent));
  }

  // Evaluate V4 for track pt/centrality
  Double_t v4 = 0;
  if(fDistributionV4)
  {
    ptBin = fDistributionV4->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV4->GetNbinsX())
      v4 = fDistributionV4->GetBinContent(fDistributionV4->GetNbinsX(), fDistributionV4->GetYaxis()->FindBin(fToyCent));
    else if(ptBin>0)
      v4 = fDistributionV4->GetBinContent(ptBin, fDistributionV4->GetYaxis()->FindBin(fToyCent));
  }

  // Evaluate V5 for track pt/centrality
  Double_t v5 = 0;
  if(fDistributionV5)
  {
    ptBin = fDistributionV5->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV5->GetNbinsX())
      v5 = fDistributionV5->GetBinContent(fDistributionV5->GetNbinsX(), fDistributionV5->GetYaxis()->FindBin(fToyCent));
    else if(ptBin>0)
      v5 = fDistributionV5->GetBinContent(ptBin, fDistributionV5->GetYaxis()->FindBin(fToyCent));
  }

  // Add all v's
  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=phi; //store last value for comparison
    f =  phi-phi0
        +      v2*TMath::Sin(2.*(phi-(fEPV0+(TMath::Pi()/2.))))
        +2./3.*v3*TMath::Sin(3.*(phi-fRandomPsi3))
        +0.5  *v4*TMath::Sin(4.*(phi-fRandomPsi4))
        +0.4  *v5*TMath::Sin(5.*(phi-fRandomPsi5));

    fp =  1.0+2.0*(
           +v2*TMath::Cos(2.*(phi-(fEPV0+(TMath::Pi()/2.))))
           +v3*TMath::Cos(3.*(phi-fRandomPsi3))
           +v4*TMath::Cos(4.*(phi-fRandomPsi4))
           +v5*TMath::Cos(5.*(phi-fRandomPsi5))); //first derivative

    phi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,phi,precisionPhi)) break;
  }

  return phi;
}


//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsHadronToy::FillHistogram(const char * key, Double_t x)
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
inline void AliAnalysisTaskChargedJetsHadronToy::FillHistogram(const char * key, Double_t x, Double_t y)
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
inline void AliAnalysisTaskChargedJetsHadronToy::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
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
template <class T> T* AliAnalysisTaskChargedJetsHadronToy::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
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
template <class T> T* AliAnalysisTaskChargedJetsHadronToy::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
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


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
