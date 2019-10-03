#include "AliAnalysisMuMuMCGene.h"

/**
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuMCGene
 *
 * Very simple histogramming of basic distributions (eta, phi, y, pt)
 * for MC generated events for a few selected particle types
 *
 *
 */

#include "TH1.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include <set>
#include "AliMergeableCollection.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisMuonUtility.h"
#include "AliESDEvent.h"
#include "TDatabasePDG.h"
#include "AliLog.h"

ClassImp(AliAnalysisMuMuMCGene)

//_____________________________________________________________________________
AliAnalysisMuMuMCGene::AliAnalysisMuMuMCGene() : AliAnalysisMuMuBase(),
fParticlesOfInterest(),
fPDGCodeOfInterest()
{
  /// ctor

  fParticlesOfInterest.push_back("pi");
  fParticlesOfInterest.push_back("K");
  fParticlesOfInterest.push_back("phi");

  for ( std::vector<std::string>::size_type i = 0; i < fParticlesOfInterest.size(); ++i )
  {
    TString part = fParticlesOfInterest[i].c_str();

    TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(part.Data());

    if (!p)
    {
      part += "+";
      p = TDatabasePDG::Instance()->GetParticle(part.Data());

      if (!p)
      {
        AliError(Form("Don't know what to do with particle = %s. Skipping.",part.Data()));
        continue;
      }

      fPDGCodeOfInterest.insert(p->PdgCode());

      part.ReplaceAll("+","-");
      p = TDatabasePDG::Instance()->GetParticle(part.Data());
      fPDGCodeOfInterest.insert(p->PdgCode());
      part.ReplaceAll("+","-");
      continue;
    }

    fPDGCodeOfInterest.insert(p->PdgCode());
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMCGene::SelectAnyTriggerClass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses) const
{
  if ( firedTriggerClasses.Length()>0)
  {
    acceptedTriggerClasses = "NOTRIGGERSELECTION";
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
AliGenEventHeader* AliAnalysisMuMuMCGene::GetGenEventHeader(const AliVEvent& event) const
{
  if  ( AliAnalysisMuonUtility::IsAODEvent(&event) )
  {
    const AliAODEvent* eventAOD = static_cast<const AliAODEvent*>(&event);

    AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*>(eventAOD->FindListObject(AliAODMCHeader::StdBranchName()));

    if ( !aodMCHeader ) return 0x0;

    TList* lheaders = aodMCHeader->GetCocktailHeaders();

    if ( lheaders->GetEntries() > 1 )
    {
      AliWarning("There is more than one header: The simulation is a cocktail. Returning first one.");
    }
    return static_cast<AliGenEventHeader*>(lheaders->First());
  }
  else if ( event.IsA() == AliESDEvent::Class() )
  {
    return MCEvent()->GenEventHeader();
  }
  else
  {
    AliWarning(Form("event is of class %s",event.ClassName()));
  }
  return 0x0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMCGene::ParticleOfInterest(const AliVParticle& part) const
{
  return fPDGCodeOfInterest.find(part.PdgCode()) != fPDGCodeOfInterest.end();
}

//_____________________________________________________________________________
void AliAnalysisMuMuMCGene::FillHistosForMCEvent(const char* eventSelection,
                                                 const char* triggerClassName,
                                                 const char* centrality)
{
  // Fill MCEvent-wise histograms

  Int_t nMCTracks = MCEvent()->GetNumberOfTracks(); // MC number of MC tracks

  AliVParticle* p(0x0);

  for ( Int_t i = 0; i < nMCTracks ; ++i ) //Loop over generated tracks
  {
    p = MCEvent()->GetTrack(i);

    if (!p->IsPrimary()) continue;

    if ( ParticleOfInterest(*p) )
    {
      TParticlePDG* pdg = TDatabasePDG::Instance()->GetParticle(p->PdgCode());

      TString pname(pdg->GetName());
      pname.ReplaceAll("+","");
      pname.ReplaceAll("-","");

      MCHisto(eventSelection,triggerClassName,centrality,Form("Y%s",pname.Data()))->Fill(p->Y());
      MCHisto(eventSelection,triggerClassName,centrality,Form("Eta%s",pname.Data()))->Fill(p->Eta());

      MCHisto(eventSelection,triggerClassName,centrality,Form("PtY%s",pname.Data()))->Fill(p->Y(),p->Pt());

      if ( TMath::Abs(p->Y()) > 2.0 )
      {
        MCHisto(eventSelection,triggerClassName,centrality,Form("Pt%sYgt2",pname.Data()))->Fill(p->Pt());
      }
      if ( TMath::Abs(p->Y()) < 1.0 )
      {
        MCHisto(eventSelection,triggerClassName,centrality,Form("Pt%sYle1",pname.Data()))->Fill(p->Pt());
      }

      MCHisto(eventSelection,triggerClassName,centrality,Form("Pt%s",pname.Data()))->Fill(p->Pt());

      if ( p->Pt() > 1 )
      {
        MCHisto(eventSelection,triggerClassName,centrality,Form("Y%s1GeV",pname.Data()))->Fill(p->Y());
        MCHisto(eventSelection,triggerClassName,centrality,Form("Eta%s1GeV",pname.Data()))->Fill(p->Eta());
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuMCGene::DefineHistogramCollection(const char* eventSelection,
                                                      const char* triggerClassName,
                                                      const char* centrality,
                                                      Bool_t mix)
{
  /// Actually create the histograms for phyics/triggerClassName

//  AliInfo(Form("%s %s %s %d",eventSelection,triggerClassName,centrality,hasMC));

  if (HistogramCollection()->Histo(Form("/%s/%s/%s/Semaphore",eventSelection,triggerClassName,centrality)))
  {
    return;
  }

  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Semaphore","Semaphore",1,0.0,1.0);


  Double_t xmin = -6;
  Double_t xmax = 6;
  Int_t nbins = GetNbins(xmin,xmax,0.1);

  Double_t ptmin = 0;
  Double_t ptmax = 7;
  Int_t nptbins = GetNbins(ptmin,ptmax,0.1);

  for ( std::vector<std::string>::size_type i = 0; i < fParticlesOfInterest.size(); ++i )
  {
    const std::string& part = fParticlesOfInterest[i];

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Y%s",part.c_str()),Form("Rapidity of %s",part.c_str()),nbins,xmin,xmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Y%s1GeV",part.c_str()),Form("Rapidity of %s with pT > 1 GeV/c",part.c_str()),nbins,xmin,xmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Eta%s",part.c_str()),Form("Pseudo-rapidity of %s",part.c_str()),nbins,xmin,xmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Eta%s1GeV",part.c_str()),Form("Pseudo-rapidity of %s with pT > 1 GeV/c",part.c_str()),nbins,xmin,xmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Pt%s",part.c_str()),Form("pT distribution of %s",part.c_str()),nptbins,ptmin,ptmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Pt%sYgt2",part.c_str()),Form("pT distribution of %s in |Y|>2",part.c_str()),nptbins,ptmin,ptmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("Pt%sYle1",part.c_str()),Form("pT distribution of %s in |Y|<1.0",part.c_str()),nptbins,ptmin,ptmax);

    CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                      Form("PtY%s",part.c_str()),Form("pT vs Y distribution of %s",part.c_str()),nbins,xmin,xmax,nptbins,ptmin,ptmax);

  }
}

