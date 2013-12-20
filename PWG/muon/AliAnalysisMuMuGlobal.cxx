#include "AliAnalysisMuMuGlobal.h"

/**
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuGlobal
 *
 * Very simple histogramming analysis for global event properties, like vertex, bunch-crossing,
 * background, etc...
 *
 * Warning : not really ready for prime-time
 *
 */

#include "TH1.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAODTZERO.h"
#include "AliESDEvent.h"
#include "AliAnalysisMuMuBinning.h"
#include <set>
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"

ClassImp(AliAnalysisMuMuGlobal)

//_____________________________________________________________________________
AliAnalysisMuMuGlobal::AliAnalysisMuMuGlobal() : AliAnalysisMuMuBase()
{
  /// ctor
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::DefineHistogramCollection(const char* eventSelection,
                                                      const char* triggerClassName,
                                                      const char* centrality)
{
  /// Actually create the histograms for phyics/triggerClassName
  
  //  AliInfo(Form("%s %s %s %d",eventSelection,triggerClassName,centrality,hasMC));
  
  if (HistogramCollection()->Histo(Form("/%s/%s/%s/Zvertex",eventSelection,triggerClassName,centrality)))
  {
    return;
  }
  
  Double_t xmin = -40;
  Double_t xmax = +40;
  Int_t nbins = GetNbins(xmin,xmax,0.5);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Zvertex","z vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"T0Zvertex","T0 zvertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "RecZvertexVsMCZvertex","Reconstructed vertex (w/ Ncontrib>=1) vs MC vertex",nbins,xmin,xmax,nbins,xmin,xmax);
  
  xmin = -5;
  xmax = 5;
  nbins = GetNbins(xmin,xmax,0.01);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"ZvertexMinusZvertexSPD","Primary vertex z - SPD vertex",nbins,xmin,xmax);
  
  xmin = -1;
  xmax = 50;
  nbins = GetNbins(xmin,xmax,1);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"ZvertexNContributors","z vertex nof contributors",nbins,xmin,xmax);
  
  xmin = -2;
  xmax = 2;
  nbins = GetNbins(xmin,xmax,0.01);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Xvertex","x vertex",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Yvertex","y vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Nevents","number of events",2,-0.5,1.5);
  
  CreateEventHistos(kHistoForMCInput | kHistoForData,eventSelection,triggerClassName,centrality,
                    "VertexType","Type of vertexer used",10,0,10);
  
  
  xmin = 0;
  xmax = 3564;
  nbins = GetNbins(xmin,xmax,1.0);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"BCX","bunch-crossing ids",nbins,xmin-0.5,xmax-0.5);
  
  
  xmin = -30;
  xmax = +30;
  nbins = GetNbins(xmin,xmax,0.1);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V02D","V0C+V0A versus V0A-V0C;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V02DwT0BB","V0C+V0A versus V0A-V0C with T0 BB;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V02DwT0BG","V0C+V0A versus V0A-V0C with T0 background flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V02DwT0PU","V0C+V0A versus V0A-V0C with T0 pile up flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V02DwT0SAT","V0C+V0A versus V0A-V0C with T0 satellite flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
  
  xmin = 0;
  xmax = 600;
  nbins = GetNbins(xmin,xmax,1);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0AMult","V0A multiplicity",nbins,xmin,xmax);
  
  if ( !IsHistogramDisabled("Centrality") )
  {
    TObjArray* centralities = Binning()->CreateBinObjArray("centrality");
    TIter next(centralities);
    AliAnalysisMuMuBinning::Range* r;
    std::set<std::string> estimators;
    
    while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
    {
      estimators.insert(r->Quantity().Data());
    }
    
    std::set<std::string>::const_iterator it;
    
    for ( it = estimators.begin(); it != estimators.end(); ++it )
    {
      TH1* h = new TH1F("Centrality","Centrality",12,-10,110);
      HistogramCollection()->Adopt(Form("/%s/%s/%s",eventSelection,triggerClassName,it->c_str()),h);
    }
    
    delete centralities;
  }
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"PileUpEstimators","pile up estimators",10,0,10);
}


//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::FillHistosForEvent(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality)
{
  // Fill event-wise histograms
  
  if (!IsHistogramDisabled("BCX"))
  {
    Histo(eventSelection,triggerClassName,centrality,"BCX")->Fill(1.0*Event()->GetBunchCrossNumber());
  }
  if (!IsHistogramDisabled("Nevents"))
  {
    Histo(eventSelection,triggerClassName,centrality,"Nevents")->Fill(1.0);
  }
  
  const AliVVertex* vertex = Event()->GetPrimaryVertex();
  const AliVVertex* vertexFromSPD = AliAnalysisMuonUtility::GetVertexSPD(Event());
  
  if ( vertex && vertexFromSPD )
  {
    if ( vertex->GetNContributors() > 0 )
    {
      if (!IsHistogramDisabled("Xvertex"))
      {
        Histo(eventSelection,triggerClassName,centrality,"Xvertex")->Fill(vertex->GetX());
      }
      if (!IsHistogramDisabled("Yvertex"))
      {
        Histo(eventSelection,triggerClassName,centrality,"Yvertex")->Fill(vertex->GetY());
      }
      if (!IsHistogramDisabled("Zvertex"))
      {
        Histo(eventSelection,triggerClassName,centrality,"Zvertex")->Fill(vertex->GetZ());
      }
      if (!IsHistogramDisabled("ZvertexMinusZvertexSPD"))
      {
        Histo(eventSelection,triggerClassName,centrality,"ZvertexMinusZvertexSPD")->Fill(vertexFromSPD->GetZ()-vertex->GetZ());
      }
      if (!IsHistogramDisabled("VertexType"))
      {
        Histo(eventSelection,triggerClassName,centrality,"VertexType")->Fill(vertex->GetTitle(),1.0);
      }
    }
    if (!IsHistogramDisabled("ZvertexNContributors"))
    {
      Histo(eventSelection,triggerClassName,centrality,"ZvertexNContributors")->Fill(vertex->GetNContributors());
    }
  }
  
  if ( AliAnalysisMuonUtility::IsAODEvent(Event()) )
  {
    const AliAODTZERO* tzero = static_cast<const AliAODEvent*>(Event())->GetTZEROData();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      Histo(eventSelection,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0VertexRaw());
    }
  }
  else
  {
    const AliESDTZERO* tzero = static_cast<const AliESDEvent*>(Event())->GetESDTZERO();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      Histo(eventSelection,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0zVertex());
    }
  }
  
  AliVVZERO* vzero = Event()->GetVZEROData();
  
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t x = v0a-v0c;
    Float_t y = v0a+v0c;
    
    if (!IsHistogramDisabled("V0A"))
    {
      Histo(eventSelection,triggerClassName,centrality,"V0AMult")->Fill(vzero->GetMTotV0A());
    }
    if (!IsHistogramDisabled("V02D"))
    {
      Histo(eventSelection,triggerClassName,centrality,"V02D")->Fill(x,y);
    }
    
    Bool_t background,pileup,satellite;
    
    Bool_t tzero = AliAnalysisMuonUtility::EAGetTZEROFlags(Event(),background,pileup,satellite);
    
    if (tzero)
    {
      if ( background )
      {
        if (!IsHistogramDisabled("V02DwT0BG"))
        {
          Histo(eventSelection,triggerClassName,centrality,"V02DwT0BG")->Fill(x,y);
        }
      }
      
      if ( pileup )
      {
        if (!IsHistogramDisabled("V02DwT0PU"))
        {
          Histo(eventSelection,triggerClassName,centrality,"V02DwT0PU")->Fill(x,y);
        }
        
        if ( !IsHistogramDisabled("PileUpEstimators") )
        {
          Histo(eventSelection,triggerClassName,centrality,"PileUpEstimators")->Fill("TZERO",1.0);
        }
      }
      
      if ( satellite )
      {
        if (!IsHistogramDisabled("V02DwT0SAT"))
        {
          Histo(eventSelection,triggerClassName,centrality,"V02DwT0SAT")->Fill(x,y);
        }
      }
      
      if ( !background && !pileup && !satellite )
      {
        if (!IsHistogramDisabled("V02DwT0BB"))
        {
          Histo(eventSelection,triggerClassName,centrality,"V02DwT0BB")->Fill(x,y);
        }
      }
    }
  }
  
  //  /* FIXME : how to properly get multiplicity from AOD and ESD consistently ?
  //   is is doable at all ?
  
  TH1* hpileup = Histo(eventSelection,triggerClassName,centrality,"PileUpEstimators");
  
  
  //  virtual Bool_t  IsPileupFromSPD(Int_t minContributors=3, Double_t minZdist=0.8, Double_t nSigmaZdist=3., Double_t nSigmaDiamXY=2., Double_t nSigmaDiamZ=5.) const;
  
  const Double_t nSigmaZdist=3.0;
  const Double_t nSigmaDiamXY=2.0;
  const Double_t nSigmaDiamZ=5.0;
  
  for ( Int_t minContributors = 3; minContributors <= 6; ++minContributors )
  {
    for ( double minZdist = 0.6; minZdist <= 0.8; minZdist += 0.2 )
    {
      if ( Event()->IsPileupFromSPD(minContributors,minZdist,nSigmaZdist,nSigmaDiamXY,nSigmaDiamZ) )
      {
        hpileup->Fill(Form("SPD:n%dd%d",minContributors,static_cast<Int_t>(10*minZdist)),1);
      }
    }
  }
  
  
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::FillHistosForMCEvent(const char* eventSelection,
                                                 const char* triggerClassName,
                                                 const char* centrality)
{
  // Fill MCEvent-wise histograms
  
  Double_t Zvertex = AliAnalysisMuonUtility::GetMCVertexZ(Event(),MCEvent());
  
  if (!IsHistogramDisabled("Zvertex"))
  {
    MCHisto(eventSelection,triggerClassName,centrality,"Zvertex")->Fill(Zvertex);
  }

  if (!IsHistogramDisabled("RecZvertexVsMCZvertex"))
  {
    const AliVVertex* vertex = Event()->GetPrimaryVertex();
    if  (vertex && vertex->GetNContributors()>0)
    {
      MCHisto(eventSelection,triggerClassName,centrality,"RecZvertexVsMCZvertex")->Fill(Zvertex,vertex->GetZ());
    }
  }
  
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuGlobal::SelectAnyTriggerClass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses) const
{
  /// Accept all trigger classes
  if ( firedTriggerClasses.Length()>0)
  {
    acceptedTriggerClasses = "NOTRIGGERSELECTION";
    return kTRUE;
  }
  return kFALSE;
}

