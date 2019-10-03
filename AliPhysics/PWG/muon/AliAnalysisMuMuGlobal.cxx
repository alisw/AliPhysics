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
#include "AliESDUtils.h"

ClassImp(AliAnalysisMuMuGlobal)

//_____________________________________________________________________________
AliAnalysisMuMuGlobal::AliAnalysisMuMuGlobal() : AliAnalysisMuMuBase()
{
  /// ctor
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuGlobal::SelectAnyTriggerClass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses) const
{
  if ( firedTriggerClasses.Length()>0)
  {
    acceptedTriggerClasses = "NOTRIGGERSELECTION";
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::FillHistosForEvent(AliMergeableCollectionProxy& p)
{
  // Fill event-wise histograms
  
  if (!IsHistogramDisabled("BCX"))
  {
    p.Histo("BCX")->Fill(1.0*Event()->GetBunchCrossNumber());
  }
  if (!IsHistogramDisabled("Nevents"))
  {
    p.Histo("Nevents")->Fill(1.0);
  }
  
  if (!IsHistogramDisabled("EventsWOL0inputs"))
  {
    UInt_t l0 = Event()->GetHeader()->GetL0TriggerInputs();
    
    if ( l0 == 0 ) p.Histo("EventsWOL0inputs")->Fill(1.);
  }
  
  const AliVVertex* vertex = Event()->GetPrimaryVertex();
  const AliVVertex* vertexFromSPD = Event()->GetPrimaryVertexSPD();
  
  if ( vertex )
  {
    if ( vertex->GetNContributors() > 0 )
    {
      if (!IsHistogramDisabled("Xvertex"))
      {
        p.Histo("Xvertex")->Fill(vertex->GetX());
      }
      if (!IsHistogramDisabled("Yvertex"))
      {
        p.Histo("Yvertex")->Fill(vertex->GetY());
      }
      if (!IsHistogramDisabled("Zvertex"))
      {
        p.Histo("Zvertex")->Fill(vertex->GetZ());
      }
      if ( vertexFromSPD )
      {
        if (!IsHistogramDisabled("ZvertexMinusZvertexSPD"))
        {
          p.Histo("ZvertexMinusZvertexSPD")->Fill(vertexFromSPD->GetZ()-vertex->GetZ());
        }
        if (!IsHistogramDisabled("SPDXvertex"))
        {
          p.Histo("SPDXvertex")->Fill(vertexFromSPD->GetX());
        }
        if (!IsHistogramDisabled("SPDYvertex"))
        {
          p.Histo("SPDYvertex")->Fill(vertexFromSPD->GetY());
        }
        if (!IsHistogramDisabled("SPDZvertex"))
        {
          p.Histo("SPDZvertex")->Fill(vertexFromSPD->GetZ());
        }
        if (!IsHistogramDisabled("SPDZvertexNContributors"))
        {
          p.Histo("SPDZvertexNContributors")->Fill(vertexFromSPD->GetNContributors());
        }
        if (!IsHistogramDisabled("ZvertexMinusSPDZvertexNContributors"))
        {
          p.Histo("ZvertexMinusSPDZvertexNContributors")->Fill(vertexFromSPD->GetNContributors(),vertex->GetZ() - vertexFromSPD->GetZ());
        }
        if (!IsHistogramDisabled("SPDZvertexResolutionNContributors"))
        {
          Double_t cov[6]={0};
          static_cast<const AliAODVertex*>(vertexFromSPD)->GetCovarianceMatrix(cov);
          
          p.Histo("SPDZvertexResolutionNContributors")->Fill(vertexFromSPD->GetNContributors(),TMath::Sqrt(cov[5]));
        }
        if (!IsHistogramDisabled("SPDVertexType"))
        {
          p.Histo("SPDVertexType")->Fill(vertexFromSPD->GetTitle(),1.0);
        }
        
      }
      if (!IsHistogramDisabled("VertexType"))
      {
        p.Histo("VertexType")->Fill(vertex->GetTitle(),1.0);
      }
      if (!IsHistogramDisabled("VertexClass"))
      {
        p.Histo("VertexClass")->Fill(static_cast<const AliAODVertex*>(vertex)->GetType(),1.0);
      }
    }
    if (!IsHistogramDisabled("ZvertexNContributors"))
    {
      p.Histo("ZvertexNContributors")->Fill(vertex->GetNContributors());
    }
  }
  
  if ( AliAnalysisMuonUtility::IsAODEvent(Event()) )
  {
    const AliAODTZERO* tzero = static_cast<const AliAODEvent*>(Event())->GetTZEROData();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      p.Histo("T0Zvertex")->Fill(tzero->GetT0VertexRaw());
    }
  }
  else
  {
    const AliESDTZERO* tzero = static_cast<const AliESDEvent*>(Event())->GetESDTZERO();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      p.Histo("T0Zvertex")->Fill(tzero->GetT0zVertex());
    }
  }
  
  AliVVZERO* vzero = Event()->GetVZEROData();
  
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t x = v0a-v0c;
    Float_t y = v0a+v0c;
    
    if( vertexFromSPD )
    {
      Float_t multV0A = vzero->GetMTotV0A();
      Float_t multV0C = vzero->GetMTotV0C();
      Float_t multV0 = multV0A + multV0C; //Corrected V0 tot mult?
      
      Float_t v0aMult = AliESDUtils::GetCorrV0A(multV0A,vertexFromSPD->GetZ());
      Float_t v0cMult = AliESDUtils::GetCorrV0C(multV0C,vertexFromSPD->GetZ());
      
      if (!IsHistogramDisabled("V0AMult"))
      {
        p.Histo("V0AMult")->Fill(v0aMult);
      }
      if (!IsHistogramDisabled("V0CMult"))
      {
        p.Histo("V0CMult")->Fill(v0cMult);
      }
      if (!IsHistogramDisabled("V0TotMult"))
      {
        p.Histo("V0TotMult")->Fill(multV0);
      }
    }
    
    
    if (!IsHistogramDisabled("V02D"))
    {
      p.Histo("V02D")->Fill(x,y);
    }
    
    Bool_t background,pileup,satellite;
    
    Bool_t tzero = AliAnalysisMuonUtility::EAGetTZEROFlags(Event(),background,pileup,satellite);
    
    if (tzero)
    {
      if ( background )
      {
        if (!IsHistogramDisabled("V02DwT0BG"))
        {
          p.Histo("V02DwT0BG")->Fill(x,y);
        }
      }
      
      if ( pileup )
      {
        if (!IsHistogramDisabled("V02DwT0PU"))
        {
          p.Histo("V02DwT0PU")->Fill(x,y);
        }
        
        if ( !IsHistogramDisabled("PileUpEstimators") )
        {
          p.Histo("PileUpEstimators")->Fill("TZERO",1.0);
        }
      }
      
      if ( satellite )
      {
        if (!IsHistogramDisabled("V02DwT0SAT"))
        {
          p.Histo("V02DwT0SAT")->Fill(x,y);
        }
      }
      
      if ( !background && !pileup && !satellite )
      {
        if (!IsHistogramDisabled("V02DwT0BB"))
        {
          p.Histo("V02DwT0BB")->Fill(x,y);
        }
      }
    }
  }
  
  //  /* FIXME : how to properly get multiplicity from AOD and ESD consistently ?
  //   is is doable at all ?
  
  TH1* hpileup = p.Histo("PileUpEstimators");
  
  
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
void AliAnalysisMuMuGlobal::FillHistosForMCEvent(AliMergeableCollectionProxy& p)
{
  // Fill MCEvent-wise histograms
  
  Double_t Zvertex = AliAnalysisMuonUtility::GetMCVertexZ(Event(),MCEvent());
  
  if (!IsHistogramDisabled("Zvertex"))
  {
    p.Histo("Zvertex")->Fill(Zvertex);
  }
  
  if (!IsHistogramDisabled("RecZvertexVsMCZvertex"))
  {
    const AliVVertex* vertex = Event()->GetPrimaryVertex();
    if  (vertex && vertex->GetNContributors()>0)
    {
      p.Histo("RecZvertexVsMCZvertex")->Fill(Zvertex,vertex->GetZ());
    }
    
    const AliVVertex* vertexFromSPD = Event()->GetPrimaryVertexSPD();
    if  (vertexFromSPD && vertexFromSPD->GetNContributors()>0)
    {
      p.Histo("RecSPDZvertexVsMCZvertex")->Fill(Zvertex,vertexFromSPD->GetZ());
      p.Histo("NofEvWSPDZvertexVsMCZvertex")->Fill(Zvertex,1);
      
      if ( !vertexFromSPD->IsFromVertexerZ() )
      {
        p.Histo("NofEvWSPDZvertexAndNoVtexerZVsMCZvertex")->Fill(Zvertex,1);
        
        Double_t cov[6]={0};
        vertexFromSPD->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        Double_t zvertex = vertexFromSPD->GetZ();
        if ( (zRes <= 0.25) && TMath::Abs(zvertex - vertex->GetZ()) <= 0.5 ) //These events are those passing AliAnalysisMuMuEventCutter::IsSPDzQA()
        {
          p.Histo("NofEvPassingVtxQAVsMCZvertex")->Fill(Zvertex,1);
        }
        else p.Histo("NofEvNotPassingVtxResCutVsMCZvertex")->Fill(Zvertex,1);
      }
      else p.Histo("NofEvWSPDZvertexAndVtexerZVsMCZvertex")->Fill(Zvertex,1);
    }
    else p.Histo("NofEvWOSPDZvertexVsMCZvertex")->Fill(Zvertex,1);
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::FillHistosForEvent(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality)
{
  // Fill event-wise histograms
  
  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality));

  FillHistosForEvent(*proxy);
  
  delete proxy;
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::FillHistosForMCEvent(const char* eventSelection,
                                                 const char* triggerClassName,
                                                 const char* centrality)
{
  // Fill MCEvent-wise histograms

  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality));

  FillHistosForMCEvent(*proxy);
  
  delete proxy;
}

//_____________________________________________________________________________
void AliAnalysisMuMuGlobal::DefineHistogramCollection(const char* eventSelection,
                                                      const char* triggerClassName,
                                                      const char* centrality,
                                                      Bool_t mix)
{
  /// Actually create the histograms for phyics/triggerClassName
  
  
  if ( ExistSemaphoreHistogram(eventSelection,triggerClassName,centrality) )
  {
    return;
  }
  
  CreateSemaphoreHistogram(eventSelection,triggerClassName,centrality);
  
  Double_t xmin = -40;
  Double_t xmax = +40;
  Int_t nbins = GetNbins(xmin,xmax,0.5);
  
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Zvertex","z vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"SPDZvertex","SPD z vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"T0Zvertex","T0 zvertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "RecZvertexVsMCZvertex","Reconstructed vertex (w/ Ncontrib>=1) vs MC vertex",nbins,xmin,xmax,nbins,xmin,xmax);
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "RecSPDZvertexVsMCZvertex","Reconstructed SPD vertex (w/ Ncontrib>=1) vs MC vertex",nbins,xmin,xmax,nbins,xmin,xmax);
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvWSPDZvertexVsMCZvertex","Number of events with SPD vertex (w/ Ncontrib>=1) vs MC vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvWSPDZvertexAndNoVtexerZVsMCZvertex","Number of events with SPD vertex (w/ Ncontrib>=1) and no vertexer: Z vs MC vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvPassingVtxQAVsMCZvertex","Number of events with SPD vertex passing QA vs MC vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvNotPassingVtxResCutVsMCZvertex","Number of events with SPD vertex not passing resolution cut vs MC vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvWSPDZvertexAndVtexerZVsMCZvertex","Number of events with SPD vertex and vertexer: Z vs MC vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,
                    "NofEvWOSPDZvertexVsMCZvertex","Number of events without SPD vertex vs MC vertex",nbins,xmin,xmax);
  
  
  xmin = -5;
  xmax = 5;
  nbins = GetNbins(xmin,xmax,0.01);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"ZvertexMinusZvertexSPD","Primary vertex z - SPD vertex",nbins,xmin,xmax);
  
  xmin = -1;
  xmax = 50;
  nbins = GetNbins(xmin,xmax,1);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"ZvertexNContributors","z vertex vs nof contributors",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"SPDZvertexNContributors","SPD z vertex vs nof contributors",nbins,xmin,xmax);
  
  Double_t ymin = -10;
  Double_t ymax = 10;
  Int_t nbinsy = GetNbins(ymin,ymax,0.01);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"ZvertexMinusSPDZvertexNContributors","Primary vertex z - SPD vertex vs nof contributors",nbins,xmin,xmax,nbinsy,ymin,ymax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"SPDZvertexResolutionNContributors","SPD vertex resolution vs nof contributors",nbins,xmin,xmax,nbinsy,ymin,ymax);
  
  xmin = -2;
  xmax = 2;
  nbins = GetNbins(xmin,xmax,0.01);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Xvertex","x vertex",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Yvertex","y vertex",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"SPDXvertex","SPD x vertex",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"SPDYvertex","SPD y vertex",nbins,xmin,xmax);
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"Nevents","number of events",2,-0.5,1.5);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"EventsWOL0inputs","number of events w/o L0 inputs",2,-0.5,1.5);
  
  CreateEventHistos(kHistoForMCInput | kHistoForData,eventSelection,triggerClassName,centrality,
                    "VertexType","Type of vertexer used",10,0,10);
  
  CreateEventHistos(kHistoForMCInput | kHistoForData,eventSelection,triggerClassName,centrality,
                    "VertexClass","Type of vertex used",10,0,10);
  
  CreateEventHistos(kHistoForMCInput | kHistoForData,eventSelection,triggerClassName,centrality,
                    "SPDVertexType","Type of SPD vertexer used",10,0,10);
  
  
  
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
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0AMult","V0A multiplicity (Corrected by AliESDUtils);V0A mult;N_{events}",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0CMult","V0C multiplicity (Corrected by AliESDUtils);V0C mult;N_{events}",nbins,xmin,xmax);
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"V0TotMult","V0 multiplicity;V0 mult;N_{events}",nbins,xmin,xmax);
  
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
      TString name;
      
      name.Form("/%s/%s/%s/Centrality",eventSelection,triggerClassName,it->c_str());
      
      TH1* h = HistogramCollection()->Histo(name.Data());
      
      if (!h)
      {
        h = new TH1F("Centrality","Centrality",350,-50,300);
        AliInfo(Form("Creating centrality histogram for estimator %s %s %s",it->c_str(),eventSelection,triggerClassName));
        HistogramCollection()->Adopt(Form("/%s/%s/%s",eventSelection,triggerClassName,it->c_str()),h);
      }
    }
    
    delete centralities;
  }
  
  CreateEventHistos(kHistoForData,eventSelection,triggerClassName,centrality,"PileUpEstimators","pile up estimators",10,0,10);
}

