#include "AliAnalysisMuMuSingle.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuSingle
 *
 * Histogramming of single muon tracks. Mostly to get control plots for
 * the AliAnalysisMuMuMinv sub-analysis, with respect to track cuts used,
 * like Rabs, p x DCA, etc...
 *
 */

#include "TH2F.h"
#include "AliCodeTimer.h"
#include "AliMuonTrackCuts.h"
#include "AliAnalysisMuonUtility.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "AliAnalysisMuMuCutCombination.h"
#include "AliAnalysisMuMuCutRegistry.h"
#include "AliMergeableCollection.h"
#include "AliVEvent.h"

ClassImp(AliAnalysisMuMuSingle)

//_____________________________________________________________________________
AliAnalysisMuMuSingle::AliAnalysisMuMuSingle()
: AliAnalysisMuMuBase(),
fMuonTrackCuts(0x0),
fShouldSeparatePlusAndMinus(kFALSE),
fAccEffHisto(0x0),
fPtEtaSpectraPerBCX(kFALSE),
fDCAHistos(kFALSE)
{
  /// ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuSingle::~AliAnalysisMuMuSingle()
{
  /// dtor
  delete fMuonTrackCuts;
  delete fAccEffHisto;
}


//_____________________________________________________________________________
void
AliAnalysisMuMuSingle::CreateTrackHisto(const char* eventSelection,
                                        const char* triggerClassName,
                                        const char* centrality,
                                        const char* hname, const char* htitle,
                                        Int_t nbinsx, Double_t xmin, Double_t xmax,
                                        Int_t nbinsy, Double_t ymin, Double_t ymax,
                                        Bool_t separatePlusAndMinus) const
{
  /// Append histograms for single track to our histogram collection

  if ( IsHistogramDisabled(hname) ) return;

  if ( separatePlusAndMinus )
  {
    const char* suffix[] = { "Plus", "Minus" };
    const char* symbol[] = { "+", "-" };

    for ( Int_t i = 0; i < 2; ++i )
    {
      TString shtitle(htitle);
      TString shname(hname);

      shtitle.ReplaceAll("#mu",Form("#mu^{%s}",symbol[i]));

      shname += suffix[i];

      CreateTrackHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,shname.Data(),shtitle.Data(),
                        nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    }
  }
  else
  {
    CreateTrackHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,hname,htitle,
                nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  }
}

//_____________________________________________________________________________
void
AliAnalysisMuMuSingle::CreateTrackTHnSparse(const char* eventSelection,
                                        const char* triggerClassName,
                                        const char* centrality,
                                        const char* hname, const char* htitle,
                                        Int_t nDim, Int_t* nbinsx, Double_t* xmin, Double_t* xmax,
                                        Bool_t separatePlusAndMinus) const
{
  /// Append histograms for single track to our histogram collection

  if ( IsHistogramDisabled(hname) ) return;

  if ( separatePlusAndMinus )
  {
    const char* suffix[] = { "Plus", "Minus" };
    const char* symbol[] = { "+", "-" };

    for ( Int_t i = 0; i < 2; ++i )
    {
      TString shtitle(htitle);
      TString shname(hname);

      shtitle.ReplaceAll("#mu",Form("#mu^{%s}",symbol[i]));

      shname += suffix[i];

      AliAnalysisMuMuBase::CreateTrackTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,shname.Data(),shtitle.Data(),nDim,nbinsx,xmin,xmax);
    }
  }
  else
  {
    AliAnalysisMuMuBase::CreateTrackTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,hname,htitle,nDim,nbinsx,xmin,xmax);
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSingle::IsPDCAOK(const AliVParticle& part)
{
  UInt_t selectionMask = MuonTrackCuts() ? MuonTrackCuts()->GetSelectionMask(&part) : 0;

  return ( ( selectionMask & AliMuonTrackCuts::kMuPdca ) == AliMuonTrackCuts::kMuPdca );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSingle::IsRabsOK(const AliVParticle& part) const
{
  Double_t thetaAbsEndDeg = AliAnalysisMuonUtility::GetThetaAbsDeg(&part);

  return ( thetaAbsEndDeg > 2. && thetaAbsEndDeg < 10. );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuSingle::IsEtaInRange(const AliVParticle& part) const
{
  return (part.Eta() > -4.0 && part.Eta() < -2.5);
}

//_____________________________________________________________________________
void AliAnalysisMuMuSingle::SetRun(const AliInputEventHandler* eventHandler)
{
  MuonTrackCuts()->SetRun(eventHandler);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuSingle::EAGetNumberOfMuonTracks() const
{
  // Get the number of muon tracks *that are not ghosts*

  Int_t ntracks = AliAnalysisMuonUtility::GetNTracks(Event());

  for ( Int_t i = 0; i < ntracks; ++i )
  {
    AliVParticle* track = AliAnalysisMuonUtility::GetTrack(i,Event());
    if (AliAnalysisMuonUtility::IsMuonGhost(track)) --ntracks;
  }

  return ntracks;
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuSingle::EAGetTrackDCA(const AliVParticle& track) const
{
  // Get track DCA

  Double_t xdca = AliAnalysisMuonUtility::GetXatDCA(&track);
  Double_t ydca = AliAnalysisMuonUtility::GetYatDCA(&track);

  return TMath::Sqrt(xdca*xdca+ydca*ydca);
}

//_____________________________________________________________________________
void AliAnalysisMuMuSingle::DefineHistogramCollection(const char* eventSelection,
                                                      const char* triggerClassName,
                                                      const char* centrality,
                                                      Bool_t mix)
{
  /// Actually create the histograms for phyics/triggerClassName


  if ( Histo(eventSelection,triggerClassName,centrality,"AliAnalysisMuMuSingle") )
  {
    return;
  }

  AliAnalysisMuMuBase::EDataType dt = AliAnalysisMuMuBase::kHistoForData;

  // dummy histogram to signal that we already defined all our histograms (see above)
  CreateEventHistos(dt,eventSelection,triggerClassName,centrality,"AliAnalysisMuMuSingle","Dummy semaphore",1,0,1);

  Double_t ptMin  = 0;
  Double_t ptMax  = 12*3;
  Int_t nbinsPt   = GetNbins(ptMin,ptMax,0.05);
  Double_t pMin   = 0;
  Double_t pMax   = 100*3;
  Int_t nbinsP    = GetNbins(pMin,pMax,2.0);
  Double_t etaMin = -5;
  Double_t etaMax = -2;
  Int_t nbinsEta  = GetNbins(etaMin,etaMax,0.05);

  Double_t rapidityMin = -5;
  Double_t rapidityMax = -2;
  Int_t nbinsRapidity  = GetNbins(rapidityMin,rapidityMax,0.05);

  Double_t phiMin = -TMath::Pi();
  Double_t phiMax = TMath::Pi();
  Int_t nbinsPhi  = GetNbins(phiMin,phiMax,0.05);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"Chi2MatchTrigger","Chi2 Match Trigger",72,0,72);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"EtaRapidityMu", "Eta distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax,nbinsEta,etaMin,etaMax, fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"PtEtaMu", "P_{T} distribution vs Eta for #mu", nbinsEta,etaMin,etaMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"PtRapidityMu", "P_{T} distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"PtPhiMu", "P_{T} distribution vs phi for #mu", nbinsPhi,phiMin,phiMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);


  CreateTrackHisto(eventSelection,triggerClassName,centrality,"PEtaMu", "P distribution for #mu",nbinsEta,etaMin,etaMax,nbinsP,pMin,pMax,fShouldSeparatePlusAndMinus);

  Double_t chi2min = 0;
  Double_t chi2max = 20;
  Int_t nbinchi2   = GetNbins(chi2min,chi2max,0.05);

  CreateTrackHisto(eventSelection, triggerClassName, centrality, "Chi2Mu", "chisquare per NDF #mu", nbinchi2, chi2min, chi2max,-1, 0.0, 0.0, fShouldSeparatePlusAndMinus);

  Double_t xmin = 0;
  Double_t xmax = 150;
  Int_t nbins   = GetNbins(xmin,xmax,2.0);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"dcaP23Mu","#mu DCA vs P for 2-3 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"dcaP310Mu","#mu DCA vs P for 3-10 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"dcaPwPtCut23Mu","#mu DCA vs P for 2-3 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"dcaPwPtCut310Mu","#mu DCA vs P for 3-10 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);


  // Int_t nBins[3]={234,nbinsPt,nbinsEta};
  // Double_t xminsparse[3]={0.5,ptMin,etaMin};
  // Double_t xmaxsparse[3]={234.5,ptMax,etaMax};

  // CreateTrackTHnSparse(eventSelection,triggerClassName,centrality,"HitperTriggerLocalBoardMu","#mu P_{T}-eta distribution vs Trigger Local Board",3,nBins,xminsparse,xmaxsparse,fShouldSeparatePlusAndMinus);

  xmin  = 0;
  xmax  = 3564;
  nbins = GetNbins(xmin,xmax,1.0);

  CreateTrackHisto(eventSelection,triggerClassName,centrality,"BCX","bunch-crossing ids",nbins,xmin-0.5,xmax-0.5);
}


//_____________________________________________________________________________
void AliAnalysisMuMuSingle::FillHistosForMuonTrack(AliMergeableCollectionProxy& proxy,
                                                   const AliVParticle& track)
{
  /// Fill histograms for one track

  AliCodeTimerAuto("",0);

  if ( HasMC() )
  {
    MuonTrackCuts()->SetIsMC();
  }

  TLorentzVector p(track.Px(),track.Py(),track.Pz(),
                   TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+track.P()*track.P()));


  TString charge("");

  if ( ShouldSeparatePlusAndMinus() )
  {
    if ( track.Charge() < 0 )
    {
      charge = "Minus";
    }
    else
    {
      charge = "Plus";
    }
  }

  Double_t dca = EAGetTrackDCA(track);

  Double_t theta = AliAnalysisMuonUtility::GetThetaAbsDeg(&track);

  if (!IsHistogramDisabled("BCX"))
  {
    proxy.Histo("BCX")->Fill(1.0*Event()->GetBunchCrossNumber());
  }

  if (!IsHistogramDisabled("Chi2MatchTrigger"))
  {
    proxy.Histo("Chi2MatchTrigger")->Fill(AliAnalysisMuonUtility::GetChi2MatchTrigger(&track));
  }

  if (!IsHistogramDisabled("EtaRapidityMu*"))
  {
    proxy.Histo(Form("EtaRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Eta());
  }

  if (!IsHistogramDisabled("PtEtaMu*"))
  {
    TH1* h = proxy.Histo(Form("PtEtaMu%s",charge.Data()));

    h->Fill(p.Eta(),p.Pt());

    if  ( fPtEtaSpectraPerBCX )
    {
      if (!IsHistogramDisabled("BCX"))
      {
        TH1* hbcx = proxy.Histo(Form("PtEtaMu%sBCX%d",charge.Data(),Event()->GetBunchCrossNumber()));

        if (!hbcx)
        {
          hbcx = static_cast<TH1*>(h->Clone(Form("PtEtaMu%sBCX%d",charge.Data(),Event()->GetBunchCrossNumber())));
          proxy.Adopt(hbcx);
        }
      }
    }
  }

  if (!IsHistogramDisabled("PtRapidityMu*"))
  {
    proxy.Histo(Form("PtRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Pt());
  }

  if (!IsHistogramDisabled("PEtaMu*"))
  {
    proxy.Histo(Form("PEtaMu%s",charge.Data()))->Fill(p.Eta(),p.P());
  }

  if (!IsHistogramDisabled("PtPhiMu*"))
  {
    proxy.Histo(Form("PtPhiMu%s",charge.Data()))->Fill(p.Phi(),p.Pt());
  }

  if (!IsHistogramDisabled("Chi2Mu*"))
  {
    proxy.Histo(Form("Chi2Mu%s",charge.Data()))->Fill(AliAnalysisMuonUtility::GetChi2perNDFtracker(&track));
  }

  // if (!IsHistogramDisabled("HitperTriggerLocalBoardMu*"))
  // {
  //   if(AliAnalysisMuonUtility::GetLoCircuit(&track) !=0){
  //     Double_t x[3] = {static_cast<Double_t>(AliAnalysisMuonUtility::GetLoCircuit(&track)),p.Pt(),p.Eta()};
  //     if(proxy.GetObject(Form("HitperTriggerLocalBoardMu%s",charge.Data())))
  //       static_cast<THnSparse*>(proxy.GetObject(Form("HitperTriggerLocalBoardMu%s",charge.Data())))->Fill(x);

  //   }
  // }

  if (!fDCAHistos)
  {
    return;
  }

  if ( theta >= 2.0 && theta < 3.0 )
  {

    if (!IsHistogramDisabled("dcaP23Mu*"))
    {
      proxy.Histo(Form("dcaP23Mu%s",charge.Data()))->Fill(p.P(),dca);
    }

    if ( p.Pt() > 2 )
    {
      if (!IsHistogramDisabled("dcaPwPtCut23Mu*"))
      {
        proxy.Histo(Form("dcaPwPtCut23Mu%s",charge.Data()))->Fill(p.P(),dca);
      }
    }
  }
  else if ( theta >= 3.0 && theta < 10.0 )
  {
    if (!IsHistogramDisabled("dcaP310Mu*"))
    {
      proxy.Histo(Form("dcaP310Mu%s",charge.Data()))->Fill(p.P(),dca);
    }
    if ( p.Pt() > 2 )
    {
      if (!IsHistogramDisabled("dcaPwPtCut310Mu*"))
      {
        proxy.Histo(Form("dcaPwPtCut310Mu%s",charge.Data()))->Fill(p.P(),dca);
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuSingle::FillHistosForTrack(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality,
                                               const char* trackCutName,
                                               const AliVParticle& track)
{
  /// Fill histograms for one track

  if (!AliAnalysisMuonUtility::IsMuonTrack(&track) ) return;

  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,trackCutName));

  FillHistosForMuonTrack(*proxy,track);

  delete proxy;
}

//_____________________________________________________________________________
AliMuonTrackCuts* AliAnalysisMuMuSingle::MuonTrackCuts()
{
  /// Get (and create the first time) our internal track cuts
  if (!fMuonTrackCuts)
  {
    fMuonTrackCuts = new AliMuonTrackCuts;

    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta |
                                  AliMuonTrackCuts::kMuThetaAbs |
                                  AliMuonTrackCuts::kMuPdca |
                                  AliMuonTrackCuts::kMuMatchApt |
                                  AliMuonTrackCuts::kMuMatchLpt |
                                  AliMuonTrackCuts::kMuMatchHpt |
                                  AliMuonTrackCuts::kMuTrackChiSquare);

  }

  return fMuonTrackCuts;
}

//_____________________________________________________________________________
void AliAnalysisMuMuSingle::SetMuonTrackCuts(const AliMuonTrackCuts& trackCuts)
{
  /// Set our muontrackcuts from external source
  delete fMuonTrackCuts;
  fMuonTrackCuts = static_cast<AliMuonTrackCuts*>(trackCuts.Clone());
}
