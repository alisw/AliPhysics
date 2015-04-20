#include "AliAnalysisMuMuMinv.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuMinv
 *
 * Analysis which fills a bunch of histograms for invariant mass analysis of J/psi
 *
 * Can be used on real data and/or simulated (for instance to get Acc x Eff)
 *
 * Can optionally use as input an already computed Acc x Eff matrix that will be applied
 * when filling the invariant mass histograms.
 *
 */

#include "TH2F.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "AliAnalysisMuMuBinning.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "AliMCEvent.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"
#include "TParameter.h"
#include <cassert>

ClassImp(AliAnalysisMuMuMinv)

//_____________________________________________________________________________
AliAnalysisMuMuMinv::AliAnalysisMuMuMinv(TH2* accEffHisto, Bool_t computeMeanPt, Int_t systLevel)
: AliAnalysisMuMuBase(),
fcomputeMeanPt(computeMeanPt),
fAccEffHisto(0x0),
fsystLevel(systLevel)
{
  // FIXME : find the AccxEff histogram from HistogramCollection()->Histo("/EXCHANGE/JpsiAccEff")
  
  if ( accEffHisto )
  {
    fAccEffHisto = static_cast<TH2F*>(accEffHisto->Clone());
    fAccEffHisto->SetDirectory(0);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuMinv::~AliAnalysisMuMuMinv()
{
  /// dtor
  delete fAccEffHisto;
}

//_____________________________________________________________________________
void
AliAnalysisMuMuMinv::DefineHistogramCollection(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality)
{
  /// Define the histograms this analysis will use
  
  if ( Histo(eventSelection,triggerClassName,centrality,"AliAnalysisMuMuMinv") )
  {
    return;
  }
  
  // dummy histogram to signal that we already defined all our histograms (see above)
  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"AliAnalysisMuMuMinv","Dummy semaphore",1,0,1);

  /// Create invariant mass histograms
  
  Double_t minvMin = 0;
  Double_t minvMax = 16;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,0.025);
  
  Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);
  
  Double_t rapidityMin = -5;
  Double_t rapidityMax = -2;
  Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);

  Double_t etaMin = -5;
  Double_t etaMax = -2;
  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);

  TObjArray* bins = Binning()->CreateBinObjArray("psi","integrated,ptvsy,yvspt,pt,y,phi,nch,dnchdeta,ntrcorr,ntrcorrpt,ntrcorry,relntrcorr,v0acorr,v0ccorr,v0mcorr","");//We may include ,v0a,v0acent
  
  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","#mu+#mu- Pt distribution",
                  200,0,20,-2);
  
  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","#mu+#mu- Y distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);
  
  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","#mu+#mu- Eta distribution",
                   nbinsEta,etaMin,etaMax);
  
  
  //___Histos for pure MC
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","MCINPUT #mu+#mu- Pt distribution",
                   200,0,20,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","MCINPUT #mu+#mu- Y distribution",
                   nbinsRapidity,rapidityMin,rapidityMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","MCINPUT #mu+#mu- Eta distribution",
                   nbinsEta,etaMin,etaMax);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Pt","MCINPUT #mu+#mu- Pt distribution",
                    200,0,20,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Y","MCINPUT #mu+#mu- Y distribution",
                    nbinsRapidity,rapidityMin,rapidityMax,-2);
  
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Eta","MCINPUT #mu+#mu- Eta distribution",
                    nbinsEta,etaMin,etaMax);
  //____
  
  //  CreatePairHistos(eventSelection,triggerClassName,centrality,"BinFlowPt","#mu+#mu- BinFlowPt distribution",
  //                  200,0,20);
  
  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"PtRecVsSim","#mu+#mu- Pt distribution rec vs sim",
                  200,0,20,200,0,20);
  
  //________________
  Double_t multMin = -0.5;  //Tracklets multiplicity range
  Double_t multMax = 500.5;
  Int_t nbinsMult = GetNbins(multMin,multMax,1.);

  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"NchForJpsi","Corrected multiplicity distribution for 2.9 < m_{#mu^{+}#mu^{-}} < 3.3",
                   nbinsMult,multMin,multMax);
  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"NchForPsiP","Corrected multiplicity distribution for 3.6 < m_{#mu^{+}#mu^{-}} < 3.9",
                   nbinsMult,multMin,multMax);
  //________________

  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  Int_t nb(0);
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    TString minvName(GetMinvHistoName(*r,kFALSE));
    
    ++nb;
    
    if ( !IsHistogramDisabled(minvName.Data()) )
    {
      
      AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
      
      CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                       Form("#mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                            r->AsString().Data()),
                       nMinvBins,minvMin,minvMax,-2);
      
      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                       Form("MCINPUT #mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                            r->AsString().Data()),
                       nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo
      
      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),minvName.Data(),
                        Form("MCINPUT #mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                             r->AsString().Data()),
                        nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo 

      
      if ( fcomputeMeanPt )
      {
        TString mPtName(Form("MeanPtVs%s",minvName.Data()));
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                         Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",
                              r->AsString().Data()),
                         nMinvBins,minvMin,minvMax,0);
        
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                         Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",
                              r->AsString().Data()),
                         nMinvBins,minvMin,minvMax,0); //Pure MC Histo
        
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),mPtName.Data(),
                          Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",
                               r->AsString().Data()),
                          nMinvBins,minvMin,minvMax,0); //Pure MC Histo
      }
      
//      if ( HasMC() )
//      {
//        TH1* h = new TH1F(minvName.Data(),Form("MC #mu+#mu- inv. mass %s",r->AsString().Data()),
//                          nMCMinvBins,minvMin,minvMax);
//        
//        HistogramCollection()->Adopt(Form("/%s/ALL",MCInputPrefix()),h);
//        
//        HistogramCollection()->Adopt(Form("/%s/INYRANGE",MCInputPrefix()),static_cast<TH1*>(h->Clone()));
//      }
    }
    
    if ( ShouldCorrectDimuonForAccEff() )
    {
      minvName = GetMinvHistoName(*r,kTRUE);
      
      if ( !IsHistogramDisabled(minvName.Data()) )
      {
        
        AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
        
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                         Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}}(GeV/c^{2});Counts",
                              r->AsString().Data()),
                         nMinvBins,minvMin,minvMax,-2);
        
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                         Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                              r->AsString().Data()),
                         nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo
        
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),minvName.Data(),
                          Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",
                               r->AsString().Data()),
                          nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo
        
        if ( fcomputeMeanPt )
        {
          TString mPtName(Form("MeanPtVs%s",minvName.Data()));
          CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                           Form("#mu+#mu- mean p_{T} %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});<p_{T}^{#mu^{+}#mu^{-}}>",
                                r->AsString().Data()),
                           nMinvBins,minvMin,minvMax,0);
        }
        
//        if ( HasMC() )
//        {
//          TH1*  h = new TH1F(minvName.Data(),Form("MC #mu+#mu- inv. mass %s",r->AsString().Data()),
//                             nMCMinvBins,minvMin,minvMax);
//          
//          h->Sumw2();
//          
//          HistogramCollection()->Adopt(Form("/%s/ALL",MCInputPrefix()),h);
//          
//          HistogramCollection()->Adopt(Form("/%s/INYRANGE",MCInputPrefix()),static_cast<TH1*>(h->Clone()));
//        }
      }
    }
  }

  delete bins;
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::FillHistosForPair(const char* eventSelection, const char* triggerClassName,
                                            const char* centrality, const char* pairCutName,
                                            const AliVParticle& tracki,
                                            const AliVParticle& trackj)
{
  
  TLorentzVector pi(tracki.Px(),tracki.Py(),tracki.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki.P()*tracki.P()));
  
  if (!AliAnalysisMuonUtility::IsMuonTrack(&tracki) ) return;
  if (!AliAnalysisMuonUtility::IsMuonTrack(&trackj) ) return;
  
  TLorentzVector pair4Momentum(trackj.Px(),trackj.Py(),trackj.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj.P()*trackj.P()));
  
  pair4Momentum += pi;
  
  
////  if (!IsHistogramDisabled("Chi12"))
////  {
////    Histo(eventSelection,triggerClassName,centrality,pairCutName,"Chi12")
////    ->Fill(
////           AliAnalysisMuonUtility::GetChi2perNDFtracker(&tracki),
////           AliAnalysisMuonUtility::GetChi2perNDFtracker(&trackj));
////  }
////  
////  if (!IsHistogramDisabled("Rabs12"))
////  {
////    Histo(eventSelection,triggerClassName,centrality,pairCutName,"Rabs12")
////    ->Fill(AliAnalysisMuonUtility::GetRabs(&tracki),
////           AliAnalysisMuonUtility::GetRabs(&trackj));
////  }
  
  if ( ( tracki.Charge() != trackj.Charge() ) )
  {
    Double_t inputWeight = WeightDistribution(pair4Momentum.Pt(),pair4Momentum.Rapidity());
    
    if ( !IsHistogramDisabled("Pt") )
    {
      Histo(eventSelection,triggerClassName,centrality,pairCutName,"Pt")->Fill(pair4Momentum.Pt(),inputWeight);
    }
    if ( !IsHistogramDisabled("Y") )
    {
      Histo(eventSelection,triggerClassName,centrality,pairCutName,"Y")->Fill(pair4Momentum.Rapidity(),inputWeight);
    }
    if ( !IsHistogramDisabled("Eta") )
    {
      Histo(eventSelection,triggerClassName,centrality,pairCutName,"Eta")->Fill(pair4Momentum.Eta());
    }
    
    TLorentzVector* pair4MomentumMC(0x0);
    Double_t inputWeightMC(1.);
    if ( HasMC() )
    {
      Int_t labeli = tracki.GetLabel();
      Int_t labelj = trackj.GetLabel();
      
      if ( labeli < 0 || labelj < 0 )
      {
        AliError("Got negative labels!");
      }
      else
      {
        AliVParticle* mcTracki = MCEvent()->GetTrack(labeli);
        AliVParticle* mcTrackj = MCEvent()->GetTrack(labelj);
        
        TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),
                            TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTracki->P()*mcTracki->P()));
        TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),
                            TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTrackj->P()*mcTrackj->P()));
        
        mcpj += mcpi;
        
        inputWeightMC = WeightDistribution(mcpj.Pt(),mcpj.Rapidity());
        
        Histo(eventSelection,triggerClassName,centrality,pairCutName,"PtRecVsSim")->Fill(mcpj.Pt(),pair4Momentum.Pt());
        
        if ( !IsHistogramDisabled("Pt") )
        {
          MCHisto(eventSelection,triggerClassName,centrality,pairCutName,"Pt")->Fill(mcpj.Pt(),inputWeightMC);
        }
        if ( !IsHistogramDisabled("Y") )
        {
          MCHisto(eventSelection,triggerClassName,centrality,pairCutName,"Y")->Fill(mcpj.Rapidity(),inputWeightMC);
        }
        if ( !IsHistogramDisabled("Eta") )
        {
          MCHisto(eventSelection,triggerClassName,centrality,pairCutName,"Eta")->Fill(mcpj.Eta());
        }
        pair4MomentumMC = &mcpj;
        
      }
    }
    
    
    TObjArray* bins = Binning()->CreateBinObjArray("psi","integrated,ptvsy,yvspt,pt,y,phi,nch,dnchdeta,ntrcorr,ntrcorrpt,ntrcorry,relntrcorr,v0acorr,v0ccorr,v0mcorr",""); // We may include: ,v0a,v0acent
    TIter nextBin(bins);
    AliAnalysisMuMuBinning::Range* r;
    
    while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
    {
      Bool_t ok(kFALSE);
      Bool_t okMC(kFALSE);
      
      if ( r->IsIntegrated() )
      {
        ok = kTRUE;
        if ( pair4MomentumMC ) okMC = kTRUE;

        //_________________________
        TH1* h(0x0);
        if ( pair4Momentum.M() >= 2.9 && pair4Momentum.M() <= 3.3 )
        {
          h = Histo(eventSelection,triggerClassName,centrality,pairCutName,"NchForJpsi");

          Double_t ntrcorr = (-1.);
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }

              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

              if ( TString(p->GetName()).Contains("NtrCorr") )
              {
                parFound = kTRUE;
                ntrcorr = p->GetVal();
              }
            }
          }

          h->Fill(ntrcorr);
        }
        else if ( pair4Momentum.M() >= 3.6 && pair4Momentum.M() <= 3.9)
        {
          h = Histo(eventSelection,triggerClassName,centrality,pairCutName,"NchForPsiP");

          Double_t ntrcorr = (-1.);
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }

              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

              if ( TString(p->GetName()).Contains("NtrCorr") )
              {
                parFound = kTRUE;
                ntrcorr = p->GetVal();
              }
            }
          }
          h->Fill(ntrcorr);
        }
        //_________________________

      }
      else if ( r->Is2D() )
      {
        if ( r->AsString().BeginsWith("PTVSY") )
        {
          ok = r->IsInRange(pair4Momentum.Rapidity(),pair4Momentum.Pt());
          if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity(),pair4MomentumMC->Pt());
        }
        else if ( r->AsString().BeginsWith("YVSPT") )
        {
          ok = r->IsInRange(pair4Momentum.Pt(),pair4Momentum.Rapidity());
          if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity());
        }
        else if ( r->Quantity() == "NTRCORRPT" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }
              
              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));
              
              if ( TString(p->GetName()).Contains("NtrCorr") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal(),pair4Momentum.Pt());
              }
            }
          }
          
        }
        else if ( r->Quantity() == "NTRCORRY" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }
              
              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));
              
              if ( TString(p->GetName()).Contains("NtrCorr") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal(),pair4Momentum.Rapidity());
              }
            }
          }
          
        }

        else
        {
          AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
        }
      }
      else
      {
        if ( r->Quantity() == "PT" )
        {
          ok = r->IsInRange(pair4Momentum.Pt());
          if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt());
        }
        else if ( r->Quantity() == "Y" )
        {
          ok = r->IsInRange(pair4Momentum.Rapidity());
          if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity());
        }
        else if ( r->Quantity() == "PHI" )
        {
          ok = r->IsInRange(pair4Momentum.Phi());
          if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Phi());
        }
        else if ( r->Quantity() == "DNCHDETA" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }
              
              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));
              
              if ( TString(p->GetName()).Contains("dNchdEta") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal());
              }
            }
          }
          
        }
        else if ( r->Quantity() == "NTRCORR" || r->Quantity() == "RELNTRCORR" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }
              
              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));
              
              if ( TString(p->GetName()).Contains("NtrCorr") )
              {
                parFound = kTRUE;
                if ( r->Quantity() == "NTRCORR" ) ok = r->IsInRange(p->GetVal());
                else ok = r->IsInRange(p->GetVal()/5.97);
              }
            }
          }
          
        }
        else if ( r->Quantity() == "V0ACORR" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }

              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

              if ( TString(p->GetName()).Contains("V0ACorr") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal());
              }
            }
          }

        }
        else if ( r->Quantity() == "V0CCORR" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }

              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

              if ( TString(p->GetName()).Contains("V0CCorr") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal());
              }
            }
          }
          //          else AliFatal("No ntrcorr info on Event");

        }
        else if ( r->Quantity() == "V0MCORR" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            Int_t i(-1);
            Bool_t parFound(kFALSE);
            while ( i < list->GetEntries() - 1 && !parFound )
            {
              i++;
              while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) // In case there is a diferent object, just to skip it
              {
                i++;
              }

              TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

              if ( TString(p->GetName()).Contains("V0MCorr") )
              {
                parFound = kTRUE;
                ok = r->IsInRange(p->GetVal());
              }
            }
          }

        }


      }
      
      if ( ok || okMC )
      {
        TString minvName = GetMinvHistoName(*r,kFALSE);
        
        if (!IsHistogramDisabled(minvName.Data()))
        {
          TH1* h(0x0);
          if ( ok )
          {
            h = Histo(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
            
            if (!h)
            {
              AliError(Form("Could not get %s",minvName.Data()));
              //continue;
            }
            else h->Fill(pair4Momentum.M(),inputWeight);
          }
          if( okMC )
          {
            h = MCHisto(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
            
            if (!h)
            {
              AliError(Form("Could not get MC %s",minvName.Data()));
              //continue;
            }
            else h->Fill(pair4MomentumMC->M(),inputWeightMC);
          }
          
          if ( fcomputeMeanPt )
          {
            TString hprofName("");
            
            if ( ok )
            {
              hprofName= Form("MeanPtVs%s",minvName.Data());
              
              TProfile* hprof = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
              
              if ( !hprof )
              {
                AliError(Form("Could not get %s",hprofName.Data()));
              }
              
              else
              {
                //              hprof->Approximate(); //I dont think its necessary here
                hprof->Fill(pair4Momentum.M(),pair4Momentum.Pt(),inputWeight);
              }
            }
            if ( okMC )
            {
              hprofName= Form("MeanPtVs%s",minvName.Data());
              
              TProfile* hprof = MCProf(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
              
              if ( !hprof )
              {
                AliError(Form("Could not get MC %s",hprofName.Data()));
              }
              
              else
              {
                //              hprof->Approximate(); //I dont think its necessary here
                hprof->Fill(pair4MomentumMC->M(),pair4MomentumMC->Pt(),inputWeightMC);
              }
            }

          }
          
        }
        
        if ( ShouldCorrectDimuonForAccEff() )
        {
          
          Double_t AccxEff(0);
          Bool_t okAccEff(kFALSE);
          if ( ok )
          {
            AccxEff = GetAccxEff(pair4Momentum.Pt(),pair4Momentum.Rapidity());
            if ( AccxEff <= 0.0 )
            {
              AliError(Form("AccxEff < 0 for pt = %f & y = %f ",pair4Momentum.Pt(),pair4Momentum.Rapidity()));
              //            continue;
            }
            else okAccEff = kTRUE;
          }
          
          Double_t AccxEffMC(0);
          Bool_t okAccEffMC(kFALSE);
          if ( okMC )
          {
             AccxEffMC= GetAccxEff(pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity());
            if ( AccxEffMC <= 0.0 )
            {
              AliError(Form("AccxEff < 0 for MC pair with pt = %f & y = %f ",pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity()));
              //            continue;
            }
            else okAccEffMC = kTRUE;
          }
          
          minvName = GetMinvHistoName(*r,kTRUE);
          
          if (!IsHistogramDisabled(minvName.Data()))
          {
            TH1* hCorr = Histo(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
            
            if (!hCorr)
            {
              AliError(Form("Could not get %sr",minvName.Data()));
            }
            
            else  if ( okAccEff ) hCorr->Fill(pair4Momentum.M(),inputWeight/AccxEff);
            
            if( okAccEffMC )
            {
              hCorr = MCHisto(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
              
              if (!hCorr)
              {
                AliError(Form("Could not get MC %s",minvName.Data()));
                //continue;
              }
              else hCorr->Fill(pair4MomentumMC->M(),inputWeightMC/AccxEffMC);
            }
            
            if ( fcomputeMeanPt )
            {
              TString hprofCorrName("");
              if( ok )
              {
                hprofCorrName = Form("MeanPtVs%s",minvName.Data());
                
                TProfile* hprofCorr = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofCorrName.Data());
                
                if ( !hprofCorr )
                {
                  AliError(Form("Could not get %s",hprofCorrName.Data()));
                }
                else if ( okAccEff )
                {
                  //                hprofCorr->Approximate(); //I dont know if its necessary here
                  hprofCorr->Fill(pair4Momentum.M(),pair4Momentum.Pt(),inputWeight/AccxEff);
                }
              }
              if( okMC )
              {
                hprofCorrName = Form("MeanPtVs%s",minvName.Data());
                
                TProfile* hprofCorr = MCProf(eventSelection,triggerClassName,centrality,pairCutName,hprofCorrName.Data());
                
                if ( !hprofCorr )
                {
                  AliError(Form("Could not get MC %s",hprofCorrName.Data()));
                }
                else if ( okAccEffMC )
                {
                  //                hprofCorr->Approximate(); //I dont know if its necessary here
                  hprofCorr->Fill(pair4MomentumMC->M(),pair4MomentumMC->Pt(),inputWeightMC/AccxEffMC);
                }
              }

            }
            
          }
          
        }
      }
    }
    
    delete bins;
  }
  
}


//_____________________________________________________________________________
void AliAnalysisMuMuMinv::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  // Fill the input Monte-Carlo histograms related to muons. Intended to be used on pure simulations so we wont use eventSelection, triggerClassName and centrality variables.
  
  if ( !HasMC() ) return;
  
  // Specific things for MC // These histos should go in the AliAnalysisMuMuGlobal task, but then we have to loop 2 times on the MCTracks...
//  if (!Histo(MCInputPrefix(),"ALL","Pt"))
//  {
//    HistogramCollection()->Adopt(Form("/%s/ALL",MCInputPrefix()),new TH1F("Pt","Pt",200,0,25));
//    HistogramCollection()->Adopt(Form("/%s/INYRANGE",MCInputPrefix()),new TH1F("Pt","Pt",200,0,25));
//    
//    Double_t rapidityMin = -5;
//    Double_t rapidityMax = -2;
//    Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);
//    
//    HistogramCollection()->Adopt(Form("/%s/ALL",MCInputPrefix()),new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
//    HistogramCollection()->Adopt(Form("/%s/INYRANGE",MCInputPrefix()),new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
//    
//    Double_t etaMin = -5;
//    Double_t etaMax = -2;
//    Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);
//    
//    HistogramCollection()->Adopt(Form("/%s/ALL",MCInputPrefix()),new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
//    HistogramCollection()->Adopt(Form("/%s/INYRANGE",MCInputPrefix()),new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
//  }
  
  Int_t nMCTracks = MCEvent()->GetNumberOfTracks();
  
  TObjArray* bins = Binning()->CreateBinObjArray("psi","integrated,ptvsy,yvspt,pt,y,phi,nch,dnchdeta,ntrcorr,ntrcorrpt,ntrcorry,relntrcorr,v0acorr,v0ccorr,v0mcorr","");//We may include: ,v0a,v0acent.
  TIter nextBin(bins);
  AliAnalysisMuMuBinning::Range* r;
  
  for ( Int_t i = 0; i < nMCTracks; ++i )
  {
    AliVParticle* part = MCEvent()->GetTrack(i);
    
    if  (AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) &&
         AliAnalysisMuonUtility::GetMotherIndex(part)==-1)
    {
      Double_t inputWeight = WeightDistribution(part->Pt(),part->Y());
      
      MCHisto(eventSelection,triggerClassName,centrality,"Pt")->Fill(part->Pt(),inputWeight);
      MCHisto(eventSelection,triggerClassName,centrality,"Y")->Fill(part->Y(),inputWeight);
      MCHisto(eventSelection,triggerClassName,centrality,"Eta")->Fill(part->Eta());
      
      if ( part->Y() < -2.5 && part->Y() > -4.0 )
      {
        MCHisto(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Pt")->Fill(part->Pt(),inputWeight);
        MCHisto(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Y")->Fill(part->Y(),inputWeight);
        MCHisto(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Eta")->Fill(part->Eta());
      }
      
      nextBin.Reset();
      
      while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
      {
        Bool_t ok(kFALSE);
        
        if ( r->IsIntegrated() )
        {
          ok = kTRUE;
        }
        else if ( r->Is2D() )
        {
          if ( r->AsString().BeginsWith("PTVSY") )
          {
            ok = r->IsInRange(part->Y(),part->Pt());
          }
          else if ( r->AsString().BeginsWith("YVSPT") )
          {
            ok = r->IsInRange(part->Pt(),part->Y());
          }
          else
          {
            AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
          }
        }
        else
        {
          if ( r->Quantity() == "PT" )
          {
            ok = r->IsInRange(part->Pt());
          }
          else if ( r->Quantity() == "Y" )
          {
            ok = r->IsInRange(part->Y());
          }
          else if ( r->Quantity() == "PHI" )
          {
            ok = r->IsInRange(part->Phi());
          }
        }
        
        if ( ok )
        {
          TString hname = GetMinvHistoName(*r,kFALSE);
          
          if (!IsHistogramDisabled(hname.Data()))
          {
            
            TH1* h = MCHisto(eventSelection,triggerClassName,centrality,hname.Data());
            
            if (!h)
            {
              AliError(Form("Could not get /%s/%s/%s/%s/ %s",MCInputPrefix(),eventSelection,triggerClassName,centrality,hname.Data()));
              continue;
            }
            
            h->Fill(part->M(),inputWeight);
            
            if ( part->Y() < -2.5 && part->Y() > -4.0 )
            {
              h = MCHisto(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),hname.Data());
              if (!h)
              {
                AliError(Form("Could not get /%s/%s/%s/%s/INYRANGE %s",MCInputPrefix(),eventSelection,triggerClassName,centrality,hname.Data()));
                continue;
              }
              h->Fill(part->M(),inputWeight);
            }
            
          }
          
          if ( fcomputeMeanPt )
          {
            TString hprofName= Form("MeanPtVs%s",hname.Data());
            
            TProfile* hprof = MCProf(eventSelection,triggerClassName,centrality,hprofName.Data());
            
            if ( !hprof )
            {
              AliError(Form("Could not get %s",hprofName.Data()));
            }
            else
            {
              //              hprof->Approximate(); //I dont think its necessary here
              hprof->Fill(part->M(),part->Pt(),inputWeight);
            }
            
            if ( part->Y() < -2.5 && part->Y() > -4.0 )
            {
              hprof = MCProf(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),hprofName.Data());
              
              if ( !hprof )
              {
                AliError(Form("Could not get %s",hprofName.Data()));
              }
              else
              {
                //              hprof->Approximate(); //I dont think its necessary here
                hprof->Fill(part->M(),part->Pt(),inputWeight);
              }

            }

            
          }
          
//          if ( ShouldCorrectDimuonForAccEff() ) //What is the sense of this?? We should not correct the input
//          {
//            Double_t AccxEff = GetAccxEff(part->Pt(),part->Y());
//            if ( AccxEff <= 0.0 )
//            {
//              AliError(Form("AccxEff < 0 for pt = %f & y = %f ",part->Pt(),part->Y()));
//              continue;
//            }
//            hname = GetMinvHistoName(*r,kTRUE);
//            
//            if (!IsHistogramDisabled(hname.Data()))
//            {
//              
//              TH1* h = MCHisto(eventSelection,triggerClassName,centrality,hname.Data());
//              
//              if (!h)
//              {
//                AliError(Form("Could not get /%s/%s/%s/%s/ %s",MCInputPrefix(),eventSelection,triggerClassName,centrality,hname.Data()));
//                continue;
//              }
//              
//              h->Fill(part->M(),inputWeight/AccxEff);
//              
//              if ( part->Y() < -2.5 && part->Y() > -4.0 )
//              {
//                h = MCHisto(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),hname.Data());
//                if (!h)
//                {
//                  AliError(Form("Could not get /%s/%s/%s/%s/INYRANGE %s",MCInputPrefix(),eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),hname.Data()));
//                  continue;
//                }
//                h->Fill(part->M(),inputWeight./AccxEff);
//              }
//              
//            }
//            
//          }
          
        }
      }
    }
  }
  
  delete bins;

  
//  for ( Int_t i = 0; i < nMCTracks; ++i )
//  {
//    AliVParticle* part = MCEvent()->GetTrack(i);
//    
//    if  (AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) &&
//         AliAnalysisMuonUtility::GetMotherIndex(part)==-1)
//    {
//      
//      Histo(MCInputPrefix(),"ALL","Pt")->Fill(part->Pt());
//      Histo(MCInputPrefix(),"ALL","Y")->Fill(part->Y());
//      Histo(MCInputPrefix(),"ALL","Eta")->Fill(part->Eta());
//      
//      if ( part->Y() < -2.5 && part->Y() > -4.0 )
//      {
//        Histo(MCInputPrefix(),"INYRANGE","Pt")->Fill(part->Pt());
//        Histo(MCInputPrefix(),"INYRANGE","Y")->Fill(part->Y());
//        Histo(MCInputPrefix(),"INYRANGE","Eta")->Fill(part->Eta());
//      }
//      
//      nextBin.Reset();
//      
//      while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
//      {
//        Bool_t ok(kFALSE);
//        
//        if ( r->IsIntegrated() )
//        {
//          ok = kTRUE;
//        }
//        else if ( r->Is2D() )
//        {
//          if ( r->AsString().BeginsWith("PTVSY") )
//          {
//            ok = r->IsInRange(part->Y(),part->Pt());
//          }
//          else if ( r->AsString().BeginsWith("YVSPT") )
//          {
//            ok = r->IsInRange(part->Pt(),part->Y());
//          }
//          else
//          {
//            AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
//          }
//        }
//        else
//        {
//          if ( r->Quantity() == "PT" )
//          {
//            ok = r->IsInRange(part->Pt());
//          }
//          else if ( r->Quantity() == "Y" )
//          {
//            ok = r->IsInRange(part->Y());
//          }
//          else if ( r->Quantity() == "PHI" )
//          {
//            ok = r->IsInRange(part->Phi());
//          }
//        }
//        
//        if ( ok )
//        {
//          TString hname = GetMinvHistoName(*r,kFALSE);
//          
//          if (!IsHistogramDisabled(hname.Data()))
//          {
//            
//            TH1* h = Histo(MCInputPrefix(),"ALL",hname.Data());
//            
//            if (!h)
//            {
//              AliError(Form("Could not get ALL %s",hname.Data()));
//              continue;
//            }
//            
//            h->Fill(part->M());
//            
//            if ( part->Y() < -2.5 && part->Y() > -4.0 )
//            {
//              h = Histo(MCInputPrefix(),"INYRANGE",hname.Data());
//              if (!h)
//              {
//                AliError(Form("Could not get INYRANGE %s",hname.Data()));
//                continue;
//              }
//              h->Fill(part->M());
//            }
//            
//          }
//          
//          if ( ShouldCorrectDimuonForAccEff() )
//          {
//            Double_t AccxEff = GetAccxEff(part->Pt(),part->Y());
//            if ( AccxEff <= 0.0 )
//            {
//              AliError(Form("AccxEff < 0 for pt = %f & y = %f ",part->Pt(),part->Y()));
//              continue;
//            }
//            hname = GetMinvHistoName(*r,kTRUE);
//            
//            if (!IsHistogramDisabled(hname.Data()))
//            {
//              
//              TH1* h = Histo(MCInputPrefix(),"ALL",hname.Data());
//              
//              if (!h)
//              {
//                AliError(Form("Could not get ALL %s",hname.Data()));
//                continue;
//              }
//              
//              h->Fill(part->M(),1./AccxEff);
//              
//              if ( part->Y() < -2.5 && part->Y() > -4.0 )
//              {
//                h = Histo(MCInputPrefix(),"INYRANGE",hname.Data());
//                if (!h)
//                {
//                  AliError(Form("Could not get INYRANGE %s",hname.Data()));
//                  continue;
//                }
//                h->Fill(part->M(),1./AccxEff);
//              }
//              
//            }
//            
//          }
//          
//        }
//      }
//    }
//    
//    delete bins;
//  }
}

//_____________________________________________________________________________
TString AliAnalysisMuMuMinv::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const
{
  return TString::Format("MinvUS%s%s",r.AsString().Data(),
                         accEffCorrected ? "_AccEffCorr" : "");
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuMinv::GetAccxEff(Double_t pt,Double_t rapidity)
{
  if (!fAccEffHisto)
  {
    AliError("ERROR: No AccxEff histo");
    return 0;
  }
  Int_t bin = fAccEffHisto->FindBin(pt,rapidity);
  Double_t accXeff = fAccEffHisto->GetBinContent(bin);
  
  return accXeff;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuMinv::WeightDistribution(Double_t pt,Double_t rapidity)
{
  //Return a weight for a dimuon pt and y, which depend on the varied distributions.
  // FIXME: hard coded, find a clean way to fix the distribution parameters from outside
  
  if (!HasMC() ) return 1.;
  
  //================ p-Pb ==============//
  //value for input distribution: this is the nominal pt and y distribution
  Double_t paryPPB[2] = {1.0,0.174};
  Double_t parptPPB[3] = {1.0,0.0557,3.52};
  
  Double_t paryHardPPB = 0.1344, parySoftPPB = 0.1971;
  Double_t par1ptHardPPB = 5.51e-2, par2ptHardPPB = 3.47,
  par1ptSoftPPB = 5.67e-2, par2ptSoftPPB = 3.68;
  //systematic 1: hardest in y  x softest in pt
  Double_t pary1[2] = {1.0,paryHardPPB};
  Double_t parpt1[3] = {1.0,par1ptSoftPPB,par2ptSoftPPB};
  //systematic 2: hardest in y x hardest in pt
  Double_t pary2[2] = {1.0,paryHardPPB};
  Double_t parpt2[3] = {1.0,par1ptHardPPB,par2ptHardPPB};
  //systematic 3: softest in y  x softest in pt
  Double_t pary3[2] = {1.0, parySoftPPB};
  Double_t parpt3[3] = {1.0,par1ptSoftPPB,par2ptSoftPPB};
  //systematic 4: softest in y  x hardest in pt
  Double_t pary4[2] = {1.0,parySoftPPB};
  Double_t parpt4[3] = {1.0,par1ptHardPPB,par2ptHardPPB};
  
  Double_t funcptvalPPB = powerLaw3Par(&pt,parptPPB);
  Double_t funcyvalPPB = normPol12Par(&rapidity,paryPPB);
  
  //================ Pb-p ==============//
  //value for input distribution
  Double_t paryPBP[2] = {1.0,0.189};
  Double_t parptPBP[3] = {1.0,0.0592,3.92};
  
  Double_t paryHardPBP = 0.1517, parySoftPBP = 0.2191;
  Double_t par1ptHardPBP = 5.58e-2, par2ptHardPBP = 3.83,
  par1ptSoftPBP = 5.59e-2, par2ptSoftPBP = 4.31;
  //systematic 5: hardest in y  x softest in pt
  Double_t pary5[2] = {1.0,paryHardPBP};
  Double_t parpt5[3] = {1.0,par1ptSoftPBP,par2ptSoftPBP};
  //systematic 6: hardest in y x hardest in pt
  Double_t pary6[2] = {1.0,paryHardPBP};
  Double_t parpt6[3] = {1.0,par1ptHardPBP,par2ptHardPBP};
  //systematic 7: softest in y  x softest in pt
  Double_t pary7[2] = {1.0, parySoftPBP};
  Double_t parpt7[3] = {1.0,par1ptSoftPBP,par2ptSoftPBP};
  //systematic 8: softest in y  x hardest in pt
  Double_t pary8[2] = {1.0,parySoftPBP};
  Double_t parpt8[3] = {1.0,par1ptHardPBP,par2ptHardPBP};
  
  Double_t funcptvalPBP = powerLaw3Par(&pt,parptPBP);
  Double_t funcyvalPBP = normPol12Par(&rapidity,paryPBP);
  
  //================ pp ==============//
  //value for input distribution
  Double_t paryPP[2] = {3.0,0.514/3.};
  Double_t parptPP[3] = {1.0,0.0546,3.90};
  
  Double_t paryHardPP = 0.4125/3., parySoftPP = 0.5958/3.;
  Double_t par1ptHardPP = 4.78e-2, par2ptHardPP = 3.65,//4.84e-2/3.45
  par1ptSoftPP = 6.12e-2, par2ptSoftPP = 4.31;//5.47e-2//4.29
  //systematic 9: hardest in y  x softest in pt
  Double_t pary9[2] = {3.0,paryHardPP};
  Double_t parpt9[3] = {1.0,par1ptSoftPP,par2ptSoftPP};
  //systematic 10: hardest in y x hardest in pt
  Double_t pary10[2] = {3.0,paryHardPP};
  Double_t parpt10[3] = {1.0,par1ptHardPP,par2ptHardPP};
  //systematic 11: softest in y  x softest in pt
  Double_t pary11[2] = {3.0, parySoftPP};
  Double_t parpt11[3] = {1.0,par1ptSoftPP,par2ptSoftPP};
  //systematic 12: softest in y  x hardest in pt
  Double_t pary12[2] = {3.0,parySoftPP};
  Double_t parpt12[3] = {1.0,par1ptHardPP,par2ptHardPP};
  
  Double_t funcptvalPP = powerLaw3Par(&pt,parptPP);
  Double_t funcyvalPP = normPol12Par(&rapidity,paryPP);

  //______
  Double_t weight(1.),funcptsyst(0.),funcysyst(0.);
  switch ( fsystLevel )
  {
    case 0:
      weight = 1;
      break;
    case 1:
      funcptsyst = powerLaw3Par(&pt,parpt1);
      if ( funcptvalPPB > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPPB;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary1);
      if ( funcyvalPPB > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPPB;
      break;
    case 2:
      funcptsyst = powerLaw3Par(&pt,parpt2);
      if ( funcptvalPPB > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPPB;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary2);
      if ( funcyvalPPB > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPPB;
      break;
    case 3:
      funcptsyst = powerLaw3Par(&pt,parpt3);
      if ( funcptvalPPB > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPPB;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary3);
      if ( funcyvalPPB > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPPB;
      break;
    case 4:
      funcptsyst = powerLaw3Par(&pt,parpt4);
      if ( funcptvalPPB > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPPB;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary4);
      if ( funcyvalPPB > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPPB;
      break;
    case 5:
      funcptsyst = powerLaw3Par(&pt,parpt5);
      if ( funcptvalPBP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPBP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary5);
      if ( funcyvalPBP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPBP;
      break;
    case 6:
      funcptsyst = powerLaw3Par(&pt,parpt6);
      if ( funcptvalPBP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPBP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary6);
      if ( funcyvalPBP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPBP;
      break;
    case 7:
      funcptsyst = powerLaw3Par(&pt,parpt7);
      if ( funcptvalPBP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPBP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary7);
      if ( funcyvalPBP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPBP;
      break;
    case 8:
      funcptsyst = powerLaw3Par(&pt,parpt8);
      if ( funcptvalPBP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPBP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary8);
      if ( funcyvalPBP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPBP;
      break;
    case 9:
      funcptsyst = powerLaw3Par(&pt,parpt9);
      if ( funcptvalPP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary9);
      if ( funcyvalPP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPP;
      break;
    case 10:
      funcptsyst = powerLaw3Par(&pt,parpt10);
      if ( funcptvalPP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary10);
      if ( funcyvalPP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPP;
      break;
    case 11:
      funcptsyst = powerLaw3Par(&pt,parpt11);
      if ( funcptvalPP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary11);
      if ( funcyvalPP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPP;
      break;
    case 12:
      funcptsyst = powerLaw3Par(&pt,parpt12);
      if ( funcptvalPP > 0 && funcptsyst > 0 ) weight = funcptsyst/funcptvalPP;
      else  weight = 1;
      funcysyst = normPol12Par(&rapidity,pary12);
      if ( funcyvalPP > 0 && funcysyst > 0 ) weight *= funcysyst/funcyvalPP;
      break;
  }
  
  return weight;
}

//________________________________________________________________________
Double_t AliAnalysisMuMuMinv::powerLaw3Par(Double_t *x, Double_t *par)
{
  //3 parameters
  Double_t arg = 0;
  
  arg = par[0]*x[0] / TMath::Power( 1 + par[1]*x[0]*x[0], par[2]);
  
  return arg;
}


//________________________________________________________________________
Double_t AliAnalysisMuMuMinv::normPol12Par(Double_t *x, Double_t *par)
{
  //2 parameters
  Double_t arg1 = 0;
  
  arg1 = par[0] * ( 1 + par[1]*x[0] );
  
  
  return arg1;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMinv::IsPtInRange(const AliVParticle& t1, const AliVParticle& t2, Double_t& ptmin, Double_t& ptmax) const
{
  /// Whether the pair passes the rapidity cut
  
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t pt = total.Pt();
  
  return  ( pt < ptmax && pt > ptmin );
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::NameOfIsPtInRange(TString& name, Double_t& ptmin, Double_t& ptmax) const
{
  name.Form("PAIRPTIN%2.1f-%2.1f",ptmin,ptmax);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMinv::IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2) const
{
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t y = total.Rapidity();

  return  ( y < -2.5 && y > -4.0 );
}
