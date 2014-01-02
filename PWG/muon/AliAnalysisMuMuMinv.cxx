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
AliAnalysisMuMuMinv::AliAnalysisMuMuMinv()
: AliAnalysisMuMuBase(),
fAccEffHisto(0x0)
{
  // FIXME : find the AccxEff histogram from HistogramCollection()->Histo("/EXCHANGE/JpsiAccEff")
  
//  if ( accEff )
//  {
//    fAccEffHisto = static_cast<TH2F*>(accEff->Clone());
//    fAccEffHisto->SetDirectory(0);
//  }
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
  
  //  TObjArray* bins = fBinning->CreateBinObjArray("psi","y vs pt,integrated,pt,y,phi","");
  TObjArray* bins = Binning()->CreateBinObjArray("psi");
  
  CreatePairHistos(eventSelection,triggerClassName,centrality,"Pt","#mu+#mu- Pt distribution",
                  200,0,20);

  Int_t nbinsy = 6;
  Double_t ymin = -4.0;
  Double_t ymax = -2.5;
  
  CreatePairHistos(eventSelection,triggerClassName,centrality,"y","#mu+#mu- y distribution",nbinsy,ymin,ymax);

  //  CreatePairHistos(eventSelection,triggerClassName,centrality,"BinFlowPt","#mu+#mu- BinFlowPt distribution",
  //                  200,0,20);
  
  CreatePairHistos(eventSelection,triggerClassName,centrality,"PtRecVsSim","#mu+#mu- Pt distribution rec vs sim",
                  200,0,20,200,0,20);
  
  if (!Histo("INPUT","ALL","Pt"))
  {
    HistogramCollection()->Adopt("/INPUT/ALL",new TH1F("Pt","Pt",200,0,20));
    HistogramCollection()->Adopt("/INPUT/INYRANGE",new TH1F("Pt","Pt",200,0,20));
    HistogramCollection()->Adopt("/INPUT/ALL",new TH1F("Y","Y",nbinsy,ymin,ymax));
    HistogramCollection()->Adopt("/INPUT/INYRANGE",new TH1F("Y","Y",nbinsy,ymin,ymax));
  }
  
  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  Int_t nb(0);
  
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    TString minvName(GetMinvHistoName(*r,ShouldCorrectDimuonForAccEff()));
    
    if ( IsHistogramDisabled(minvName.Data()) ) continue;
    
    ++nb;
    
    AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
    
    CreatePairHistos(eventSelection,triggerClassName,centrality,minvName.Data(),
                     Form("#mu+#mu- inv. mass %s %s;M_{#mu^+#mu^-} (GeV/c^2)",
                          ShouldCorrectDimuonForAccEff() ? "(AccxEff corrected)":"",
                          r->AsString().Data()),
                     nMinvBins,minvMin,minvMax);
    
    TString hname(GetMinvHistoName(*r,kFALSE));
    
    TH1* h = HistogramCollection()->Histo("/INPUT/ALL",hname.Data());
    if (!h)
    {
      h = new TH1F(hname.Data(),Form("MC #mu+#mu- inv. mass %s",r->AsString().Data()),
                   nMCMinvBins,minvMin,minvMax);
      
      HistogramCollection()->Adopt("/INPUT/ALL",h);
      
      HistogramCollection()->Adopt("/INPUT/INYRANGE",static_cast<TH1*>(h->Clone()));
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
  /** Fill histograms for muon track pairs
   */
  
  if (!AliAnalysisMuonUtility::IsMuonTrack(&tracki) ) return;
  if (!AliAnalysisMuonUtility::IsMuonTrack(&trackj) ) return;

  TLorentzVector pi(tracki.Px(),tracki.Py(),tracki.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki.P()*tracki.P()));
  
  
  TLorentzVector pair4Momentum(trackj.Px(),trackj.Py(),trackj.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj.P()*trackj.P()));
  
  pair4Momentum += pi;
  
  
  if (!IsHistogramDisabled("Chi12"))
  {
    Histo(eventSelection,triggerClassName,centrality,pairCutName,"Chi12")
    ->Fill(
           AliAnalysisMuonUtility::GetChi2perNDFtracker(&tracki),
           AliAnalysisMuonUtility::GetChi2perNDFtracker(&trackj));
  }
  
  if (!IsHistogramDisabled("Rabs12"))
  {
    Histo(eventSelection,triggerClassName,centrality,pairCutName,"Rabs12")
    ->Fill(AliAnalysisMuonUtility::GetRabs(&tracki),
           AliAnalysisMuonUtility::GetRabs(&trackj));
  }
  
  if ( ( tracki.Charge() != trackj.Charge() ) )
  {
    if ( !IsHistogramDisabled("Pt") )
    {
      Histo(eventSelection,triggerClassName,centrality,pairCutName,"Pt")->Fill(pair4Momentum.Pt());
    }

    if ( !IsHistogramDisabled("y") )
    {
      Histo(eventSelection,triggerClassName,centrality,pairCutName,"y")->Fill(pair4Momentum.Y());
    }

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
        
        Histo(eventSelection,triggerClassName,centrality,pairCutName,"PtRecVsSim")->Fill(mcpj.Pt(),pair4Momentum.Pt());
        
      }
    }
    

    TObjArray* bins = Binning()->CreateBinObjArray("psi","ptvsy,yvspt,pt,y,phi","");
    TIter nextBin(bins);
    AliAnalysisMuMuBinning::Range* r;
    
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
          ok = r->IsInRange(pair4Momentum.Rapidity(),pair4Momentum.Pt());
        }
        else if ( r->AsString().BeginsWith("YVSPT") )
        {
          ok = r->IsInRange(pair4Momentum.Pt(),pair4Momentum.Rapidity());
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
        }
        else if ( r->Quantity() == "Y" )
        {
          ok = r->IsInRange(pair4Momentum.Rapidity());
        }
        else if ( r->Quantity() == "PHI" )
        {
          ok = r->IsInRange(pair4Momentum.Phi());
        }
        else if ( r->Quantity() == "NCH" )
        {
          TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
          if (list)
          {
            TIter next(list);
            TParameter<Double_t>* p;
            
            while ( ( p = static_cast<TParameter<Double_t>*>(next())) )
            {
              if (TString(eventSelection).Contains(p->GetName()))
              {
                ok = r->IsInRange(p->GetVal());
                break;
              }
            }
          }
        }
      }
      
      if ( ok )
      {
        TString minvName = GetMinvHistoName(*r,kFALSE);
        
        if (!IsHistogramDisabled(minvName.Data()))
        {
          TH1* h = Histo(eventSelection,triggerClassName,centrality,pairCutName,minvName.Data());
          
          if (!h)
          {
            AliError(Form("Could not get %s",minvName.Data()));
            continue;
          }
          h->Fill(pair4Momentum.M());
          
          if ( fAccEffHisto )
          {
            // FIXME : fill Minv with weight = 1/AccxEff
          }
        }
      }
    }
    
    delete bins;
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::FillHistosForMCEvent(const char* /*eventSelection*/,
                                               const char* /*triggerClassName*/,
                                               const char* /*centrality*/)
{
  /** Fill histograms for MC event
   *
   * FIXME: this is here we should streamline a bit the code to carefully
   * compute the measured and truth values for Pt, Y, and (Pt,Y), in order
   * to be able to investigate unfolding techniques later on...
   */
  
  Int_t nMCTracks = MCEvent()->GetNumberOfTracks();

  for ( Int_t i = 0; i < nMCTracks; ++i )
  {
    AliVParticle* part = MCEvent()->GetTrack(i);
    
    if  ( AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) &&
          AliAnalysisMuonUtility::GetMotherIndex(part)==-1 )
    {
            Histo("INPUT","ALL","Pt")->Fill(part->Pt());
            Histo("INPUT","ALL","Y")->Fill(part->Y());
            if ( part->Y() < -2.5 && part->Y() > -4.0 )
            {
              Histo("INPUT","INYRANGE","Pt")->Fill(part->Pt());
              Histo("INPUT","INYRANGE","Y")->Fill(part->Y());
            }

    }
    
  }
  
//  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(Event());
//
//  Double_t measuredY(0.0);
//  
//  for (Int_t i = 0; i < nTracks; ++i)
//  {
//    AliVParticle* tracki = AliAnalysisMuonUtility::GetTrack(i,Event());
//
//    if (!AliAnalysisMuonUtility::IsMuonTrack(tracki) ) continue;
//
//    for (Int_t j = i+1; j < nTracks; ++j )
//    {
//      AliVParticle* trackj = AliAnalysisMuonUtility::GetTrack(j,Event());
//      
//      if (!AliAnalysisMuonUtility::IsMuonTrack(trackj) ) continue;
//
//      TLorentzVector pi(tracki->Px(),tracki->Py(),tracki->Pz(),
//                        TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki->P()*tracki->P()));
//      
//      TLorentzVector pair4Momentum(trackj->Px(),trackj->Py(),trackj->Pz(),
//                                   TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj->P()*trackj->P()));
//      
//      pair4Momentum += pi;
//      
//      measuredY = pi.Y();
//      
//      rur->Fill(measuredY,trueY); // FIXME : is this working if more than one pair is found in the reconstructed tracks ???
//      // should we try to find the closest one in minv and assign the other(s)
//      // to rur->Fake() ?
//    }
//    
//  }
//
//  if ( ! ( measuredY < 0.0 ) )
//  {
//    rur->Miss(trueY);
//  }
}

////_____________________________________________________________________________
//void AliAnalysisMuMuMinv::FillMC()
//{
//  // Fill the input Monte-Carlo histograms related to muons
//  
//  // Specific things for MC
//  if (!Histo("INPUT","ALL","Pt"))
//  {
//    HistogramCollection()->Adopt("/INPUT/ALL",new TH1F("Pt","Pt",200,0,20));
//    HistogramCollection()->Adopt("/INPUT/INYRANGE",new TH1F("Pt","Pt",200,0,20));
//    
//    Double_t rapidityMin = -5;
//    Double_t rapidityMax = -2;
//    Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);
//    
//    HistogramCollection()->Adopt("/INPUT/ALL",new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
//    HistogramCollection()->Adopt("/INPUT/INYRANGE",new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
//    
//    Double_t etaMin = -5;
//    Double_t etaMax = -2;
//    Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);
//    
//    HistogramCollection()->Adopt("/INPUT/ALL",new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
//    HistogramCollection()->Adopt("/INPUT/INYRANGE",new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
//  }
//  
//  Int_t nMCTracks = MCEvent()->GetNumberOfTracks();
//  
//  TObjArray* bins = Binning()->CreateBinObjArray("psi","ptvsy,yvspt,pt,y,phi","");
//  TIter nextBin(bins);
//  AliAnalysisMuMuBinning::Range* r;
//  
//  for ( Int_t i = 0; i < nMCTracks; ++i )
//  {
//    AliVParticle* part = MCEvent()->GetTrack(i);
//    
//    //    std::cout << "part " << i << " isprimary=" << AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) << " motherindex=" << AliAnalysisMuonUtility::GetMotherIndex(part) << std::endl;
//    //
//    //    part->Print();
//    
//    if  (AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) &&
//         AliAnalysisMuonUtility::GetMotherIndex(part)==-1)
//    {
//      
//      Histo("INPUT","ALL","Pt")->Fill(part->Pt());
//      Histo("INPUT","ALL","Y")->Fill(part->Y());
//      Histo("INPUT","ALL","Eta")->Fill(part->Eta());
//      
//      if ( part->Y() < -2.5 && part->Y() > -4.0 )
//      {
//        Histo("INPUT","INYRANGE","Pt")->Fill(part->Pt());
//        Histo("INPUT","INYRANGE","Y")->Fill(part->Y());
//        Histo("INPUT","INYRANGE","Eta")->Fill(part->Eta());
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
//            TH1* h = Histo("INPUT","ALL",hname.Data());
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
//              h = Histo("INPUT","INYRANGE",hname.Data());
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
//        }
//      }
//    }
//  }
//  
//  delete bins;
//}
//

//_____________________________________________________________________________
TString AliAnalysisMuMuMinv::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const
{
  /// Get the invariant mass histogram name
  return TString::Format("MinvUS%s%s",r.AsString().Data(),
                         accEffCorrected ? "_AccEffCorr" : "");
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMinv::IsRapidityInRange(const AliVParticle& t1, const AliVParticle& t2, Double_t& ymin, Double_t& ymax) const
{
  /// Whether the particle pair has its rapidity within [ymin,ymax[
  
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t y = total.Rapidity();

  return  ( y < ymax && y > ymin );
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::NameOfIsRapidityInRange(TString& name, Double_t& ymin, Double_t& ymax) const
{
  /** Get the name of the rapidity range (making the IsRapidityInRange method useable as an
   * AliAnalysisMuMuCutElement
   */
  
  name.Form("PAIRYIN%2.1f-%2.1f",ymin,ymax);
}

