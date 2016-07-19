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
#include "TF1.h"
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
AliAnalysisMuMuMinv::AliAnalysisMuMuMinv(TH2* accEffHisto, Int_t systLevel)
: AliAnalysisMuMuBase(),
fcomputeMeanPt(kFALSE),
fWeightMuon(kFALSE),
fAccEffHisto(0x0),
fMinvBinSeparator("+"),
fsystLevel(systLevel),
fBinsToFill(0x0),
fMinvBinSize(0.025),
fMinvMin(0.0),
fMinvMax(16.0),
fmcptcutmin(0.0),
fmcptcutmax(12.0),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0)
{
  // FIXME ? find the AccxEff histogram from HistogramCollection()->Histo("/EXCHANGE/JpsiAccEff")

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
  delete fBinsToFill;
}

//_____________________________________________________________________________
void
AliAnalysisMuMuMinv::DefineHistogramCollection(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality)
{
  /// Define the histograms this analysis will use

  // Check if histo is not already here
  if ( ExistSemaphoreHistogram(eventSelection,triggerClassName,centrality) )
  {
    return;
  }

  CreateSemaphoreHistogram(eventSelection,triggerClassName,centrality);

  // no bins defined by the external steering macro, use our own defaults
  if (!fBinsToFill) SetBinsToFill("psi","integrated,ptvsy,yvspt,pt,y,phi");

  // mass range
  Double_t minvMin = fMinvMin;
  Double_t minvMax = fMinvMax;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,fMinvBinSize);

  Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);

  // Rapidity range
  Double_t rapidityMin = -5;
  Double_t rapidityMax = -2;
  Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);

  // eta range
  Double_t etaMin = -5;
  Double_t etaMax = -2;
  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);

 //Multiplicity range
  Double_t multMin = -0.5;
  Double_t multMax = 500.5;
  Int_t nbinsMult = GetNbins(multMin,multMax,1.);


  // Reconstructed pair distribution Histo
  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","#mu+#mu- Pt distribution",200,0,20,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","#mu+#mu- Y distribution",nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","#mu+#mu- Eta distribution",nbinsEta,etaMin,etaMax);

  CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"PtPaireVsPtTrack","single #mu Pt distribution vs #mu+#mu- Pt distribution ",200,0,20,200,0,20);


  // Histos for pure MC
  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Pt","MCINPUT #mu+#mu- Pt distribution",200,0,20,-2);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Y","MCINPUT #mu+#mu- Y distribution",nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,"Eta","MCINPUT #mu+#mu- Eta distribution",nbinsEta,etaMin,etaMax);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Pt","MCINPUT #mu+#mu- Pt distribution",nMinvBins,minvMax,minvMax,-2);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Y","MCINPUT #mu+#mu- Y distribution",nbinsRapidity,rapidityMin,rapidityMax,-2);

  CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),"Eta","MCINPUT #mu+#mu- Eta distribution",nbinsEta,etaMin,etaMax);

  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"PtRecVsSim","#mu+#mu- Pt distribution rec vs sim",200,0,20,200,0,20);

  //  CreatePairHistos(eventSelection,triggerClassName,centrality,"BinFlowPt","#mu+#mu- BinFlowPt distribution",200,0,20);

  // Multiplicity Histo
  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"NchForJpsi","Corrected multiplicity distribution for 2.9 < m_{#mu^{+}#mu^{-}} < 3.3",nbinsMult,multMin,multMax);

  CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"NchForPsiP","Corrected multiplicity distribution for 3.6 < m_{#mu^{+}#mu^{-}} < 3.9",nbinsMult,multMin,multMax);

  TIter next(fBinsToFill);
  AliAnalysisMuMuBinning::Range* r;
  Int_t nb(0);

  // Minv Histos for each bin
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) ){
    // Histos name
    TString minvName(GetMinvHistoName(*r,kFALSE));
    ++nb;

    // Make sure histo is wanted
    if ( !IsHistogramDisabled(minvName.Data()) ){
      AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));

      // Reconstructed pair histo
      CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                       Form("#mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",r->AsString().Data()),nMinvBins,minvMin,minvMax,-2);
      // Generated J/psi histo
      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                        Form("MCINPUT #mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",r->AsString().Data()),nMCMinvBins,minvMin,minvMax,-2);

      CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),minvName.Data(),
                        Form("MCINPUT #mu+#mu- inv. mass %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",r->AsString().Data()),nMCMinvBins,minvMin,minvMax,-2); // Pure MC histo

      // Mean pt minv histo
      if ( fcomputeMeanPt ){
        TString mPtName(Form("MeanPtVs%s",minvName.Data()));
        // Reconstructed pair histo
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                         Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        // Generated J/psi histo
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                          Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);

        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),mPtName.Data(),
                          Form("#mu+#mu- mean p_{T} %s;M_{#mu^{+}#mu^{-}} (GeV/c^2);<p_{T}^{#mu^{+}#mu^{-} (GeV/c^2)}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
      }
    }

    // Create corrected histo
    if ( ShouldCorrectDimuonForAccEff() )
    {
      minvName = GetMinvHistoName(*r,kTRUE);

      if ( !IsHistogramDisabled(minvName.Data()) ){

        AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
        // Reconstructed pair histo
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                         Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}}(GeV/c^{2});Counts",r->AsString().Data()),nMinvBins,minvMin,minvMax,-2);
        // Generated J/psi histo
        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,centrality,minvName.Data(),
                          Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",r->AsString().Data()),nMCMinvBins,minvMin,minvMax,-2);

        CreateEventHistos(kHistoForMCInput,eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),minvName.Data(),
                          Form("#mu+#mu- inv. mass %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});Counts",r->AsString().Data()),nMCMinvBins,minvMin,minvMax,-2);
        // Mean pt accxeff corrected
        if ( fcomputeMeanPt ){
          TString mPtName(Form("MeanPtVs%s",minvName.Data()));

          CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mPtName.Data(),
                           Form("#mu+#mu- mean p_{T} %s (Acc #times Eff Corrected);M_{#mu^{+}#mu^{-}} (GeV/c^{2});<p_{T}^{#mu^{+}#mu^{-}}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::DefineMinvRange(Double_t minvMin, Double_t minvMax, Double_t minvBinSize)
{
  /// Define the Minv histogram range

  fMinvMin     = minvMin;
  fMinvMax     = minvMax;
  fMinvBinSize = minvBinSize;
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::FillHistosForPair(const char* eventSelection,
                                            const char* triggerClassName,
                                            const char* centrality,
                                            const char* pairCutName,
                                            const AliVParticle& tracki,
                                            const AliVParticle& trackj)
{
  /// Fill histograms for unlike-sign reconstructed  muon pairs.
  /// For the MC case, we check that only tracks with an associated MC label are selected (usefull when running on embedding).
  /// A weight is also applied for MC case at the pair or the muon track level according to SetMuonWeight() and systLevel.

  // Usual cuts
  if ( ( tracki.Charge() == trackj.Charge() ) ) return;
  if (!AliAnalysisMuonUtility::IsMuonTrack(&tracki) ) return;
  if (!AliAnalysisMuonUtility::IsMuonTrack(&trackj) ) return;

  // Pointers in case running on MC
  Int_t labeli               = 0;
  Int_t labelj               = 0;
  AliVParticle               * mcTracki(0x0);
  AliVParticle               * mcTrackj(0x0);
  TLorentzVector             * pair4MomentumMC(0x0);
  Double_t inputWeightMC(1.);
  AliMergeableCollectionProxy* mcProxy(0x0);

  // Make sure we have an associated tracks in simulation stack if running on MC
  if(HasMC()){
    // Get label
    labeli = tracki.GetLabel();
    labelj = trackj.GetLabel();
    if ( labeli < 0 || labelj < 0 ){
      AliError("Got negative labels!");
      return;
    }

    //Check if first track is a muon
    mcTracki = MCEvent()->GetTrack(labeli);
    if ( TMath::Abs(mcTracki->PdgCode()) != 13 ) return;

    //Check if second track is a muon
    mcTrackj = MCEvent()->GetTrack(labelj);
    if ( TMath::Abs(mcTrackj->PdgCode()) != 13 ) return;

    //Check if tracks has the same mother
    Int_t currMotheri = mcTracki->GetMother();
    Int_t currMotherj = mcTrackj->GetMother();
    if(currMotheri!=currMotherj) return;
    if(currMotheri<0) return;

    // Check if mother is J/psi
    AliMCParticle* mother = static_cast<AliMCParticle*>(MCEvent()->GetTrack(currMotheri));
    if(mother->PdgCode() !=443) return;

    // Create proxy for MC
    mcProxy = HistogramCollection()->CreateProxy(BuildMCPath(eventSelection,triggerClassName,centrality,pairCutName));
    TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTracki->P()*mcTracki->P()));
    TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTrackj->P()*mcTrackj->P()));
    mcpj+=mcpi;

    // Weight tracks if specified
    if(!fWeightMuon)      inputWeightMC = WeightPairDistribution(mother->Pt(),mother->Y());
    else if(fWeightMuon)  inputWeightMC = WeightMuonDistribution(mcTracki->Px()) * WeightMuonDistribution(mcTrackj->Px());
  }

  // Construct dimuons vector
  TLorentzVector pi(tracki.Px(),tracki.Py(),tracki.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki.P()*tracki.P()));
  TLorentzVector pair4Momentum(trackj.Px(),trackj.Py(),trackj.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj.P()*trackj.P()));
  pair4Momentum += pi;

  // Create proxy in AliMergeableCollection
  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,pairCutName));

  // Weight tracks if specified
  Double_t inputWeight=0.;
  if(!fWeightMuon)      inputWeight = WeightPairDistribution(pair4Momentum.Pt(),pair4Momentum.Rapidity());
  else if(fWeightMuon)  inputWeight = WeightMuonDistribution(tracki.Pt()) * WeightMuonDistribution(trackj.Pt());

  // Fill some distribution histos
  if ( !IsHistogramDisabled("Pt") )  proxy->Histo("Pt")->Fill(pair4Momentum.Pt(),inputWeight);
  if ( !IsHistogramDisabled("Y") )   proxy->Histo("Y")->Fill(pair4Momentum.Rapidity(),inputWeight);
  if ( !IsHistogramDisabled("Eta") ) proxy->Histo("Eta")->Fill(pair4Momentum.Eta());
  if ( !IsHistogramDisabled("PtPaireVsPtTrack") ) {
    static_cast<TH2*>( proxy->Histo("PtPaireVsPtTrack"))->Fill(pair4Momentum.Pt(),tracki.Pt(),inputWeight);
    static_cast<TH2*>( proxy->Histo("PtPaireVsPtTrack"))->Fill(pair4Momentum.Pt(),trackj.Pt(),inputWeight);
  }

  // Fill histos with MC stack info
  if ( HasMC() ){
    // Get 4-vector pairs from MC stack
    TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTracki->P()*mcTracki->P()));
    TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTrackj->P()*mcTrackj->P()));
    mcpj+=mcpi;

    // Fill histo
    proxy->Histo("PtRecVsSim")->Fill(mcpj.Pt(),pair4Momentum.Pt());
    if ( !IsHistogramDisabled("Pt") )  mcProxy->Histo("Pt")->Fill(mcpj.Pt(),inputWeightMC);
    if ( !IsHistogramDisabled("Y") )   mcProxy->Histo("Y")->Fill(mcpj.Rapidity(),inputWeightMC);
    if ( !IsHistogramDisabled("Eta") ) mcProxy->Histo("Eta")->Fill(mcpj.Eta());

    // set pair4MomentumMC for the rest of the function
    pair4MomentumMC = &mcpj;
  }

  TIter nextBin(fBinsToFill);
  AliAnalysisMuMuBinning::Range* r;

  // Loop over all bin ranges
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ){

    //In this loop we first check if the pairs pass some tests and we fill histo accordingly.

    // Flag for cuts and ranges
    Bool_t ok(kFALSE);
    Bool_t okMC(kFALSE);

    //Fully integrated case
    if ( r->IsIntegrated() ){

      ok = kTRUE;
      if ( pair4MomentumMC ) okMC = kTRUE;

      TH1* h(0x0);

      // Fill NchForJpsi histo
      if ( pair4Momentum.M() >= 2.9 && pair4Momentum.M() <= 3.3 ){

        h = proxy->Histo("NchForJpsi");

        Double_t ntrcorr = (-1.);
        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));

        if (list){
          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it
            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("NtrCorr") ){
              parFound = kTRUE;
              ntrcorr = p->GetVal();
            }
          }
        }
        h->Fill(ntrcorr);
      }
      else if ( pair4Momentum.M() >= 3.6 && pair4Momentum.M() <= 3.9){

        h = proxy->Histo("NchForPsiP");
        Double_t ntrcorr = (-1.);

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){

          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){
            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++; // In case there is a different object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("NtrCorr") ){
              parFound = kTRUE;
              ntrcorr = p->GetVal();
            }
          }
        }
        h->Fill(ntrcorr);
      }
    }
    // 2D Binning
    else if ( r->Is2D() ){

      if ( r->AsString().BeginsWith("PTVSY") ){
        ok                          = r->IsInRange(pair4Momentum.Rapidity(),pair4Momentum.Pt());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity(),pair4MomentumMC->Pt());
      }
      else if ( r->AsString().BeginsWith("YVSPT") ){
        ok                          = r->IsInRange(pair4Momentum.Pt(),pair4Momentum.Rapidity());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity());
      }
      else if ( r->Quantity() == "NTRCORRPT" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){
          Int_t i(-1);
          Bool_t parFound(kFALSE);

          while ( i < list->GetEntries() - 1 && !parFound ){
            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("NtrCorr") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal(),pair4Momentum.Pt());
            }
          }
        }
      }
      else if ( r->Quantity() == "NTRCORRY" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){
          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("NtrCorr") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal(),pair4Momentum.Rapidity());
            }
          }
        }
      }
      else AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
    }
    // The rest
    else{

      if ( r->Quantity() == "PT" ){
        ok                          = r->IsInRange(pair4Momentum.Pt());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt());
      }
      else if ( r->Quantity() == "Y" ){
        ok                          = r->IsInRange(pair4Momentum.Rapidity());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity());
      }
      else if ( r->Quantity() == "PHI" ){
        ok                          = r->IsInRange(pair4Momentum.Phi());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Phi());
      }
      else if ( r->Quantity() == "DNCHDETA" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){
          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("dNchdEta") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal());
            }
          }
        }
      }
      else if ( r->Quantity() == "NTRCORR" || r->Quantity() == "RELNTRCORR" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){

          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("NtrCorr") ){
              parFound = kTRUE;
              if ( r->Quantity() == "NTRCORR" ) ok = r->IsInRange(p->GetVal());
              else ok = r->IsInRange(p->GetVal()/5.97);
            }
          }
        }
      }
      else if ( r->Quantity() == "V0ACORR" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){

          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("V0ACorr") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal());
            }
          }
        }
      }
      else if ( r->Quantity() == "V0CCORR" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){

          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 )  i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("V0CCorr") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal());
            }
          }
        }
      }
      else if ( r->Quantity() == "V0MCORR" ){

        TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
        if (list){

          Int_t i(-1);
          Bool_t parFound(kFALSE);
          while ( i < list->GetEntries() - 1 && !parFound ){

            i++;
            while ( list->At(i)->IsA() != TParameter<Double_t>::Class() && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it

            TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(list->At(i));

            if ( TString(p->GetName()).Contains("V0MCorr") ){
              parFound = kTRUE;
              ok = r->IsInRange(p->GetVal());
            }
          }
        }
      }
    }

    // Check if pair pass all conditions, either MC or not, and fill Minv Histogrames
    if ( ok || okMC ){

      // Get Minv histo name associated to the bin
      TString minvName = GetMinvHistoName(*r,kFALSE);

      //Create, fill and store Minv histo
      if (!IsHistogramDisabled(minvName.Data())){

        TH1* h(0x0);

        if ( ok ){
          h = proxy->Histo(minvName.Data());
          if (!h) AliError(Form("Could not get %s",minvName.Data()));
          else h->Fill(pair4Momentum.M(),inputWeight);
        }

        if( okMC ){
          h = mcProxy->Histo(minvName.Data());
          if (!h) AliError(Form("Could not get MC %s",minvName.Data()));
          else h->Fill(pair4MomentumMC->M(),inputWeightMC);
        }

        // Fill Mean pT
        if ( fcomputeMeanPt ){
          TString hprofName("");

          if ( ok ){
            hprofName= Form("MeanPtVs%s",minvName.Data());
            TProfile* hprof = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
            if ( !hprof )AliError(Form("Could not get %s",hprofName.Data()));
            else hprof->Fill(pair4Momentum.M(),pair4Momentum.Pt(),inputWeight);
          }

          if ( okMC ){
            hprofName= Form("MeanPtVs%s",minvName.Data());
            TProfile* hprof = MCProf(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
            if ( !hprof )AliError(Form("Could not get MC %s",hprofName.Data()));
            else hprof->Fill(pair4MomentumMC->M(),pair4MomentumMC->Pt(),inputWeightMC);
          }
        }
      }

      // Create, fill and store Minv histo already corrected with accxeff
      if ( ShouldCorrectDimuonForAccEff() ){

        Double_t AccxEff(0);
        Bool_t okAccEff(kFALSE);

        // Protection
        if ( ok ){
          AccxEff = GetAccxEff(pair4Momentum.Pt(),pair4Momentum.Rapidity());
          if ( AccxEff <= 0.0 ) AliError(Form("AccxEff < 0 for pt = %f & y = %f ",pair4Momentum.Pt(),pair4Momentum.Rapidity()));
          else okAccEff = kTRUE;
        }

        Double_t AccxEffMC(0);
        Bool_t okAccEffMC(kFALSE);

        // Protection
        if ( okMC ){
          AccxEffMC= GetAccxEff(pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity());
          if ( AccxEffMC <= 0.0 ) AliError(Form("AccxEff < 0 for MC pair with pt = %f & y = %f ",pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity()));
          else okAccEffMC = kTRUE;
        }

        // Get histo name
        minvName = GetMinvHistoName(*r,kTRUE);

        // fill histo
        if (!IsHistogramDisabled(minvName.Data())){

          TH1* hCorr = proxy->Histo(minvName.Data());

          if (!hCorr) AliError(Form("Could not get %sr",minvName.Data()));
          else if ( okAccEff ) hCorr->Fill(pair4Momentum.M(),inputWeight/AccxEff);

          if( okAccEffMC ){
            hCorr = mcProxy->Histo(minvName.Data());
            if (!hCorr) AliError(Form("Could not get MC %s",minvName.Data()));
            else hCorr->Fill(pair4MomentumMC->M(),inputWeightMC/AccxEffMC);
          }

          if ( fcomputeMeanPt ){

            TString hprofCorrName("");

            if( ok ){
              hprofCorrName = Form("MeanPtVs%s",minvName.Data());
              TProfile* hprofCorr = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofCorrName.Data());
              if ( !hprofCorr ) AliError(Form("Could not get %s",hprofCorrName.Data()));
              else if ( okAccEff ) hprofCorr->Fill(pair4Momentum.M(),pair4Momentum.Pt(),inputWeight/AccxEff);
            }

            if( okMC ){
              hprofCorrName = Form("MeanPtVs%s",minvName.Data());
              TProfile* hprofCorr = MCProf(eventSelection,triggerClassName,centrality,pairCutName,hprofCorrName.Data());
              if ( !hprofCorr ) AliError(Form("Could not get MC %s",hprofCorrName.Data()));
              else if ( okAccEffMC )hprofCorr->Fill(pair4MomentumMC->M(),pair4MomentumMC->Pt(),inputWeightMC/AccxEffMC);
            }
          }
        }
      }
    }
  }
  delete proxy;
  delete mcProxy;
}


//_____________________________________________________________________________
void AliAnalysisMuMuMinv::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  ///
  /// Fill MC inputs histograms.
  ///

  if ( !HasMC() ) return;

  // Create general proxies to the Histogram Collection
  TString mcPath = BuildMCPath(eventSelection,triggerClassName,centrality);
  AliMergeableCollectionProxy* mcProxy = HistogramCollection()->CreateProxy(mcPath);

  // Create proxy to the Histogram Collection for input particles satisfying Y cut
  TString mcInYRangeProxyPath = mcPath;
  mcInYRangeProxyPath += "INYRANGE/";
  AliMergeableCollectionProxy* mcInYRangeProxy = HistogramCollection()->CreateProxy(mcInYRangeProxyPath);
  if(!mcInYRangeProxy) printf("Warning :  unable to create proxy for mcInYRangeProxy, will not be filled \n");

  // number of tracks in Event
  Int_t nMCTracks = MCEvent()->GetNumberOfTracks();

  TIter nextBin(fBinsToFill);
  AliAnalysisMuMuBinning::Range* r;

  // Loop over all events
  for ( Int_t i = 0; i < nMCTracks; ++i ){
    // Get particle
    AliVParticle* part = MCEvent()->GetTrack(i);

    // Select only primary particles
    if  (AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) && part->GetMother()==-1){

      //pt cut. 0-12 by default
      if(part->Pt()>fmcptcutmax || part->Pt()<fmcptcutmin  ) return;

      // Get the default WeightPairDistribution
      Double_t inputWeight = WeightPairDistribution(part->Pt(),part->Y());

      // Fill Pt, Y, Eta histos
      mcProxy->Histo("Pt")->Fill(part->Pt(),inputWeight);
      mcProxy->Histo("Y")->Fill(part->Y(),inputWeight);
      mcProxy->Histo("Eta")->Fill(part->Eta());

      // Fill Pt, Y, Eta histos if tracks rapidity in range
      if ( part->Y() < -2.5 && part->Y() > -4.0 && mcInYRangeProxy ){
        mcInYRangeProxy->Histo("Pt")->Fill(part->Pt(),inputWeight);
        mcInYRangeProxy->Histo("Y")->Fill(part->Y(),inputWeight);
        mcInYRangeProxy->Histo("Eta")->Fill(part->Eta());
      }

      nextBin.Reset();

      // Loop on all range in order to fill Histo
      while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ){

        // Check if particles pass all the cuts for different bins
        Bool_t ok(kFALSE);

        // always true for integrated
        if ( r->IsIntegrated() )ok = kTRUE;
        // 2D Binning
        else if ( r->Is2D() ){
          if ( r->AsString().BeginsWith("PTVSY") )     ok = r->IsInRange(part->Y(),part->Pt());
          else if ( r->AsString().BeginsWith("YVSPT") )ok = r->IsInRange(part->Pt(),part->Y());
          else AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
        }
        // 1D Binning
        else {
          if ( r->Quantity()      == "PT" ) ok = r->IsInRange(part->Pt());
          else if ( r->Quantity() == "Y" )  ok = r->IsInRange(part->Y());
          else if ( r->Quantity() == "PHI" )ok = r->IsInRange(part->Phi());
        }

        // Fill Minv histo if bin is in range
        if ( ok ){

          // Get histo name
          TString hname = GetMinvHistoName(*r,kFALSE);

          // Chek if histo disabled
          if (!IsHistogramDisabled(hname.Data())){
            TH1* h = mcProxy->Histo(hname.Data());
            if (!h) {
              AliError(Form("Could not get /%s/%s/%s/%s/ %s",MCInputPrefix(),eventSelection,triggerClassName,centrality,hname.Data()));
              continue;
            }
            h->Fill(part->M(),inputWeight);

            if ( part->Y() < -2.5 && part->Y() > -4.0 && mcInYRangeProxy  ){
              h = mcInYRangeProxy->Histo(hname.Data());
              if (!h){
                AliError(Form("Could not get /%s/%s/%s/%s/INYRANGE %s",MCInputPrefix(),eventSelection,triggerClassName,centrality,hname.Data()));
                continue;
              }
              h->Fill(part->M(),inputWeight);
            }
          }

          // Fill compute mean pt histo
          if ( fcomputeMeanPt ){

            TString hprofName = Form("MeanPtVs%s",hname.Data());
            TProfile* hprof   = MCProf(eventSelection,triggerClassName,centrality,hprofName.Data());

            if ( !hprof )AliError(Form("Could not get %s",hprofName.Data()));
            else hprof->Fill(part->M(),part->Pt(),inputWeight);

            if ( part->Y() < -2.5 && part->Y() > -4.0 ){
              hprof = MCProf(eventSelection,triggerClassName,Form("%s/INYRANGE",centrality),hprofName.Data());
              if ( !hprof )AliError(Form("Could not get %s",hprofName.Data()));
              else hprof->Fill(part->M(),part->Pt(),inputWeight);
            }
          }
        }
      }
    }
  }

  delete mcProxy;
  delete mcInYRangeProxy;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuMinv::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const
{
  return TString::Format("MinvUS%s%s%s",
                         accEffCorrected ? "_AccEffCorr" : "",fMinvBinSeparator.Data(),r.AsString().Data());
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuMinv::GetAccxEff(Double_t pt,Double_t rapidity)
{
  if (!fAccEffHisto){
    AliError("ERROR: No AccxEff histo");
    return 0;
  }
  Int_t bin        = fAccEffHisto->FindBin(pt,rapidity);
  Double_t accXeff = fAccEffHisto->GetBinContent(bin);

  return accXeff;
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuMinv::WeightMuonDistribution(Double_t pt)
{
  ///Return a weight for a single pt and y, which depend on the varied distributions.
  // FIXME: hard coded, find a clean way to fix the distribution parameters from outside

  if (!HasMC() ) return 1.;
  if (!fWeightMuon) return 1.;


  //================ Trigger Efficiency distribution pp@5TeV ==============//
  //value for input distribution
  Double_t parptpp5 [8]     = {6.24341e-01,2.68680e-01,6.57929e-01,1.00000e+00,1.00000e+00,1.00000e+00,-1.00000e+00,9.67072e-01}; //Initial dunction
  Double_t parptpp5corr [8] = {6.24341e-01,2.70085e-01,5.85241e-01,1.00000e+00,1.00000e+00,1.00000e+00,-1.00000e+00,9.67072e-01}; //Corrected one

  //______
  Double_t weight(1.),funcptsyst(0.);
  switch ( fsystLevel )
  {
    case 0:
      weight = 1;
      break;
    case 1:
      funcptsyst = TriggerLptApt(&pt,parptpp5);
      if ( funcptsyst > 0 ) weight = funcptsyst;
      else  weight = 1;
      break;
    case 2:
      funcptsyst = TriggerLptApt(&pt,parptpp5corr);
      if ( funcptsyst > 0 ) weight = funcptsyst;
      else  weight = 1;
      break;
  }

  return weight;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuMinv::WeightPairDistribution(Double_t pt,Double_t rapidity)
{
  //Return a weight for a dimuon pt and y, which depend on the varied distributions.
  // FIXME: hard coded, find a clean way to fix the distribution parameters from outside

  if ( !HasMC() )    return 1.;
  if ( fWeightMuon ) return 1.;
  if ( !fPtFuncOld ) return 1.;
  if ( !fPtFuncNew ) return 1.;
  if ( !fYFuncOld )  return 1.;
  if ( !fYFuncNew )  return 1.;

  Double_t weight = fPtFuncNew->Eval(pt) / fPtFuncOld->Eval(pt) * fYFuncNew->Eval(rapidity) / fYFuncOld->Eval(rapidity);
  return weight;
}

//______________________________________________
Double_t AliAnalysisMuMuMinv::TriggerLptApt ( Double_t* xVal, Double_t* par )
{
  // trigger response function
  Double_t xx             = xVal[0];
  Double_t currX          = TMath::Max(xx,par[6]);
  Double_t sqrtTwo        = TMath::Sqrt(2.);
  Double_t yVal           = par[7]+par[0]*(TMath::Erf((currX-par[1])/par[2]/sqrtTwo)-1.);
  if ( xx < par[6] ) yVal += par[3]*(TMath::Erf((-xx-par[4])/par[5]/sqrtTwo) - TMath::Erf((-par[6]-par[4])/par[5]/sqrtTwo));

  return yVal;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuMinv::IsPtInRange(const AliVParticle& t1, const AliVParticle& t2, Double_t& ptmin, Double_t& ptmax) const
{
  /// Whether the pair passes the pT cut

  TLorentzVector total(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));

  total += p2;

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
  /// Whether the pair passes the rapidity cut
  TLorentzVector total(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));

  total += p2;

  Double_t y = total.Rapidity();

  return  ( y < -2.5 && y > -4.0 );
}

//_____________________________________________________________________________
void AliAnalysisMuMuMinv::SetBinsToFill(const char* particle, const char* bins)
{
  delete fBinsToFill;
  fBinsToFill = Binning()->CreateBinObjArray(particle,bins,"");
}

//________________________________________________________________________
void AliAnalysisMuMuMinv::SetOriginPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.

  assert(param);

  delete fPtFuncOld;
  fPtFuncOld = new TF1("fPtFuncOld", formula.Data(), xMin, xMax);
  fPtFuncOld->SetParameters(param);
  NormFunc(fPtFuncOld, xMin, xMax);

}

//________________________________________________________________________
void AliAnalysisMuMuMinv::SetNewPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted pT distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.

  assert(param);

  delete fPtFuncNew;
  fPtFuncNew = new TF1("fPtFuncNew", formula.Data(), xMin, xMax);
  fPtFuncNew->SetParameters(param);
  NormFunc(fPtFuncNew, xMin, xMax);

}

//________________________________________________________________________
void AliAnalysisMuMuMinv::SetOriginYFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
{
  /// Create the original function with the parameters used in simulation to generate the y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution for cross-check.

  assert(param);

  delete fYFuncOld;
  fYFuncOld = new TF1("fYFuncOld", formula.Data(), xMin, xMax);
  fYFuncOld->SetParameters(param);
  NormFunc(fYFuncOld, xMin, xMax);

}

//________________________________________________________________________
void AliAnalysisMuMuMinv::SetNewYFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax)
{
  /// Create the new function with its initial parameters to fit the generated/weighted y distribution.
  /// It must be in the form [0]*(...) to allow for global normalization.
  /// The [xMin,xMax] range is used to normalized the function.
  /// Some parameters can be fixed when fitting the generated distribution.

  assert(param);

  delete fYFuncNew;
  fYFuncNew = new TF1("fYFuncNew", formula.Data(), xMin, xMax);
  fYFuncNew->SetParameters(param);
  NormFunc(fYFuncNew, xMin, xMax);

}

//________________________________________________________________________
void AliAnalysisMuMuMinv::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
   f->SetNpx(100.*(max-min));
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}
