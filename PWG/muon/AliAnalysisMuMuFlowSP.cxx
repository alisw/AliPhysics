#include "AliAnalysisMuMuFlowSP.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuFlowSP
 *
 * Analysis which fills a bunch of histograms for invariant mass analysis of J/psi
 * The flow vectors are recovered from Qn correction framework (PWGPP/EVCHAR/FlowVectorCorrections)
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
#include "TVector2.h"
#include "TString.h"
#include "AliMCEvent.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"
#include "TParameter.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include <cassert>

ClassImp(AliAnalysisMuMuFlowSP)

//_____________________________________________________________________________
AliAnalysisMuMuFlowSP::AliAnalysisMuMuFlowSP(TH2* accEffHisto, Int_t systLevel)
: AliAnalysisMuMuBase(),
fWeightMuon(kFALSE),
fAccEffHisto(0x0),
fMinvBinSeparator("+"),
fsystLevel(systLevel),
fPtFuncOld(0x0),
fPtFuncNew(0x0),
fYFuncOld(0x0),
fYFuncNew(0x0),
fBinsToFill(0x0),
fMinvBinSize(0.025),
fMinvMin(0.0),
fMinvMax(16.0),
// fmcptcutmin(0.0),
// fmcptcutmax(12.0),
fNDetectors(3)
{
  // FIXME ? find the AccxEff histogram from HistogramCollection()->Histo("/EXCHANGE/JpsiAccEff")

  fHar = 2;
  if ( accEffHisto )
  {
    fAccEffHisto = static_cast<TH2F*>(accEffHisto->Clone());
    fAccEffHisto->SetDirectory(0);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuFlowSP::~AliAnalysisMuMuFlowSP()
{
  /// dtor
  delete fAccEffHisto;
  delete fBinsToFill;
}

//_____________________________________________________________________________
void
AliAnalysisMuMuFlowSP::DefineHistogramCollection(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality,
                                               Bool_t mix)
{
  /// Define the histograms this analysis will use

  // Check if histo is not already here
  if ( ExistSemaphoreHistogram(eventSelection,triggerClassName,centrality) )
  {
    return;
  }

  CreateSemaphoreHistogram(eventSelection,triggerClassName,centrality);

  // no bins defined by the external steering macro, use our own defaults
  //if (!fBinsToFill) SetBinsToFill("psi","integrated,ptvsy,yvspt,pt,y,phi");
  if (!fBinsToFill) SetBinsToFill("psi","integrated,pt,y,ptvsy,yvspt");

  // mass range
  Double_t minvMin = fMinvMin;
  Double_t minvMax = fMinvMax;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,fMinvBinSize);

  // Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);

 //  // Rapidity range
 //  Double_t rapidityMin = -5;
 //  Double_t rapidityMax = -2;
 //  Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);

 //  // eta range
 //  Double_t etaMin = -5;
 //  Double_t etaMax = -2;
 //  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);

 // //Multiplicity range
 //  Double_t multMin = -0.5;
 //  Double_t multMax = 500.5;
 //  Int_t nbinsMult = GetNbins(multMin,multMax,1.);


  TIter next(fBinsToFill);
  AliAnalysisMuMuBinning::Range* r;
  Int_t nb(0);


  for(Int_t i=0; i<fNDetectors; i++){
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qn_%s",fDetectors[i].Data()),Form("Qn vector from %s; Qn_{x};Qn_{y}",fDetectors[i].Data()),
                     600, -1., 1., 600, -1., 1.);
    for(Int_t j=i+1; j<fNDetectors;j++){
      CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("R_%svs%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- resolution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                         600, -2., 2.);//dphi corrected to be in [O,pi]
    }
  }
  // if( ShouldCorrectDimuonForAccEff())
  // {
  //   CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"U_AccEffCorr",
  //                          "#mu+#mu- Q vector (Acc #times Eff Corrected);u_{x};u_{y}",600, -0.01, 3.2, 600, -0.01, 3.2);
  // }

  // Minv Histos for each bin
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    // Histos name
    TString minvName(GetMinvHistoName(*r,kFALSE));
    ++nb;
    // Make sure histo is wanted
    if ( !IsHistogramDisabled(minvName.Data()) ){
      AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));

      TString mUName(Form("U_%s",minvName.Data()));
      CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mUName.Data(),
                         "#mu+#mu- Q vector ;u_{x};u_{y}",600, -1., 1., 600, -1., 1.);

      for(Int_t i=0; i<fNDetectors; i++){
        TString mSPName(Form("SP_%s_vs%s",fDetectors[i].Data(),minvName.Data()));
        // Reconstructed pair histo
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mSPName.Data(),
                         Form("#mu+#mu- v_{2} {SP} with %s vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);

        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPName.Data()),
                         Form("#mu+#mu- v_{2} {SP} (corrected with |Qn|) with %s vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
      }

    }
    // Create corrected histo
    if ( ShouldCorrectDimuonForAccEff() )
    {
      minvName = GetMinvHistoName(*r,kTRUE);

      if ( !IsHistogramDisabled(minvName.Data()) ){
        AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));

        TString mUAEName(Form("U_%s",minvName.Data()));
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mUAEName.Data(),
                         "#mu+#mu- Q vector (Acc #times Eff Corrected) ;u_{x};u_{y}",600, -1., 1., 600, -1., 1.);
        //scalar product
         for(Int_t i=0; i<fNDetectors; i++){
          TString mSPAEName(Form("SP_%s_vs%s",fDetectors[i].Data(),minvName.Data()));
          CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,mSPAEName.Data(),
            Form("#mu+#mu- v_{2} {SP} with %s (Acc #times Eff Corrected) vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
          CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPAEName.Data()),
            Form("#mu+#mu- v_{2} {SP} (corrected with |Qn|) with %s (Acc #times Eff Corrected) vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlowSP::DefineMinvRange(Double_t minvMin, Double_t minvMax, Double_t minvBinSize)
{
  /// Define the Minv histogram range

  fMinvMin     = minvMin;
  fMinvMax     = minvMax;
  fMinvBinSize = minvBinSize;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlowSP::FillHistosForPair(const char* eventSelection,
                                            const char* triggerClassName,
                                            const char* centrality,
                                            const char* pairCutName,
                                            const AliVParticle& tracki,
                                            const AliVParticle& trackj,
                                            const Bool_t IsMixedHisto)
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
  TVector2 UMC(0.,0.);
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

  TVector2 U(cos(fHar*pair4Momentum.Phi()),sin(fHar*pair4Momentum.Phi()));//Unitary Q vector of the dimuon

  // Create proxy in AliMergeableCollection
  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,pairCutName));

  // // Weight tracks if specified
  Double_t inputWeight=0.;
  if(!fWeightMuon)      inputWeight = WeightPairDistribution(pair4Momentum.Pt(),pair4Momentum.Rapidity());
  else if(fWeightMuon)  inputWeight = WeightMuonDistribution(tracki.Pt()) * WeightMuonDistribution(trackj.Pt());

  // Fill some distribution histos
  TVector2 Qn[3];//Q vectors (2nd harmonic) for each detector
  Double_t SP[3];//Scalar product

  // if ( !IsHistogramDisabled("U") ) proxy->Histo("U")->Fill(U.X(),U.Y());
  for(Int_t i=0; i<fNDetectors; i++){
    if(i==0) Qn[i]= GetQn(fDetectors[i].Data(),4); //twist for SPD
    else Qn[i]= GetQn(fDetectors[i].Data(),3); //twist

    SP[i] = U*Qn[i];
  }


  // Fill histos with MC stack info
  if ( HasMC() ){
    AliWarning("MC is not implemented for flow analysis");
    // Get 4-vector pairs from MC stack
    TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTracki->P()*mcTracki->P()));
    TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTrackj->P()*mcTrackj->P()));
    mcpj+=mcpi;


    // set pair4MomentumMC for the rest of the function
    pair4MomentumMC = &mcpj;
    UMC.Set(cos(fHar*pair4MomentumMC->Phi()),sin(fHar*pair4MomentumMC->Phi()));//Unitary Q vector of the dimuon
    Double_t SPMC[3];//Scalar product

    if ( !IsHistogramDisabled("U") )mcProxy->Histo("U")->Fill(U.X(),U.Y());

    for(Int_t i=0; i<fNDetectors; i++){
      SPMC[i] = UMC*Qn[i];
    }
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
    }
    // 2D Binning
    else if ( r->Is2D() ){
      if (r->AsString().BeginsWith("PTVSY") )
      {
        ok = (r->IsInRange(pair4Momentum.Pt(),pair4Momentum.Rapidity()));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),pair4MomentumMC->Rapidity());
      }
      if (r->AsString().BeginsWith("YVSPT") )
      {
        ok = (r->IsInRange(pair4Momentum.Rapidity(),pair4Momentum.Pt()));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity(),pair4MomentumMC->Pt());
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

    }

    // Check if pair pass all conditions, either MC or not, and fill Minv Histogrames
    if ( ok || okMC ){

      // Get Minv histo name associated to the bin
      TString minvName = GetMinvHistoName(*r,kFALSE);
      //Create, fill and store Minv histo
      if (!IsHistogramDisabled(minvName.Data())){

        if ( ok ){

          TString hprofName("");
          for(Int_t i=0; i<fNDetectors;i++){
            hprofName = Form("SP_%s_vs%s",fDetectors[i].Data(),minvName.Data());
            TProfile* hprof = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
            if ( !hprof)AliError(Form("Could not get %s",hprofName.Data()));
            else hprof->Fill(pair4Momentum.M(),SP[i],inputWeight);
            TProfile* hprofcorr = Prof(eventSelection,triggerClassName,centrality,pairCutName,Form("%s_corr",hprofName.Data()));
            if ( !hprofcorr)AliError(Form("Could not get %s",Form("%s_corr",hprofName.Data())));
            else hprofcorr->Fill(pair4Momentum.M(),SP[i]-Qn[i].X()*Qn[i].X()+Qn[i].Y()*Qn[i].Y(),inputWeight);
          }
          TString mUName(Form("U_%s",minvName.Data()));
          if ( !IsHistogramDisabled(mUName.Data() )) proxy->Histo(mUName.Data())->Fill(U.X(),U.Y());
        }

        if( okMC ){

          TString hprofMCName("");
          for(Int_t i=0; i<fNDetectors;i++){
            hprofMCName= Form("SP_MC_%s_vs%s",fDetectors[i].Data(),minvName.Data());
            TProfile* hprofMC = Prof(eventSelection,triggerClassName,centrality,hprofMCName.Data());
            if ( !hprofMC)AliError(Form("Could not get %s",hprofMCName.Data()));
            else hprofMC->Fill(pair4Momentum.M(),SP[i],inputWeight);
            TProfile* hprofMCcorr = Prof(eventSelection,triggerClassName,centrality,Form("%s_corr",hprofMCName.Data()));
            if ( !hprofMCcorr)AliError(Form("Could not get %s",Form("%s_corr",hprofMCName.Data())));
            else hprofMCcorr->Fill(pair4MomentumMC->M(),sqrt(SP[i]-Qn[i].X()*Qn[i].X()+Qn[i].Y()*Qn[i].Y()),inputWeight);//Correct
          }

          TString mUMCName(Form("UMC_%s",minvName.Data()));
          if ( !IsHistogramDisabled(mUMCName.Data() )) proxy->Histo(mUMCName.Data())->Fill(UMC.X(),UMC.Y());
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

          if ( ok ){
            TString hprofNameAE("");
            for(Int_t i=0; i<fNDetectors;i++){
              hprofNameAE= Form("SP_%s_vs%s",fDetectors[i].Data(),minvName.Data());
              TProfile* hprof = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofNameAE.Data());
              if ( !hprof)AliError(Form("Could not get SP_%s_vs%s",fDetectors[i].Data(),minvName.Data()));
              else hprof->Fill(pair4Momentum.M(),SP[i],inputWeight/AccxEff);
            }
            TString mUName(Form("U_%s",minvName.Data()));
            if ( !IsHistogramDisabled(mUName.Data() )) proxy->Histo(mUName.Data())->Fill(U.X(),U.Y());
          }

          if( okAccEffMC ){
            TString hprofMCNameAE("");
            for(Int_t i=0; i<fNDetectors;i++){
              hprofMCNameAE= Form("SP_%s_vs%s",fDetectors[i].Data(),minvName.Data());
              TProfile* hprofMC = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofMCNameAE.Data());
              if ( !hprofMC)AliError(Form("Could not get SP_%s_vs%s",fDetectors[i].Data(),minvName.Data()));
              else hprofMC->Fill(pair4MomentumMC->M(),SP[i],inputWeight/AccxEff);
            }
            TString mUAEMCName(Form("UMC_%s",minvName.Data()));
            if ( !IsHistogramDisabled(mUAEMCName.Data() )) proxy->Histo(mUAEMCName.Data())->Fill(UMC.X(),UMC.Y());
          }
        }
      }
    }
  }
  delete proxy;
  delete mcProxy;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlowSP::FillHistosForEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
    // Create proxy in AliMergeableCollection
  AliMergeableCollectionProxy* proxyEv = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality));

  TVector2 Qn[3];//Q vectors (2nd harmonic) for each detector
  for(Int_t i=0; i<fNDetectors; i++){
    if(i==0) Qn[i]= GetQn(fDetectors[i].Data(),4); //twist for SPD
    else Qn[i]= GetQn(fDetectors[i].Data(),3); //twist

    if ( !IsHistogramDisabled(Form("Qn_%s",fDetectors[i].Data())) ) proxyEv->Histo(Form("Qn_%s",fDetectors[i].Data()))->Fill(Qn[i].X(),Qn[i].Y());

    for(Int_t j=i+1; j<fNDetectors;j++){
      if ( !IsHistogramDisabled(Form("R_%svs%s",fDetectors[i].Data(),fDetectors[j].Data())))  proxyEv->Histo(Form("R_%svs%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(Qn[i]*Qn[j]);
    }
  }

  delete proxyEv;
}
//_____________________________________________________________________________
void AliAnalysisMuMuFlowSP::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  ///
  /// Fill MC inputs histograms.
  ///

  if ( !HasMC() ) return;
  AliError("AliMuMuFlow does not deal with MC event yet, implement it !");
}

//_____________________________________________________________________________
TString AliAnalysisMuMuFlowSP::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const
{
  return TString::Format("MinvUS%s%s%s",
                         accEffCorrected ? "_AccEffCorr" : "",fMinvBinSeparator.Data(),r.AsString().Data());
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuFlowSP::GetAccxEff(Double_t pt,Double_t rapidity)
{
  if (!fAccEffHisto){
    AliError("ERROR: No AccxEff histo");
    return 0;
  }
  //valid for Mohammad's map
  Double_t xhisto = 4*rapidity+16;
  Double_t pthisto = -1.;
  if(pt < 0.3) {pthisto = 0.3*pt;}
  else if(pt < 1.) {pthisto = 2*pt - 0.3;}
  else if(pt < 6.) {pthisto = pt+1;}
  else if(pt < 8.) {pthisto = pt/2. + 4.;}
  else if(pt < 12.) {pthisto = pt/4. + 6.;}

  Int_t bin        = fAccEffHisto->FindBin(xhisto, pthisto);//pt,-rapidity);
  Double_t accXeff = fAccEffHisto->GetBinContent(bin);

  return accXeff;
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuFlowSP::WeightMuonDistribution(Double_t pt)
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
Double_t AliAnalysisMuMuFlowSP::WeightPairDistribution(Double_t pt,Double_t rapidity)
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
Double_t AliAnalysisMuMuFlowSP::TriggerLptApt ( Double_t* xVal, Double_t* par )
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
void AliAnalysisMuMuFlowSP::SetBinsToFill(const char* particle, const char* bins)
{
  delete fBinsToFill;
  fBinsToFill = Binning()->CreateBinObjArray(particle,bins,"");
}

//_____________________________________________________________________________
TVector2 AliAnalysisMuMuFlowSP::GetQn(const char* detector, Int_t step)
{
  // The function access the corrected Qn vector from the Qn correction framework (PWGPP/EVCHAR/FlowVectorCorrections)
  // Check the documentation at https://twiki.cern.ch/twiki/bin/view/ALICE/StartUsingR2FlowVectorCorrections
  //
  Double_t Qx=0.;
  Double_t Qy=0.;
  AliQnCorrectionsManager *flowQnMgr;
  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
      static_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != NULL) {
    flowQnMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliFatal("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!");
    return NULL;
  }

  //Get a TList containing Tlist of Qnvectors for each detector
  //Detector > correction step (cf QnCorr : raw, plain, rec, align + info quality) > harmonic (QX, QY, EP)
  TList* detectorlist = flowQnMgr->GetQnVectorList();
  if (!detectorlist) AliError("detectorlist is missing... You should check what happened");
  else{
    TList* qnlist = static_cast<TList*> (detectorlist->FindObject(detector));
    if (!qnlist) AliError("Detectorlist was found but there is no entry for your detector");
    else{
      AliQnCorrectionsQnVector* qn = static_cast<AliQnCorrectionsQnVector*> (qnlist->FindObject(fEqSteps[step].Data())); //last interesting step for us

      if (qn == NULL) {
      // align step was not found, trying to get something else
        AliError(Form("%s step was not found for detector %s",fEqSteps[step].Data(),detector));
      }
      else {
        Qx = static_cast<Double_t> (qn->Qx(fHar));
        Qy = static_cast<Double_t> (qn->Qy(fHar));
      }
      if(Qx == 0.) AliError(Form(" Qx=0 but qn vector is not null for detector %s at step %s",detector,fEqSteps[step].Data()));
    }
  }
  TVector2 Q(Qx,Qy);
  return Q;
}

//________________________________________________________________________
void AliAnalysisMuMuFlowSP::SetOriginPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlowSP::SetNewPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlowSP::SetOriginYFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlowSP::SetNewYFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlowSP::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
   f->SetNpx(100.*(max-min));
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}
