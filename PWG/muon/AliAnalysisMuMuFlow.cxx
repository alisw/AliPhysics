#include "AliAnalysisMuMuFlow.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuFlow
 *
 * Analysis which fills a bunch of histograms for elliptic flow analysis of J/psi
 * The flow vectors and event plane angles are recovered from Qn correction framework (PWGPP/EVCHAR/FlowVectorCorrections)
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
#include "THnSparse.h"
#include "AliMCEvent.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuonUtility.h"
#include "TParameter.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include <cassert>

ClassImp(AliAnalysisMuMuFlow)

//_____________________________________________________________________________
AliAnalysisMuMuFlow::AliAnalysisMuMuFlow(TH2* accEffHisto, TList *q2Map, Int_t systLevel)
: AliAnalysisMuMuBase(),
fcomputeMeanV2(kTRUE),
fcomputeEP(kTRUE),
fcomputeSP(kTRUE),
fESE(kTRUE),
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
fmcptcutmin(0.0),
fmcptcutmax(12.0),
fNDetectors(3),
fHar(2),
EP{0.,0.,0.},
Q2{{0,0},{0,0},{0,0}},
fq2Map{0x0,0x0}
{
  // FIXME ? find the AccxEff histogram from HistogramCollection()->Histo("/EXCHANGE/JpsiAccEff")
  if ( accEffHisto )
  {
    fAccEffHisto = static_cast<TH2F*>(accEffHisto->Clone());
    fAccEffHisto->SetDirectory(0);
  }
  if(q2Map){
    fESE = kTRUE;
    fq2Map[0] = static_cast<TH1F*>(q2Map->FindObject("smallq2")->Clone());
    fq2Map[1] = static_cast<TH1F*>(q2Map->FindObject("largeq2")->Clone());
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuFlow::~AliAnalysisMuMuFlow()
{
  /// dtor
  delete fAccEffHisto;
  delete fBinsToFill;
}

//_____________________________________________________________________________
void
AliAnalysisMuMuFlow::DefineHistogramCollection(const char* eventSelection,
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
  if (!fBinsToFill) SetBinsToFill("psi","integrated,pt,y,dphiSPD, dphiV0A, dphiV0C,dphivsptSPD,dphivsptV0A,dphivsptV0C");

  // mass range
  Double_t minvMin = fMinvMin;
  Double_t minvMax = fMinvMax;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,fMinvBinSize);
  // cent range
  Double_t centMin = 20.;
  Double_t centMax = 40.;
  Int_t nCentBins = GetNbins(centMin,centMax,1.);

  // Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);
  //Event shape : ThnSparses
  // const Int_t nDimThNS = 8;
  // const Int_t nDimThNS_SP = 6;
  // //                   minv   pt dphi cos2dp q2 cent EPp  EPev
  // Int_t nBins[nDimThNS]={nMinvBins,200,200,500,500,200,100,100};
  // Double_t xMin[nDimThNS]={minvMin,0.,0.,-1.,0.,20.,-1.,-1.};
  // Double_t xMax[nDimThNS]={minvMax,12.,3.2,1.,10.,40.,1.,1.};
  // // //                   minv   pt SP SP corr cent
  // Int_t nBinsSP[nDimThNS_SP]={nMinvBins,200,500,500,100,300};
  // Double_t xMinSP[nDimThNS_SP]={minvMin,0.,0.,0.,20.};
  // Double_t xMaxSP[nDimThNS_SP]={minvMax,12.,10.,10.,40.};

  // if(fESE){
  //   CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SPD","#mu+#mu+ v2 distribution",nDimThNS,nBins,xMin,xMax);
  //   // CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_V0A","#mu+#mu+ v2 distribution",nDimThNS,nBins,xMin,xMax);

  //   CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SP_SPD","#mu+#mu+ v2 distribution",nDimThNS_SP,nBinsSP,xMinSP,xMaxSP);
  //   // CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SP_V0A","#mu+#mu+ v2 distribution",nDimThNS,nBins,xMin,xMax);
  //   if(ShouldCorrectDimuonForAccEff()){
  //     CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SPD_AE","#mu+#mu+ v2 distribution (Acc #times Eff Corrected)",nDimThNS,nBins,xMin,xMax);
  //     // CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_V0A","#mu+#mu+ v2 distribution",nDimThNS,nBins,xMin,xMax);

  //     CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SP_SPD_AE","#mu+#mu+ v2 distribution (Acc #times Eff Corrected)",nDimThNS_SP,nBinsSP,xMinSP,xMaxSP);
  //     // CreatePairTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"ESE_SP_V0A","#mu+#mu+ v2 distribution",nDimThNS,nBins,xMin,xMax);
  //   }
  // }
  for(Int_t i=0; i<fNDetectors;i++){
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("EVENTPLANE_%s",fDetectors[i].Data()),Form("#mu+#mu- event plane distributionwith %s",fDetectors[i].Data()),
                     500, -1.6, 1.6,-2);
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Cos2EP_%s",fDetectors[i].Data()),Form("#mu+#mu- cos2EP distributionwith %s",fDetectors[i].Data()),
                     100, -1., 1.,-2);
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Sin2EP_%s",fDetectors[i].Data()),Form("#mu+#mu- Sin2EP distributionwith %s",fDetectors[i].Data()),
                     100, -1., 1.,-2);
    CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("DPHI_%s",fDetectors[i].Data()),Form("#mu+#mu- Dphi distribution with %s",fDetectors[i].Data()),
                     500, -0.01, 3.2,-2);//dphi corrected to be in [O,pi]
    CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("SP_%s",fDetectors[i].Data()),Form("#mu+#mu- Dphi distribution with %s",fDetectors[i].Data()),
                     500, -0.01, 3.2,-2);//dphi corrected to be in [O,pi]

    //Form("EPforpairs_%s",fDetectors[i].Data()
    CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("EPforpairs_%s",fDetectors[i].Data()),Form("EP distribution for each pair for%s",fDetectors[i].Data()),
                     500, -3.2, 3.2,-2);//dphi corrected to be in [O,pi]
    CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qnforpairs_%s",fDetectors[i].Data()),Form("Qn distribution for each pair for%s",fDetectors[i].Data()),
                     500, 0., 2., -2);//dphi corrected to be in [O,pi]
    CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qnforpairsvscent_%s",fDetectors[i].Data()),Form("Qn vector for pairs from %s; centrality (%%);q_{2}^{%s}",fDetectors[i].Data(),fDetectors[i].Data()),
                     50, 0., 90., 500, 0., 2.);
    ///
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qn_%s",fDetectors[i].Data()),Form("Qn vector from %s; q_{2}^{%s};N_{entries}",fDetectors[i].Data(),fDetectors[i].Data()),
                     500, 0., 2., -2);
    CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qnvscent_%s",fDetectors[i].Data()),Form("Qn vector from %s; centrality (%%);q_{2}^{%s}",fDetectors[i].Data(),fDetectors[i].Data()),
                     50, 0., 90., 500, 0., 2.);
    for(Int_t j=i+1; j<fNDetectors;j++){
      if(fcomputeEP){
        CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("DPHI%svsDPHI%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- Dphi distribution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                         500, -0.01, 3.2,500, -0.01, 3.2);//dphi corrected to be in [O,pi]
        CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("EP%svsEP%sforpairs",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- event plane distribution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                         500, -1.6, 1.6,500, -1.6, 1.6);
        CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("EP%svsEP%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- event plane distribution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                         500, -1.6, 1.6,500, -1.6, 1.6);
        CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("hEvPlaneReso%s_%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("Cos2(Psi_%s - Psi_%s)",fDetectors[i].Data(),fDetectors[j].Data()),
                          500, -1.2, 1.2, -2);
      }

      if(fcomputeSP) {
        CreateEventHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Qn%svsQn%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- resolution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                       500, 0., 2., 500,0.,2.);
        CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("SP%svsSP%s",fDetectors[i].Data(),fDetectors[j].Data()),Form("#mu+#mu- resolution : %s vs %s",fDetectors[i].Data(),fDetectors[j].Data()),
                       500, -1., 1.,500, -1., 1.);
        CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("Q%s*Q%s",fDetectors[i].Data(),fDetectors[j].Data()),
                  Form("SP Q%s*Q%s ;Centrality",fDetectors[i].Data(),fDetectors[j].Data()),nCentBins,centMin,centMax,0);

      }
    }
  }


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

      if(fcomputeSP){
        TString mUName(Form("U_%s",minvName.Data()));
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mUName.Data(),
                         "#mu+#mu- Q vector ;u_{x};u_{y}",600, -1., 1., 600, -1., 1.);
      }
      if ( fcomputeMeanV2 && !minvName.Contains("PHI")){
        TString mYName(Form("MeanYVs%s",minvName.Data()));
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mYName.Data(),
                  Form("#mu+#mu- mean y %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});<y^{#mu^{+}#mu^{-} (GeV/c^{2})}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
      }

      for(Int_t i=0; i<3;i++){
        if(fcomputeMeanV2){
        TString mV2Name = Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
        // Reconstructed pair histo
          CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mV2Name.Data(),
                         Form("#mu+#mu- mean v_{2}^{obs} %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2}^{obs} =< cos {2(#varphi_{#mu^{+}#mu^{-}}- #Psi_{EP,2})} > with %s",r->AsString().Data(),fDetectors[i].Data()),nMinvBins,minvMin,minvMax,0);
        }
        if(fcomputeSP){
        TString mSPName(Form("SPVs%s_EP_%s",minvName.Data(),fDetectors[i].Data()));
        // Reconstructed pair histo
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mSPName.Data(),
                         Form("#mu+#mu- v_{2} {SP} with %s vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPName.Data()),
                         Form("#mu+#mu- v_{2} {SP} with %s vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        // CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPName.Data()),
        //                  Form("#mu+#mu- v_{2} {SP} (corrected with |Qn|) with %s vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        }
      }
    }

    // Create corrected histo
    if ( ShouldCorrectDimuonForAccEff() )
    {
      minvName = GetMinvHistoName(*r,kTRUE);

      if ( !IsHistogramDisabled(minvName.Data()) ){

        AliDebug(1,Form("bin %d %s histoname = %s",nb,r->AsString().Data(),minvName.Data()));
        if ( fcomputeMeanV2 && !minvName.Contains("PHI") ){
          TString mYName(Form("MeanYVs%s",minvName.Data()));
          CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mYName.Data(),
                           Form("#mu+#mu- mean (AE corrected) y %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});<y^{#mu^{+}#mu^{-} (GeV/c^{2})}>",r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
        }
        if(fcomputeSP){
        TString mUAEName(Form("U_%s",minvName.Data()));
        CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mUAEName.Data(),
                         "#mu+#mu- Q vector (Acc #times Eff Corrected) ;u_{x};u_{y}",600, -1., 1., 600, -1., 1.);
        }
        for(Int_t i=0; i<fNDetectors;i++){
          if(fcomputeMeanV2){
            TString mV2Name = Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
          // Reconstructed pair histo
            CreatePairHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,mV2Name.Data(),
                           Form("#mu+#mu- mean (AE corrected) v_{2}^{obs} %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2}^{obs} =< cos {2(#varphi_{#mu^{+}#mu^{-}}- #Psi_{EP,2})} > with %s",r->AsString().Data(),fDetectors[i].Data()),nMinvBins,minvMin,minvMax,0);
          }
          if(fcomputeSP){
            TString mSPAEName(Form("SPVs%s_EP_%s",minvName.Data(),fDetectors[i].Data()));
            CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,mSPAEName.Data(),
              Form("#mu+#mu- v_{2} {SP} with %s (Acc #times Eff Corrected) vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
            CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPAEName.Data()),
              Form("#mu+#mu- v_{2} {SP} with %s (Acc #times Eff Corrected) vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);
            // CreatePairHistos(kHistoForData| kHistoForMCInput,eventSelection,triggerClassName,centrality,Form("%s_corr",mSPAEName.Data()),
            //   Form("#mu+#mu- v_{2} {SP} (corrected with |Qn|) with %s (Acc #times Eff Corrected) vs %s;M_{#mu^{+}#mu^{-}} (GeV/c^{2});v_{2} {SP}",fDetectors[i].Data(),r->AsString().Data()),nMinvBins,minvMin,minvMax,0);

          }
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlow::DefineMinvRange(Double_t minvMin, Double_t minvMax, Double_t minvBinSize)
{
  /// Define the Minv histogram range

  fMinvMin     = minvMin;
  fMinvMax     = minvMax;
  fMinvBinSize = minvBinSize;
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlow::FillHistosForPair(const char* eventSelection,
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

  // Weight tracks if specified
  Double_t inputWeight=0.;
  if(!fWeightMuon)      inputWeight = WeightPairDistribution(pair4Momentum.Pt(),pair4Momentum.Rapidity());
  else if(fWeightMuon)  inputWeight = WeightMuonDistribution(tracki.Pt()) * WeightMuonDistribution(trackj.Pt());

  //********************
  // Get QN and EP from QnCorrections Framework
  //********************
  TVector2 Qn[3];
  Double_t phiEP[3];//PHIEP from 2nd harmonic
  // Fill some distribution histos
  // Double_t phiEP[3];//PHIEP from 2nd harmonic
  Double_t dphi[3];//relative angle for each EP detector (V0A, SPD, V0C)
  // TVector2 Qn[3];//Q vectors (2nd harmonic) for each detector
  Double_t SP[3];//Scalar product

  for(Int_t i=0; i<3; i++){
    // if(i==0) {
      phiEP[i]= EP[i]; //twist for SPD
      Qn[i].Set(Q2[i][0],Q2[i][1]); //twist for SPD HERE
    // }
    // else {
    //   phiEP[i]= GetEventPlane(fDetectors[i].Data(),3); //twist
    //   Qn[i]= GetQn(fDetectors[i].Data(),3); //twist
    // }

    dphi[i] = phiEP[i]-pair4Momentum.Phi();
    if( dphi[i] <  0 ) dphi[i]+=2*TMath::Pi();
    if( dphi[i] >=TMath::Pi()) dphi[i] -= TMath::Pi();

    SP[i] = U.X()*Qn[i].X()+U.Y()*Qn[i].Y();
    if ( !IsHistogramDisabled(Form("EPforpairs_%s",fDetectors[i].Data())) ) proxy->Histo(Form("EPforpairs_%s",fDetectors[i].Data()))->Fill(phiEP[i]);
    if ( !IsHistogramDisabled(Form("Qnforpairs_%s",fDetectors[i].Data())) ) proxy->Histo(Form("Qnforpairs_%s",fDetectors[i].Data()))->Fill(sqrt(Q2[i][0]*Q2[i][0]+Q2[i][1]*Q2[i][1]));
    if ( !IsHistogramDisabled(Form("Qnforpairsvscent_%s",fDetectors[i].Data())) ) proxy->Histo(Form("Qnforpairsvscent_%s",fDetectors[i].Data()))->Fill(GetCentrality(),sqrt(Q2[i][0]*Q2[i][0]+Q2[i][1]*Q2[i][1]));
    if ( !IsHistogramDisabled(Form("DPHI_%s",fDetectors[i].Data())) ) proxy->Histo(Form("DPHI_%s",fDetectors[i].Data()))->Fill(dphi[i]);
    if (fcomputeSP && !IsHistogramDisabled(Form("SP_%s",fDetectors[i].Data())) ) proxy->Histo(Form("SP_%s",fDetectors[i].Data()))->Fill(SP[i]);
  }
  // Double_t x_ESE[8]={pair4Momentum.M(),pair4Momentum.Pt(),dphi[0],cos(2*dphi[0]),sqrt(Qn[0]*Qn[0]),GetCentrality(),phiEP[0],phiEP[1]}; //minv      pt     dphi     q2     cent    EPp  EPev
  // if( fESE && !IsHistogramDisabled("ESE_SPD"))static_cast<THnSparse*>(proxy->GetObject("ESE_SPD"))->Fill(x_ESE,inputWeight);
  // x_ESE[2] = dphi[1];
  // x_ESE[3] = cos(2*dphi[1]);
  // x_ESE[4] = sqrt(Qn[1]*Qn[1]);
  // x_ESE[6] = phiEP[1];
  // x_ESE[7] = phiEP[0];
  // if( !IsHistogramDisabled("ESE_V0A"))static_cast<THnSparse*>(proxy->GetObject("ESE_V0A"))->Fill(x_ESE,inputWeight);

  // Double_t x_ESE_SP[6]={pair4Momentum.M(),pair4Momentum.Pt(),SP[0]/sqrt(Qn[1]*Qn[2]),(SP[0]-cos(Qn[0].X())*cos(Qn[0].X())-sin(Qn[0].X())*sin(Qn[0].X()))/sqrt(Qn[1]*Qn[2]-Qn[1].X()*Qn[2].X()-Qn[1].Y()*Qn[2].Y())/sqrt(Qn[0]*Qn[0]),GetCentrality(),sqrt(Qn[0].X()*Qn[0].X()+Qn[0].X()*Qn[0].X())}; //minv      pt     SP     SPcorr     cent   q2

  // if( fESE && !IsHistogramDisabled("ESE_SP_SPD"))static_cast<THnSparse*>(proxy->GetObject("ESE_SP_SPD"))->Fill(x_ESE_SP,inputWeight);


  //AE
  // if(ShouldCorrectDimuonForAccEff()){
  //   Double_t AccxEff = GetAccxEff(pair4Momentum.Pt(),pair4Momentum.Rapidity());
  //   if ( AccxEff <= 0.0 ) AliError(Form("AccxEff < 0 for pt = %f & y = %f ",pair4Momentum.Pt(),pair4Momentum.Rapidity()));
  //   else {
  //     if( !IsHistogramDisabled("ESE_SPD_AE"))static_cast<THnSparse*>(proxy->GetObject("ESE_SPD_AE"))->Fill(x_ESE,inputWeight/AccxEff);
  //     if( !IsHistogramDisabled("ESE_SP_SPD_AE"))static_cast<THnSparse*>(proxy->GetObject("ESE_SP_SPD_AE"))->Fill(x_ESE_SP,inputWeight/AccxEff);
  //   }
  // }
  Double_t detSP[3]={0.,0.,0.};
  for(Int_t i=0; i<3; i++){
    for(Int_t j=i+1; j<fNDetectors;j++){
      detSP[i+j-1]= Qn[i].X()*Qn[j].X()+Qn[i].Y()*Qn[j].Y();
      if ( !IsHistogramDisabled(Form("DPHI%svsDPHI%s",fDetectors[i].Data(),fDetectors[j].Data())) ) proxy->Histo(Form("DPHI%svsDPHI%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(dphi[i],dphi[j]);
      if ( !IsHistogramDisabled(Form("SP%svsSP%s",fDetectors[i].Data(),fDetectors[j].Data())) ) proxy->Histo(Form("SP%svsSP%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(SP[i],SP[j]);
      if ( !IsHistogramDisabled(Form("EP%svsEP%sforpairs",fDetectors[i].Data(),fDetectors[j].Data()) )) proxy->Histo(Form("EP%svsEP%sforpairs",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(phiEP[i],phiEP[j]);

    }
  }
  // Fill histos with MC stack info
  if ( HasMC() ){
    AliWarning("MC is not implemented for flow analysis");

    // Get 4-vector pairs from MC stack
    // TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTracki->P()*mcTracki->P()));
    // TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+mcTrackj->P()*mcTrackj->P()));
    // mcpj+=mcpi;

    // TVector2 UMC(0.,0.);

    // set pair4MomentumMC for the rest of the function
    // pair4MomentumMC = &mcpj;
    // UMC.Set(cos(fHar*pair4MomentumMC->Phi()),sin(fHar*pair4MomentumMC->Phi()));//Unitary Q vector of the dimuon
    // Double_t SPMC[3];//Scalar product

    // if ( !IsHistogramDisabled("UMC") )mcProxy->Histo("UMC")->Fill(UMC.X(),UMC.Y());

    // for(Int_t i=0; i<fNDetectors; i++){
    //   SPMC[i] = UMC*Qn[i];
    // }
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
      if (r->AsString().BeginsWith("DPHIVSPTSPD") )
      {
        ok = (r->IsInRange(pair4Momentum.Pt(),dphi[0]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),dphi[0]);
      }
      else if (r->AsString().BeginsWith("DPHIVSPTV0A") )
      {
        ok = (r->IsInRange(pair4Momentum.Pt(),dphi[1]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),dphi[1]);
      }
      else if (r->AsString().BeginsWith("DPHIVSPTV0C") )
      {
        ok = (r->IsInRange(pair4Momentum.Pt(),dphi[2]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt(),dphi[2]);
      }
      else if (r->AsString().BeginsWith("DPHIVSY") )
      {
        ok = (r->IsInRange(pair4Momentum.Rapidity(),dphi[0]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity(),dphi[0]);
      }
      else AliError(Form("Don't know how to deal with 2D bin %s",r->AsString().Data()));
    }
    else{
      if ( r->Quantity() == "PT" ){
        ok                          = r->IsInRange(pair4Momentum.Pt());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Pt());
      }
      else if ( r->Quantity() == "Y" ){
        ok                          = r->IsInRange(pair4Momentum.Rapidity());
        if ( pair4MomentumMC ) okMC = r->IsInRange(pair4MomentumMC->Rapidity());
      }
      else if ( r->Quantity() == "DPHISPD" )
      {
        ok = (r->IsInRange(dphi[0]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(dphi[0]);
      }
       else if ( r->Quantity() == "DPHIV0A" )
      {
        ok = (r->IsInRange(dphi[1]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(dphi[1]);
      }
       else if ( r->Quantity() == "DPHIV0C" )
      {
        ok = (r->IsInRange(dphi[2]));
        if ( pair4MomentumMC ) okMC = r->IsInRange(dphi[2]);
      }
    }

    // Check if pair pass all conditions, either MC or not, and fill Minv Histogrames
    if ( ok || okMC ){

      // Get Minv histo name associated to the bin
      TString minvName = GetMinvHistoName(*r,kFALSE);

      //Create, fill and store Minv histo
      if (!IsHistogramDisabled(minvName.Data())){

        if ( (r->Quantity() == "PT"||r->IsIntegrated())){
          if(fcomputeMeanV2){
            TString hprofPtName("");
            TString hprofYName("");
            TString hprofmV2Name("");

            if ( ok ){
              // rapidity
              hprofYName= Form("MeanYVs%s",minvName.Data());
              TProfile* hprofY = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofYName.Data());
              if ( !hprofY )AliError(Form("Could not get %s",hprofYName.Data()));
              else hprofY->Fill(pair4Momentum.M(),pair4Momentum.Rapidity(),inputWeight);
              // Costheta
              for(Int_t i=0; i<3;i++){
                hprofmV2Name= Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
                TProfile* hprofmV2 = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofmV2Name.Data());
                if ( !hprofmV2)AliError(Form("Could not get %s",hprofmV2Name.Data()));
                else hprofmV2->Fill(pair4Momentum.M(),cos(2*dphi[i]),inputWeight);
              }
            }

            if ( okMC ){
              // rapidity
              hprofYName= Form("MeanYVs%s",minvName.Data());
              TProfile* hprofY = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofYName.Data());
              if ( !hprofY )AliError(Form("Could not get %s",hprofYName.Data()));
              else hprofY->Fill(pair4MomentumMC->M(),pair4MomentumMC->Rapidity(),inputWeightMC);
              // Costheta
              for(Int_t i=0; i<3;i++){
                hprofmV2Name= Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
                TProfile* hprofmV2 = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofmV2Name.Data());
                if ( !hprofmV2)AliError(Form("Could not get %s",hprofmV2Name.Data()));
                else hprofmV2->Fill(pair4MomentumMC->M(),cos(2*dphi[i]),inputWeightMC);
              }
            }
          }
          if(fcomputeSP){
            if ( ok ){
              TString hprofName("");
              for(Int_t i=0; i<fNDetectors;i++){
                hprofName = Form("SPVs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
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

            // if( okMC ){
            //   TString hprofMCName("");
            //   for(Int_t i=0; i<fNDetectors;i++){
            //     hprofMCName= Form("SPMCvs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
            //     TProfile* hprofMC = Prof(eventSelection,triggerClassName,centrality,hprofMCName.Data());
            //     if ( !hprofMC)AliError(Form("Could not get %s",hprofMCName.Data()));
            //     else hprofMC->Fill(pair4Momentum.M(),SPMC[i],inputWeight);
            //     TProfile* hprofMCcorr = Prof(eventSelection,triggerClassName,centrality,Form("%s_corr",hprofMCName.Data()));
            //     if ( !hprofMCcorr)AliError(Form("Could not get %s",Form("%s_corr",hprofMCName.Data())));
            //     else hprofMCcorr->Fill(pair4MomentumMC->M(),sqrt(SPMC[i]-Qn[i].X()*Qn[i].X()+Qn[i].Y()*Qn[i].Y()),inputWeight);//Correct
            //   }
              // TString mUMCName(Form("UMC_%s",minvName.Data()));
              // if ( !IsHistogramDisabled(mUMCName.Data() )) proxy->Histo(mUMCName.Data())->Fill(UMC.X(),UMC.Y());
            // }
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

          if (r->Quantity() == "PT"||r->IsIntegrated()){
            if(fcomputeMeanV2){
              TString hprofCorrName("");
              TString hprofYName("");
              TString hprofmV2Name("");

              // rapidity
              hprofYName= Form("MeanYVs%s",minvName.Data());
              TProfile* hprofY = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofYName.Data());
              if ( !hprofY )AliError(Form("Could not get %s",hprofYName.Data()));
              else if ( okAccEff )hprofY->Fill(pair4Momentum.M(),pair4Momentum.Rapidity(),inputWeight/AccxEff);
              // Costheta
              for(Int_t i=0; i<fNDetectors;i++){
              hprofmV2Name= Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
              TProfile* hprofmV2 = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofmV2Name.Data());
              if ( !hprofmV2)AliError(Form("Could not get %s",hprofmV2Name.Data()));
              else if ( okAccEff )hprofmV2->Fill(pair4Momentum.M(),cos(2*dphi[i]),inputWeight/AccxEff);
              }

              if( okMC ){
                // rapidity
                hprofYName= Form("MeanYVs%s",minvName.Data());
                TProfile* hprofY = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofYName.Data());
                if ( !hprofY )AliError(Form("Could not get %s",hprofYName.Data()));
                else if ( okAccEffMC )hprofY->Fill(pair4MomentumMC->M(),pair4MomentumMC->Rapidity(),inputWeight/AccxEff);
                // Costheta
                for(Int_t i=0; i<3;i++){
                  hprofmV2Name= Form("MeanV2Vs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
                  TProfile* hprofmV2 = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofmV2Name.Data());
                  if ( !hprofmV2)AliError(Form("Could not get %s",hprofmV2Name.Data()));
                  else if ( okAccEffMC )hprofmV2->Fill(pair4MomentumMC->M(),cos(2*dphi[i]),inputWeight/AccxEff);
                }
              }
            }
            if(fcomputeSP){
              if ( ok ){
                TString hprofName("");
                for(Int_t i=0; i<fNDetectors;i++){
                  hprofName = Form("SPVs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
                  TProfile* hprof = Prof(eventSelection,triggerClassName,centrality,pairCutName,hprofName.Data());
                  if ( !hprof)AliError(Form("Could not get %s",hprofName.Data()));
                  else hprof->Fill(pair4Momentum.M(),SP[i],inputWeight/AccxEff);
                  TProfile* hprofcorr = Prof(eventSelection,triggerClassName,centrality,pairCutName,Form("%s_corr",hprofName.Data()));
                  if ( !hprofcorr)AliError(Form("Could not get %s",Form("%s_corr",hprofName.Data())));
                  else hprofcorr->Fill(pair4Momentum.M(),(SP[i]-Qn[i].X()*Qn[i].X()+Qn[i].Y()*Qn[i].Y()),inputWeight/AccxEff);
                }
              }

              if( okMC ){
              TString hprofMCName("");
              for(Int_t i=0; i<fNDetectors;i++){
                hprofMCName= Form("SPMCvs%s_EP_%s",minvName.Data(),fDetectors[i].Data());
                TProfile* hprofMC = Prof(eventSelection,triggerClassName,centrality,hprofMCName.Data());
                if ( !hprofMC)AliError(Form("Could not get %s",hprofMCName.Data()));
                else hprofMC->Fill(pair4Momentum.M(),SP[i],inputWeight/AccxEff);
                TProfile* hprofMCcorr = Prof(eventSelection,triggerClassName,centrality,Form("%s_corr",hprofMCName.Data()));
                if ( !hprofMCcorr)AliError(Form("Could not get %s",Form("%s_corr",hprofMCName.Data())));
                else hprofMCcorr->Fill(pair4MomentumMC->M(),sqrt(SP[i]-Qn[i].X()*Qn[i].X()+Qn[i].Y()*Qn[i].Y()),inputWeight/AccxEff);//Correct
              }
              }
            }
          }
        }
      }
    }
  }
  delete proxy;
  delete mcProxy;
}

//________________________________________________________________________
void AliAnalysisMuMuFlow::FillHistosForEvent(const char* eventSelection,
                                            const char* triggerClassName,
                                            const char* centrality)
{
  // Fill histos with event planes and Qn vectors + compute the resolution with the 3 sub-event method
  // The function access the corrected Qn vector from the Qn correction framework (PWGPP/EVCHAR/FlowVectorCorrections)
  // Check the documentation at https://twiki.cern.ch/twiki/bin/view/ALICE/StartUsingR2FlowVectorCorrections
  AliQnCorrectionsManager *flowQnMgr;
  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
      static_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != NULL) {
    flowQnMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliFatal("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!");
    return;
  }
  Double_t phiEP[fNDetectors];
  TVector2 Qn[3];
  Int_t step =3;
  //Get a TList containing Tlist of Qnvectors for each detector
  //Detector > correction step (cf QnCorr : raw, plain, rec, align + info quality) > harmonic (QX, QY, EP)
  TList* detectorlist = flowQnMgr->GetQnVectorList();
  if (!detectorlist) AliError("detectorlist is missing... You should check what happened");
  //here we get and fill
  for(Int_t i=0; i<fNDetectors; i++){
    if(i==0)step=4;
    else step=3;
    TList* qnlist = static_cast<TList*> (detectorlist->FindObject(fDetectors[i].Data()));
    if (!qnlist) AliError("Detectorlist was found but there is no entry for your detector");
    else{
      AliQnCorrectionsQnVector* qn = static_cast<AliQnCorrectionsQnVector*> (qnlist->FindObject(fEqSteps[step].Data())); //last step

      if (qn == NULL) {
      // align step was not found, trying to get something else
        AliError(Form("%s step was not found for detector %s",fEqSteps[step].Data(),fDetectors[i].Data()));
      }
      else {
        //EP
        phiEP[i] = static_cast<Double_t> (qn->EventPlane(2)); //2nd harmonic
        if(phiEP[i] == 0.) AliError(Form("EP=0 but qn vector is not null for detector %s at step %s",fDetectors[i].Data(),fEqSteps[step].Data()));
        else EP[i]=phiEP[i];
        //Qn
        Qn[i].Set(static_cast<Double_t> (qn->Qx(fHar)),static_cast<Double_t> (qn->Qy(fHar)));
        if(Qn[i].X() == 0.|| Qn[i].Y() == 0.) AliError(Form(" Qx=0 but qn vector is not null for detector %s at step %s",fDetectors[i].Data(),fEqSteps[step].Data()));
        else {
          Q2[i][0]=Qn[i].X();
          Q2[i][1]=Qn[i].Y();
        }
      }
    }
  }
  //Filling the histos
  for(Int_t i=0;i<3;i++){
    if( !IsHistogramDisabled(Form("EVENTPLANE_%s",fDetectors[i].Data())) ) Histo(eventSelection,triggerClassName,centrality,Form("EVENTPLANE_%s",fDetectors[i].Data()))->Fill(phiEP[i]);
    if( !IsHistogramDisabled(Form("Cos2EP_%s",fDetectors[i].Data())) ) Histo(eventSelection,triggerClassName,centrality,Form("Cos2EP_%s",fDetectors[i].Data()))->Fill(TMath::Cos(2*phiEP[i]));
    if( !IsHistogramDisabled(Form("Sin2EP_%s",fDetectors[i].Data())) ) Histo(eventSelection,triggerClassName,centrality,Form("Sin2EP_%s",fDetectors[i].Data()))->Fill(TMath::Sin(2*phiEP[i]));
    if( !IsHistogramDisabled(Form("Qn_%s",fDetectors[i].Data()))) Histo(eventSelection,triggerClassName,centrality,Form("Qn_%s",fDetectors[i].Data()))->Fill(sqrt(Qn[i]*Qn[i]));
    if( !IsHistogramDisabled(Form("Qnvscent_%s",fDetectors[i].Data()))) Histo(eventSelection,triggerClassName,centrality,Form("Qnvscent_%s",fDetectors[i].Data()))->Fill(GetCentrality(),sqrt(Qn[i]*Qn[i]));
    for(Int_t j=i+1; j<3;j++){
      //EP
      Double_t deltaEP =phiEP[i]-phiEP[j];
      if(TMath::Abs(deltaEP)>TMath::Pi()/fHar){
        if(deltaEP>0.) deltaEP-=2.*TMath::Pi()/fHar;
        else deltaEP+=2.*TMath::Pi()/fHar;
      }
      if(!IsHistogramDisabled(Form("hEvPlaneReso%s_%s",fDetectors[i].Data(),fDetectors[j].Data())))
        Histo(eventSelection,triggerClassName,centrality,Form("hEvPlaneReso%s_%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(TMath::Cos(fHar*deltaEP));

      //Fill Qn vector histos
      if ( !IsHistogramDisabled(Form("EP%svsEP%s",fDetectors[i].Data(),fDetectors[j].Data()) )) Histo(eventSelection,triggerClassName,centrality,Form("EP%svsEP%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(phiEP[i],phiEP[j]);
      if ( !IsHistogramDisabled(Form("Qn%svsQn%s",fDetectors[i].Data(),fDetectors[j].Data())) ) Histo(eventSelection,triggerClassName,centrality,Form("Qn%svsQn%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(sqrt(Qn[i]*Qn[i]),sqrt(Qn[j]*Qn[j]));
      if ( fcomputeSP && !IsHistogramDisabled(Form("Q%s*Q%s",fDetectors[i].Data(),fDetectors[j].Data())) )Prof(eventSelection,triggerClassName,centrality,Form("Q%s*Q%s",fDetectors[i].Data(),fDetectors[j].Data()))->Fill(GetCentrality(),Qn[i].X()*Qn[j].X()+Qn[i].Y()*Qn[j].Y());
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisMuMuFlow::FillHistosForMCEvent(const char* eventSelection,const char* triggerClassName,const char* centrality)
{
  ///
  /// Fill MC inputs histograms.
  ///

  if ( !HasMC() ) return;
  AliError("AliMuMuFlow does not deal with MC event, please implement it !");
}
//_____________________________________________________________________________
TString AliAnalysisMuMuFlow::GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r, Bool_t accEffCorrected) const
{
  return TString::Format("MinvUS%s%s%s",
                         accEffCorrected ? "_AccEffCorr" : "",fMinvBinSeparator.Data(),r.AsString().Data());
}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuFlow::GetAccxEff(Double_t pt,Double_t rapidity)
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
Double_t AliAnalysisMuMuFlow::WeightMuonDistribution(Double_t pt)
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
Double_t AliAnalysisMuMuFlow::WeightPairDistribution(Double_t pt,Double_t rapidity)
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
Double_t AliAnalysisMuMuFlow::TriggerLptApt ( Double_t* xVal, Double_t* par )
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
Bool_t AliAnalysisMuMuFlow::IsDPhiInPlane(const AliVParticle& t1, const AliVParticle& t2) const
{
  /// Whether the pair passes the dphi cut

  TLorentzVector pi(t1.Px(),t1.Py(),t1.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector pair4Momentum(t2.Px(),t2.Py(),t2.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  pair4Momentum += pi;

  Double_t dphi = EP[0] - pair4Momentum.Phi();
  if( dphi <  0 ) dphi+=2*TMath::Pi();
  if( dphi >=TMath::Pi()) dphi-= TMath::Pi();

  return  (( dphi < 3.142 && dphi > 2.356  )||( dphi < 0.785 && dphi > 0 ));
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuFlow::IsDPhiOutOfPlane(const AliVParticle& t1, const AliVParticle& t2) const
{
  /// Whether the pair passes the dphi cut

  TLorentzVector pi(t1.Px(),t1.Py(),t1.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t1.P()*t1.P()));
  TLorentzVector pair4Momentum(t2.Px(),t2.Py(),t2.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+t2.P()*t2.P()));
  pair4Momentum += pi;

  Double_t dphi = EP[0] - pair4Momentum.Phi();
  if( dphi <  0 ) dphi+=2*TMath::Pi();
  if( dphi >=TMath::Pi()) dphi-= TMath::Pi();

  return  ( dphi < 2.356 && dphi > 0.785 ); // dphi in [-pi/4,pi/4]
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlow::NameOfIsDPhiInPlane(TString& name) const
{
  name.Form("INPLANE");
}
//_____________________________________________________________________________
void AliAnalysisMuMuFlow::NameOfIsDPhiOutOfPlane(TString& name) const
{
  name.Form("OUTOFPLANE");
}
//________________________________________________________________________
Bool_t AliAnalysisMuMuFlow::Isq2InSmallRange(const AliVEvent& event) const
{
  if(!fq2Map[0]) {
    AliWarning("ERROR : no q2SmallMap provided");
    return kFALSE;
  }
  Double_t q2 = sqrt(Q2[0][0]*Q2[0][0]+Q2[0][1]*Q2[0][1]);
  Int_t bin    = fq2Map[0]->FindBin(GetCentrality());
  return q2 < fq2Map[0]->GetBinContent(bin);
}
//________________________________________________________________________
Bool_t AliAnalysisMuMuFlow::Isq2InLargeRange(const AliVEvent& event) const
{
  if(!fq2Map[1]) {
    AliWarning("ERROR : no q2SmallMap provided");
    return kFALSE;
  }
  Double_t q2 = sqrt(Q2[0][0]*Q2[0][0]+Q2[0][1]*Q2[0][1]);
  Int_t bin    = fq2Map[1]->FindBin(GetCentrality());
  return q2 > fq2Map[1]->GetBinContent(bin);
}

//_____________________________________________________________________________
void AliAnalysisMuMuFlow::NameOfIsq2InSmallRange(TString& name) const
{
  name.Form("smallq2");
}
//_____________________________________________________________________________
void AliAnalysisMuMuFlow::NameOfIsq2InLargeRange(TString& name) const
{
  name.Form("Largeq2");
}
//_____________________________________________________________________________
void AliAnalysisMuMuFlow::SetBinsToFill(const char* particle, const char* bins)
{
  delete fBinsToFill;
  fBinsToFill = Binning()->CreateBinObjArray(particle,bins,"");
}

//________________________________________________________________________
Double_t AliAnalysisMuMuFlow::GetCentrality() const{
  AliMultSelection* multSelection = static_cast<AliMultSelection*>(Event()->FindListObject("MultSelection"));
  if ( multSelection ) return multSelection->GetMultiplicityPercentile("V0M");
  else return 0.;
}
//________________________________________________________________________
void AliAnalysisMuMuFlow::SetOriginPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlow::SetNewPtFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlow::SetOriginYFunc(TString formula, const Double_t *param,Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlow::SetNewYFunc(TString formula, const Double_t *param, Double_t xMin, Double_t xMax)
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
void AliAnalysisMuMuFlow::NormFunc(TF1 *f, Double_t min, Double_t max)
{
  /// normalize the function to its integral in the given range
   f->SetNpx(100.*(max-min));
  Double_t integral = f->Integral(min, max);
  if (integral != 0.) f->SetParameter(0, f->GetParameter(0)/integral);
}
