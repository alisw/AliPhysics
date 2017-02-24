#include "AliAnalysisMuMuTriggerResponse.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuTriggerResponse
 *
 *   Class to study the trigger response per group of local board.
 *   This class can run either on MC and Data. One should run it first on Data to get and histogram collection with the lpt/apt distributions from data/MC from single muons.
 *   Then, one should use the macro ComputeLptAptRatio.C on the output to add histograms to the Mergeable collection
 *   Finally, add this new mergeable collection as data member and re-run it on MC to get the weighted trigger response.
 *   For details : https://indico.cern.ch/event/247499/contributions/537130/attachments/425512/590690/trigSystMethod.pdf
 *
 *   WARNING : TO BE USED WITH ALIANALYSISMUMUSINGLE
 *
 */

#include "AliLog.h"
#include "TString.h"
#include "THnSparse.h"
#include "AliMCEvent.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TCanvas.h"
#include "AliMergeableCollection.h"
#include "AliVParticle.h"
#include "AliAnalysisMuonUtility.h"
#include "TParameter.h"
#include "TGraphAsymmErrors.h"

ClassImp(AliAnalysisMuMuTriggerResponse)

//_____________________________________________________________________________
AliAnalysisMuMuTriggerResponse::AliAnalysisMuMuTriggerResponse(AliMergeableCollection* OC)
: AliAnalysisMuMuBase(),
fOC(0x0),
fLptAptEtaRanges(0x0),
fNEtaBin(0x0)
{

  if(fOC) delete fOC;
  if(OC)
  {
    fOC = static_cast<AliMergeableCollection*>(OC->Clone());
    fOC->SetOwner(kTRUE);
  }

}

//_____________________________________________________________________________
AliAnalysisMuMuTriggerResponse::~AliAnalysisMuMuTriggerResponse()
{
  /// dtor
  delete fOC;
  delete[] fLptAptEtaRanges;

}

//_____________________________________________________________________________
void
AliAnalysisMuMuTriggerResponse::DefineHistogramCollection(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality,
                                               Bool_t mix)
{
  /// Actually create the histograms for phyics/triggerClassName


  // Check if histo is not already here
  if ( Histo(eventSelection,triggerClassName,centrality,"AliAnalysisMuMuTriggerResponse") )
  {
    return;
  }

  AliAnalysisMuMuBase::EDataType dt = AliAnalysisMuMuBase::kHistoForData;

  // dummy histogram to signal that we already defined all our histograms (see above)
  CreateEventHistos(dt,eventSelection,triggerClassName,centrality,"AliAnalysisMuMuTriggerResponse","Dummy semaphore",1,0,1);

  // pt range
  Double_t ptMin  = 0;
  Double_t ptMax  = 12.;
  Int_t nbinsPt   = GetNbins(ptMin,ptMax,0.05);

  // eta range
  Double_t etaMin = -4;
  Double_t etaMax = -2.5;
  Int_t nbinsEta  = GetNbins(etaMin,etaMax,0.05);

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_lb")  )
  {
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMu_lb","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMin_lb","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMax_lb","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
  }

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_lbg")  )
  {
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMu_lbg","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMin_lbg","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMax_lbg","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
  }

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_y")  )
  {
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMu_y","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMin_y","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
    CreatePairHistos(kHistoForData,eventSelection,triggerClassName,centrality,"TriggerHitperLocalBoardMuMuErrorMax_y","#mu+#mu- P_{T}-eta distribution vs Trigger Local Board",nbinsPt,ptMin,ptMax,nbinsEta,etaMin,etaMax);
  }

  Int_t nBins[3]={233,nbinsPt,nbinsEta};
  Double_t xminsparse[3]={0.5,ptMin,etaMin};
  Double_t xmaxsparse[3]={234.5,ptMax,etaMax};

  if ( !IsHistogramDisabled("HitperTriggerLocalBoardMu")  )
  CreateTrackTHnSparse(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,"HitperTriggerLocalBoardMu","#mu P_{T}-eta distribution vs Trigger Local Board",3,nBins,xminsparse,xmaxsparse);


}


//_____________________________________________________________________________
void AliAnalysisMuMuTriggerResponse::FillHistosForTrack(const char* eventSelection,
                                               const char* triggerClassName,
                                               const char* centrality,
                                               const char* trackCutName,
                                               const AliVParticle& track)
{
  /// Fill histograms for one track

  if (!AliAnalysisMuonUtility::IsMuonTrack(&track) ) return;

  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,trackCutName));

  if ( HasMC() ) {
    // Select muons
    if ( track.GetLabel() < 0 ) return;
    AliVParticle* mcTrack = MCEvent()->GetTrack(track.GetLabel());
    if ( TMath::Abs(mcTrack->PdgCode()) != 13 ) return;
  }

  TLorentzVector p(track.Px(),track.Py(),track.Pz(),
                   TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+track.P()*track.P()));

  if (!IsHistogramDisabled("HitperTriggerLocalBoardMu"))
  {
    if(AliAnalysisMuonUtility::GetLoCircuit(&track) !=0){
      Double_t x[3] = {static_cast<Double_t>(AliAnalysisMuonUtility::GetLoCircuit(&track)),p.Pt(),p.Eta()};
      if(proxy->GetObject(Form("HitperTriggerLocalBoardMu")))
        static_cast<THnSparse*>(proxy->GetObject(Form("HitperTriggerLocalBoardMu")))->Fill(x);

    }
  }

  delete proxy;
}



//_____________________________________________________________________________
void AliAnalysisMuMuTriggerResponse::FillHistosForPair(const char* eventSelection,
                                            const char* triggerClassName,
                                            const char* centrality,
                                            const char* pairCutName,
                                            const AliVParticle& tracki,
                                            const AliVParticle& trackj,
                                            const Bool_t IsMixedHisto)
{

  /// Fill histograms for unlike-sign reconstructed  muon pairs.
  /// For the MC case, we check that only tracks with an associated MC label are selected (usefull when running on embedding).
  /// A weight is also applied for MC case at the the muon track level according to WeightFromLocalBoard() and systLevel.

  // Usual cuts
  if (!AliAnalysisMuonUtility::IsMuonTrack(&tracki) || !AliAnalysisMuonUtility::IsMuonTrack(&trackj) ) return;

  // Get total charge in order to get the correct histo name
  Double_t PairCharge = tracki.Charge() + trackj.Charge();
  if(PairCharge !=0) return;
  if(IsMixedHisto)   return;
  if(!fOC)           return;
  if(!HasMC())       return;

  // Pointers in case running on MC
  Int_t labeli               = 0;
  Int_t labelj               = 0;
  AliVParticle               * mcTracki(0x0);
  AliVParticle               * mcTrackj(0x0);

  // Create proxy in AliMergeableCollection
  AliMergeableCollectionProxy* proxy = HistogramCollection()->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,pairCutName));

  // Construct dimuons vector
  TLorentzVector pi(tracki.Px(),tracki.Py(),tracki.Pz(),
                    TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+tracki.P()*tracki.P()));
  TLorentzVector pair4Momentum(trackj.Px(),trackj.Py(),trackj.Pz(),
                               TMath::Sqrt(AliAnalysisMuonUtility::MuonMass2()+trackj.P()*trackj.P()));
  pair4Momentum += pi;

  // Get label
  labeli = tracki.GetLabel();
  labelj = trackj.GetLabel();
  if ( labeli < 0 || labelj < 0 ){
    AliError("Got negative labels!");
    return;
  }

  // Check if first track is a muon
  mcTracki = MCEvent()->GetTrack(labeli);
  if ( TMath::Abs(mcTracki->PdgCode()) != 13 ) {
    delete proxy;
    return;
  }

  // Check if second track is a muon
  mcTrackj = MCEvent()->GetTrack(labelj);
  if ( TMath::Abs(mcTrackj->PdgCode()) != 13 ) {
    delete proxy;
    return;
  }

  // Check if tracks has the same mother
  Int_t currMotheri = mcTracki->GetMother();
  Int_t currMotherj = mcTrackj->GetMother();
  if( currMotheri!=currMotherj ) {
    delete proxy;
    return;
  }
  if( currMotheri<0 ) {
    delete proxy;
    return;
  }

  // Check if mother is J/psi
  AliMCParticle* mother = static_cast<AliMCParticle*>(MCEvent()->GetTrack(currMotheri));
  if(!mother){
    delete proxy;
    return;
  }
  if(mother->PdgCode() !=443) {
    delete proxy;
    return;
  }

  // Get the local board of each muons
  Double_t locbi = static_cast<Double_t>(AliAnalysisMuonUtility::GetLoCircuit(&tracki));
  Double_t locbj = static_cast<Double_t>(AliAnalysisMuonUtility::GetLoCircuit(&trackj));

  Double_t inputWeight=1.;
  Double_t x[2] = {pair4Momentum.Pt(),pair4Momentum.Rapidity()};

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_lbg") )
  {
    AliMergeableCollectionProxy* proxyoclptaptHisto = fOC->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,"LocalBoardGroup"));

    // Weight tracks if specified
    inputWeight = WeightFromLocalBoardGroup(tracki,locbi,*proxyoclptaptHisto,0) * WeightFromLocalBoardGroup(trackj,locbj,*proxyoclptaptHisto,0);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMu_lbg") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMu_lbg"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromLocalBoardGroup(tracki,locbi,*proxyoclptaptHisto,1) * WeightFromLocalBoardGroup(trackj,locbj,*proxyoclptaptHisto,1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_lbg") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_lbg"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromLocalBoardGroup(tracki,locbi,*proxyoclptaptHisto,-1) * WeightFromLocalBoardGroup(trackj,locbj,*proxyoclptaptHisto,-1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_lbg") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_lbg"))->Fill(x[0],x[1],inputWeight);

    delete proxyoclptaptHisto;

  }

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_lb") )
  {
    AliMergeableCollectionProxy* proxyoclptaptHisto = fOC->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,"LocalBoard"));

    // Weight tracks if specified
    inputWeight = WeightFromLocalBoard(tracki,locbi,*proxyoclptaptHisto,0) * WeightFromLocalBoard(trackj,locbj,*proxyoclptaptHisto,0);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMu_lb") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMu_lb"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromLocalBoard(tracki,locbi,*proxyoclptaptHisto,1) * WeightFromLocalBoard(trackj,locbj,*proxyoclptaptHisto,1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_lb") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_lb"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromLocalBoard(tracki,locbi,*proxyoclptaptHisto,-1) * WeightFromLocalBoard(trackj,locbj,*proxyoclptaptHisto,-1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_lb") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_lb"))->Fill(x[0],x[1],inputWeight);

    delete proxyoclptaptHisto;

  }

  if ( !IsHistogramDisabled("TriggerHitperLocalBoardMuMu_y") )
  {
    AliMergeableCollectionProxy* proxyoclptaptHisto = fOC->CreateProxy(BuildPath(eventSelection,triggerClassName,centrality,"Range"));

    // Weight tracks if specified
    inputWeight = WeightFromEta(tracki,*proxyoclptaptHisto,0) * WeightFromEta(trackj,*proxyoclptaptHisto,0);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMu_y") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMu_y"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromEta(tracki,*proxyoclptaptHisto,1) * WeightFromEta(trackj,*proxyoclptaptHisto,1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_y") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMax_y"))->Fill(x[0],x[1],inputWeight);

    inputWeight = WeightFromEta(tracki,*proxyoclptaptHisto,-1) * WeightFromEta(trackj,*proxyoclptaptHisto,-1);
    if (proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_y") && inputWeight > 0)
      static_cast<TH2*>(proxy->GetObject("TriggerHitperLocalBoardMuMuErrorMin_y"))->Fill(x[0],x[1],inputWeight);

    delete proxyoclptaptHisto;

  }

  delete proxy;

}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuTriggerResponse::WeightFromEta(
  const AliVParticle& track,
  AliMergeableCollectionProxy& proxylptaptHisto,
  int value)
{
  /// return weight according to eta.
  /// value = 0  : mean value
  /// value = 1  : + error stat
  /// value = -1 : - error stat

  if(fNEtaBin==0) return 0;

  // Get Histo
  TString lptaptHistoName = "";
  for (int i = 0; i < fNEtaBin; ++i)
  {
    if(fLptAptEtaRanges[i] < track.Eta() && track.Eta() < fLptAptEtaRanges[i+1])
      lptaptHistoName = Form("lpt_apt_pt_eta_%.2f_%.2f",fLptAptEtaRanges[i],fLptAptEtaRanges[i+1]);
  }

  //Get lpt TGraphAsymmErrors
  TGraphAsymmErrors* h_lptapt=0x0;
  if(proxylptaptHisto.GetObject(lptaptHistoName.Data()))
      h_lptapt = static_cast<TGraphAsymmErrors*>(proxylptaptHisto.GetObject(lptaptHistoName.Data()));
  if(!h_lptapt){
      printf("cannot get %s/%s \n",proxylptaptHisto.GetName(),lptaptHistoName.Data());
      return 0;
  }

  //Get axis
  TAxis* axis=0x0;
  if(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())))
      axis = static_cast<TAxis*>(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())));
  if(!axis){
      printf("cannot get %s/xaxis \n",proxylptaptHisto.GetName());
      return 0;
  }

  // Weight according to pT
  Double_t weight =0.;

  double xy[2]     ={0.,0.};
  double xyhigh[2] ={0.,0.};
  double xylow[2]  ={0.,0.};
  int bin          = 0;

  if(axis && axis->FindBin(track.Pt())) bin = axis->FindBin(track.Pt());
  if(bin!=0){
    h_lptapt->GetPoint(bin-1,xy[0],xy[1]);
    xylow[1] = h_lptapt->GetErrorYlow(bin-1);
    xyhigh[1] = h_lptapt->GetErrorYhigh(bin-1);
  }
  if(value      == 0) weight = xy[1];
  else if(value == 1) weight = xy[1] + xyhigh[1];
  else if(value == -1)weight = xy[1] - xylow[1];
  return weight;

}

//_____________________________________________________________________________
Double_t AliAnalysisMuMuTriggerResponse::WeightFromLocalBoard(
  const AliVParticle& track,
  Double_t LocalBoardNumber,
  AliMergeableCollectionProxy& proxylptaptHisto,
  int value)
{
  /// return weight according to the group of local board touched.
  /// value = 0  : mean value
  /// value = 1  : + error stat
  /// value = -1 : - error stat

  if(fNEtaBin==0) return 0;

  // Get sum of lpt/apt histo for a group of local boards
  TString lptaptHistoName = Form("lpt_apt_pt_eta_%.2f_%.2f_localboard_%.0f",fLptAptEtaRanges[0],fLptAptEtaRanges[fNEtaBin-1],LocalBoardNumber);

  //Get lpt TGraphAsymmErrors
  TGraphAsymmErrors* h_lptapt=0x0;
  if(proxylptaptHisto.GetObject(lptaptHistoName.Data()))
      h_lptapt = static_cast<TGraphAsymmErrors*>(proxylptaptHisto.GetObject(lptaptHistoName.Data()));
  if(!h_lptapt){
      printf("cannot get %s/%s \n",proxylptaptHisto.GetName(),lptaptHistoName.Data());
      return 0;
  }

  //Get axis
  TAxis* axis=0x0;
  if(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())))
      axis = static_cast<TAxis*>(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())));
  if(!axis){
      printf("cannot get %s/xaxis \n",proxylptaptHisto.GetName());
      return 0;
  }

  // Weight according to pT
  Double_t weight =0.;

  int bin = 0;

  if(axis && axis->FindBin(track.Pt())) bin= axis->FindBin(track.Pt());

  double xy[2]     ={0.,0.};
  double xyhigh[2] ={0.,0.};
  double xylow[2]  ={0.,0.};
  if(bin!=0) {
    h_lptapt->GetPoint(bin-1,xy[0],xy[1]);
    xylow[1] = h_lptapt->GetErrorYlow(bin-1);
    xyhigh[1] = h_lptapt->GetErrorYhigh(bin-1);
  }
  if(value      == 0) weight = xy[1];
  else if(value == 1) weight = xy[1] + xyhigh[1];
  else if(value == -1)weight = xy[1] - xylow[1];
  return weight;

}


//_____________________________________________________________________________
Double_t AliAnalysisMuMuTriggerResponse::WeightFromLocalBoardGroup(
  const AliVParticle& track,
  Double_t LocalBoardNumber,
  AliMergeableCollectionProxy& proxylptaptHisto,
  int value)
{
  /// return weight according to the group of local board touched.
  /// value = 0  : mean value
  /// value = 1  : + error stat
  /// value = -1 : - error stat

  if(!fLptAptEtaRanges || !fNEtaBin ) return 0;

  // the local board group (yes it is ugly but I was on a hurry...)
  int group[6][48] = {
    {26,27,28,29,48,49,50,51,165,166,167,168,143,144,145,146,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {6,7,8,22,23,24,25,44,45,46,47,9,10,11,30,31,32,33,52,53,54,55,123,124,125,139,140,141,142,161,162,163,164,169,170,171,172,147,148,149,150,126,127,128,0,0,0,0},
    {215,216,217,218,219,220,199,200,201,202,203,204,183,184,185,186,187,188,66,67,68,69,70,82,83,84,85,86,87,98,99,100,101,102,103,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {213,214,197,198,181,182,159,160,137,138,121,122,4,5,20,21,42,43,64,65,80,81,96,97,221,222,205,206,189,190,173,174,151,152,129,130,12,13,34,35,56,57,72,73,88,89,104,105},
    {211,212,195,196,179,180,157,158,135,136,119,120,2,3,18,19,40,41,62,63,78,79,94,95,223,224,207,208,191,192,175,176,153,154,131,132,14,15,36,37,58,59,74,75,90,91,106,107},
    {1,17,39,61,77,93,109,109,110,111,112,113,114,115,116,117,108,92,76,60,38,16,133,155,177,193,209,225,234,233,232,231,230,229,228,227,226,0,0,0,0,0,0}
  };
  int Group =0;

  // Find the groups
  for (int i = 0; i < 6; ++i)
  {
      for (int j = 0; j < 48; ++j)
      {
          if(group[i][j]==0) continue;
          if(group[i][j]==LocalBoardNumber) Group = i+1;
          else continue;
      }
  }

  // Get sum of lpt/apt histo for a group of local boards
  TString lptaptHistoName = "";
  for (int i = 0; i < fNEtaBin; ++i)
  {
    if(fLptAptEtaRanges[i] < track.Eta() && track.Eta() < fLptAptEtaRanges[i+1])
      lptaptHistoName = Form("lpt_apt_pt_eta_%.2f_%.2f_group_%d",fLptAptEtaRanges[i],fLptAptEtaRanges[i+1],Group);
  }

  //Get lpt TGraphAsymmErrors
  TGraphAsymmErrors* h_lptapt=0x0;
  if(proxylptaptHisto.GetObject(lptaptHistoName.Data()))
      h_lptapt = static_cast<TGraphAsymmErrors*>(proxylptaptHisto.GetObject(lptaptHistoName.Data()));
  if(!h_lptapt){
      printf("cannot get %s/%s \n",proxylptaptHisto.GetName(),lptaptHistoName.Data());
      return 0;
  }

  //Get axis
  TAxis* axis=0x0;
  if(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())))
      axis = static_cast<TAxis*>(proxylptaptHisto.GetObject(Form("%s_axis",lptaptHistoName.Data())));
  if(!axis){
      printf("cannot get %s/xaxis \n",proxylptaptHisto.GetName());
      return 0;
  }

  // Weight according to pT
  Double_t weight =0.;

  int bin = 0;

  if(axis && axis->FindBin(track.Pt())) bin= axis->FindBin(track.Pt());

  double xy[2]     ={0.,0.};
  double xyhigh[2] ={0.,0.};
  double xylow[2]  ={0.,0.};
  if(bin!=0) {
    h_lptapt->GetPoint(bin-1,xy[0],xy[1]);
    xylow[1]  = h_lptapt->GetErrorYlow(bin-1);
    xyhigh[1] = h_lptapt->GetErrorYhigh(bin-1);
  }
  if(value == 0)      weight = xy[1];
  else if(value == 1) weight = xy[1] + xyhigh[1];
  else if(value == -1)weight = xy[1] - xylow[1];
  return weight;

}

//_____________________________________________________________________________
void AliAnalysisMuMuTriggerResponse::SetLptAptEtaRange(Double_t* x, int xsize )
{

  if(!x) return;
  fNEtaBin = xsize;
  if(fNEtaBin==0)
  {
    fNEtaBin         = 0;
    return;
  }

  fLptAptEtaRanges = new Double_t[fNEtaBin];

  for (int i = 0; i < fNEtaBin ; ++i)
  {
    fLptAptEtaRanges[i]=x[i];
  }

 return;
}

