// $Id$
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//
//  Analysis Taks for reconstructed particle correlation 
//  (first implementation done for D mesons) with jets 
//  (use the so called Emcal framework)
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include "TROOT.h"
#include <TH3F.h>

#include "AliAnalysisTaskFlavourJetCorrelations.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEmcalJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliPicoTrack.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"

ClassImp(AliAnalysisTaskFlavourJetCorrelations)

//__________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations() :
  AliAnalysisTaskEmcalJet(),
  fUseMCInfo(kTRUE), 
  fCandidateType(),
  fPDGmother(),
  fNProngs(),
  fPDGdaughters(),
  fBranchName(),
  fOutput(0),
  fCuts(0),
  fMinMass(),
  fMaxMass(),  
  fJetArrName(0),
  fCandArrName(0),
  fLeadingJetOnly(kFALSE)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype, TString jetArrName) :
  AliAnalysisTaskEmcalJet(name),
  fUseMCInfo(kTRUE),
  fCandidateType(),
  fPDGmother(),
  fNProngs(),
  fPDGdaughters(),
  fBranchName(),
  fOutput(0),
  fCuts(0),
  fMinMass(),
  fMaxMass(),  
  fJetArrName(0),
  fCandArrName(0),
  fLeadingJetOnly(kFALSE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
 
  Info("AliAnalysisTaskFlavourJetCorrelations","Calling Constructor");
  fCuts=cuts;
  fCandidateType=candtype;
  const Int_t nptbins=fCuts->GetNPtBins();
  Float_t defaultSigmaD013[13]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040};

  switch(fCandidateType){
  case 0:
    fPDGmother=421;
    fNProngs=2;
    fPDGdaughters[0]=211;//pi 
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=0; //empty
    fPDGdaughters[3]=0; //empty
    fBranchName="D0toKpi";
    fCandArrName="D0";
    break;
  case 1: 
    fPDGmother=413;
    fNProngs=3;
    fPDGdaughters[1]=211;//pi soft
    fPDGdaughters[0]=421;//D0
    fPDGdaughters[2]=211;//pi fromD0
    fPDGdaughters[3]=321; // K from D0
    fBranchName="Dstar";
    fCandArrName="Dstar";
    //fSigmaD0=new Float_t[nptbins];
    if(nptbins<=13){
      for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]= defaultSigmaD013[ipt];
    } else {
      AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
    }
    break;
  default:
    printf("%d not accepted!!\n",fCandidateType);
    break;
  }

  if(fCandidateType==kD0toKpi)SetMassLimits(0.15,fPDGmother);
  if(fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
  fJetArrName=jetArrName;
  Printf("Jet read is %s",fJetArrName.Data());



  DefineOutput(1,TList::Class()); // histos
  DefineOutput(2,AliRDHFCuts::Class()); // my cuts
}
//___________________________________________________________________________
AliAnalysisTaskFlavourJetCorrelations::~AliAnalysisTaskFlavourJetCorrelations() {
  //
  // destructor
  //

  Info("~AliAnalysisTaskFlavourJetCorrelations","Calling Destructor");  
 
    delete fOutput;
    delete fCuts;
    
}

//___________________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::Init(){
  //
  // Initialization
  //
  if(fDebug > 1) printf("AnalysisTaskRecoJetCorrelations::Init() \n");
  switch(fCandidateType){
  case 0:
    {
      AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      copyfCuts->SetName("AnalysisCutsDzero");
      // Post the data
      PostData(2,copyfCuts);
    }
    break;
  case 1:
    {
      AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      copyfCuts->SetName("AnalysisCutsDStar");
      // Post the cuts
      PostData(2,copyfCuts);
    }
    break;
  default:
    return;
  }

  return;
}

//_________________________________________________
void AliAnalysisTaskFlavourJetCorrelations::UserExec(Option_t *)
{
  // user exec

  // Load the event
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
 
  TClonesArray *arrayDStartoD0pi=0;
  TClonesArray *trackArr = 0;
  TClonesArray *clusArr = 0;
  TClonesArray *jetArr = 0;
  TClonesArray *candidatesArr = 0;
//TClonesArray *isselArr = 0;

  if (!aodEvent && AODEvent() && IsStandardAOD()) {

    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());

    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
    }
  } else {
    arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
  }

  if (!arrayDStartoD0pi) {
    AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
//  return;
  } else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));   

  TClonesArray* mcArray = 0x0;
  if (fUseMCInfo) {
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      return;
    }
  }

  //retrieve jets
  trackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("PicoTracks"));
  clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClustersCorr"));
  jetArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrName));

  if(!trackArr){
    AliInfo(Form("Could not find the track array\n"));
    return;
  }

  if(!jetArr){
    Printf("JET array not found");
    return;
  }

  //retrieve reconstructed particles selected
  /*
  TString listname="AliAnalysisTaskSEDmesonsForJetCorrelations";
  TList *l=dynamic_cast<TList*>(InputEvent()->FindListObject(listname));
  TClonesArray *cla=dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(listname));
  // if(l){
  //   l->ls();
  // } else{
  //   Printf("list not found!!!!!!!!!!!");
  //   return;
  // } 
  if(!cla){
    Printf("cla not found!!!!!!!!!!!");
    return;
  } else {
    cla->ls();
  }
  */

  TString arrname="fCandidateArray";
  candidatesArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s%s",arrname.Data(),fCandArrName.Data())));
  if (!candidatesArr) {
    Printf("%s%s array not found",arrname.Data(),fCandArrName.Data());
    InputEvent()->GetList()->ls();
    return;
  }

  /*
  arrname="fIsSelectedArray";
  isselArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s%s",arrname.Data(),fCandArrName.Data())));
  if(!isselArr){
    Printf("%s%s array not found",arrname.Data(),fCandArrName.Data());
    InputEvent()->ls();
    return;
  }
  */

  //Histograms
  TH1I* hstat=(TH1I*)fOutput->FindObject("hstat");
  TH1F* hPtJetTrks=(TH1F*)fOutput->FindObject("hPtJetTrks");
  TH1F* hPhiJetTrks=(TH1F*)fOutput->FindObject("hPhiJetTrks");
  TH1F* hEtaJetTrks=(TH1F*)fOutput->FindObject("hEtaJetTrks");
  TH1F* hEjetTrks=(TH1F*)fOutput->FindObject("hEjetTrks");
  TH1F* hPtJet=(TH1F*)fOutput->FindObject("hPtJet");
  TH1F* hPhiJet=(TH1F*)fOutput->FindObject("hPhiJet");
  TH1F* hEtaJet=(TH1F*)fOutput->FindObject("hEtaJet");
  TH1F* hEjet=(TH1F*)fOutput->FindObject("hEjet");
  TH1F* hNtrArr=(TH1F*)fOutput->FindObject("hNtrArr");
  TH1F* hNJetPerEv=(TH1F*)fOutput->FindObject("hNJetPerEv");
  TH1F* hdeltaRJetTracks=((TH1F*)fOutput->FindObject("hdeltaRJetTracks"));

  hstat->Fill(0);

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;

  //Event selection
  Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
  TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
  if(!iseventselected) return;
  
  hstat->Fill(1);
  
  //trigger on jets
  Int_t njets=jetArr->GetEntriesFast();
  if(njets == 0) return;

  const Int_t nD=arrayDStartoD0pi->GetEntriesFast();
  hstat->Fill(2,nD);
  
  // counters for efficiencies
  Int_t icountReco = 0;
  
  //D* and D0 prongs needed to MatchToMC method
  // Int_t pdgDgDStartoD0pi[2]={421,211};
  // Int_t pdgDgD0toKpi[2]={321,211};
  //D0 from D0 bar
 
  // we start with jets
  Double_t ejet   = 0;
  Double_t phiJet = 0;
  Double_t etaJet = 0;
  Double_t ptjet = 0;
  Double_t leadingJet =0;

  Int_t ntrarr=trackArr->GetEntriesFast();
  hNtrArr->Fill(ntrarr);
  
  for(Int_t i=0;i<ntrarr;i++){
    AliVTrack *jtrack=static_cast<AliVTrack*>(trackArr->At(i));
    //reusing histograms
    hPtJetTrks->Fill(jtrack->Pt());
    hPhiJetTrks->Fill(jtrack->Phi());
    hEtaJetTrks->Fill(jtrack->Eta());
    hEjetTrks->Fill(jtrack->E());
  }

  
  //option to use only the leading jet
  if(fLeadingJetOnly){
    for (Int_t iJetsL = 0; iJetsL<njets; iJetsL++) {    
      AliEmcalJet* jetL = (AliEmcalJet*)jetArr->At(iJetsL);
      ptjet   = jetL->Pt();
      if(ptjet>leadingJet ) leadingJet = ptjet;

    }
  }

  Int_t cntjet=0;
  //loop over jets and consider only the leading jet in the event
  for (Int_t iJets = 0; iJets<njets; iJets++) {    
    AliEmcalJet* jet = (AliEmcalJet*)jetArr->At(iJets);
    if(!AcceptJet(jet)) continue;

    vector<int> DmesonInJetLabels(10);
    //Int_t iD=0;
    //jets variables
    ejet   = jet->E();
    phiJet = jet->Phi();
    etaJet = jet->Eta();
    ptjet = jet->Pt();
    
    // choose the leading jet
    if(fLeadingJetOnly && (ejet<leadingJet)) continue;
    //used for normalization
    hstat->Fill(3);
    cntjet++;
    // fill energy, eta and phi of the jet
    hEjet   ->Fill(ejet);
    hPhiJet ->Fill(phiJet);
    hEtaJet ->Fill(etaJet);
    hPtJet  ->Fill(ptjet);
    
    //loop on jet particles
    Int_t ntrjet=  jet->GetNumberOfTracks();    
    for(Int_t itrk=0;itrk<ntrjet;itrk++){
      
      AliPicoTrack* jetTrk=(AliPicoTrack*)jet->TrackAt(itrk,trackArr);     
      hdeltaRJetTracks->Fill(DeltaR(jet,jetTrk));
      
      // check MC in for the traks within the jet, look at D mesons
      // within the jet area
      //if(fUseMCInfo) FillMCDJetInfo(jetTrk,jet,mcArray,ptjet);
      
    }//end loop on jet tracks
    
    //retrieve charm candidates selected
    Int_t candidates = candidatesArr->GetEntriesFast();
    for(Int_t ic = 0; ic < candidates; ic++) {
      // D* candidates
      AliAODRecoDecayHF* charm = 0x0;
      charm=(AliAODRecoDecayHF*)candidatesArr->At(ic);
      
      FillHistogramsRecoJetCorr(charm, jet);      

    }
  } // end of jet loop 

  hNJetPerEv->Fill(cntjet);

  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));
  
  PostData(1,fOutput);

}

//________________________________________ terminate ___________________________
void AliAnalysisTaskFlavourJetCorrelations::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  Info("Terminate"," terminate");
  AliAnalysisTaskSE::Terminate();

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
}
//___________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::UserCreateOutputObjects() { 
 // output 
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  // define histograms
  DefineHistoForAnalysis();

  PostData(1,fOutput);
  return;
}
//_________________________________________________________________
void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t range, Int_t pdg){
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  // compute the Delta mass for the D*
   if(fCandidateType==kDstartoKpipi){
     Float_t mass1=0;
     mass1=TDatabasePDG::Instance()->GetParticle(421)->Mass();
     mass = mass-mass1;
  }
   cout<<"mass ---------------"<<endl;
  fMinMass = mass-range;
  fMaxMass = mass+range;
  
  AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
  if (fMinMass<0 || fMaxMass<=0 || fMaxMass<fMinMass) AliFatal("Wrong mass limits!\n");
}
//_________________________________________________________________
void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t lowlimit, Double_t uplimit){
  if(uplimit>lowlimit)
    {
      fMinMass = lowlimit;
      fMaxMass = uplimit;
    }
  else{
    printf("Error! Lower limit larger than upper limit!\n");
    fMinMass = uplimit - uplimit*0.2;
    fMaxMass = uplimit;
  }
}

//__________________________________________________________________
Bool_t AliAnalysisTaskFlavourJetCorrelations::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
  if(nptbins>30) {
    AliInfo("Maximum number of bins allowed is 30!");
    return kFALSE;
  }
  if(!width) return kFALSE;
  for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
  return kTRUE;
}

//__________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelations::Z(AliVParticle* part,AliEmcalJet* jet) const{
  if(!part || !jet){
    printf("Particle or jet do not exist!\n");
    return -99;
  }
  Double_t p[3],pj[3];
  Bool_t okpp=part->PxPyPz(p);
  Bool_t okpj=jet->PxPyPz(pj);
  if(!okpp || !okpj){
    printf("Problems getting momenta\n");
    return -999;
  }
  Double_t pjet=jet->P();
  Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet*pjet);
  return z;
}
//___________________________________ histograms _______________________________________

Bool_t  AliAnalysisTaskFlavourJetCorrelations::DefineHistoForAnalysis(){
 
  // Statistics 
  TH1I* hstat=new TH1I("hstat","Statistics",5,-0.5,4.5);
  hstat->GetXaxis()->SetBinLabel(1,"N ev anal");
  hstat->GetXaxis()->SetBinLabel(2,"N ev sel");
  hstat->GetXaxis()->SetBinLabel(3,"N cand sel cuts");
  hstat->GetXaxis()->SetBinLabel(4,"N jets");
  hstat->GetXaxis()->SetBinLabel(5,"N cand in jet");
  // if(fUseMCInfo) {
  //   hstat->GetXaxis()->SetBinLabel(7,"N D");
  //   hstat->GetXaxis()->SetBinLabel(8,"N D in jet");

  // }

  hstat->SetNdivisions(1);
  fOutput->Add(hstat);

  const Int_t nbinsmass=200;

  
  if(fCandidateType==kDstartoKpipi) 
    {
      // TH2F *hDiff = new TH2F("hDiff","M(kpipi)-M(kpi)",500,fMinMass,fMaxMass,100, 0.,30.);
      // hDiff->SetStats(kTRUE);
      // hDiff->GetXaxis()->SetTitle("M(kpipi)-M(kpi) GeV/c^{2}");
      // hDiff->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
      // fOutput->Add(hDiff);
  
      TH2F* hDiffSideBand = new TH2F("hDiffSideBand","M(kpipi)-M(kpi) Side Band Background",nbinsmass,fMinMass,fMaxMass,100, 0.,30.);
      hDiffSideBand->SetStats(kTRUE);
      hDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^{2}");
      hDiffSideBand->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
      fOutput->Add(hDiffSideBand); 
 
      //correlation histograms
      // fPhi = new TH1F("phi","Delta phi between Jet axis and DStar ",25,-1.57,4.72);
      // fPhi->SetStats(kTRUE);
      // fPhi->GetXaxis()->SetTitle("#Delta #phi (rad)");
      // fPhi->GetYaxis()->SetTitle("Entries");
      // fOutput->Add(fPhi);

      // fDphiD0Dstar = new TH1F("phiD0Dstar","Delta phi between D0 and DStar ",1000,-6.5,6.5);
      // fOutput->Add(fDphiD0Dstar);

      // fPhiBkg = new TH1F("phiBkg","Delta phi between Jet axis and DStar background ",25,-1.57,4.72);
      // fPhiBkg->SetStats(kTRUE);
      // fPhiBkg->GetXaxis()->SetTitle("#Delta #phi (rad)");
      // fPhiBkg->GetYaxis()->SetTitle("Entries");
      // fOutput->Add(fPhiBkg);

      // TH1F* hRECOPtDStar = new TH1F("hRECODStar","RECO DStar pt distribution",30,0,30);
      // hRECOPtDStar->SetStats(kTRUE);
      // hRECOPtDStar->SetLineColor(2);
      // hRECOPtDStar->GetXaxis()->SetTitle("GeV/c");
      // hRECOPtDStar->GetYaxis()->SetTitle("Entries");
      // fOutput->Add(hRECOPtDStar);
  
      // TH1F* hRECOPtBkg = new TH1F("hRECOptBkg","RECO pt distribution side bands",30,0,30);
      // hRECOPtBkg->SetStats(kTRUE);
      // hRECOPtBkg->SetLineColor(2);
      // hRECOPtBkg->GetXaxis()->SetTitle("GeV/c");
      // hRECOPtBkg->GetYaxis()->SetTitle("Entries");
      // fOutput->Add(hRECOPtBkg);

      TH1F* hPtPion = new TH1F("hPtPion","Primary pions candidates pt ",500,0,10);
      hPtPion->SetStats(kTRUE);
      hPtPion->GetXaxis()->SetTitle("GeV/c");
      hPtPion->GetYaxis()->SetTitle("Entries");
      fOutput->Add(hPtPion);

    }

  // jet related fistograms
  
  TH1F* hEjetTrks      = new TH1F("hEjetTrks",  "Jet tracks energy distribution;Energy (GeV)",500,0,200);
  TH1F* hPhiJetTrks    = new TH1F("hPhiJetTrks","Jet tracks #phi distribution; #phi (rad)",  200,0,6.30);
  TH1F* hEtaJetTrks    = new TH1F("hEtaJetTrks","Jet tracks #eta distribution; #eta",  100,-1.5,1.5);
  TH1F* hPtJetTrks     = new TH1F("hPtJetTrks",  "Jet tracks Pt distribution; p_{T} (GeV/c^{2})",500,0,200);
  
  TH1F* hEjet      = new TH1F("hEjet",  "Jet energy distribution;Energy (GeV)",500,0,200);
  TH1F* hPhiJet    = new TH1F("hPhiJet","Jet #phi distribution; #phi (rad)",  200,0,6.30);
  TH1F* hEtaJet    = new TH1F("hEtaJet","Jet #eta distribution; #eta",  100,-1.5,1.5);
  TH1F* hPtJet      = new TH1F("hPtJet",  "Jet Pt distribution; p_{T} (GeV/c^{2})",500,0,200);

  TH1F* hPtJetWithD=new TH1F("hPtJetWithD","D-Jet Pt distribution; p_{T} (GeV/c^{2})",500,0,200);
  TH1F* hdeltaRJetTracks=new TH1F("hdeltaRJetTracks","Delta R of tracks in the jets",200, 0.,10.);

  TH1F* hNtrArr= new TH1F("hNtrArr", "Number of tracks in the array of jets; number of tracks",500,0,1000);
  TH1F *hNJetPerEv=new TH1F("hNJetPerEv","Number of jets used per event; number of jets/ev",10,-0.5,9.5);

  fOutput->Add(hEjetTrks);
  fOutput->Add(hPhiJetTrks);
  fOutput->Add(hEtaJetTrks);
  fOutput->Add(hPtJetTrks);
  fOutput->Add(hEjet);
  fOutput->Add(hPhiJet);
  fOutput->Add(hEtaJet);
  fOutput->Add(hPtJet);
  fOutput->Add(hPtJetWithD);
  fOutput->Add(hdeltaRJetTracks);
  fOutput->Add(hNtrArr);
  fOutput->Add(hNJetPerEv);

  /*
  if(fUseMCInfo){
    fhMomjetpartPdg=new TH1F("fhMomjetpartPdg","Jet particles' mothers distribution;PDG code;Counts;",1100,-550,550);
    fOutput->Add(fhMomjetpartPdg);
    fptDinjetallvsptjMC=new TH2F("fptDinjetallvsptjMC","p_{t} distribution of D within a jet (all DeltaR) vs p_{t}^{jet};p_{t}^{D};p_{t}^{jet}",100, 0.,30.,500,0.,200.);
    fOutput->Add(fptDinjetallvsptjMC);
    fptJwithDMC=new TH1F("fptJwithDMC","p_{t}^{jet} with a D meson (all DeltaR);p_{t}^{jet};Counts",500,0.,200.);
    fOutput->Add(fptJwithDMC);

    fptDinjet04vsptjMC=new TH2F("fptDinjet04vsptjMC","p_{t} distribution of D within a jet (DeltaR 0.4) vs p_{t}^{jet};p_{t}^{D};p_{t}^{jet}",100, 0.,30.,500,0.,200.);
    fOutput->Add(fptDinjet04vsptjMC);
    TH1F* hdeltaRDMC=new TH1F("hdeltaRDMC","#Delta R for MC tagged D mesons; #Delta R_{D}^{MC}",200, 0.,10.);
    fOutput->Add(hdeltaRDMC);
  }
  */
  TH1F* hDeltaRD=new TH1F("hDeltaRD","#Delta R distribution of D candidates selected;#Delta R",200, 0.,10.);
  fOutput->Add(hDeltaRD);

  TH3F* hdeltaPhiDja=new TH3F("hdeltaPhiDja", "#Delta#phi D-jet (jet p_{T} > threshold)",nbinsmass,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  hdeltaPhiDja->GetXaxis()->SetTitle("mass (GeV/c)");
  hdeltaPhiDja->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
  hdeltaPhiDja->GetYaxis()->SetTitle("#Delta#phi (rad)");
  // TH3F* hdeltaPhiDjl=new TH3F("hdeltaPhiDjl", Form("#Delta#phi D-jet (jet p_{T} < %.0f (GeV/c^{2}))",fJetPtThr),500,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  // hdeltaPhiDjl->GetXaxis()->SetTitle("mass (GeV/c)");
  // hdeltaPhiDjl->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
  // hdeltaPhiDjl->GetYaxis()->SetTitle("#Delta#phi (rad)");
  // TH3F* hdeltaPhiDjh=new TH3F("hdeltaPhiDjh", Form("#Delta#phi D-jet (jet p_{T} > %.0f (GeV/c^{2}))",fJetPtThr),500,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  // hdeltaPhiDjh->GetXaxis()->SetTitle("mass (GeV/c)");
  // hdeltaPhiDjh->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
  // hdeltaPhiDjh->GetYaxis()->SetTitle("#Delta#phi (rad)");
  fOutput->Add(hdeltaPhiDja);
  // fOutput->Add(hdeltaPhiDjl);
  // fOutput->Add(hdeltaPhiDjh);

  //background (side bands for the Dstar and like sign for D0)

  TH3F* hdeltaPhiDjaB=new TH3F("hdeltaPhiDjaB", "#Delta#phi D-jet (all jet p_{T})",nbinsmass,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  hdeltaPhiDjaB->GetXaxis()->SetTitle("mass (GeV/c)");
  hdeltaPhiDjaB->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
  hdeltaPhiDjaB->GetYaxis()->SetTitle("#Delta#phi (rad)");
  // TH3F* hdeltaPhiDjlB=new TH3F("hdeltaPhiDjlB", Form("#Delta#phi D-jet (jet p_{T} < %.0f (GeV/c^{2}))",fJetPtThr),1500,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  // hdeltaPhiDjlB->GetXaxis()->SetTitle("mass (GeV/c)");
  // hdeltaPhiDjlB->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
  // hdeltaPhiDjlB->GetYaxis()->SetTitle("#Delta#phi (rad)");
  // TH3F* hdeltaPhiDjhB=new TH3F("hdeltaPhiDjhB", Form("#Delta#phi D-jet (jet p_{T} > %.0f (GeV/c^{2}))",fJetPtThr),1500,fMinMass,fMaxMass,100, 0.,30.,50,(-1)*TMath::Pi()/2.,3./2.*TMath::Pi());
  // hdeltaPhiDjhB->GetXaxis()->SetTitle("mass (GeV/c)");
  // hdeltaPhiDjhB->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c^{2})");
  // hdeltaPhiDjhB->GetYaxis()->SetTitle("#Delta#phi (rad)");
  fOutput->Add(hdeltaPhiDjaB);
  // fOutput->Add(hdeltaPhiDjlB);
  // fOutput->Add(hdeltaPhiDjhB);

  TH2F* hInvMassptD = new TH2F("hInvMassptD","D (Delta R < 0.4) invariant mass distribution p_{T}^{j} > threshold",nbinsmass,fMinMass,fMaxMass,100,0.,50.);
  hInvMassptD->SetStats(kTRUE);
  hInvMassptD->GetXaxis()->SetTitle("mass (GeV/c)");
  hInvMassptD->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");

  fOutput->Add(hInvMassptD);
  // fMasspjDeltaR=new TH3F("fMasspjDeltaR","Mass vs p^{jet} vs #Delta R;Mass (Gev/c);p^{jet}(Gev/c^{2});#Delta R",1500,fMinMass,fMaxMass,100,0.,50.,100,0.,1.);
  // fMasspjDeltaR->SetStats(kFALSE);
  // fOutput->Add(fMasspjDeltaR);

  TH3F* hzptD=new TH3F("hzptD","Fragmentation function (DeltaR < 0.4)",100,0.,1.2,nbinsmass,fMinMass,fMaxMass,100,0,50);
  hzptD->SetStats(kTRUE);
  hzptD->GetXaxis()->SetTitle("z=p_{D} #cdot p_{j}/|p_{j}|^{2}");
  hzptD->GetYaxis()->SetTitle("mass (GeV/c)");
  hzptD->GetZaxis()->SetTitle("p_{t}^{D} (GeV/c)");
  fOutput->Add(hzptD);

  TH3F* hzptDB=new TH3F("hzptDB","Fragmentation function (DeltaR < 0.4) - Side Bands",100,0.,1.2,nbinsmass,fMinMass,fMaxMass,100,0.,50.);
  hzptDB->SetStats(kTRUE);
  hzptDB->GetXaxis()->SetTitle("z=p_{D} #cdot p_{j}/|p_{j}|^{2}");
  hzptDB->GetYaxis()->SetTitle("mass (GeV/c)");
  hzptDB->GetZaxis()->SetTitle("p_{t}^{D} (GeV/c)");
  fOutput->Add(hzptDB);
  //TH1F* hResZ      = new TH1F("hResZ","Fragmentation function ",50,0,1);
  //  TH1F* hResZBkg   = new TH1F("hResZBkg","Fragmentation function background",50,0,1);  

  //fOutput->Add(hResZ);
  //fOutput->Add(hResZBkg);


  return kTRUE; 
}

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsRecoJetCorr(AliAODRecoDecayHF* candidate, AliEmcalJet *jet){

  Double_t ptD=candidate->Pt();
  Double_t ptjet=jet->Pt();
  Double_t deltaR=DeltaR(candidate,jet);
  Double_t deltaphi = jet->Phi()-candidate->Phi();
  if(deltaphi<=-(TMath::Pi())/2) deltaphi = deltaphi+2*(TMath::Pi());
  if(deltaphi>(3*(TMath::Pi()))/2) deltaphi = deltaphi-2*(TMath::Pi());
  Double_t z=Z(candidate,jet);

  TH1F* hDeltaRD=(TH1F*)fOutput->FindObject("hDeltaRD");
  hDeltaRD->Fill(deltaR); 
  TH1I* hstat=(TH1I*)fOutput->FindObject("hstat");
  hstat->Fill(4);
  TH1F* hPtJetWithD=(TH1F*)fOutput->FindObject("hPtJetWithD");
  hPtJetWithD->Fill(ptjet);

  if(fCandidateType==kD0toKpi) {

    FillHistogramsD0JetCorr(candidate,deltaphi,z,ptD,deltaR, AODEvent());

  }

  if(fCandidateType==kDstartoKpipi) {
    AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)candidate;
    FillHistogramsDstarJetCorr(dstar,deltaphi,z,ptD,deltaR);

  }

}

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t deltaR, AliAODEvent* aodEvent){
  Double_t masses[2]={0.,0.};
  Int_t pdgdaughtersD0[2]={211,321};//pi,K 
  Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 

  masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
  masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
  
 
  Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);
  if(isselected==1 || isselected==3) {

    FillHistogramsD(masses[0],dPhi,z, ptD, deltaR);
  }
  if(isselected>=2) {

    FillHistogramsD(masses[1],dPhi,z, ptD, deltaR);
  }

}

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi, Double_t z, Double_t ptD, Double_t deltaR){
  AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
  Double_t deltamass= dstar->DeltaInvMass(); 
  Double_t massD0= dstar->InvMassD0();

  TH1F* hPtPion=(TH1F*)fOutput->FindObject("hPtPion");
  hPtPion->Fill(softpi->Pt());
  FillHistogramsD(deltamass,dPhi,z, ptD, deltaR);
  // evaluate side band background
  TH2F* hDiffSideBand=(TH2F*)fOutput->FindObject("hDiffSideBand");

  TH3F* hdeltaPhiDjaB=(TH3F*)fOutput->FindObject("hdeltaPhiDjaB");

  TH3F* hzptDB=(TH3F*)fOutput->FindObject("hzptDB");

  Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  Int_t bin = fCuts->PtBin(ptD);
  Float_t fourSigmal= mPDGD0-3.5*fSigmaD0[bin] , sixSigmal= mPDGD0-5.*fSigmaD0[bin];
  Float_t fourSigmar= mPDGD0+3.5*fSigmaD0[bin] , sixSigmar= mPDGD0+5.*fSigmaD0[bin];

  if((massD0>sixSigmal && massD0<fourSigmal) || (massD0>fourSigmar && massD0<=sixSigmar)){  
    hDiffSideBand->Fill(deltamass,ptD); // M(Kpipi)-M(Kpi) side band background    
    hdeltaPhiDjaB->Fill(deltamass,ptD,dPhi);

    if(deltaR<0.4){  // evaluate in the near side	
      hzptDB->Fill(z,deltamass,ptD);	
    }

  }  //  SideBandBackground(dstar, dPhi, z, ptD, deltaR);
 
}

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsD(Double_t mass,Double_t dphi, Double_t z,Double_t ptD, Double_t deltaR){
  TH3F* hdeltaPhiDja=((TH3F*)fOutput->FindObject("hdeltaPhiDja"));
  hdeltaPhiDja->Fill(mass,ptD,dphi);

  if(deltaR<0.4) {
    TH3F* hzptD=(TH3F*)fOutput->FindObject("hzptD");
    hzptD->Fill(z,mass,ptD);

    TH2F* hInvMassptD=(TH2F*)fOutput->FindObject("hInvMassptD");
    hInvMassptD->Fill(mass,ptD);
  }
}
//______________________________ side band background for D*___________________________________

// void AliAnalysisTaskFlavourJetCorrelations::SideBandBackground(AliAODRecoCascadeHF *candDstar, Double_t dPhi, Double_t z, Double_t ptD, Double_t deltaR){

//   //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
//   // (expected detector resolution) on the left and right frm the D0 mass. Each band
//   //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
//   TH2F* hDiffSideBand=(TH2F*)fOutput->FindObject("hDiffSideBand");

//   TH3F* hdeltaPhiDjaB=(TH3F*)fOutput->FindObject("hdeltaPhiDjaB");

//   TH3F* hzptDB=(TH3F*)fOutput->FindObject("hzptDB");

//   Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();

//   Int_t bin = fCuts->PtBin(candDstar->Pt());
//   Float_t fourSigmal= mPDGD0-3.5*fSigmaD0[bin] , sixSigmal= mPDGD0-5.*fSigmaD0[bin];
//   Float_t fourSigmar= mPDGD0+3.5*fSigmaD0[bin] , sixSigmar= mPDGD0+5.*fSigmaD0[bin];

//   Double_t invM=candDstar->InvMassD0(),  deltaM=candDstar->DeltaInvMass(); //invMDstar=candDstar->InvMassDstarKpipi(),
//   Printf("Inv mass = %f between %f and %f or %f and %f?",invM, sixSigmal,fourSigmal,fourSigmar,sixSigmar);
//   Double_t ptD=candDstar->Pt();
//   //Double_t ptjet=jet->Pt();
//   Double_t dPhi=jet->Phi()-candDstar->Phi();
//   Double_t deltaR=DeltaR(candDstar,jet);
//   if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
//   if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();

//   if((invM>sixSigmal && invM<fourSigmal) || (invM>fourSigmar && invM<=sixSigmar)){  
//     hDiffSideBand->Fill(deltaM,ptD); // M(Kpipi)-M(Kpi) side band background    
//     hdeltaPhiDjaB->Fill(deltaM,ptD,dPhi);

//     if(deltaR<0.4){  // evaluate in the near side	
//       hzptDB->Fill(Z(candDstar,jet),deltaM,ptD);	
//     }

//   }
// }

Float_t AliAnalysisTaskFlavourJetCorrelations::DeltaR(AliVParticle *p1, AliVParticle *p2) const {
  //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)

  if(!p1 || !p2) return -1;
  Double_t phi1=p1->Phi(),eta1=p1->Eta();
  Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;

  Double_t dPhi=phi1-phi2;
  if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
  if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());

  Double_t dEta=eta1-eta2;
  Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
  return deltaR;
 
}

/*
//_____________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::IsD(Int_t pdg) const{
  Int_t abspdg=TMath::Abs(pdg);
  if(abspdg>400 && abspdg<500) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::IsD(Int_t pdg,Int_t abspdgD) const{
  Int_t abspdg=TMath::Abs(pdg);
  if(abspdg==abspdgD) return kTRUE;
  else return kFALSE;
}

//_____________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::PartFromC(AliMCParticle* mother) const{
  Int_t pdgmoth=mother->PdgCode();
  if(TMath::Abs(pdgmoth)==4) return kTRUE;
  else return kFALSE;
}

Int_t AliAnalysisTaskFlavourJetCorrelations::GetFirstMother(Int_t labpart, TClonesArray *mcArray) const{
  AliAODMCParticle* partMC=(AliAODMCParticle*)mcArray->UncheckedAt(labpart);
  if (!partMC) return -2;
  Int_t labmom=labpart;
  Printf("Starting from %d",labmom);
  while(labmom>-1){
    labmom=labpart;
    partMC=(AliAODMCParticle*)mcArray->UncheckedAt(labmom);
    if (!partMC) return -2;
    labmom= partMC->GetMother();
    Printf("Lab mom %d",labmom);
  }
  Printf("Return labpart %d", labpart);
  return labpart;
}


// -------------------------------------- check the PDG -----------------------------------------

Int_t AliAnalysisTaskFlavourJetCorrelations::FindPDGInFamily(Int_t labpart,Int_t pdgcode, TClonesArray *mcArray) const{

  //return the label of the particle which is a "pdgcode" in the family
  Printf("FindPDGInFamily label %d, pdg %d, mcarray %p",labpart,pdgcode,mcArray);
  AliAODMCParticle* partMC=(AliAODMCParticle*)mcArray->UncheckedAt(labpart);
  if (!partMC) return -2;
  Int_t labmom=labpart;
  Printf("Starting from %d",labmom);
  while(labmom>-1){

    partMC=(AliAODMCParticle*)mcArray->UncheckedAt(labmom);
    if (!partMC) return -2;
    Int_t pdgmom=partMC->GetPdgCode();
    if(pdgmom==pdgcode) return labmom;
    labmom= partMC->GetMother();
    Printf("Lab mom %d",labmom);
  }

  return -1;

}

// ------------------------- check on MC the distance between D meson and jet ----------------------

Bool_t AliAnalysisTaskFlavourJetCorrelations::FillMCDJetInfo(AliPicoTrack *jetTrk, AliEmcalJet* jet, TClonesArray *mcArray, Double_t ptjet){
  
  Bool_t foundD = kFALSE;
  vector<int> DmesonInJetLabels(10);
  
  Int_t jtlabel=jetTrk->GetLabel();
  if(jtlabel<=0) return foundD;
  AliAODMCParticle* jetMCpart=(AliAODMCParticle*)mcArray->UncheckedAt(jtlabel);
  if(!jetMCpart) return foundD;
  printf("AliMCParticle %d, %p\n",1,jetMCpart);
  
  Int_t labDmeson=FindPDGInFamily(jtlabel,fPDGmother, mcArray);
  if(labDmeson>0){
    AliAODMCParticle *partDmeson=(AliAODMCParticle*)mcArray->UncheckedAt(labDmeson);
    fhMomjetpartPdg->Fill(partDmeson->GetPdgCode());
    
    //tmp
    Int_t momjetpartlabel=labDmeson;
    
    Int_t iD=5;
    Bool_t exists=kFALSE;
    for(Int_t k=0;k<iD;k++){
      if(momjetpartlabel==DmesonInJetLabels[k]) {//mother already found
	exists=kTRUE;
	break;
      }
    }
    if(!exists) {
      DmesonInJetLabels[iD]=momjetpartlabel;
      AliDebug(2,Form("D meson number %d found: label %d\n",iD,DmesonInJetLabels[iD]));
      hstat->Fill(6);
      iD++;
      foundD=kTRUE;
      
    }
  }
  
  if(fUseMCInfo && foundD) {
    fptJwithDMC->Fill(ptjet); //filled only once per jet, not per each D meson
    Int_t iD=5;
    // loop over the D within the jet 
    for(Int_t kD=0;kD<iD;kD++){
      
      AliAODMCParticle* momjetpart=(AliAODMCParticle*)mcArray->At(DmesonInJetLabels[kD]);
      Double_t ptD=momjetpart->Pt(),etaD=momjetpart->Eta(), phiD=momjetpart->Phi();
      fptDinjetallvsptjMC->Fill(ptD,ptjet);
      
      Double_t deltaRD=DeltaR(jet,momjetpart);
      
      ((TH1F*)fOutput->FindObject("hdeltaRDMC"))->Fill(deltaRD); //Delta R of D mesons (MC)
      
      Double_t z=Z(momjetpart,jet);
      
      if(deltaRD<0.4) {
	hstat->Fill(7);
	//comment if you prefer to ask for DeltaR of the daughters < 0.4 and uncomment below
	fptDinjet04vsptjMC->Fill(ptD,ptjet);
	fzptDptj->Fill(z,ptD,ptjet);
      }
      
      
    }//end loop on MC D
    
    return foundD;
    
  }
}
*/  
