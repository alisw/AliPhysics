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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison with the MC truth and cut variables distributions.
//
// Authors: A.Dainese, andrea.dainese@lnl.infn.it
// Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Carmelo Di Giglio, carmelo.digiglio@ba.infn.it (like sign)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Mass.h"


ClassImp(AliAnalysisTaskSED0Mass)


//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass():
AliAnalysisTaskSE(),
fOutputMass(0),
fDistr(0),
fNentries(0), 
fChecks(0),
//fCutList(0),
fCuts(0),
fArray(0),
fReadMC(0),
fCutOnDistr(0),
fNPtBins(1),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)


{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass(const char *name,AliRDHFCutsD0toKpi* cuts):
AliAnalysisTaskSE(name),
fOutputMass(0), 
fDistr(0),
fNentries(0),
fChecks(0),
//fCutList(0),
fCuts(0),
fArray(0),
fReadMC(0),
fCutOnDistr(0),
fNPtBins(1),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)


{
  // Default constructor

  fNPtBins=cuts->GetNPtBins();
  fTotPosPairs=new Int_t[fNPtBins];
  fTotNegPairs=new Int_t[fNPtBins];
  for(Int_t i=0;i<fNPtBins;i++) {fTotPosPairs[i]=0; fTotNegPairs[i]=0;}
    
  //fCuts=new AliRDHFCutsD0toKpi(*cuts);
  fCuts=cuts;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TList container (distributions)
  DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TH1F container (number of events)
  DefineOutput(3,TH1F::Class());  //My private output
  // Output slot #4 writes into a TList container (quality check)
  DefineOutput(4,TList::Class());  //My private output
  // Output slot #5 writes into a TList container (cuts)
  DefineOutput(5,AliRDHFCutsD0toKpi::Class());  //My private output
  //DefineOutput(5,TList::Class());  //My private output

}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
   if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fDistr) {
    delete fDistr;
    fDistr = 0;
  }
  if (fChecks) {
    delete fChecks;
    fChecks = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
 
//   if(fCutList){
//     delete fCutList;
//     fCutList=0;
//   }
 
}  

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");

  
//   fCutList=new TList();
//   fCutList->SetOwner();
//   fCutList->SetName("listofCutsObj");
//   fCutList->Add(fCuts);
  
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
  // Post the data
  PostData(5,copyfCuts);
  //PostData(5,fCutList);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fDistr = new TList();
  fDistr->SetOwner();
  fDistr->SetName("distributionslist");

  fChecks = new TList();
  fChecks->SetOwner();
  fChecks->SetName("checkHistograms");

  TString nameMass=" ",nameSgn27=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassNocutsS =" ",nameMassNocutsB =" ", namedistr=" ";

  for(Int_t i=0;i<fCuts->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn27="histSgn27_";
    nameSgn27+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;
    nameMassNocutsS="hMassS_";
    nameMassNocutsS+=i;
    nameMassNocutsB="hMassB_";
    nameMassNocutsB+=i;

    //histograms of cut variable distributions

    //  pT
    namedistr="hptpiS_";
    namedistr+=i;
    TH1F *hptpiS = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptKS_";
    namedistr+=i;
    TH1F *hptKS = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB_";
    namedistr+=i;
    TH1F *hptB = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    //  pT no mass cut
    namedistr="hptpiSnoMcut_";
    namedistr+=i;
    TH1F *hptpiSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptKSnoMcut_";
    namedistr+=i;
    TH1F *hptKSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB1prongnoMcut_";
    namedistr+=i;
    TH1F *hptB1pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    namedistr="hptB2prongsnoMcut_";
    namedistr+=i;
    TH1F *hptB2pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

    
    //  dca
    namedistr="hdcaS_";
    namedistr+=i;
    TH1F *hdcaS = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
    namedistr="hdcaB_";
    namedistr+=i;
    TH1F *hdcaB = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);

    //  costhetastar
    namedistr="hcosthetastarS_";
    namedistr+=i;
    TH1F *hcosthetastarS = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);
    namedistr="hcosthetastarB_";
    namedistr+=i;
    TH1F *hcosthetastarB = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

    // impact parameter
    namedistr="hd0piS_";
    namedistr+=i;
    TH1F *hd0piS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,-0.1,0.1);

    namedistr="hd0KS_";
    namedistr+=i;
    TH1F *hd0KS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,-0.1,0.1);
    namedistr="hd0B_";
    namedistr+=i;
    TH1F *hd0B = new TH1F(namedistr.Data(), "Impact parameter distribution;d0 [cm]",200,-0.1,0.1);

    namedistr="hd0d0S_";
    namedistr+=i;
    TH1F *hd0d0S = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
    namedistr="hd0d0B_";
    namedistr+=i;
    TH1F *hd0d0B = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

    //  costhetapoint
    namedistr="hcosthetapointS_";
    namedistr+=i;
    TH1F *hcosthetapointS = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);
    namedistr="hcosthetapointB_";
    namedistr+=i;
    TH1F *hcosthetapointB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

    namedistr="hcosthpointd0d0S_";
    namedistr+=i;
    TH2F *hcosthpointd0d0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);
    namedistr="hcosthpointd0d0B_";
    namedistr+=i;
    TH2F *hcosthpointd0d0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

    fDistr->Add(hptpiS);
    fDistr->Add(hptKS);
    fDistr->Add(hptB);

    fDistr->Add(hptpiSnoMcut);
    fDistr->Add(hptKSnoMcut);
    fDistr->Add(hptB1pnoMcut);
    fDistr->Add(hptB2pnoMcut);

    fDistr->Add(hdcaS);
    fDistr->Add(hdcaB);

    fDistr->Add(hd0piS);
    fDistr->Add(hd0KS);
    fDistr->Add(hd0B);

    fDistr->Add(hd0d0S);
    fDistr->Add(hd0d0B);

    fDistr->Add(hcosthetastarS);
    fDistr->Add(hcosthetastarB);

    fDistr->Add(hcosthetapointS);
    fDistr->Add(hcosthetapointB);

    fDistr->Add(hcosthpointd0d0S);
    fDistr->Add(hcosthpointd0d0B);


    //histograms of invariant mass distributions

    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpMl=(TH1F*)tmpMt->Clone();
    tmpMt->Sumw2();
    tmpMl->Sumw2();

    //to compare with AliAnalysisTaskCharmFraction
    TH1F* tmpS27t = new TH1F(nameSgn27.Data(),"D^{0} invariant mass in M(D^{0}) +/- 27 MeV - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpS27l=(TH1F*)tmpS27t->Clone();
    tmpS27t->Sumw2();
    tmpS27l->Sumw2();
 
    //distribution w/o cuts
    //   TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,0.7,3.);
    TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,1.56484,2.16484); //range (MD0-300MeV, mD0 + 300MeV)
    TH1F *tmpMB=(TH1F*)tmpMS->Clone();
    tmpMB->SetName(nameMassNocutsB.Data());
    tmpMS->Sumw2();
    tmpMB->Sumw2();

    //MC signal and background
    TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpSl=(TH1F*)tmpSt->Clone();
    tmpSt->Sumw2();
    tmpSl->Sumw2();

    TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpBl=(TH1F*)tmpBt->Clone();
    tmpBt->Sumw2();
    tmpBl->Sumw2();

    //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
    TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpRl=(TH1F*)tmpRt->Clone();
    tmpRt->Sumw2();
    tmpRl->Sumw2();
  //  printf("Created histograms %s\t%s\t%s\t%s\n",tmpM->GetName(),tmpS->GetName(),tmpB->GetName(),tmpR->GetName());

    fOutputMass->Add(tmpMt);
    fOutputMass->Add(tmpSt);
    fOutputMass->Add(tmpS27t);
    fOutputMass->Add(tmpBt);
    fOutputMass->Add(tmpRt);

    fDistr->Add(tmpMS);
    fDistr->Add(tmpMB);


  }

  //histograms for vertex checking
  TString checkname="hptGoodTr";

  TH1F* hptGoodTr=new TH1F(checkname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Number",200,0.,8.);
  hptGoodTr->SetTitleOffset(1.3,"Y");
  checkname="hdistrGoodTr";

  TH1F* hdistrGoodTr=new TH1F(checkname.Data(),"Distribution of number of good tracks per event;no.good-tracks/ev;Entries",31,0,31);
  hdistrGoodTr->SetTitleOffset(1.3,"Y");

  //conta gli eventi con vertice buoni e almeno due tracce utilizzabili 
  fChecks->Add(hptGoodTr);
  fChecks->Add(hdistrGoodTr);

  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();
  cout<<nameoutput<<endl;
  fNentries=new TH1F(nameoutput, "nentriesD0->Integral(1,2) = number of AODs *** nentriesD0->Integral(3,4) = number of candidates selected with cuts *** nentriesD0->Integral(5,6) = number of D0 selected with cuts *** nentriesD0->Integral(7,8) = events with good vertex", 5,0.,5.);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandidatesSelected");
  fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtx");
  fNentries->GetXaxis()->SetBinLabel(5,"nEventsGoodVtx+>2tracks");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

 
  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(4,fChecks);
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TString bname;
  if(fArray==0){ //D0 candidates
    // load D0->Kpi candidates
    //cout<<"D0 candidates"<<endl;
    bname="D0toKpi";
  } else { //LikeSign candidates
    // load like sign candidates
    //cout<<"LS candidates"<<endl;
    bname="LikeSign2Prong";
  }

  TClonesArray *inputArray=0;
 
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }


  if(!inputArray) {
    printf("AliAnalysisTaskSED0Mass::UserExec: input branch not found!\n");
    return;
  }
  
  
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());
  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
    

  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Bool_t isGoodVtx=kFALSE;

   //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }

  //cout<<"Start checks"<<endl;
  Int_t ntracks=0,isGoodTrack=0;

  if(aod) ntracks=aod->GetNTracks();

  //cout<<"ntracks = "<<ntracks<<endl;
 
  //loop on tracks in the event
  for (Int_t k=0;k<ntracks;k++){
    AliAODTrack* track=aod->GetTrack(k);
    //cout<<"in loop"<<endl;
    //check clusters of the tracks
    Int_t nclsTot=0,nclsSPD=0;
    
    for(Int_t l=0;l<6;l++) {
      if(TESTBIT(track->GetITSClusterMap(),l)) {
	nclsTot++; if(l<2) nclsSPD++;
      }
    }
    
    if (track->Pt()>0.3 &&
	track->GetStatus()&AliESDtrack::kTPCrefit &&
	track->GetStatus()&AliESDtrack::kITSrefit &&
	nclsTot>3 &&
	nclsSPD>0) {//fill hist good tracks
      //cout<<"in if"<<endl;
      ((TH1F*)fChecks->FindObject("hptGoodTr"))->Fill(track->Pt());
      isGoodTrack++;
    }
    //cout<<"isGoodTrack = "<<isGoodTrack<<endl;
    ((TH1F*)fChecks->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
  }
  //number of events with good vertex and at least 2 good tracks
  if (isGoodTrack>=2 && isGoodVtx) fNentries->Fill(4);


  // loop over candidates
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>1) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    //Int_t nPosPairs=0, nNegPairs=0;
    //cout<<"inside the loop"<<endl;
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
  
    //check reco daughter in acceptance
    Double_t eta0=d->EtaProng(0);
    Double_t eta1=d->EtaProng(1);
   
    if (TMath::Abs(eta0) < 0.9 && TMath::Abs(eta1) < 0.9) {
      FillVarHists(d,mcArray,fCuts,fDistr);
      FillMassHists(d,mcArray,fCuts,fOutputMass);
    }
  
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  } //end for prongs


  // Post the data
  PostData(1,fOutputMass);
  PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(4,fChecks);

  cout<<"Other PostData"<<endl;
  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Mass::FillVarHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout){
  //
  // function used in UserExec to fill variable histograms:
  //

  //add distr here
  UInt_t pdgs[2];
    
  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  pdgs[0]=211;
  pdgs[1]=321;
  Double_t minvD0 = part->InvMassD0();
  pdgs[1]=211;
  pdgs[0]=321;
  Double_t minvD0bar = part->InvMassD0bar();
  //cout<<"inside mass cut"<<endl;

  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t lab=-9999;
  if(fReadMC) lab=part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
  //Double_t pt = d->Pt(); //mother pt
  Bool_t isSelected=kFALSE;

  if(fCutOnDistr){
    isSelected = cuts->IsSelected(part,AliRDHFCuts::kAll);
    if (!isSelected){
      //cout<<"Not Selected"<<endl;
      return;
    }
  }

  TString fillthispi="",fillthisK="",fillthis="";

  Int_t ptbin=cuts->PtBin(part->Pt());
  if(!fCutOnDistr || (fCutOnDistr && isSelected)){ //if no cuts or cuts passed 
    //printf("\nif no cuts or cuts passed\n");
    if(lab>=0 && fReadMC){ //signal

      //check pdg of the prongs
      AliAODTrack *prong0=(AliAODTrack*)part->GetDaughter(0);
      AliAODTrack *prong1=(AliAODTrack*)part->GetDaughter(1);
      Int_t labprong[2];
      labprong[0]=prong0->GetLabel();
      labprong[1]=prong1->GetLabel();
      AliAODMCParticle *mcprong=0;
      Int_t pdgProng[2]={0,0};
      for (Int_t iprong=0;iprong<2;iprong++){
	if(labprong[iprong]>=0)  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	pdgProng[iprong]=mcprong->GetPdgCode();
      }

      //no mass cut ditributions: ptbis
	
      fillthispi="hptpiSnoMcut_";
      fillthispi+=ptbin;

      fillthisK="hptKSnoMcut_";
      fillthisK+=ptbin;

      if (TMath::Abs(pdgProng[0]) == 211 && TMath::Abs(pdgProng[1]) == 321){
	((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(0));
	((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(1));
      }else {
	if (TMath::Abs(pdgProng[0]) == 321 && TMath::Abs(pdgProng[1]) == 211){
	  ((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(0));
	  ((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(1));
	}
      }
      
      //no mass cut ditributions: mass
      fillthis="hMassS_";
      fillthis+=ptbin;
      
      if (((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421){//D0
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      }
      else { //D0bar
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
      }
      
      //apply cut on invariant mass on the pair
      if(TMath::Abs(minvD0-mPDG)<0.03 || TMath::Abs(minvD0bar-mPDG)<0.03){
	
	if(fArray==1) cout<<"LS signal: ERROR"<<endl;
	for (Int_t iprong=0; iprong<2; iprong++){
	  AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(iprong);
	  labprong[iprong]=prong->GetLabel();
	  
	  //cout<<"prong name = "<<prong->GetName()<<" label = "<<prong->GetLabel()<<endl;
	  if(labprong[iprong]>=0)  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	  Int_t pdgprong=mcprong->GetPdgCode();
	  
	  if(TMath::Abs(pdgprong)==211) {
	    //cout<<"pi"<<endl;
	    
	    fillthispi="hptpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(part->PtProng(iprong));
	    
	    fillthisK="hd0KS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(part->Getd0Prong(iprong));
	  }
	  
	  if(TMath::Abs(pdgprong)==321) {
	    //cout<<"kappa"<<endl;
	    
	    fillthisK="hptKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(part->PtProng(iprong));
	    
	    fillthisK="hd0KS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(part->Getd0Prong(iprong));
	  }

	  fillthis="hdcaS_";
	  fillthis+=ptbin;	  
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

	  fillthis="hcosthetapointS_";
	  fillthis+=ptbin;	  
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());

	  fillthis="hcosthpointd0d0S_";
	  fillthis+=ptbin;	  
	  ((TH2F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Prodd0d0());
	  
	  fillthis="hcosthetastarS_";
	  fillthis+=ptbin;
	  if (((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0());
	  else ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0bar());

	  fillthis="hd0d0S_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());
 
	}
      }
    } else{ //Background or LS
      //cout<<"is background"<<endl;
     
      //no mass cut distributions: mass, ptbis
      fillthis="hMassB_";
      fillthis+=ptbin;
      
      if (!fCutOnDistr || (fCutOnDistr && (isSelected==1 || isSelected==3))) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
      if (!fCutOnDistr || (fCutOnDistr && isSelected>1)) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);

      fillthis="hptB1prongnoMcut_";
      fillthis+=ptbin;
      
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));
      
      fillthis="hptB2prongsnoMcut_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(1));

      //apply cut on invariant mass on the pair
      if(TMath::Abs(minvD0-mPDG)<0.03 || TMath::Abs(minvD0bar-mPDG)<0.03){


	AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(0);
	if(!prong) cout<<"No daughter found";
	else{
	  if(prong->Charge()==1) {fTotPosPairs[ptbin]++;} else {fTotNegPairs[ptbin]++;}
	}
	
	//normalise pt distr to half afterwards
	fillthis="hptB_";
	fillthis+=ptbin;

	((TH1F*)listout->FindObject(fillthis))->Fill(part->PtProng(0));((TH1F*)listout->FindObject("hptB_1"))->Fill(part->PtProng(1));

	fillthis="hd0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Getd0Prong(0));

	fillthis="hdcaB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

	fillthis="hcosthetastarB_";
	fillthis+=ptbin;
	if (!fCutOnDistr || (fCutOnDistr && (isSelected==1 || isSelected==3)))((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0());
	    if (!fCutOnDistr || (fCutOnDistr && isSelected>1))((TH1F*)listout->FindObject(fillthis))->Fill(part->CosThetaStarD0bar());	

	fillthis="hd0d0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());

	fillthis="hcosthetapointB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle());

	fillthis="hcosthpointd0d0B_";
	fillthis+=ptbin;
	((TH2F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngle(),part->Prodd0d0());
	
      }	
    }
  }
  return;
}
//____________________________________________________________________________

void AliAnalysisTaskSED0Mass::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  // function used in UserExec to fill mass histograms:
  //


  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  Int_t isSelected=cuts->IsSelected(part,AliRDHFCuts::kAll); //selected
  //cout<<"is selected = "<<isSelected<<endl;

  //cout<<"check cuts = "<<endl;
  //cuts->PrintAll();
  if (!isSelected){
    //cout<<"Not Selected"<<endl;
    return;
  }

  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  //printf("SELECTED\n");
  Int_t ptbin=cuts->PtBin(part->Pt());

  AliAODTrack *prong=(AliAODTrack*)part->GetDaughter(0);
  if(!prong) cout<<"No daughter found";
  else{
    if(prong->Charge()==1) {fTotPosPairs[ptbin]++;} else {fTotNegPairs[ptbin]++;}
  }

  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);
  //PostData(3,fNentries);

  if (isSelected==1 || isSelected==3) { //D0
    fillthis="histMass_";
    fillthis+=ptbin;
    //cout<<"Filling "<<fillthis<<endl;

    //printf("Fill mass with D0");
    ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
 
    if(labD0>=0) {
      if(fArray==1) cout<<"LS signal ERROR"<<endl;

      AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
      Int_t pdgD0 = partD0->GetPdgCode();
      //cout<<"pdg = "<<pdgD0<<endl;
      if (pdgD0==421){ //D0
	//cout<<"Fill S with D0"<<endl;
	fillthis="histSgn_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	if(TMath::Abs(invmassD0 - mPDG) < 0.027){
	  fillthis="histSgn27_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	}
      } else{ //it was a D0bar
	fillthis="histRfl_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      }
    } else {//background
      fillthis="histBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
    }
      
  }
  if (isSelected>1) { //D0bar
    fillthis="histMass_";
    fillthis+=ptbin;
    //printf("Fill mass with D0bar");
    ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);
 
    if(labD0>=0) {
      if(fArray==1) cout<<"LS signal ERROR"<<endl;
      AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
      Int_t pdgD0 = partD0->GetPdgCode();
      //cout<<" pdg = "<<pdgD0<<endl;
      if (pdgD0==-421){ //D0bar
	fillthis="histSgn_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	if (TMath::Abs(invmassD0bar - mPDG) < 0.027){
	  fillthis="histSgn27_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}

	  
      } else{
	fillthis="histRfl_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    } else {//background or LS
      fillthis="histBkg_";
      fillthis+=ptbin;
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
    }
  }

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");


  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }
  fDistr = dynamic_cast<TList*> (GetOutputData(2));
  if (!fDistr) {
    printf("ERROR: fDistr not available\n");
    return;
  }
 
//   fNentries = dynamic_cast<TH1F*>(GetOutputData(3));
//   if(!fNentries){
//     printf("ERROR: fNEntries not available\n");
//     return;
//   }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(3));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }

  fChecks = dynamic_cast<TList*> (GetOutputData(4));
  if (!fChecks) {
    printf("ERROR: fChecks not available\n");
    return;
  }
  
  if(fArray==1){
    for(Int_t ipt=0;ipt<5;ipt++){ //change 5 in GetNPtBins when sure it is written and check
      fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[ipt]*fTotNegPairs[ipt]); 


      if(fLsNormalization>1e-6) {
	
	TString massName="histMass_";
	massName+=ipt;
	((TH1F*)fOutputMass->FindObject(massName))->Scale((1/fLsNormalization)*((TH1F*)fOutputMass->FindObject(massName))->GetEntries());

      }
    

      fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[4]*fTotNegPairs[4]);

      if(fLsNormalization>1e-6) {

	TString nameDistr="hptB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hdcaB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetastarB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0d0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetapointB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthpointd0d0B_";
	nameDistr+=ipt;
	((TH2F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH2F*)fDistr->FindObject(nameDistr))->GetEntries());

      }
    }
  }
  TString cvname;

  if (fArray==0){
    cvname="D0invmass";
  } else cvname="LSinvmass";

  TCanvas *cMass=new TCanvas(cvname,cvname);
  cMass->cd();
  ((TH1F*)fOutputMass->FindObject("histMass_3"))->Draw();

  TCanvas* cStat=new TCanvas("cstat","Stat");
  cStat->cd();
  cStat->SetGridy();
  //((TH1F*)fDistr->FindObject("nEntriesD0"))->Draw("htext0");
  fNentries->Draw("htext0");

  return;
}


