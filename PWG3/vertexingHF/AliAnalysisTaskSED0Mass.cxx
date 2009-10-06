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
// and Chiara Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

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
fOutputtight(0),
fOutputloose(0),
fDistr(0), 
fVHFtight(0),
fVHFloose(0),
fNentries(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass(const char *name):
AliAnalysisTaskSE(name),
fOutputtight(0), 
fOutputloose(0), 
fDistr(0),
fVHFtight(0),
fVHFloose(0),
fNentries(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TList container
  DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TH1F container
  DefineOutput(3,TH1F::Class());  //My private output
  // Output slot #4 writes into a TList container
  DefineOutput(4,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
   if (fOutputtight) {
    delete fOutputtight;
    fOutputtight = 0;
  }
  if (fVHFtight) {
    delete fVHFtight;
    fVHFtight = 0;
  }
  if (fOutputloose) {
    delete fOutputloose;
    fOutputloose = 0;
  }

  if (fDistr) {
    delete fDistr;
    fDistr = 0;
  }

  if (fVHFloose) {
    delete fVHFloose;
    fVHFloose = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }


}  

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");

  // 2 sets of dedidcated cuts -- defined in UserExec
  fVHFtight = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHFloose = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputtight = new TList();
  fOutputtight->SetOwner();
  fOutputtight->SetName("listtight");

  fOutputloose = new TList();
  fOutputloose->SetOwner();
  fOutputloose->SetName("listloose");

  fDistr = new TList();
  fDistr->SetOwner();
  fDistr->SetName("distributionslist");

  const Int_t nhist=4;

  TString nameMass=" ", nameSgn=" ", nameBkg=" ", nameRfl=" ", namedistr=" ";

  for(Int_t i=0;i<nhist;i++){
    nameMass="histMass_";
    nameMass+=i+1;
    nameSgn="histSgn_";
    nameSgn+=i+1;
    nameBkg="histBkg_";
    nameBkg+=i+1;
    nameRfl="histRfl_";
    nameRfl+=i+1;

    //histograms of invariant mass distributions

    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpMl=(TH1F*)tmpMt->Clone();
    tmpMt->Sumw2();
    tmpMl->Sumw2();

    TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpSl=(TH1F*)tmpSt->Clone();
    tmpSt->Sumw2();
    tmpSl->Sumw2();

    TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpBl=(TH1F*)tmpBt->Clone();
    tmpBt->Sumw2();
    tmpBl->Sumw2();

    //Reflection: histo filled with D0 which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
    TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    TH1F *tmpRl=(TH1F*)tmpRt->Clone();
    tmpRt->Sumw2();
    tmpRl->Sumw2();
  //  printf("Created histograms %s\t%s\t%s\t%s\n",tmpM->GetName(),tmpS->GetName(),tmpB->GetName(),tmpR->GetName());

    fOutputtight->Add(tmpMt);
    fOutputtight->Add(tmpSt);
    fOutputtight->Add(tmpBt);
    fOutputtight->Add(tmpRt);

    fOutputloose->Add(tmpMl);
    fOutputloose->Add(tmpSl);
    fOutputloose->Add(tmpBl);
    fOutputloose->Add(tmpRl);
    
  }

  //histograms of cut variable distributions
  //  pT
  namedistr="hptpiS";
  TH1F *hptpiS = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);
//   namedistr="hptpiR";
//   TH1F *hptpiR = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

  namedistr="hptKS";
  TH1F *hptKS = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

  namedistr="hptB";
  TH1F *hptB = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);
//   namedistr="hptKR";
//   TH1F *hptKR = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

  //  dca
  namedistr="hdcaS";
  TH1F *hdcaS = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
  namedistr="hdcaB";
  TH1F *hdcaB = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
//   namedistr="hdcaR";
//   TH1F *hdcaR = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);

  //  costhetastar
  namedistr="costhetastarS";
  TH1F *hcosthetastarS = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);
  namedistr="costhetastarB";
  TH1F *hcosthetastarB = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);
//   namedistr="costhetastarR";
//   TH1F *hcosthetastarR = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

  // impact parameter
  namedistr="hd0piS";
  TH1F *hd0piS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,-0.1,0.1);
//   namedistr="hd0piR";
//   TH1F *hd0piR = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,0.,1.);

  namedistr="hd0KS";
  TH1F *hd0KS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,-0.1,0.1);
  namedistr="hd0B";
  TH1F *hd0B = new TH1F(namedistr.Data(), "Impact parameter distribution;d0 [cm]",200,-0.1,0.1);
//   namedistr="hd0KR";
//   TH1F *hd0KR = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,0.,0.1);

  namedistr="hd0d0S";
  TH1F *hd0d0S = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
  namedistr="hd0d0B";
  TH1F *hd0d0B = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
//   namedistr="hd0d0R";
//   TH1F *hd0d0R = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (pions);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

  //  costhetapoint
  namedistr="hcosthetapointS";
  TH1F *hcosthetapointS = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);
  namedistr="hcosthetapointB";
  TH1F *hcosthetapointB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);
//   namedistr="costhetapointR";
//   TH1F *hcosthetapointR = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,-1,1.);

  namedistr="hcosthpointd0d0S";
  TH2F *hcosthpointd0d0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);
  namedistr="hcosthpointd0d0B";
  TH2F *hcosthpointd0d0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

  fDistr->Add(hptpiS);
  //fDistr->Add(hptpiR);
  fDistr->Add(hptKS);
  fDistr->Add(hptB);
  //fDistr->Add(hptKR);

  fDistr->Add(hdcaS);
  fDistr->Add(hdcaB);
  //fDistr->Add(hdcaR);

  fDistr->Add(hd0piS);
  //fDistr->Add(hd0piR);
  fDistr->Add(hd0KS);
  fDistr->Add(hd0B);
  //fDistr->Add(hd0KR);

  fDistr->Add(hd0d0S);
  fDistr->Add(hd0d0B);
  //fDistr->Add(hd0d0R);

  fDistr->Add(hcosthetastarS);
  fDistr->Add(hcosthetastarB);
  //fDistr->Add(hcosthetastarR);

  fDistr->Add(hcosthetapointS);
  fDistr->Add(hcosthetapointB);
  //fDistr->Add(hcosthetapointR);

  fDistr->Add(hcosthpointd0d0S);
  fDistr->Add(hcosthpointd0d0B);

  fNentries=new TH1F("nentriesD0", "Look at the number of entries! it is = to the  number of AODs", 2,1.,2.);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
 
  // load D0->Kpi candidates                                                   
  TClonesArray *inputArrayD0toKpi =
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi");
  if(!inputArrayD0toKpi) {
    printf("AliAnalysisTaskSECompareHFpt::UserExec: D0toKpi branch not found!\n");
    return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //vtx1->Print();
  
  // load MC particles
  TClonesArray *mcArray = 
    (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!mcArray) {
    printf("AliAnalysisTaskSED0Mass::UserExec: MC particles branch not found!\n");
    return;
  }
  
  // load MC header
  AliAODMCHeader *mcHeader = 
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSED0Mass::UserExec: MC header branch not found!\n");
    return;
  }
  

  //histogram filled with 1 for every AOD
  fNentries->Fill(1);
  PostData(3,fNentries);
  // loop over D0->Kpi candidates
  Int_t nInD0toKpi = inputArrayD0toKpi->GetEntriesFast();
  printf("Number of D0->Kpi: %d\n",nInD0toKpi);
  
  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    //cout<<"inside the loop"<<endl;
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    
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

    //add distr here
    UInt_t pdgs[2];
    
    Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    pdgs[0]=211;
    pdgs[1]=321;
    Double_t minvD0 = d->InvMassD0();
    pdgs[1]=211;
    pdgs[0]=321;
    Double_t minvD0bar = d->InvMassD0bar();
    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvD0-mPDG)<0.03 || TMath::Abs(minvD0bar-mPDG)<0.03){
      Int_t pdgDgD0toKpi[2]={321,211};
      Int_t lab=d->MatchToMC(421,mcArray,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

      //      Int_t lab=d->MatchToMC(421,mcArray); //|pdg| requested
//       AliAODMCParticle *mcmother = (AliAODMCParticle*)mcArray->At(lab);
//       cout<<"mcmother name = "<<mcmother->GetName()<<" label = "<<mcmother->GetLabel()<<endl;
      if(lab>=0){ //signal
	//cout<<"is signal"<<endl;
	for (Int_t iprong=0; iprong<2; iprong++){
	  AliAODTrack *prong=(AliAODTrack*)d->GetDaughter(iprong);
	  Int_t labprong=prong->GetLabel();

	  //cout<<"prong name = "<<prong->GetName()<<" label = "<<prong->GetLabel()<<endl;
	  AliAODMCParticle *mcprong=0;
	  if(labprong>=0)  mcprong= (AliAODMCParticle*)mcArray->At(labprong);
	  Int_t pdgprong=mcprong->GetPdgCode();
	  if(TMath::Abs(pdgprong)==211) {
	    //cout<<"pi"<<endl;
	    ((TH1F*)fDistr->FindObject("hptpiS"))->Fill(d->PtProng(iprong));
	    ((TH1F*)fDistr->FindObject("hd0piS"))->Fill(d->Getd0Prong(iprong));
	  }

	  if(TMath::Abs(pdgprong)==321) {
	    //cout<<"kappa"<<endl;
	    ((TH1F*)fDistr->FindObject("hptKS"))->Fill(d->PtProng(iprong));
	    ((TH1F*)fDistr->FindObject("hd0KS"))->Fill(d->Getd0Prong(iprong));
	  }
	  ((TH1F*)fDistr->FindObject("hdcaS"))->Fill(d->GetDCA());

	}

	if (((AliAODMCParticle*)mcArray->At(lab))->GetPdgCode() == 421)
	((TH1F*)fDistr->FindObject("costhetastarS"))->Fill(d->CosThetaStarD0());
	else ((TH1F*)fDistr->FindObject("costhetastarS"))->Fill(d->CosThetaStarD0bar());

	((TH1F*)fDistr->FindObject("hd0d0S"))->Fill(d->Prodd0d0());

	((TH1F*)fDistr->FindObject("hcosthetapointS"))->Fill(d->CosPointingAngle());
	((TH1F*)fDistr->FindObject("hcosthpointd0d0S"))->Fill(d->CosPointingAngle(),d->Prodd0d0());

	//cout<<"ok point"<<endl;

      } else{ //Background
	//cout<<"is background"<<endl;
	//	AliAODTrack *prong=(AliAODTrack*)d->GetDaughter(0);
	((TH1F*)fDistr->FindObject("hptB"))->Fill(d->PtProng(0));
	//cout<<"ptok"<<endl;
	((TH1F*)fDistr->FindObject("hd0B"))->Fill(d->Getd0Prong(0));
	//cout<<"d0ok"<<endl;
	((TH1F*)fDistr->FindObject("hdcaB"))->Fill(d->GetDCA());
	//cout<<"dcaok"<<endl;
	((TH1F*)fDistr->FindObject("costhetastarB"))->Fill(d->CosThetaStarD0());
	((TH1F*)fDistr->FindObject("costhetastarB"))->Fill(d->CosThetaStarD0bar());	
	((TH1F*)fDistr->FindObject("hd0d0B"))->Fill(d->Prodd0d0());
	//cout<<"d0d0ok"<<endl;
	((TH1F*)fDistr->FindObject("hcosthetapointB"))->Fill(d->CosPointingAngle());
	((TH1F*)fDistr->FindObject("hcosthpointd0d0B"))->Fill(d->CosPointingAngle(),d->Prodd0d0());

	//cout<<"pointok"<<endl;

      }

    } //inv mass cut
  

    Double_t pt = d->Pt();

    Int_t ptbin=0;
    
    //cout<<"P_t = "<<pt<<endl;
    if (pt>0. && pt<=1.) {
      ptbin=1;
      fVHFtight->SetD0toKpiCuts(0.7,0.04,0.8,0.5,0.5,0.05,0.05,-0.0003,0.7);
      fVHFloose->SetD0toKpiCuts(0.7,0.04,0.8,0.5,0.5,0.05,0.05,-0.00025,0.7);
//       fVHFtight->SetD0toKpiCuts(0.7,0.04,0.8,0.5,0.5,0.05,0.05,-0.0002,0.7);
//       fVHFloose->SetD0toKpiCuts(0.7,0.04,0.8,0.5,0.5,1,1,-0.00015,0.5);
      //printf("I'm in the bin %d\n",ptbin);
      
    }
    
    if(pt>1. && pt<=3.) {
      ptbin=2;  
      fVHFtight->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,0.05,0.05,-0.0003,0.9);
      fVHFloose->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
      //printf("I'm in the bin %d\n",ptbin);
    }
 
    if(pt>3. && pt<=5.){
	ptbin=3;  
	fVHFtight->SetD0toKpiCuts(0.7,0.015,0.8,0.7,0.7,0.05,0.05,-0.0002,0.9);
	fVHFloose->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.8);
	//printf("I'm in the bin %d\n",ptbin);
    }
    if(pt>5.){
      ptbin=4;
      fVHFtight->SetD0toKpiCuts(0.7,0.015,0.8,0.7,0.7,0.05,0.05,-0.0002,0.95);
      fVHFloose->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
    }//if(pt>5)
  
    //printf("I'm in the bin %d\n",ptbin);
    //old
    //fVHF->SetD0toKpiCuts(0.7,0.03,0.8,0.06,0.06,0.05,0.05,-0.0002,0.6); //2.p-p vertex reconstructed    

    FillHists(ptbin,d,mcArray,fVHFtight,fOutputtight);
    FillHists(ptbin,d,mcArray,fVHFloose,fOutputloose);
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  
     
  
   
  // Post the data
  PostData(1,fOutputtight);
  PostData(2,fOutputloose);
  PostData(4,fDistr);
  return;
}
//____________________________________________________________________________*
void AliAnalysisTaskSED0Mass::FillHists(Int_t ptbin, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAnalysisVertexingHF *vhf, TList *listout){
  //
  // function used in UserExec:
  //
  Int_t okD0=0,okD0bar=0;

  if(part->SelectD0(vhf->GetD0toKpiCuts(),okD0,okD0bar)) {//selected
    Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
    //printf("SELECTED\n");
    Int_t pdgDgD0toKpi[2]={321,211};
    Int_t labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

    //Int_t labD0 = part->MatchToMC(421,arrMC); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
    //printf("labD0 %d",labD0);
    
    TString fillthis="";
    
    if (okD0==1) {
      fillthis="histMass_";
      fillthis+=ptbin;
      //cout<<"Filling "<<fillthis<<endl;

      //printf("Fill mass with D0");
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//cout<<"pdg = "<<pdgD0<<endl;
	if (pdgD0==421){ //D0
	  //cout<<"Fill S with D0"<<endl;
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);

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
    if (okD0bar==1) {
      fillthis="histMass_";
      fillthis+=ptbin;
      //printf("Fill mass with D0bar");
      ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);
      
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//cout<<" pdg = "<<pdgD0<<endl;
	if (pdgD0==-421){ //D0bar
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	  
	} else{
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}
      } else {//background
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    }


  } //else cout<<"NOT SELECTED"<<endl;

}
//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");

  fOutputtight = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputtight) {     
    printf("ERROR: fOutputthight not available\n");
    return;
  }
  fOutputloose = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputloose) {     
    printf("ERROR: fOutputloose not available\n");
    return;
  }
  fDistr = dynamic_cast<TList*> (GetOutputData(3));
  if (!fDistr) {
    printf("ERROR: fDistr not available\n");
    return;
  }


  return;
}

