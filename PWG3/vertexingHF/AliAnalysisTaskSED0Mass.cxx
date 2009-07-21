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
// and comparison with the MC truth.
//
// Authors: A.Dainese, andrea.dainese@lnl.infn.it
// and Chiara Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>


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
//fNtupleD0Cmp(0),
/*
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
*/
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
//fNtupleD0Cmp(0),
/*
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
*/
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
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
  // Destructor
  /*
  if (fhistMass){
    delete fhistMass;
    fhistMass=0;
  }

  if (fhistSgn){
    delete fhistSgn;
    fhistSgn=0;
  }

  if (fhistBkg){
    delete fhistBkg;
    fhistBkg=0;
  }
  */
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
  if (fVHFloose) {
    delete fVHFloose;
    fVHFloose = 0;
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

  const Int_t nhist=4;

  TString nameMass=" ", nameSgn=" ", nameBkg=" ", nameRfl=" ";

  for(Int_t i=0;i<nhist;i++){
    nameMass="histMass_";
    nameMass+=i+1;
    nameSgn="histSgn_";
    nameSgn+=i+1;
    nameBkg="histBkg_";
    nameBkg+=i+1;
    nameRfl="histRfl_";
    nameRfl+=i+1;

    TH1F* tmpM = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
    tmpM->Sumw2();

    TH1F* tmpS = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    tmpS->Sumw2();

    TH1F* tmpB = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    tmpB->Sumw2();

    //Reflection: histo filled with D0 which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
    TH1F* tmpR = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    tmpR->Sumw2();

  /*
    tmpM->SetName(nameMass);
    tmpS->SetName(nameSgn);
    tmpB->SetName(nameBkg);
    tmpR->SetName(nameRfl);
  */
  //  printf("Created histograms %s\t%s\t%s\t%s\n",tmpM->GetName(),tmpS->GetName(),tmpB->GetName(),tmpR->GetName());

    fOutputtight->Add(tmpM);
    fOutputtight->Add(tmpS);
    fOutputtight->Add(tmpB);
    fOutputtight->Add(tmpR);

    fOutputloose->Add(tmpM);
    fOutputloose->Add(tmpS);
    fOutputloose->Add(tmpB);
    fOutputloose->Add(tmpR);
    
//   delete tmpM;
//   delete tmpS;
//   delete tmpB;
//   delete tmpR;

  }


  //fOutputloose->ls();
  //fOutputtight->ls();

  fNentries=new TH1F("nentries", "Look at the number of entries! = number of AODs", 2,1.,2.);

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
    
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    
    //  cout<<"Cuts applied: ";
    //     Double_t *cutsapplied=(Double_t*)fVHF->GetD0toKpiCuts();
    //     for (Int_t i=0;i<9;i++){
    //       printf(cutsapplied[i]\t);
    //     }
    //    cout<<endl;

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
    else {
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
	
      } else{
	ptbin=4;
	fVHFtight->SetD0toKpiCuts(0.7,0.015,0.8,0.7,0.7,0.05,0.05,-0.0002,0.95);
	fVHFloose->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,0.05,0.05,-0.00015,0.9);
      }//if(pt>5)
    }
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
    Int_t labD0 = part->MatchToMC(421,arrMC); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
    //printf("labD0 %d",labD0);
    
    TString fillthis="";
    TString fillrfl= "";

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
	fillthis="histSgn_";
	fillrfl="histRfl_";
	fillthis+=ptbin;
	fillrfl+=ptbin;
	if (pdgD0==421){ //D0
	  //cout<<"Fill S with D0"<<endl;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);

	} else{ //it was a D0bar
	  ((TH1F*)(listout->FindObject(fillrfl)))->Fill(invmassD0);
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
	  ((TH1F*)(listout->FindObject(fillrfl)))->Fill(invmassD0bar);
	  
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
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutputloose = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputloose) {     
    printf("ERROR: fOutput not available\n");
    return;
  }


  return;
}

