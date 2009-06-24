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
fOutput(0), 
//fNtupleD0Cmp(0),
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass(const char *name):
AliAnalysisTaskSE(name),
fOutput(0), 
//fNtupleD0Cmp(0),
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
  // Destructor

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

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }

}  

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  // set dedidcated cuts
  fVHF->SetD0toKpiCuts(0.7,0.03,0.8,0.06,0.06,0.05,0.05,-0.0002,0.6); //2.p-p vertex reconstructed
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

  const Int_t nhist=3;
  fhistMass=new TClonesArray("TH1F",nhist);
  fhistMass->SetName("fhistMass");
  fhistSgn =new TClonesArray("TH1F",nhist);
  fhistSgn->SetName("fhistSgn");
  fhistBkg =new TClonesArray("TH1F",nhist);
  fhistBkg->SetName("fhistBkg");

  TString nameMass, nameSgn, nameBkg;

  for(Int_t i=0;i<nhist;i++){
    nameMass="fhistMass_";
    nameMass+=i+1;
    cout<<"Creating TH1F with name "<<nameMass<<endl;
    new ((*fhistMass)[i]) TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
    nameSgn="fhistSgn_";
    nameSgn+=i+1;
    new ((*fhistSgn)[i]) TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    nameBkg="fhistBkg_";
    nameBkg+=i+1;
    new ((*fhistBkg)[i]) TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.765,1.965);
    printf("Created histograms %s\t%s\t%s\n",fhistMass->UncheckedAt(i)->GetName(),fhistSgn->UncheckedAt(i)->GetName(),fhistBkg->UncheckedAt(i)->GetName());
  }
 
  fOutput->Add(fhistMass);
  fOutput->Add(fhistSgn);
  fOutput->Add(fhistBkg);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
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
    
    
    Double_t pt = d->Pt();
    Int_t ptbin=0;
    
    //cout<<"P_t = "<<pt<<endl;
    if (pt>0. && pt<=1.) {
      ptbin=1; 
      //printf("I'm in the bin %d\n",ptbin);
    }
    else {
      if(pt>1. && pt<=3.) {
	ptbin=2;  
	//printf("I'm in the bin %d\n",ptbin);
      }
      else {
	ptbin=3;  
	//printf("I'm in the bin %d\n",ptbin);
      }//if(pt>3)
    }
    
    FillHists(ptbin,d,mcArray);
   
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
 
  }
     
  
   
  // Post the data
  PostData(1,fOutput);


  return;
}
//____________________________________________________________________________*
void AliAnalysisTaskSED0Mass::FillHists(Int_t ptbin, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC){
  //
  // function used in UserExec:
  //
  Int_t okD0=0,okD0bar=0;

  if(part->SelectD0(fVHF->GetD0toKpiCuts(),okD0,okD0bar)) {//selected
    Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
    //printf("SELECTED\n");
    Int_t labD0 = part->MatchToMC(421,arrMC); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
    //printf("labD0 %d",labD0);

    if(labD0>=0) {
      
      AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);

      Int_t pdgD0 = partD0->GetPdgCode();
      
      //printf(" pdgD0 %d\n",pdgD0);
      
      if (pdgD0==421){ //D0
	//cout<<"Fill S with D0"<<endl;
	((TH1F*)(fhistSgn->UncheckedAt(ptbin-1)))->Fill(invmassD0);
      }
      else {//D0bar  
	//printf("Fill S with D0bar");
	((TH1F*)(fhistSgn->UncheckedAt(ptbin-1)))->Fill(invmassD0bar);
      }

    }
    else {//background
	  //printf("Fill background");
      ((TH1F*)(fhistBkg->UncheckedAt(ptbin-1)))->Fill(invmassD0);
      ((TH1F*)(fhistBkg->UncheckedAt(ptbin-1)))->Fill(invmassD0bar);
    }
    
    //no MC info, just cut selection
    if (okD0==1) {
      //printf("Fill mass with D0");
      ((TH1F*)(fhistMass->UncheckedAt(ptbin-1)))->Fill(invmassD0);
      
    }
    if (okD0bar==1) {
      ((TH1F*)fhistMass->UncheckedAt(ptbin-1))->Fill(invmassD0bar);
      //printf("Fill mass with D0bar");
    }
      
  } //else cout<<"NOT SELECTED"<<endl;

}
//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  
  fhistMass = dynamic_cast<TClonesArray*>(fOutput->FindObject("fhistMass"));
  fhistSgn = dynamic_cast<TClonesArray*>(fOutput->FindObject("fhistSgn"));
  fhistBkg = dynamic_cast<TClonesArray*>(fOutput->FindObject("fhistBkg"));
 
  return;
}

