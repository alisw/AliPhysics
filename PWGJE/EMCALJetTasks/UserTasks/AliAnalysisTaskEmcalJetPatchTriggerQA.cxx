// ChristineQA task
//
// Author: J Mazer

#include "AliAnalysisTaskEmcalJetPatchTriggerQA.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"
#include "AliRhoParameter.h"
#include "AliEmcalParticle.h"
#include "AliLocalRhoParameter.h"
#include "AliAnalysisTaskLocalRho.h"
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>

ClassImp(AliAnalysisTaskEmcalJetPatchTriggerQA)

//________________________________________________________________________
AliAnalysisTaskEmcalJetPatchTriggerQA::AliAnalysisTaskEmcalJetPatchTriggerQA() : 
  AliAnalysisTaskEmcalJet("ChristineQA",kFALSE), 
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0), 
  fLocalRhoVal(0),
  fHistNjetvsCent(0),
  fhnJetTriggerQA(0x0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetPatchTriggerQA::AliAnalysisTaskEmcalJetPatchTriggerQA(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0), 
  fLocalRhoVal(0),
  fHistNjetvsCent(0),
  fhnJetTriggerQA(0x0)
{ 

   SetMakeGeneralHistograms(kTRUE);
 
   // DefineInput(0,TChain::Class());
   // DefineOutput(1, TList::Class());
}

//_______________________________________________________________________
AliAnalysisTaskEmcalJetPatchTriggerQA::~AliAnalysisTaskEmcalJetPatchTriggerQA()
{
  // destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetPatchTriggerQA::UserCreateOutputObjects()
{
  if (! fCreateHisto)
    return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistNjetvsCent            = new TH2F("NjetvsCent", "NjetvsCent", 100, 0.0, 100.0, 100, 0, 100);
  fOutput->Add(fHistNjetvsCent);

  UInt_t bitcode = 0; // bit codes, see GetDimParams() below
  //bitcode = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9;
  bitcode = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 ;
  fhnJetTriggerQA = NewTHnSparseF("fhnJetTriggerQA", bitcode);
  fhnJetTriggerQA->Sumw2();
  fOutput->Add(fhnJetTriggerQA);

  PostData(1, fOutput);
}

//________________________________________________________
void AliAnalysisTaskEmcalJetPatchTriggerQA::ExecOnce()
{
  //  Initialize the analysis 
  AliAnalysisTaskEmcalJet::ExecOnce();

} // end of ExecOnce

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetPatchTriggerQA::Run()
{
  // TEST TEST TEST TEST for OBJECTS!!
  if(!fLocalRho) {
    AliWarning(Form("%s: Could not retrieve LocalRho object!", GetName()));
    fLocalRho = GetLocalRhoFromEvent(fLocalRhoName);
  }

  // check to see if we have jet object
  if (!fJets)  return kTRUE;

  // find NUMBER of jets
  const Int_t Njets = fJets->GetEntries();
  Int_t NjetAcc = 0;

  // loop over jets in the event and make appropriate cuts
  //cout<<"I at least get in the event and I have "<<Njets<<" jets"<<endl;
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
     AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
     if (!jet)  // see if we have a jet
       continue; 
     if (!AcceptMyJet(jet)) 
       continue;

     //cout<<"jet pt "<<jet->Pt()<<" area "<<jet->Area()<<" maxtrackpt "<<jet->MaxTrackPt()<<endl;
     //This somehow needs to be fixed but I'm not sure what it does yet.  It seems the defaults are wacky.
//     if (! AcceptJet(jet)) // sees if jet is accepted
//     continue;
     //  jets.push_back(jet);

     NjetAcc++;
     
     // Initializations and Calculations
     // Double_t jetphi = jet->Phi();
     //Double_t jeteta = jet->Eta();    			        // ETA of jet  
    
     Double_t jetPtraw = jet->Pt();    				// raw pT of jet
     Double_t jetPtLocal = -500;                    // initialize corr jet pt LOCAL
     Double_t jetPtGlobal = -500;					// initialize corr jet pt GLOBAL
     Double_t jetarea = -500;					// initialize jet area
     jetarea = jet->Area();		           		// jet area

     Double_t dEP = -500;			                // initialize angle between jet and event plane
     dEP = RelativeEPJET(jet->Phi(),fEPV0);

     // get LOCAL rho from event and fill histo's
     fRhoVal = fRho->GetVal();
     fLocalRhoVal = fLocalRho->GetLocalVal(jet->Phi(), 0.2);  // jet angle, then jet radius  

     jetPtLocal = jet->Pt() - jetarea*fLocalRhoVal;              // corrected pT of jet from LOCAL rho value     
     jetPtGlobal = jet->Pt() - jetarea*fRhoVal;                  // corrected pT of jet from GLOBAL rho value

     // need to update entries soon
     //Double_t Entries[10] = {centbin, jetPtLocal, jetPtGlobal, jetPtraw, jet->Eta(), jet->Phi(), dEP, fLocalRhoVal, fRhoVal, fEPV0};
     Double_t Entries[9] = {fCent, jetPtLocal, jetPtGlobal, jetPtraw, jet->Phi(), dEP, fLocalRhoVal, fRhoVal, fEPV0};
     fhnJetTriggerQA->Fill(Entries);        // fill Sparse Histo with entries

     // in plane and out of plane histo's
     if( dEP>0 && dEP<=(TMath::Pi()/6) ){
        // we are IN plane, lets fill some histo's
     }else if( dEP>(TMath::Pi()/3) && dEP<=(TMath::Pi()/2) ){
        // we are OUT of PLANE, lets fill some histo's
	 }

    fHistNjetvsCent->Fill(fCent,NjetAcc);

  } // LOOP over JETS in event
  
  return kTRUE;
}   

//________________________________________________________________________
void AliAnalysisTaskEmcalJetPatchTriggerQA::Terminate(Option_t *)
{

} // end of terminate

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetPatchTriggerQA::GetCentBin(Double_t cent) const
{  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0; 
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
} 

//_________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetPatchTriggerQA:: RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{ // function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi = (EPAng - jetAng);

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  }

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){	
    dphi = dphi + 1*TMath::Pi();
  }
  
  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 )
    AliWarning(Form("%s: dPHI not constrained from [0, Pi/2]", GetName()));

  return dphi;   // dphi in [0, Pi/2]
}

//______________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetPatchTriggerQA::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){

         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF

void AliAnalysisTaskEmcalJetPatchTriggerQA::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "V0 centrality (%)";
      nbins = 20;
      xmin = 0.;
      xmax = 100.;
      break;

   case 1:
      label = "Corrected jet p_{T}: Local #rho";
      nbins = 500;
      xmin = -250.;
      xmax = 250.;
      break;

   case 2:
      label = "Corrected jet p_{T}: Global #rho";
      nbins = 500;
      xmin = -250.;
      xmax = 250.;
      break;

   case 3:
      label = "Raw jet p_{T}";
      nbins = 250;
      xmin = 0;
      xmax = 250.;
      break;

   case 4:
      label = "Jet Phi";
      nbins = 144;
      xmin = 1.4;//-1.0*pi;
      xmax = 3.4;//2.0*pi;
      break;

  case 5:
      label = "Relative angle between jet and Reaction plane";
      nbins = 72;
      xmin = 0;
      xmax = 0.5*pi;
      break;

  case 6:
      label = "Local #rho value";
      nbins = 120;
      xmin = 0.0;
      xmax = 300.0;
      break;

  case 7: 
      label = "Global #rho value";
      nbins = 120;
      xmin = 0.0;
      xmax = 300.0;
      break;

  case 8: 
      label = "fEPV0 event plane";
      nbins = 72;
      xmin = -TMath::Pi()/2.0;
      xmax = TMath::Pi()/2.0;
      break;

   } // end of switch
} // end of getting dim-params

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetPatchTriggerQA::AcceptMyJet(AliEmcalJet *jet) {
  //applies all jet cuts except pt
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax)) return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax)) return 0;
  if (jet->Area()<fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0) return 0;
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100) return 0;

  //passed all above cuts
  return 1;
}
