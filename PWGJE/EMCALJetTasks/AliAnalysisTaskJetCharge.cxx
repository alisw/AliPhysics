//
// Basic analysis task.
//
// Basic analysis task template for analysis jets storing information in both tree
// branches and histograms

#include <TH1F.h>
#include <TTree.h>
#include <TList.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include <TProfile.h>
#include <TChain.h>

// aliroot Headers
#include "AliAODEvent.h"
#include "AliParticleContainer.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliEmcalJet.h"

//My Header
#include "AliAnalysisTaskJetCharge.h"




ClassImp(AliAnalysisTaskJetCharge)

//________________________________________________________________________
AliAnalysisTaskJetCharge::AliAnalysisTaskJetCharge() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetCharge", kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.2),
  JetChargeK(0.5),
  JetMidPt(40),
  JetHighPt(80),

  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  fhJetCharge(0x0),
  fhJetChargeLow(0x0),
  fhJetChargeMid(0x0),
  fhJetChargeHigh(0x0),


  fTreeJets(0)
  {

  }



//________________________________________________________________________
AliAnalysisTaskJetCharge::AliAnalysisTaskJetCharge(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.2),
  JetChargeK(0.5),
  JetMidPt(40),
  JetHighPt(80),

  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  fhJetCharge(0x0),
  fhJetChargeLow(0x0),
  fhJetChargeMid(0x0),
  fhJetChargeHigh(0x0),



  fTreeJets(0)
{
  // Standard constructor.
  for(Int_t i=0;i<nBranchesJetCharge;i++){
    fTreeBranch[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetCharge::~AliAnalysisTaskJetCharge()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskJetCharge::UserCreateOutputObjects()
{
  // Echo jet radius
  Info("TaskJets","Using jet radius R=%f",fJetRadius);

  //fEventCuts.SetManualMode(); ///Enable manual mode
  //fEventCuts.SetupRun2pp();

  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create output TTree
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJets = new TTree(nameoutput, nameoutput);
  // Names for the branches
  TString *fTreeBranchName = new TString [nBranchesJetCharge];

  // Name the branches of your TTree here
  fTreeBranchName[0]  = "Pt";
  fTreeBranchName[1]  = "Phi";
  fTreeBranchName[2]  = "Eta";

  fTreeBranchName[3]  = "JetCharge";

  fTreeBranchName[4]  = "JetChargeLow";
  fTreeBranchName[5]  = "JetChargeMid";
  fTreeBranchName[6]  = "JetChargeHigh";




  // Associate the branches
  for(Int_t iBranch=0; iBranch < nBranchesJetCharge; iBranch++){
    std::cout<<"looping over variables"<<std::endl;
    fTreeJets->Branch(fTreeBranchName[iBranch].Data(), &fTreeBranch[iBranch], Form("%s/D", fTreeBranchName[iBranch].Data()));
  }

  // Define histograms
  fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );
  fOutput->Add(fhJetPt);
  fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
  fOutput->Add(fhJetPhi);
  fhJetEta= new TH1F("fhJetEta", "Jet Eta",100,-2,2);
  fOutput->Add(fhJetEta);

  fhJetCharge= new TH1F("fhJetCharge", "Jet Charge", 25, -3, 3);
  fOutput->Add(fhJetCharge);

  fhJetChargeLow= new TH1F("fhJetChargeLow", "Jet Charge", 25, -3, 3);
  fOutput->Add(fhJetChargeLow);

  fhJetChargeMid= new TH1F("fhJetChargeMid", "Jet Charge", 25, -3, 3);
  fOutput->Add(fhJetChargeMid);

  fhJetChargeHigh= new TH1F("fhJetChargeHigh", "Jet Charge", 25, -3, 3);
  fOutput->Add(fhJetChargeHigh);


  // Make sure that the outputs get written out
  PostData(1,fOutput);
  PostData(2,fTreeJets);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetCharge::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().



  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetCharge::FillHistograms()
{
  // Centrality selection, if enabled
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin))
      return 0;
  }


  AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
  Double_t JetPhi=0;
  Double_t JetPt_ForThreshold=0;
  if(JetCont) {
    // Technical detail; fix possibly corrupted jet container ID
    JetCont->ResetCurrentID();
    // Jet is acceptable?
    while((Jet1=JetCont->GetNextAcceptJet())) {
      if(!Jet1)
        continue;


      // Get the jet constituents



      AliParticleContainer *fTrackCont = JetCont->GetParticleContainer();
      UInt_t nJetConstituents = Jet1->GetNumberOfTracks();

      //cout << nTest << "::::" << nMCConstituents << "::::" << nJetConstituents << endl;



      // Must have at least two constituents
      if( nJetConstituents < 2 )
      {
        continue;
      }
      // Check if Pt is above the Momentum
        JetPt_ForThreshold = Jet1->Pt();

      if(JetPt_ForThreshold<fPtThreshold)
      {
        continue;
      }
        //Check for leading track Pt to reduce cominatorial jets.

        if(Jet1->GetLeadingTrack()->Pt() < 5)
        {

            continue;
            //std::cout << "LEADING TRACK TO SMALL!!!" << std::endl;
        }

        // Filling the histograms here
        Double_t JetPt = 0;
          JetPt = Jet1->Pt();
          fTreeBranch[0]= Jet1->Pt();
          fhJetPt->Fill(JetPt);

          JetPhi = Jet1->Phi();
          if(JetPhi < -1*TMath::Pi())
            JetPhi += (2*TMath::Pi());
          else if (JetPhi > TMath::Pi())
            JetPhi -= (2*TMath::Pi());

          fTreeBranch[1]=JetPhi;
          fhJetPhi->Fill(JetPhi);
          fTreeBranch[2]=Jet1->Eta();
          fhJetEta->Fill(Jet1->Eta());

          Double_t jetCharge = 0;

          // Loop over the consituents
          for (UInt_t iJetConst = 0; iJetConst < nJetConstituents; iJetConst++ )
          {
            AliVParticle *fJetConst = Jet1->Track(iJetConst);

            //Looping over the Mc Particles to find the matching particle.

            jetCharge += fJetConst->Charge()*pow(fJetConst->Pt(),JetChargeK);

          }

          jetCharge/=pow(Jet1->Pt(),JetChargeK);

          fTreeBranch[3] = jetCharge;
          fhJetCharge->Fill(jetCharge);

          if(JetPt < JetMidPt)
          {
            fTreeBranch[4] = jetCharge;
            fhJetChargeLow->Fill(jetCharge);
          }
          else if( JetPt > JetMidPt && JetPt < JetHighPt)
          {
            fTreeBranch[5] = jetCharge;
            fhJetChargeMid->Fill(jetCharge);
          }
          else
          {
            fTreeBranch[6] = jetCharge;
            fhJetChargeHigh->Fill(jetCharge);
          }





          fTreeJets->Fill();



        //cout << "End of Jet" << endl;

    }
  }
  return kTRUE;
}



//________________________________________________________________________
Bool_t AliAnalysisTaskJetCharge::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetCharge::Terminate(Option_t *)
{
  // Called once at the end of the analysis.

  // Normalise historgrams over number of Jets considered
  //fhJetCharge->Scale(1/fhJetCharge->GetEntries()*fhJetCharge->GetXaxis()->GetBinWidth(1));

}
