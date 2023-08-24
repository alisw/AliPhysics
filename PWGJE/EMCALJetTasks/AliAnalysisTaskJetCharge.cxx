//
// Basic analysis task.
//
// Basic analysis task template for analysis jets storing information in both tree
// branches and histograms


#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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

//Globals
using std::cout;
using std::endl;



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
  fJetRadius(0),
  JetChargeK(0.5),
  JetChargeMid(-999),
  JetChargeHigh(-999),


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
  fJetRadius(0),
  JetChargeK(0.5),
  JetChargeMid(-999),
  JetChargeHigh(-999),





  fTreeJets(0)
{

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
  //Info("TaskJets","Using jet radius R=%f",fJetRadius);

  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create output TTree
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJets = new TTree(nameoutput, nameoutput);
  // Names for the branches

  // Name the branches of your TTree here

  fTreeJets->Branch("Pt",&Pt,"Pt/F");
  fTreeJets->Branch("Phi",&Phi,"Phi/F");
  fTreeJets->Branch("Eta",&Eta,"Eta/F");
  fTreeJets->Branch("JetCharge",&JetCharge,"JetCharge/F");
  fTreeJets->Branch("JetChargeMid",&JetChargeMid,"JetChargeMid/F");
  fTreeJets->Branch("JetChargeHigh",&JetChargeHigh,"JetChargeHigh/F");

  fTreeJets->Branch("LeadingTrackPt",&LeadingTrackPt,"LeadingTrackPt/F");




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

  // Reset the Tree Parameters
  Pt = 0;
  Phi = 0;
  Eta = 0;
  JetCharge = 0;
  JetChargeMid = 0;
  JetChargeHigh = 0;
  LeadingTrackPt = 0;

  // Initialise jet pointer
  //cout << "Running Fill Histograms" << endl;
  AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
  Int_t nAcceptedJets = JetCont->GetNAcceptedJets();

  //TClonesArray *trackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("HybridTracks"));
  Float_t JetPhi=0;
  Float_t JetPt_ForThreshold=0;

  if(JetCont) {
    // Technical detail; fix possibly corrupted jet container ID
    // Jet is acceptable?
    for (auto Jet1 : JetCont->accepted())
    {

      if(!Jet1)
      {

        continue;
      }


      AliParticleContainer *fTrackCont = JetCont->GetParticleContainer();
      UInt_t nJetConstituents = Jet1->GetNumberOfTracks();



      //cout << nTest << "::::" << nMCConstituents << "::::" << nJetConstituents << endl;


      // Must have at least two constituents
      if( nJetConstituents < 2 )
      {
        continue;
      }


      JetPt_ForThreshold = Jet1->Pt();

      if(JetPt_ForThreshold<fPtThreshold)
      {

        continue;
      }
      else
      {

        //Check the jet leading track is Greater than 5 GeV to reduce cominatorial jets

        if(Jet1->GetLeadingTrack()->Pt() < 5)
        {
            continue;
        }

        LeadingTrackPt = Jet1->GetLeadingTrack()->Pt();


          // Filling the Tree here

          JetPhi=Jet1->Phi();
          if(JetPhi < -1*TMath::Pi())
            JetPhi += (2*TMath::Pi());
          else if (JetPhi > TMath::Pi())
            JetPhi -= (2*TMath::Pi());

          // Filling the TTree branches here
          Pt          = Jet1->Pt();
          Phi         = JetPhi;
          Eta         = Jet1->Eta();


          // Now Caluclate the jet Charge
          Float_t jetCharge = 0;
          Float_t MidjetCharge = 0;
          Float_t HighjetCharge = 0;

          double_t DetMidPt = 0;
          double_t DetHighPt = 0;
          Int_t MidCount = 0;
          Int_t HighCount = 0;

          // Loop over the consituents
          for (UInt_t iJetConst = 0; iJetConst < nJetConstituents; iJetConst++ )
          {
            AliVParticle *JetParticle = Jet1->Track(iJetConst);
            jetCharge += JetParticle->Charge()*pow(JetParticle->Pt(),JetChargeK);
            if(JetParticle->Pt() > 0.8)
            {
              MidjetCharge += JetParticle->Charge()*pow(abs(JetParticle->Pt()),JetChargeK);
              DetMidPt += JetParticle->Pt();
              MidCount += 1;
            }
            if(JetParticle->Pt() > 3.0)
            {
              HighjetCharge += JetParticle->Charge()*pow(abs(JetParticle->Pt()),JetChargeK);
              DetHighPt += JetParticle->Pt();
              HighCount += 1;
            }
          }




        // Normalise the Non Flavoured Jet CHarge
        jetCharge/=pow(Jet1->Pt(),0.5);
        MidjetCharge/=pow(DetMidPt,JetChargeK);
        HighjetCharge/=pow(DetHighPt,JetChargeK);

        //insure that the particle count is greater than 1 otherwise sets the final value to 0
        if(MidCount < 2)
        {
          MidjetCharge = 0;
        }
        if(HighCount < 2)
        {
          HighjetCharge = 0;
        }     

        //Put The Jet Charge in the right place
        JetCharge = jetCharge;
        JetChargeMid = MidjetCharge;
        JetChargeHigh = HighjetCharge;



        fTreeJets->Fill();

      }

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

}
