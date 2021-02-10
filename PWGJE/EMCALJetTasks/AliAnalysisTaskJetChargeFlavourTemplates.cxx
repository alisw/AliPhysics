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
#include "AliAnalysisTaskJetChargeFlavourTemplates.h"

//Globals
using std::cout;
using std::endl;



ClassImp(AliAnalysisTaskJetChargeFlavourTemplates)

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::AliAnalysisTaskJetChargeFlavourTemplates() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetChargeFlavourTemplates", kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.2),


  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  JC(0x0),

  JCUp(0x0),
  JCDown(0x0),
  JCGluon(0x0),
  JCOther(0x0),
  JCUnmatched(0x0),

  JCLow(0x0),

  JCUpLow(0x0),
  JCDownLow(0x0),
  JCGluonLow(0x0),
  JCOtherLow(0x0),
  JCUnmatchedLow(0x0),


  JCMid(0x0),

  JCUpMid(0x0),
  JCDownMid(0x0),
  JCGluonMid(0x0),
  JCOtherMid(0x0),
  JCUnmatchedMid(0x0),

  JCHigh(0x0),

  JCUpHigh(0x0),
  JCDownHigh(0x0),
  JCGluonHigh(0x0),
  JCOtherHigh(0x0),
  JCUnmatchedHigh(0x0),



  fTreeJets(0)
{

}

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::AliAnalysisTaskJetChargeFlavourTemplates(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fContainer(0),
  pChain(0),
  fPtThreshold(-9999.),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fJetRadius(0.2),


  fhJetPt(0x0),
  fhJetPhi(0x0),
  fhJetEta(0x0),

  JC(0x0),

  JCUp(0x0),
  JCDown(0x0),
  JCGluon(0x0),
  JCOther(0x0),
  JCUnmatched(0x0),


  JCLow(0x0),

  JCUpLow(0x0),
  JCDownLow(0x0),
  JCGluonLow(0x0),
  JCOtherLow(0x0),
  JCUnmatchedLow(0x0),


  JCMid(0x0),

  JCUpMid(0x0),
  JCDownMid(0x0),
  JCGluonMid(0x0),
  JCOtherMid(0x0),
  JCUnmatchedMid(0x0),

  JCHigh(0x0),

  JCUpHigh(0x0),
  JCDownHigh(0x0),
  JCGluonHigh(0x0),
  JCOtherHigh(0x0),
  JCUnmatchedHigh(0x0),



  fTreeJets(0)
{
  // Standard constructor.
  for(Int_t i=0;i<nBranchesJetChargeFlavourTemplates;i++){
    fTreeBranch[i]=0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetChargeFlavourTemplates::~AliAnalysisTaskJetChargeFlavourTemplates()
{
  // Destructor.
}

//________________________________________________________________________
 void AliAnalysisTaskJetChargeFlavourTemplates::UserCreateOutputObjects()
{
  // Echo jet radius
  Info("TaskJets","Using jet radius R=%f",fJetRadius);

  // Create user output.
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  //create output TTree
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJets = new TTree(nameoutput, nameoutput);
  // Names for the branches
  TString *fTreeBranchName = new TString [nBranchesJetChargeFlavourTemplates];

  // Name the branches of your TTree here
  fTreeBranchName[0]  = "Pt";
  fTreeBranchName[1]  = "Phi";
  fTreeBranchName[2]  = "Eta";

  fTreeBranchName[3]  = "JetCharge";

  fTreeBranchName[4]  = "LowJetCharge";

  fTreeBranchName[5]  = "MidJetCharge";

  fTreeBranchName[6]  = "HighJetCharge";

  fTreeBranchName[7] = "JCUp";
  fTreeBranchName[8] = "Low_JCUp";
  fTreeBranchName[9] = "Mid_JCUp";
  fTreeBranchName[10] = "High_JCUp";

  fTreeBranchName[11] = "JCDown";
  fTreeBranchName[12] = "Low_JCDown";
  fTreeBranchName[13] = "Mid_JCDown";
  fTreeBranchName[14] = "High_JCDown";


  fTreeBranchName[15] = "JCGluon";
  fTreeBranchName[16] = "Low_JCGluon";
  fTreeBranchName[17] = "Mid_JCGluon";
  fTreeBranchName[18] = "High_JCGluon";

  fTreeBranchName[19] = "JCOther";
  fTreeBranchName[20] = "Low_JCOther";
  fTreeBranchName[21] = "Mid_JCOther";
  fTreeBranchName[22] = "High_JCOther";

  fTreeBranchName[23] = "JCUnmatched";
  fTreeBranchName[24] = "Low_JCUnmatched";
  fTreeBranchName[25] = "Mid_JCUnmatched";
  fTreeBranchName[26] = "High_JCUnmatched";

  // Associate the branches
  for(Int_t iBranch=0; iBranch < nBranchesJetChargeFlavourTemplates; iBranch++){
    cout<<"looping over variables"<<endl;
    fTreeJets->Branch(fTreeBranchName[iBranch].Data(), &fTreeBranch[iBranch], Form("%s/D", fTreeBranchName[iBranch].Data()));
  }

  // Define histograms
  fhJetPt= new TH1F("fhJetPt", "Jet Pt",1500,-0.5,149.5 );
  fOutput->Add(fhJetPt);
  fhJetPhi= new TH1F("fhJetPhi", "Jet Phi",360 , -1.5*(TMath::Pi()), 1.5*(TMath::Pi()));
  fOutput->Add(fhJetPhi);
  fhJetEta= new TH1F("fhJetEta", "Jet Eta",100,-2,2);
  fOutput->Add(fhJetEta);

  /*
  fhEventCounter= new TH1F("fhEventCounter", "Event Counter",10,10,10);
  fOutput->Add(fhEventCounter);

  fhRunNumberCounter= new TH1F("fhRunNumberCounter", "Event Counter",10,10,10);
  fOutput->Add(fhRunNumberCounter);
  */

  JC= new TH1F("JC", "Jet Charge", 25, -3, 3);
  fOutput->Add(JC);

  JCUp= new TH1F("JCUp", "Jet Charge Up", 25, -3, 3);
  fOutput->Add(JCUp);
  JCDown= new TH1F("JCDown", "Jet Charge Down", 25, -3, 3);
  fOutput->Add(JCDown);
  JCGluon= new TH1F("JCGluon", "Jet Charge Gluon", 25, -3, 3);
  fOutput->Add(JCGluon);
  JCOther= new TH1F("JCOther", "Jet Charge Other", 25, -3, 3);
  fOutput->Add(JCOther);
  JCUnmatched= new TH1F("JCUnmatched", "Jet Charge Unmatched", 25, -3, 3);
  fOutput->Add(JCUnmatched);




  JCLow= new TH1F("JCLow", "Jet Charge Low Pt ", 25, -3, 3);
  fOutput->Add(JCLow);

  JCUpLow= new TH1F("JCUpLow", "Jet Charge Up Low Pt ", 25, -3, 3);
  fOutput->Add(JCUpLow);
  JCDownLow= new TH1F("JCDownLow", "Jet Charge Down Low Pt", 25, -3, 3);
  fOutput->Add(JCDownLow);
  JCGluonLow= new TH1F("JCGluonLow", "Jet Charge Gluon Low Pt", 25, -3, 3);
  fOutput->Add(JCGluonLow);
  JCOtherLow= new TH1F("JCOtherLow", "Jet Charge Other Low Pt", 25, -3, 3);
  fOutput->Add(JCOtherLow);
  JCUnmatchedLow= new TH1F("JCUnmatchedLow", "Jet Charge Unmatched Low Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedLow);



  JCMid= new TH1F("JCMid", "Jet Charge Mid Pt ", 25, -3, 3);
  fOutput->Add(JCMid);

  JCUpMid= new TH1F("JCUpMid", "Jet Charge Up Mid Pt ", 25, -3, 3);
  fOutput->Add(JCUpMid);
  JCDownMid= new TH1F("JCDownMid", "Jet Charge Down Mid Pt", 25, -3, 3);
  fOutput->Add(JCDownMid);
  JCGluonMid= new TH1F("JCGluonMid", "Jet Charge Gluon Mid Pt", 25, -3, 3);
  fOutput->Add(JCGluonMid);
  JCOtherMid= new TH1F("JCOtherMid", "Jet Charge Other Mid Pt", 25, -3, 3);
  fOutput->Add(JCOtherMid);
  JCUnmatchedMid= new TH1F("JCUnmatchedMid", "Jet Charge Unmatched Mid Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedMid);

  JCHigh= new TH1F("JCHigh", "Jet Charge High Pt ", 25, -3, 3);
  fOutput->Add(JCHigh);

  JCUpHigh= new TH1F("JCUpHigh", "Jet Charge Up High Pt ", 25, -3, 3);
  fOutput->Add(JCUpHigh);
  JCDownHigh= new TH1F("JCDownHigh", "Jet Charge Down High Pt", 25, -3, 3);
  fOutput->Add(JCDownHigh);
  JCGluonHigh= new TH1F("JCGluonHigh", "Jet Charge Gluon High Pt", 25, -3, 3);
  fOutput->Add(JCGluonHigh);
  JCOtherHigh= new TH1F("JCOtherHigh", "Jet Charge Other High Pt", 25, -3, 3);
  fOutput->Add(JCOtherHigh);
  JCUnmatchedHigh= new TH1F("JCUnmatchedHigh", "Jet Charge Unmatched High Pt", 25, -3, 3);
  fOutput->Add(JCUnmatchedHigh);



  // Make sure that the outputs get written out
  PostData(1,fOutput);
  PostData(2,fTreeJets);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().




  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::FillHistograms()
{
  // Centrality selection, if enabled
  if (fCentSelectOn){
    if ((fCent>fCentMax) || (fCent<fCentMin))
      return 0;
  }

  // Initialise jet pointer
  //cout << "Running Fill Histograms" << endl;
  AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
  AliParticleContainer *MCParticleContainer = GetParticleContainer("mcparticles");
  Int_t nAcceptedJets = JetCont->GetNAcceptedJets();
  //TClonesArray *trackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("HybridTracks"));
  Double_t JetPhi=0;
  Double_t JetPt_ForThreshold=0;

  if(JetCont) {
    // Technical detail; fix possibly corrupted jet container ID
    JetCont->ResetCurrentID();
    // Jet is acceptable?
    while((Jet1=JetCont->GetNextAcceptJet())) {
      if(!Jet1)
      {

        continue;
      }

      // Jet is above threshold?

      //cout << "Running A Jet" << endl;

      // Get the jet constituents

      //AliParticleContainer *fMCContainer = GetParticleContainer(1);
      //UInt_t nMCConstituents = fMCContainer->GetNParticles();

      //AliParticleContainer * MCTracks  = dynamic_cast<AliParticleContainer *> (GetParticleContainer("embeddedTracks"));
      //UInt_t nMCTracks = MCTracks->GetNParticles();

      /*
      TIter nextPartColl(&fParticleCollArray);
      AliParticleContainer* tracks = 0;
      while ((tracks = static_cast<AliParticleContainer*>(nextPartColl())))
      {
      AliParticleContainer* fMCContainer = tracks;
      }
      */

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
      else {
      	// Filling the histograms here
      	fhJetPt->Fill(Jet1->Pt());
      	JetPhi=Jet1->Phi();
      	if(JetPhi < -1*TMath::Pi())
      	  JetPhi += (2*TMath::Pi());
      	else if (JetPhi > TMath::Pi())
      	  JetPhi -= (2*TMath::Pi());
      	fhJetPhi->Fill(JetPhi);
      	fhJetEta->Fill(Jet1->Eta());
      	// Filling the TTree branch(es) here
        Double_t JetPt;

        JetPt = Jet1->Pt();
        fTreeBranch[0]=Jet1->Pt();

      	fTreeBranch[1]=JetPhi;
      	fTreeBranch[2]=Jet1->Eta();


        // Initialising the Tagged PYTHIA Jet.
        Bool_t kHasTruthJet = kFALSE;
        AliEmcalJet* TruthJet;
        UInt_t nTruthConstituents;

        if(Jet1->ClosestJet() != NULL)
        {
          kHasTruthJet = kTRUE;
          TruthJet = Jet1->ClosestJet();
          nTruthConstituents = TruthJet->GetNumberOfTracks();

        }


      	// Add other branches below this line ...
        // Initialise Jet shapes
        Int_t fPdgCodes[100] = {};
        Int_t fCurrentPdg = 0;
        Int_t nMothers = 0;

        Double_t jetCharge = 0;

        //Finding Pdg Code using TRUTH Jet.


        if(kHasTruthJet)
        {
          //cout << "Has Matched Jet " << endl;

          for (UInt_t iTruthConst = 0; iTruthConst < nTruthConstituents; iTruthConst++ )
          {
            AliMCParticle* TruthParticle = (AliMCParticle*) TruthJet->Track(iTruthConst);

            AliMCParticle *MotherParticle = (AliMCParticle*)  MCParticleContainer->GetParticle(TruthParticle->GetMother());

            while(MotherParticle->GetMother() > 0)
            {
              MotherParticle = (AliMCParticle*) MCParticleContainer->GetParticle(MotherParticle->GetMother());
              //cout <<"Mother Get Result: " << MotherParticle->GetMother() << endl;
              //cout << "Particle Label: " << MotherParticle->Label() << endl;
              //cout << "MParticle Pdg: " << MotherParticle->PdgCode() << endl;

            }

            if(MotherParticle->E() < 3400)
            {
            fCurrentPdg = MotherParticle->PdgCode();
            fPdgCodes[nMothers] = fCurrentPdg;
            nMothers++;
            }
          }


            // New Techniqe for determing Jet PDG Code

            Int_t UniquePdgCodes[20] = {};            //To be filled, maximium is that there are  20 uniques
            Double_t UniquePdgFrequency[20] = {};
            Int_t nUniques = 0;

            //Loop of PDG code found
            for(int i = 0; i < nMothers; i++)
            {
              // consider one pdgCode at the time.
              fCurrentPdg = fPdgCodes[i];
              // Loop over unique arry to be filled
              for(int j = 0; j < nMothers; j++)
              {
                // If it hasnt matched and the current unique value is empty list the new value and increment frequncy by 1 and then break
                if(UniquePdgCodes[j] != fCurrentPdg && UniquePdgCodes[j] == 0)
                {
                  UniquePdgCodes[j] = fCurrentPdg;
                  UniquePdgFrequency[j] = UniquePdgFrequency[j] + 1.;
                  nUniques ++;
                  break;
                }
                //Check if the PDG is already lsited if it matched increase the frequency counter by 1
                else if(UniquePdgCodes[j] == fCurrentPdg)
                {
                  UniquePdgFrequency[j] = UniquePdgFrequency[j] + 1.;
                  break;
                }

                // Otherwise the PDG hasnt matched and the value isnt zero check the next value

              }
            }

            // Setting Final Pdg Code
            // normalising frequency array

            for(unsigned int i = 0; i < nUniques; i++)
            {
              UniquePdgFrequency[i] = UniquePdgFrequency[i]/nMothers;
            }

            //Find the index of the max
            int IndexOfMaximum = -1;

            // assigens index of maximum if the over limit factation of mother particles agree
            Double_t limitFraction = 0.80;
            for(unsigned int i = 0; i < nUniques; i++)
            {
              if(UniquePdgFrequency[i] > limitFraction)
              {
                IndexOfMaximum = i;
              }
            }

            // Set final PDG Code to be used
            fCurrentPdg = fPdgCodes[IndexOfMaximum];

            if(IndexOfMaximum < 0)
            {
              fCurrentPdg = 0;
            }

                  //Outputs for Checking
          /*
                  if(nUniques > 1)
                  {
                    // output PDG Codes for testing

                    cout << "PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << fPdgCodes[i] << ",";
                    }
                    cout << "]" << endl;

                    // output Unique PDG Codes for testing

                    cout << "Unique PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << UniquePdgCodes[i] << ",";
                    }
                    cout << "]" << endl;

                    // ouput PDG Codes fort testing

                    cout << "Frequency PdgCodes : [" ;
                    for(unsigned int i = 0; i < nMothers; i++)
                    {
                      cout << UniquePdgFrequency[i] << ",";
                    }
                    cout << "]" << endl;

                    //cout << "Number of Mothers: " << nMothers << endl;
                    //cout << "Number of Uniques: " << nUniques << endl;
                    cout << "Final PDG Choice: " << fCurrentPdg << endl << endl;
                  }
          */

          // Loop over the consituents
          for (UInt_t iJetConst = 0; iJetConst < nJetConstituents; iJetConst++ )
          {
            AliVParticle *JetParticle = Jet1->Track(iJetConst);
            jetCharge += JetParticle->Charge()*pow(JetParticle->Pt(),0.5);
          }






        // Create the Non Flavoured Jet CHarge

        jetCharge/=pow(Jet1->Pt(),0.5);

        fTreeBranch[3] = jetCharge;
        JC->Fill(jetCharge);

        //Split the Jet in to apporpate momentum bin.

        if(JetPt < 40.)
        {
          fTreeBranch[4] = jetCharge;
          JCLow->Fill(jetCharge);
        }
        else if( JetPt > 40. && JetPt < 80.)
        {
          fTreeBranch[5] = jetCharge;
          JCMid->Fill(jetCharge);
        }
        else
        {
          fTreeBranch[6] = jetCharge;
          JCHigh->Fill(jetCharge);
        }



        //Add Up JetCharge
        if(fCurrentPdg == 2)
        {
          fTreeBranch[7] = jetCharge;
          JCUp->Fill(jetCharge);


            if(JetPt < 40.)
            {
              fTreeBranch[8] = jetCharge;
              JCUpLow->Fill(jetCharge);
            }
            else if( JetPt > 40. && JetPt < 80.)
            {
              fTreeBranch[9] = jetCharge;
              JCUpMid->Fill(jetCharge);
            }
            else
            {
              fTreeBranch[10] = jetCharge;
              JCUpHigh->Fill(jetCharge);
            }



        }
        //Add Down JetCharge
        else if(fCurrentPdg == 1)
        {
          fTreeBranch[11] = jetCharge;
          JCDown->Fill(jetCharge);


            if(JetPt < 40.)
            {
              fTreeBranch[12] = jetCharge;
              JCDownLow->Fill(jetCharge);
            }
            else if( JetPt > 40. && JetPt < 80.)
            {
              fTreeBranch[13] = jetCharge;
              JCDownMid->Fill(jetCharge);
            }
            else
            {
              fTreeBranch[14] = jetCharge;
              JCDownHigh->Fill(jetCharge);
            }
        }


        //Add Gluon JetCharge
        else if(fCurrentPdg == 21)
        {
          fTreeBranch[15] = jetCharge;
          JCGluon->Fill(jetCharge);


            if(JetPt < 40.)
            {
              fTreeBranch[16] = jetCharge;
              JCGluonLow->Fill(jetCharge);
            }
            else if( JetPt > 40. && JetPt < 80.)
            {
              fTreeBranch[17] = jetCharge;
              JCGluonMid->Fill(jetCharge);
            }
            else
            {
              fTreeBranch[18] = jetCharge;
              JCGluonHigh->Fill(jetCharge);
            }


        }

        //Add Unmatched JetCharge Catagory
        else if(IndexOfMaximum == -1)
        {

          fTreeBranch[23] = jetCharge;
          JCUnmatched->Fill(jetCharge);


            if(JetPt < 40.)
            {

              fTreeBranch[24] = jetCharge;
              JCUnmatchedLow->Fill(jetCharge);
            }
            else if( JetPt > 40. && JetPt < 80.)
            {

              fTreeBranch[25] = jetCharge;
              JCUnmatchedMid->Fill(jetCharge);
            }
            else
            {

              fTreeBranch[26] = jetCharge;
              JCUnmatchedHigh->Fill(jetCharge);
            }


        }

        //Adding Other Flavour JetCharge
        else
        {
          fTreeBranch[19] = jetCharge;
          JCOther->Fill(jetCharge);


            if(JetPt < 40.)
            {
              fTreeBranch[20] = jetCharge;
              JCOtherLow->Fill(jetCharge);
            }
            else if( JetPt > 40. && JetPt < 80.)
            {
              fTreeBranch[21] = jetCharge;
              JCOtherMid->Fill(jetCharge);
            }
            else
            {
              fTreeBranch[22] = jetCharge;
              JCOtherHigh->Fill(jetCharge);
            }

        }


        fTreeJets->Fill();



      }

      //cout << "End of Jet" << endl;
      }
    }
  }
  return kTRUE;
}



//________________________________________________________________________
Bool_t AliAnalysisTaskJetChargeFlavourTemplates::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetChargeFlavourTemplates::Terminate(Option_t *)
{
  // Called once at the end of the analysis.


  // Normalise historgrams over number of Jets considered

}
