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
  fJetRadius(0),
  JetChargeK(0.5),
  Pt(-999),
  Phi(-999),
  Eta(-999),
  JetCharge(-999),
  ParticlePt(-999),
  ParticlePhi(-999),
  ParticleEta(-999),
  ParticleJetCharge(-999),
  LeadingTrackPt(-999),
  PdgCode(0),
  PtMatchedPdgCode(0),
  GeoMatchedPdgCode(0),
  ProgenetorFraction(-1),


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
  fJetRadius(0),
  JetChargeK(0.5),
  Pt(-999),
  Phi(-999),
  Eta(-999),
  JetCharge(-999),
  ParticlePt(-999),
  ParticlePhi(-999),
  ParticleEta(-999),
  ParticleJetCharge(-999),
  LeadingTrackPt(-999),
  PdgCode(0),
  PtMatchedPdgCode(0),
  GeoMatchedPdgCode(0),
  ProgenetorFraction(-1),




  fTreeJets(0)
{

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
  fTreeJets->Branch("ParticlePt",&ParticlePt,"ParticlePt/F");
  fTreeJets->Branch("ParticlePhi",&ParticlePhi,"ParticlePhi/F");
  fTreeJets->Branch("ParticleEta",&ParticleEta,"ParticleEta/F");
  fTreeJets->Branch("ParticleJetCharge",&ParticleJetCharge,"ParticleJetCharge/F");
  fTreeJets->Branch("LeadingTrackPt",&LeadingTrackPt,"LeadingTrackPt/F");
  fTreeJets->Branch("PdgCode",&PdgCode,"pdgcode/I");
  fTreeJets->Branch("PtMatchedPdgCode",&PtMatchedPdgCode,"PtMatchedPdgCode/I");
  fTreeJets->Branch("GeoMatchedPdgCode",&GeoMatchedPdgCode,"GeoMatchedPdgCode/I");
  
  fTreeJets->Branch("ProgenetorFraction",&ProgenetorFraction,"ProgenetorFraction/F");







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

  // Reset the Tree Parameters
  Pt = -999;
  Phi = -999;
  Eta = -999;
  JetCharge = -999;
  ParticlePt  = -999;
  ParticlePhi = -999;
  ParticleEta = -999;
  ParticleJetCharge = -999;
  LeadingTrackPt = -999;
  PdgCode = 0;
  PtMatchedPdgCode = 0;
  GeoMatchedPdgCode = 0;
  ProgenetorFraction = -1;

  // Initialise Jet shapes
  Int_t fPdgCodes[100] = {};
  Int_t fCurrentPdg = 0;
  Int_t nMothers = 0;
  Float_t jetCharge = 0;
  Float_t jetChargeParticle = 0;

  Int_t fParticleUniqueID[100] = {};
  Int_t fCurrentParticleUniqueID = 0;

  Double_t fMotherParticlePt[100] = {};
  Double_t fCurrentMotherParticlePt = {};


  // Initialise jet pointer
  //cout << "Running Fill Histograms" << endl;
  AliEmcalJet *Jet1 = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetCont= GetJetContainer(0); //Jet Container for event
  AliEmcalJet *TruthJet = NULL; //Original Jet in the event                                                                                                         // Get jet container (0 = ?)
  AliJetContainer *JetGen= GetJetContainer(1); //Jet Container for event


  AliMCParticleContainer* MCParticleContainer = (AliMCParticleContainer*) JetGen->GetParticleContainer();

  TClonesArray*  MCParticleCloneContainer = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("mcparticles"));


  Int_t nAcceptedJets = JetCont->GetNAcceptedJets();

  Float_t JetPhi=0;
  Float_t JetParticlePhi=0;
  Float_t JetPt_ForThreshold=0;


  if(JetCont)
  {

    for (auto Jet1 : JetCont->accepted()) 
    {
    
      if(!Jet1)
      {

        continue;
      }



      AliParticleContainer *fTrackCont = JetCont->GetParticleContainer();
      UInt_t nJetConstituents = Jet1->GetNumberOfTracks();

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

        //Check the jet leading track is Greater than 5 GeV to reduce cominatorial jets

        if(Jet1->GetLeadingTrack()->Pt() < 5)
        {
            continue;
        }

        LeadingTrackPt = Jet1->GetLeadingTrack()->Pt();

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

        //Finding Pdg Code using TRUTH Jet.


        if(kHasTruthJet)
        {

          //Initalise the variables for the new technique.
          Int_t mcEntries=MCParticleCloneContainer->GetEntriesFast();
          Double_t ptpart=-1;
          Double_t dR=-99;
          const Int_t arraySize=99;
          AliAODMCParticle* CountParticle;

          Int_t countpart[arraySize] = {};
          Int_t countpartcode[arraySize] = {};
          Int_t maxInd=-1;
          Int_t count=0;
          Double_t maxPt=0;
          Int_t FoundBottomOrCharm = kFALSE;
          
          


          //Loop Through All MC Particles          
          for(Int_t i=0;i<mcEntries;i++)
          {
            AliAODMCParticle* part =  (AliAODMCParticle*)  MCParticleCloneContainer->At(i);
            // If there is no particle move to the next entry.
            if(!part)
            {
              continue;
            }
            //Gather the pdgcode from the parton
            Int_t partpdgcode=part->GetPdgCode();
            
            //Checks that the particle is a parton
            if(abs(partpdgcode)==21 || ( abs(partpdgcode)>=1 && abs(partpdgcode)<=5))
            {
              
              //Gets the Pt of the Parton and the distance from the jet radius
              ptpart=part->Pt();
              dR = TruthJet->DeltaR(part);

              // Checks if the distance between the jet is within the jet radius
              if(dR<fJetRadius)
              {
                
                //Checks if the Parton is a Bottom Quark/antiquarks
                if(abs(partpdgcode)==5)
                {
                //cout << "Parton Bottom: " <<  part->GetPdgCode() << endl;
                GeoMatchedPdgCode = part->GetPdgCode();
                FoundBottomOrCharm = kTRUE;
                //break;
                }

                else
                {

                  //This should only happen if there are too many particles such that the count falls outside of the maximise size of the array
                  if (count >arraySize-1) 
                  {
                    return 0x0; 
                  }    
                  
                  countpartcode[count]=partpdgcode;
                  countpart[count]=i;
                  
                  //Find the Index of the Parton matching the critiea with he maximium parton.
                  if(ptpart>maxPt)
                  {
                    maxPt=ptpart;
                    maxInd=i;
                  }
                  count++;
                }
              }
            }
            
          }

          //Step through the parton pdg code list and check for any charm quark/antiquarks
          for(Int_t i=0;i<count;i++)
          {
            if(abs(countpartcode[i])==4 && FoundBottomOrCharm == kFALSE )
            {
              CountParticle = (AliAODMCParticle*) MCParticleCloneContainer->At(countpart[i]);
              GeoMatchedPdgCode = CountParticle->GetPdgCode();
              //break;
            } 
          }

          //If the Maximum index has been found Collect the PDGCode for it
          if(maxInd>-1  && FoundBottomOrCharm == kFALSE)
          {
            
            AliAODMCParticle* partMax = (AliAODMCParticle*)MCParticleCloneContainer->At(maxInd);
            GeoMatchedPdgCode = partMax->GetPdgCode();
            //cout << endl << "partMax Label" << partMax->GetLabel() << endl;
            
          }

          // For Testing Perposes.
          //PDGCode = GeoMatchedPdgCode;


          if(TruthJet->Pt()<fPtThreshold)
          {
            continue;
          }



          // Filling the Tree here

          JetPhi=Jet1->Phi();
          if(JetPhi < -1*TMath::Pi())
            JetPhi += (2*TMath::Pi());
          else if (JetPhi > TMath::Pi())
            JetPhi -= (2*TMath::Pi());

          JetParticlePhi=TruthJet->Phi();
          if(JetParticlePhi < -1*TMath::Pi())
            JetParticlePhi += (2*TMath::Pi());
          else if (JetParticlePhi > TMath::Pi())
            JetParticlePhi -= (2*TMath::Pi());

          // Filling the TTree branches here
          Pt          = Jet1->Pt();
          Phi         = JetPhi;
          Eta         = Jet1->Eta();

          ParticlePt  = TruthJet->Pt();
          ParticlePhi = JetParticlePhi;
          ParticleEta = TruthJet->Eta();


          // Identifing PDG Codes of Mothers
          for (UInt_t iTruthConst = 0; iTruthConst < nTruthConstituents; iTruthConst++ )
          {
            AliAODMCParticle* TruthParticle = (AliAODMCParticle*) TruthJet->Track(iTruthConst);
            AliAODMCParticle *MotherParticle = (AliAODMCParticle*)  MCParticleContainer->GetParticle(TruthParticle->GetMother());
            AliAODMCParticle *OneSetBackParticle = MotherParticle;
            AliAODMCParticle *PtTrackerParticle = MotherParticle;
            
            //cout << "Start: " << endl;
            while(MotherParticle->GetMother() > 0)
            {
              MotherParticle = (AliAODMCParticle*) MCParticleContainer->GetParticle(MotherParticle->GetMother());
              
              //cout << "MotherParticle PDG: " << MotherParticle->PdgCode() << endl;
              //cout << "MotherParticle Label: " << MotherParticle->Label() << endl;
              //cout << "MotherParticle Pt: " << MotherParticle->Pt() << endl;
              
              if(MotherParticle->PdgCode() == 2212)
              {
                //cout << "Proton Label: " <<MotherParticle->Label() << endl;
                //cout << "One Down Label: " << OneSetBackParticle->Label() << endl;
                MotherParticle  = OneSetBackParticle;
                break;
              }
              
              if(MotherParticle->Pt() < 0.0001)
              {
                PtTrackerParticle = OneSetBackParticle;
                break; 
              }
              
              else
              {
              OneSetBackParticle = MotherParticle;           
              }
              
            }

            fCurrentParticleUniqueID = MotherParticle->GetLabel();
            fCurrentMotherParticlePt = MotherParticle->Pt();

            fCurrentPdg = MotherParticle->PdgCode();
            fPdgCodes[nMothers] = fCurrentPdg;            
            fParticleUniqueID[nMothers] = fCurrentParticleUniqueID;
            fMotherParticlePt[nMothers] = fCurrentMotherParticlePt;
            
            //cout << "Generator Index: " << MotherParticle->GetGeneratorIndex() << endl;
            //cout << fPdgCodes[nMothers] << endl;
            nMothers++;
          
          }

          
          Int_t UniquePdgCodes[20] = {};            //To be filled, maximium is that there are  20 uniques
          Float_t UniquePdgFrequency[20] = {};
          Int_t nUniques = 0;
          

          //Loop of PDG code found
          for(int i = 0; i < nMothers; i++)
          {
            
            // Consider one pdgCode at the time.
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
                //cout << UniquePdgFrequency[j] << endl;
                break;
              }
              // Otherwise the PDG hasnt matched and the value isnt zero check the next value

            }
          }

          
          // Setting Final Pdg Code
          // normalising frequency array

          for(int i = 0; i < nUniques; i++)
          {
            UniquePdgFrequency[i] = UniquePdgFrequency[i]/nMothers;
          }

          
          
         // Loop to store largest number And corresponding pdg code
         
          Float_t CurrentFraction;

          int IndexOfMaximum = -1;
          CurrentFraction = UniquePdgFrequency[0];
          for(int i = 0; i < nUniques; i++) 
          {
            if(UniquePdgFrequency[i] >= CurrentFraction)
            {
              IndexOfMaximum = i;
              CurrentFraction = UniquePdgFrequency[i];
              //cout << CurrentFraction << endl;
            }
          }
          

          ProgenetorFraction = CurrentFraction;
          fCurrentPdg = fPdgCodes[IndexOfMaximum];
          
          
          //Testing Picking the Highest momentum particle
          
        
          Int_t IndexOfMaximumPt = 0;
          Float_t CurrentHighestPt = 0;

          for(int i = 0; i < nMothers; i++) 
          {
            if(fMotherParticlePt[i] > CurrentHighestPt)
            {
              //cout << "Current Highest Pt: " << CurrentHighestPt << " - VS - " << fMotherParticlePt[i] << endl;
              IndexOfMaximumPt = i;
              CurrentHighestPt = fMotherParticlePt[i];
            }
          }

          PtMatchedPdgCode = fPdgCodes[IndexOfMaximumPt];

          

          // Now Caluclate the jet Charge



          // Loop over the consituents
          for (UInt_t iJetConst = 0; iJetConst < nJetConstituents; iJetConst++ )
          {
            AliVParticle *JetParticle = Jet1->Track(iJetConst);
            jetCharge += JetParticle->Charge()*pow(JetParticle->Pt(),JetChargeK);
          }

          for (UInt_t iTruthConst = 0; iTruthConst < nTruthConstituents; iTruthConst++ )
          {
            AliAODMCParticle* TruthParticle = (AliAODMCParticle*) TruthJet->Track(iTruthConst);
            //Divided by 3 since its parton level and in units of e/3
            jetChargeParticle += TruthParticle->Charge()/3*pow(TruthParticle->Pt(),JetChargeK);
          }



          // Normalise the Non Flavoured Jet CHarge
          jetCharge/=pow(Jet1->Pt(),0.5);

          // Normalise Particle level jet charge
          jetChargeParticle/=pow(TruthJet->Pt(),0.5);


          //Put The Jet Charge in the right place
          JetCharge = jetCharge;
          ParticleJetCharge = jetChargeParticle;

          PdgCode = fCurrentPdg;
          //PdgCode = PtMatchedPdgCode;
          //Outputs for Checking
          
          
          //cout << "Original PDG CODE: " << fCurrentPdg << endl;
          //cout << "Geo PDG CODE: " << GeoMatchedPdgCode << endl;
          //cout << "Pt PDG CODE: " << PtMatchedPdgCode << endl;
          
          //cout << "Progenetor Fraction: " << ProgenetorFraction << endl;
         /*
          //Output Frequency of Uniques
          cout << "Frequency of Uniques: [" ;
          for(int i = 0; i < nUniques; i++)
          {
            cout << UniquePdgFrequency[i] << "," ;
          }
          cout << "] " << endl;
        
          //Output for checking the list of raw pdgCodes
          cout << "Pdg Codes: [" ;
          for(int i = 0; i < nMothers; i++)
          {
            cout << fPdgCodes[i] << "," ;
          }
          cout << "] " << endl;

          
          //Outputs for checking Pts
          cout << "Particle Pt: [" ;
          for(int i = 0; i < nMothers; i++)
          {
            cout <<  fMotherParticlePt[i] << "," ;
          }
          cout << "] " << endl;
          
         
          //Outputs for checking Raw ParticleUniqueIDs
          cout << "Particle Label: [" ;
          for(int i = 0; i < nMothers; i++)
          {
            cout << fParticleUniqueID[i] << "," ;
          }
          cout << "] " << endl;
         
        
          cout << endl;
         */
         
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
