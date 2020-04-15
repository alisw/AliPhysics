// AliEmcalCorrectionClusterNonLinearityMCAfterburner.cxx
//

#include "AliEmcalCorrectionClusterNonLinearityMCAfterburner.h"
#include <TList.h>
#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterNonLinearityMCAfterburner);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearityMCAfterburner> AliEmcalCorrectionClusterNonLinearityMCAfterburner::reg("AliEmcalCorrectionClusterNonLinearityMCAfterburner");

const std::map <std::string, AliEmcalCorrectionClusterNonLinearityMCAfterburner::EMCAfterburnerMethod_t> AliEmcalCorrectionClusterNonLinearityMCAfterburner::fgkMCNonlinearityAfterburnerMethodMap = {
    { "kPCM_EMCal", AliEmcalCorrectionClusterNonLinearityMCAfterburner::kPCM_EMCal },
    { "kEMCal_EMCal", AliEmcalCorrectionClusterNonLinearityMCAfterburner::kEMCal_EMCal },
    { "kPCM_EMCalFunctionRatio", AliEmcalCorrectionClusterNonLinearityMCAfterburner::kPCM_EMCalFunctionRatio },
    { "kEMCal_EMCalFunctionRatio", AliEmcalCorrectionClusterNonLinearityMCAfterburner::kEMCal_EMCalFunctionRatio },
    { "kNoCorrection", AliEmcalCorrectionClusterNonLinearityMCAfterburner::kNoCorrection }
};

/**
 * Default constructor
 */
AliEmcalCorrectionClusterNonLinearityMCAfterburner::AliEmcalCorrectionClusterNonLinearityMCAfterburner() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterNonLinearityMCAfterburner"),
  fEnergyDistBefore(0),
  fEnergyTimeHistBefore(0),
  fEnergyDistAfter(0),
  fEnergyTimeHistAfter(0),
  fAfterburnerMethod(kNoCorrection),
  fMCPeriod(""),
  fSetForceClusterE(kFALSE)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterNonLinearityMCAfterburner::~AliEmcalCorrectionClusterNonLinearityMCAfterburner()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterNonLinearityMCAfterburner::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  std::string nonLinMethodStr = "";
  GetProperty("afterburnerVersion", nonLinMethodStr);
  EMCAfterburnerMethod_t fAfterburnerMethod = fgkMCNonlinearityAfterburnerMethodMap.at(nonLinMethodStr);
  GetProperty("setMCPeriod", fMCPeriod);
  GetProperty("setForceClusterE", fSetForceClusterE);

  //Load the right period string
  TString period = SummarizeMCProductions(fMCPeriod);
  //Load the right parameters for the MC afterburner correction
  InitNonLinearityParam(period, fAfterburnerMethod);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterNonLinearityMCAfterburner::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
  
  // Create my user objects.
  if (fCreateHisto){
    fEnergyDistBefore = new TH1F("hEnergyDistBefore","hEnergyDistBefore;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistBefore);
    fEnergyTimeHistBefore = new TH2F("hEnergyTimeDistBefore","hEnergyTimeDistBefore;E_{clus} (GeV);time (s)",1500,0,150,500,-1e-6,1e-6);
    fOutput->Add(fEnergyTimeHistBefore);
    fEnergyDistAfter = new TH1F("hEnergyDistAfter","hEnergyDistAfter;E_{clus} (GeV)",1500,0,150);
    fOutput->Add(fEnergyDistAfter);
    fEnergyTimeHistAfter = new TH2F("hEnergyTimeDistAfter","hEnergyTimeDistAfter;E_{clus} (GeV);time (s)",1500,0,150,500,-1e-6,1e-6);
    fOutput->Add(fEnergyTimeHistAfter);
    
    // Take ownership of output list
    fOutput->SetOwner(kTRUE);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterNonLinearityMCAfterburner::Run()
{
  AliEmcalCorrectionComponent::Run();

  // loop over clusters
  AliVCluster *clus = 0;
  AliTLorentzVector clusVec;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);

  //Apply the nonlinearity afterburner correction for the MC.
  //Parameters are set in InitNonLinearityParam() and copied from AliCaloPhotonCuts.cxx
  //Apply the MC nonlinearity functions with the set parameters
  if (fNonLinearityAfterburnerFunction != 0) {

    while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {

      if (!clusCont) return kFALSE;
      auto clusItCont = clusCont->all_momentum();

      for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
        clus    = static_cast<AliVCluster *>(clusIterator->second);
        clusVec = static_cast<AliTLorentzVector>(clusIterator->first);

        if (!clus->IsEMCAL()) continue;
        //currently only for EMCal clusters. Skip DCal clusters
        if (clusVec.Phi_0_2pi() > 4.) {continue;}

        Double_t oldEnergy = clus->GetNonLinCorrEnergy();
        Double_t newEnergy = -1;

        if (fCreateHisto) {
          fEnergyDistBefore->Fill(oldEnergy);
          fEnergyTimeHistBefore->Fill(oldEnergy, clus->GetTOF());
        }

        //This is function Iteration-1 FunctionNL_kSDM combined with Iteration-2 FunctionNL_kSDM from AliCaloPhotonCuts
        if (fNonLinearityAfterburnerFunction == 1) {
          newEnergy = oldEnergy;
          newEnergy/= fNLAfterburnerPara[0] + exp( fNLAfterburnerPara[1] + (fNLAfterburnerPara[2] * newEnergy));
          newEnergy/= fNLAfterburnerPara[3] + exp( fNLAfterburnerPara[4] + (fNLAfterburnerPara[5] * newEnergy));
          AliDebugStream(2) << "Using non-lin afterburner function No.1" << std::endl;
          AliDebugStream(2) << "old Energy: "<< oldEnergy<<", new Energy: "<<newEnergy<< std::endl;
         }
        //This is function FunctionNL_DPOW combined with FunctionNL_kSDM from AliCaloPhotonCuts
        if (fNonLinearityAfterburnerFunction == 2) {
          newEnergy = oldEnergy;
          newEnergy/= (fNLAfterburnerPara[0] + fNLAfterburnerPara[1] * TMath::Power(oldEnergy, fNLAfterburnerPara[2]))/(fNLAfterburnerPara[3] + fNLAfterburnerPara[4] * TMath::Power(oldEnergy, fNLAfterburnerPara[5]));
          newEnergy/= fNLAfterburnerPara[6] + exp( fNLAfterburnerPara[7] + (fNLAfterburnerPara[8] * newEnergy));
          AliDebugStream(2) << "Using non-lin afterburner function No.2" << std::endl;
          AliDebugStream(2) << "old Energy: "<< oldEnergy<<", new Energy: "<<newEnergy<< std::endl;
        }
        //This is function FunctionNL_DExp combined with FunctionNL_kSDM from AliCaloPhotonCuts
        if (fNonLinearityAfterburnerFunction == 3) {
          newEnergy = oldEnergy;
          newEnergy/= ((fNLAfterburnerPara[0] - TMath::Exp(-fNLAfterburnerPara[1]*oldEnergy +fNLAfterburnerPara[2]))/(fNLAfterburnerPara[3] - TMath::Exp(-fNLAfterburnerPara[4]*oldEnergy +fNLAfterburnerPara[5])));
          newEnergy/= fNLAfterburnerPara[6] + exp( fNLAfterburnerPara[7] + (fNLAfterburnerPara[8] * newEnergy));
          AliDebugStream(2) << "Using non-lin afterburner function No.3" << std::endl;
          AliDebugStream(2) << "old Energy: "<< oldEnergy<<", new Energy: "<<newEnergy<< std::endl;
       }

        //correction only works in a specific momentum range
        if (oldEnergy > 0.300) {
          clus->SetNonLinCorrEnergy(newEnergy);
          if (fSetForceClusterE) clus->SetE(newEnergy);
        }
        else {
          AliWarning(Form("Cluster energy out of range for correction!, E = %f \n",oldEnergy));
          clus->SetNonLinCorrEnergy(-100);
          if (fSetForceClusterE) clus->SetE(-100);
        }

        // Fill histograms only if cluster is not exotic, as in ClusterMaker (the clusters are flagged, not removed)
        if (fCreateHisto && !clus->GetIsExotic()) {
          Float_t energy = clus->GetNonLinCorrEnergy();
          if (fSetForceClusterE) energy = clus->E();

          fEnergyDistAfter->Fill(energy);
          fEnergyTimeHistAfter->Fill(energy, clus->GetTOF());
        }
      }
    }
    return kTRUE;
  }
  else {
    //No correction has been determined for this data set yet - switch this component off for this period
    AliWarning("No correction has been determined for this data set yet - switch this component off for this period");
    return kFALSE;
  }
}

// returns the period enumerator for a given period string
// can be used to obtain the enum for GetCorrectedEnergy
//________________________________________________________________________
TString AliEmcalCorrectionClusterNonLinearityMCAfterburner::SummarizeMCProductions(TString namePeriod)
{
  // Globally valid testbeam fine tuning (separate values for Run1 and Run2)
  if ( namePeriod.CompareTo("kTestBeamFinalMCRun2") == 0  )    return "kTestBeamFinalMCRun2";
  else if ( namePeriod.CompareTo("kTestBeamFinalMCRun1") == 0  )    return "kTestBeamFinalMCRun1";
  //...MC anchored to 2015 Data...
  // pp 13 TeV
  // - 2015 MB pass 2
  else if ( namePeriod.CompareTo("LHC15P2Pyt8") == 0 ||
      namePeriod.CompareTo("LHC17i4") == 0 ||
      namePeriod.CompareTo("LHC17i4_2") == 0 ||
      namePeriod.CompareTo("LHC17g7") == 0 )              return "kPP13T15P2Pyt8";
  else if ( namePeriod.CompareTo("LHC15P2EPos") == 0 ||
      namePeriod.CompareTo("LHC16d3") == 0 )              return "kPP13T15P2EPOS";
  // - 2015 HF prod
  else if ( namePeriod.CompareTo("LHC15k5a") == 0 ||
      namePeriod.CompareTo("LHC15k5b") == 0 ||
      namePeriod.CompareTo("LHC15k5c") == 0 ||
      namePeriod.CompareTo("LHC15k5a2") == 0 ||
      namePeriod.CompareTo("LHC15k5b2") == 0 ||
      namePeriod.CompareTo("LHC15k5c2") == 0 )          return "k15k5";
  //
  // pp 5 TeV
  // - 2015 MB MC pass 4
  else if ( namePeriod.CompareTo("LHC17e2") == 0  )     return "k17e2";  //Corr-1
  else if ( namePeriod.CompareTo("LHC18j3") == 0  )     return "k18j3";  //Corr-1
  // - 2015 JJ pass 3
  else if ( namePeriod.CompareTo("LHC16h3") == 0 )      return "k16h3";  //Corr-1
  //
  // PbPb 5TeV
  // - 2015 MB prods
  else if ( namePeriod.CompareTo("LHC16g1") == 0 ||
      namePeriod.CompareTo("LHC16g1a") == 0 ||
      namePeriod.CompareTo("LHC16g1b") == 0 || namePeriod.CompareTo("LHC16g1b_extra") == 0 ||
      namePeriod.CompareTo("LHC16g1c") == 0 || namePeriod.CompareTo("LHC16g1c_extra") == 0 ||
      namePeriod.CompareTo("LHC16h4") == 0)             return "kPbPb5T15HIJING";
  else if ( namePeriod.CompareTo("LHC16k3b") == 0 ||                    // special pileup prods
      namePeriod.CompareTo("LHC16k3b2") == 0 )          return "k16k3b";
  //...MC anchored to 2016 Data...
  //
  // pp 13 TeV
  // - 2016 MB prod
  else if ( namePeriod.CompareTo("LHC16P1Pyt8") == 0 ||
      namePeriod.CompareTo("LHC17f6") == 0 ||
      namePeriod.CompareTo("LHC17d17") == 0 ||
      namePeriod.CompareTo("LHC17f5") == 0 ||
      namePeriod.CompareTo("LHC17d3") == 0 ||
      namePeriod.CompareTo("LHC17e5") == 0 ||
      namePeriod.CompareTo("LHC17d20a1") == 0 ||
      namePeriod.CompareTo("LHC17d20a1_extra") == 0 ||
      namePeriod.CompareTo("LHC17d20a2") == 0 ||
      namePeriod.CompareTo("LHC17d20a2_extra") == 0 ||
      namePeriod.CompareTo("LHC17d16") == 0 ||
      namePeriod.CompareTo("LHC17d18") == 0 ||
      namePeriod.CompareTo("LHC17f9") == 0 ||
      namePeriod.CompareTo("LHC17f9_test") == 0||
      namePeriod.CompareTo("LHC18f1") == 0 ||
      namePeriod.CompareTo("LHC18d8") == 0 )            return "kPP13T16P1Pyt8";    //Corr-3
  else if ( namePeriod.CompareTo("LHC16P1Pyt8LowB") == 0 ||
      namePeriod.CompareTo("LHC17d1") == 0 )            return "kPP13T16P1Pyt8LowB"; //Corr-3
  else if ( namePeriod.CompareTo("LHC16P1EPOS") == 0 ||
      namePeriod.CompareTo("LHC17d20b1") == 0 ||
      namePeriod.CompareTo("LHC17d20b2") == 0 )         return "kPP13T16P1EPOS";
  // - 2016 JJ prod
  else if ( namePeriod.CompareTo("LHC16P1JJ") == 0 ||
      namePeriod.CompareTo("LHC17f8a") == 0 ||
      namePeriod.CompareTo("LHC17f8c") == 0 ||
      namePeriod.CompareTo("LHC17f8d") == 0 ||
      namePeriod.CompareTo("LHC17f8e") == 0 ||
      namePeriod.CompareTo("LHC17f8f") == 0 ||
      namePeriod.CompareTo("LHC17f8g") == 0 ||
      namePeriod.CompareTo("LHC17f8h") == 0 ||
      namePeriod.CompareTo("LHC17f8i") == 0 ||
      namePeriod.CompareTo("LHC17f8j") == 0 ||
      namePeriod.CompareTo("LHC17f8k") == 0 )           return "kPP13T16P1JJ";   //Corr-3
  else if ( namePeriod.CompareTo("LHC16P1JJLowB") == 0 ||
      namePeriod.CompareTo("LHC17f8b") == 0  )          return "kPP13T16P1JJLowB";
  // - 2016 HF prods
  else if ( namePeriod.CompareTo("LHC17h8a") == 0 )     return "k17h8a";
  else if ( namePeriod.CompareTo("LHC17h8b") == 0 )     return "k17h8b";
  else if ( namePeriod.CompareTo("LHC17h8c") == 0 )     return "k17h8c";
  else if ( namePeriod.CompareTo("LHC17c3b1") == 0 )    return "k17c3b1";
  else if ( namePeriod.CompareTo("LHC17c3a1") == 0 )    return "k17c3a1";
  else if ( namePeriod.CompareTo("LHC17c3b2") == 0 )    return "k17c3b2";
  else if ( namePeriod.CompareTo("LHC17c3a2") == 0 )    return "k17c3a2";
  // - 2016 GJ prod
  else if ( namePeriod.CompareTo("LHC17i3a1") == 0 )    return "k17i3a1";
  //
  // pPb 5 TeV
  // - 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f2a") == 0  )    return "kPPb5T16EPOS";
  else if ( namePeriod.CompareTo("LHC17f2b") == 0 ||
      namePeriod.CompareTo("LHC18f3")  == 0 )           return "kPPb5T16DPMJet";
  // - 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8a") == 0  )    return "k17g8a";
  // - 2016 HF prod
  else if ( namePeriod.CompareTo("LHC17d2a") == 0 )     return "k17d2a";
  else if ( namePeriod.CompareTo("LHC17d2b") == 0 )     return "k17d2b";
  //
  // pPb 8 TeV
  // - 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f3a") == 0  )    return "k17f3a";
  else if ( namePeriod.CompareTo("LHC17f3b") == 0  )    return "k17f3b";
  else if ( namePeriod.CompareTo("LHC17f4a") == 0  )    return "k17f4a";
  else if ( namePeriod.CompareTo("LHC17f4b") == 0  )    return "k17f4b";
  // - 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8b") == 0  )    return "k17g8b";
  else if ( namePeriod.CompareTo("LHC17g8c") == 0  )    return "k17g8c";
  else if ( namePeriod.CompareTo("LHC18b9b") == 0  )    return "k18b9b";
  else if ( namePeriod.CompareTo("LHC18b9c") == 0  )    return "k18b9c";
  //...MC anchored to 2017 Data...
  //
  //pp 13 TeV LHC17
  else if ( namePeriod.CompareTo("LHC17k1") ==0 )       return "k17k1"; // HF low B
  else if ( namePeriod.CompareTo("LHC17P1Pyt8NomB") ==0 ||
      namePeriod.CompareTo("LHC17k4") ==0 ||
      namePeriod.CompareTo("LHC17h11") ==0 ||
      namePeriod.CompareTo("LHC17h1") == 0 ||
      namePeriod.CompareTo("LHC17l5") == 0 ||
      namePeriod.CompareTo("LHC18c13") == 0 ||
      namePeriod.CompareTo("LHC18a8") == 0 ||
      namePeriod.CompareTo("LHC18a9") == 0 ||
      namePeriod.CompareTo("LHC18a1") == 0 ||
      namePeriod.CompareTo("LHC18c12") == 0 )           return "kPP13T17P1Pyt8";     //Corr-3
  else if ( namePeriod.CompareTo("LHC17h7b") ==0 )      return "kPP13T17P1Pho";
  else if ( namePeriod.CompareTo("LHC17h7a") ==0 )      return "kPP13T17P1Pyt6";
  else if ( namePeriod.CompareTo("LHC17j5a") ==0 ||
      namePeriod.CompareTo("LHC17j5b") ==0 ||
      namePeriod.CompareTo("LHC17j5c") ==0 ||
      namePeriod.CompareTo("LHC17j5d") ==0 ||
      namePeriod.CompareTo("LHC17j5e") ==0 )            return "kPP13T17P1Pyt8Str";
  else if ( namePeriod.CompareTo("LHC17P1Pyt8LowB") ==0 ||
      namePeriod.CompareTo("LHC17h3") == 0 )            return "kPP13T17P1Pyt8LowB"; //Corr-3
  else if ( namePeriod.CompareTo("LHC17P1JJ") == 0 ||
      namePeriod.CompareTo("LHC18f5") == 0 )            return "kPP13T17P1JJ";       //Corr-3
  // pp 5 TeV
  // - 2017 MB MC pass 1
  else if ( namePeriod.CompareTo("LHC17l4b") == 0 ||
      namePeriod.CompareTo("LHC17l4b_fast") == 0 ||
      namePeriod.CompareTo("LHC17l4b_cent") == 0 ||
      namePeriod.CompareTo("LHC17l4b_cent_woSDD") == 0) return "k17l4b";  //Corr-2
  else if ( namePeriod.CompareTo("LHC17l3b") == 0 ||
      namePeriod.CompareTo("LHC17l3b_fast") == 0 ||
      namePeriod.CompareTo("LHC17l3b_cent") == 0 ||
      namePeriod.CompareTo("LHC17l3b_cent_woSDD") == 0) return "k17l3b";  //Corr-2
  else if ( namePeriod.CompareTo("LHC18j2") == 0 ||
      namePeriod.CompareTo("LHC18j2_fast") == 0 ||
      namePeriod.CompareTo("LHC18j2_cent") == 0 ||
      namePeriod.CompareTo("LHC18j2_cent_woSDD") == 0)  return "k18j2";   //Corr-2
  // - 2017 JJ pass 3
  else if ( namePeriod.Contains("LHC18b8") )            return "k18b8";   //Corr-2
  //
  // XeXe 5.44 TeV 2017 MB MC
  else if ( namePeriod.CompareTo("LHC17j7") == 0 ||
      namePeriod.CompareTo("LHC18d2") == 0 ||
      namePeriod.CompareTo("LHC18d2_1") == 0 ||
      namePeriod.CompareTo("LHC18d2_2") == 0 ||
      namePeriod.CompareTo("LHC18d2_3") == 0 )          return "kXeXe5T17HIJING";
  //...MC anchored to 2018 Data...
  //
  //pp 13 TeV LHC18
  else if ( namePeriod.CompareTo("LHC18P1Pyt8NomB") ==0 ||
      namePeriod.CompareTo("LHC18g4") ==0 ||
      namePeriod.CompareTo("LHC18g5") ==0 ||
      namePeriod.CompareTo("LHC18g6") == 0 )            return "kPP13T18P1Pyt8";  //Corr-3
  else if ( namePeriod.CompareTo("LHC18P1Pyt8LowB") ==0 ||
      namePeriod.CompareTo("LHC18h1") ==0  )            return "kPP13T18P1Pyt8LowB";

  else return "kNoMC";
}
// This function sets the parameters for the MC nonlinearity afterburner function
// this depends on the period. Several periods share the same nonlinerarity MC correction parameters.
// Sometimes the correction was done in several iterations, in this case you will find the parameters
// for the differen iterations listed. The total correction should be a combination of both parameter sets
//
// As input you need to specify the namePeriod which summarizes several MC productions
// and you need to specify the version 0,1,2,3. 0 Is what should be used as default and 1,2,3 are
// variations for systematic uncertainty studies
//________________________________________________________________________
void AliEmcalCorrectionClusterNonLinearityMCAfterburner::InitNonLinearityParam(TString namePeriod, EMCAfterburnerMethod_t method)
{
  switch(method) {
    //. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    //These are the standard parameters extracted from the PCM method (one hit in EMCal + one converted gamma)
    case kPCM_EMCal:
      // globally valid Run2 fine tuning for the testbeam+shaper corrected clusters
      if( namePeriod=="kTestBeamFinalMCRun2") {
        fNonLinearityAfterburnerFunction = 1;
        //There are no extracted parameters yet
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0.987912;
        fNLAfterburnerPara[1] = -2.94105;
        fNLAfterburnerPara[2] = -0.273207;
        //Iteration-2 parameters
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      if( namePeriod=="kTestBeamFinalMCRun1") {
        fNonLinearityAfterburnerFunction = 1;
        //There are no extracted parameters yet
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0.987912;
        fNLAfterburnerPara[1] = -2.94105;
        fNLAfterburnerPara[2] = -0.273207;
        //Iteration-2 parameters
        fNLAfterburnerPara[3] = -0.0125; // Run1 additional correction factor of 1.25%
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      // pp 5.02 TeV (2015) - LHC15n
      if( namePeriod=="k16h3" ||  namePeriod=="k17e2" || namePeriod == "k18j3") {
        fNonLinearityAfterburnerFunction = 1;
        //There are no extracted parameters yet
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        //Iteration-2 parameters
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      // pp 5.02 TeV (2017) - LHC17pq
      else if(namePeriod=="k17l3b" || namePeriod=="k18j2" || namePeriod=="k17l4b" || namePeriod=="k18b8") {
        fNonLinearityAfterburnerFunction = 2;
        //These are parameters extracted by Nicolas Schmidt in 25.Nov 2018 (to be used after kPi0MCv3 and kBeamTestCorrectedv4)
        fNLAfterburnerPara[0] = -0.8802886739;
        fNLAfterburnerPara[1] = 1.8764944987;
        fNLAfterburnerPara[2] = -0.0020594487;
        fNLAfterburnerPara[3] = 0.9891399006;
        fNLAfterburnerPara[4] = 0.0139889085;
        fNLAfterburnerPara[5] = -2.0388063034;
        fNLAfterburnerPara[6] = 0;
        fNLAfterburnerPara[7] = 0;
        fNLAfterburnerPara[8] = 0;
      }
      //pp 13 TeV LHC16 || LHC17 || LHC18 (4)
      else if ( namePeriod=="kPP13T16P1Pyt8" || namePeriod=="kPP13T17P1Pyt8" || namePeriod=="kPP13T18P1Pyt8" || namePeriod=="kPP13T17P1JJ" || namePeriod=="kPP13T16P1JJ" || namePeriod=="kPP13T16P1Pyt8LowB" || namePeriod=="kPP13T17P1Pyt8LowB") {
        fNonLinearityAfterburnerFunction = 1;
        //There are no extracted parameters yet
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        //No Second Iteration
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      // All other periods
      else {
        //This is effectivley NO MC Nonlinearity afterburner correction
        fNonLinearityAfterburnerFunction = 0;
        AliWarning(Form("There is no defined MC Nonlinearity Afterburner correction for this period (%s)! Nothing will be corrected",namePeriod.Data()));
      }
      break;
      //. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      //These are variations of the standard parameters extracted from the EMCal method (both hits in the EMCal)
    case kEMCal_EMCal:
      // pp 5.02 TeV (2015) - LHC15n
      if( namePeriod=="k16h3" ||  namePeriod=="k17e2" || namePeriod =="k18j3") {
        fNonLinearityAfterburnerFunction = 1;
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        //No Second Iteration
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      // pp 5.02 TeV (2017) - LHC17pq
      else if(namePeriod=="k17l3b" || namePeriod=="k18j2" || namePeriod=="k17l4b" || namePeriod=="k18b8") {
        fNonLinearityAfterburnerFunction = 1;
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        //Iteration-2 paramters
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      //(4)
      //pp 13 TeV LHC16 || LHC17 || LHC18
      else if ( namePeriod=="kPP13T16P1Pyt8" || namePeriod=="kPP13T17P1Pyt8" || namePeriod=="kPP13T18P1Pyt8" || namePeriod=="kPP13T16P1Pyt8LowB" || namePeriod=="kPP13T17P1Pyt8LowB" || namePeriod=="kPP13T17P1JJ" || namePeriod=="kPP13T16P1JJ"){
        fNonLinearityAfterburnerFunction = 1;
        //Iteration-1 paramters
        fNLAfterburnerPara[0] = 0;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        //Iteration-2 paramters
        fNLAfterburnerPara[3] = 0;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
      }
      // All other periods
      else {
        //This is effectivley NO MC Nonlinearity afterburner correction
        fNonLinearityAfterburnerFunction = 0;
        AliWarning(Form("There is no defined MC Nonlinearity Afterburner correction for this period (%s)! Nothing will be corrected",namePeriod.Data()));
      }
      break;
      //. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      //These are the standard parameters extracted from the PCM method (one hit in EMCal + one converted gamma)
      //Instead of fitting the ratio of the data points in this version the masses are interpolated with a
      //fit and then the fits are divided by each other to give the right function
    case kPCM_EMCalFunctionRatio:
      // pp 5.02 TeV (2015) - LHC15n AND (2017) - LHC17pq
      if( namePeriod=="k16h3" ||  namePeriod=="k17e2" || namePeriod == "k18j3" || namePeriod=="k17l3b" || namePeriod=="k18j2" || namePeriod=="k17l4b" || namePeriod=="k18b8") {
        fNonLinearityAfterburnerFunction = 2;
        fNLAfterburnerPara[0] = 1;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 1;
        fNLAfterburnerPara[3] = 1;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 1;
        fNLAfterburnerPara[6] = 0;
        fNLAfterburnerPara[7] = 0;
        fNLAfterburnerPara[8] = 0;
      }
      //(4)
      //pp 13 TeV LHC16 || LHC17 || LHC18
      else if ( namePeriod=="kPP13T16P1Pyt8" || namePeriod=="kPP13T17P1Pyt8" || namePeriod=="kPP13T18P1Pyt8" || namePeriod=="kPP13T17P1JJ" || namePeriod=="kPP13T16P1JJ" || namePeriod=="kPP13T16P1Pyt8LowB" || namePeriod=="kPP13T17P1Pyt8LowB"){
        fNonLinearityAfterburnerFunction = 2;
        fNLAfterburnerPara[0] = 1;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 1;
        fNLAfterburnerPara[3] = 1;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 1;
        fNLAfterburnerPara[6] = 0;
        fNLAfterburnerPara[7] = 0;
        fNLAfterburnerPara[8] = 0;
      }
      //(4-2)
      //pp 13 TeV LHC16 || LHC17 || LHC18
      /*if ( DCalEnergy ){
				fNonLinearityAfterburnerFunction = 2;
				fNLAfterburnerPara[0] = 1.0167588250;
				fNLAfterburnerPara[1] = 0.0501002307;
				fNLAfterburnerPara[2] = -0.8336787497;
				fNLAfterburnerPara[3] = 0.9500009312;
				fNLAfterburnerPara[4] = 0.0944118922;
				fNLAfterburnerPara[5] = -0.1043983134;
				fNLAfterburnerPara[6] = 0;
				fNLAfterburnerPara[7] = 0;
				fNLAfterburnerPara[8] = 0;
			}*/
      // All other periods
      else {
        //This is effectivley NO MC Nonlinearity afterburner correction
        fNonLinearityAfterburnerFunction = 0;
        AliWarning(Form("There is no defined MC Nonlinearity Afterburner correction for this period (%s)! Nothing will be corrected",namePeriod.Data()));
      }
      break;
      //. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      //These are variations of the standard parameters extracted from the EMCal method (both hits in the EMCal)
      //Instead of fitting the ratio of the data points in this version the masses are interpolated with a
      //fit and then the fits are divided by each other to give the right function
    case kEMCal_EMCalFunctionRatio:
      // pp 5.02 TeV (2015) - LHC15n AND (2017) - LHC17pq
      if( namePeriod=="k16h3" ||  namePeriod=="k17e2" || namePeriod == "k18j3" || namePeriod=="k17l3b" || namePeriod=="k18j2" || namePeriod=="k17l4b" || namePeriod=="k18b8") {
        fNonLinearityAfterburnerFunction = 3;
        fNLAfterburnerPara[0] = 2;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 0;
        fNLAfterburnerPara[3] = 2;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 0;
        fNLAfterburnerPara[6] = 0;
        fNLAfterburnerPara[7] = 0;
        fNLAfterburnerPara[8] = 0;
      }
      //pp 13 TeV LHC16 || LHC17 || LHC18
      else if ( namePeriod=="kPP13T16P1Pyt8" || namePeriod=="kPP13T17P1Pyt8" || namePeriod=="kPP13T18P1Pyt8" || namePeriod=="kPP13T16P1Pyt8LowB" || namePeriod=="kPP13T17P1Pyt8LowB" || namePeriod=="kPP13T17P1JJ" || namePeriod=="kPP13T16P1JJ"){
        fNonLinearityAfterburnerFunction = 2;
        fNLAfterburnerPara[0] = 1;
        fNLAfterburnerPara[1] = 0;
        fNLAfterburnerPara[2] = 1;
        fNLAfterburnerPara[3] = 1;
        fNLAfterburnerPara[4] = 0;
        fNLAfterburnerPara[5] = 1;
        fNLAfterburnerPara[6] = 0;
        fNLAfterburnerPara[7] = 0;
        fNLAfterburnerPara[8] = 0;
      }
      // All other periods
      else {
        //This is effectivley NO MC Nonlinearity afterburner correction
        fNonLinearityAfterburnerFunction = 0;
        AliWarning(Form("There is no defined MC Nonlinearity Afterburner correction for this period (%s)! Nothing will be corrected",namePeriod.Data()));
      }
      break;
    case kNoCorrection:
      //This is effectivley NO MC Nonlinearity afterburner correction
      fNonLinearityAfterburnerFunction = 0;
      AliWarning("You have to specify which method you want to use to get the parameters for the correction. Standard is kPCM_EMCal");
    break;
  }
}
