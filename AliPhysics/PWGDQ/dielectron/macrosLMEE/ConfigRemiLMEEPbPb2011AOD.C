//#include "PWGDQ/dielectron/macrosLMEE/LMEECutLib.C"
void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void EnableMC();

TString names=("noPairing;ITSTPCTOFCentnoRej;ITSTPCTOFSemiCent1noRej;ITSTPCTOFSemiCent2noRej;ITSTPCTOFPerinoRej;ITSTPCTOFCentInvMLowRP;ITSTPCTOFSemiCent1InvMLowRP;ITSTPCTOFSemiCent2InvMLowRP;ITSTPCTOFPeriInvMLowRP;ITSTPCTOFCentInvMMiddleRP;ITSTPCTOFSemiCent1InvMMiddleRP;ITSTPCTOFSemiCent2InvMMiddleRP;ITSTPCTOFPeriInvMMiddleRP;ITSTPCTOFCentInvMHighRP;ITSTPCTOFSemiCent1InvMHighRP;ITSTPCTOFSemiCent2InvMHighRP;ITSTPCTOFPeriInvMHighRP;ITSTPCTOFCentInvMLowMag;ITSTPCTOFSemiCent1InvMLowMag;ITSTPCTOFSemiCent2InvMLowMag;ITSTPCTOFPeriInvMLowMag;ITSTPCTOFCentInvMMiddleMag;ITSTPCTOFSemiCent1InvMMiddleMag;ITSTPCTOFSemiCent2InvMMiddleMag;ITSTPCTOFPeriInvMMiddleMag;ITSTPCTOFCentInvMHighMag;ITSTPCTOFSemiCent1InvMHighMag;ITSTPCTOFSemiCent2InvMHighMag;ITSTPCTOFPeriInvMHighMag");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t MCenabled=kFALSE;


AliDielectron* ConfigRemiLMEEPbPb2011AOD(Int_t cutDefinition, Bool_t hasMC=kFALSE, Bool_t ESDanalysis=kFALSE)
{

  Int_t selectedPID=-1;
  Int_t selectedCentrality=-1;
    Int_t selectedPairInvMassCut=-1;
    Int_t selectedPairInOutCut = -1;
  Bool_t rejectionStep=kFALSE;
  Bool_t PairInvMassCut=kFALSE;
  Bool_t PairInOutCut=kFALSE;
  LMEECutLibRemi*  LMCL = new LMEECutLibRemi();

  //
  // Setup the instance of AliDielectron
  //

  MCenabled=hasMC;

  // create the actual framework object

  TString name=Form("%02d",cutDefinition);
  if ((cutDefinition)<arrNames->GetEntriesFast()){
	name=arrNames->At((cutDefinition))->GetName();
  }

  //thisCut only relevant for MC:
  AliDielectron *die =
	new AliDielectron(Form
		("%s",name.Data()),
		Form("Track cuts: %s",name.Data()));

  //  TString ZDCRecenteringfile = "alien:///alice/cern.ch/user/r/rtanizak/ZDCrpH1/ZDCRecentProf/ZDCRecenteringProfile.root";
    TString ZDCRecenteringfile = "/home/tanizaki/nfs/ZDCrpH1Recentering/ZDCRecenteringProfile.root";

  die->SetZDCRecenteringFilename(ZDCRecenteringfile);


  //Setup AnalysisSelection:
  if (cutDefinition==0) {
	//not yet implemented
  }
  else if (cutDefinition==1) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    rejectionStep = kFALSE;
    PairInvMassCut = kFALSE;

  }
  else if (cutDefinition==2) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    rejectionStep = kFALSE;
    PairInvMassCut = kFALSE;

  }
  else if (cutDefinition==3) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    rejectionStep = kFALSE;
    PairInvMassCut = kFALSE;

  }
  else if (cutDefinition==4) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    rejectionStep = kFALSE;
    PairInvMassCut = kFALSE;

  }

  //////////////////////////////////////////////////////////

  else if (cutDefinition==5) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
      }

  else if (cutDefinition==6) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==7) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==8) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==9) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==10) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==11) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==12) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==13) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==14) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==15) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==16) {
    selectedPID = LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;

  }

  //////////////////////////////////////////////////////////

  else if (cutDefinition==17) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==18) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==19) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==20) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassLow;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==21) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==22) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==23) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==24) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }

  else if (cutDefinition==25) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==26) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==27) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }
  else if (cutDefinition==28) {
    selectedPID =  LMEECutLibRemi::kPbPb2011pidITSTPCTOF;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011MassHigh;
    selectedPairInOutCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairInvMassCut = kTRUE;
    PairInOutCut = kTRUE;
  }





    /*
    ///////////////////////////////////////////////////////////
  else if (cutDefinition==21) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011RP;
    //selectedPairMCut = LMEECutLibRemi::kPbPb2011MassAll;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;
  }
  else if (cutDefinition==22) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011RP;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;
  }
  else if (cutDefinition==23) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011RP;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }
  else if (cutDefinition==24) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011RP;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }

    /////////////////////////////////////////////////////////
  else if (cutDefinition==25) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Central;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011Mag;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassAll;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }

  else if (cutDefinition==26) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral1;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011Mag;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }
  else if (cutDefinition==27) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011SemiCentral2;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011Mag;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }
  else if (cutDefinition==28) {
    selectedPID = LMEECutLibRemi::kPbPb2011TPCandTOFwide;
    selectedCentrality = LMEECutLibRemi::kPbPb2011Peripheral;
    selectedPairInvMassCut = LMEECutLibRemi::kPbPb2011Mag;
    // selectedPairMCut = LMEECutLibRemi::kPbPb2011MassMiddle;
    rejectionStep = kFALSE;
    PairCut=kTRUE;
    PolaCut=kTRUE;

  }

    */

  else {
	cout << " =============================== " << endl;
	cout << " ==== INVALID CONFIGURATION ==== " << endl;
	cout << " =============================== " << endl;
  }


  //Now configure task

  //Apply correct Pre-Filter Scheme, if necessary
  die->SetPreFilterAllSigns();

  //switch off KF PArticle:
  die->SetUseKF(kFALSE);
  /*
  if (selectedPID == LMEECutLibRemi::kPbPb2011NoPID) {
	  die->SetNoPairing();
   }
  */



  die->GetEventFilter().AddCuts( LMCL->GetCentralityCuts(selectedCentrality));
  
  //  die->GetTrackFilter().AddCuts( LMCL->GetTrackCutsAna(selectedPID) );
  die->GetTrackFilter().AddCuts( LMCL->GetPIDCutsAna(selectedPID) );
  //  die->GetPairFilter().AddCuts( LMCL->GetPairCutsAna(selectedPID,kFALSE) );

    
  if (PairInvMassCut)
    die->GetPairFilter().AddCuts( LMCL->GetPairCutsInvMass(selectedPairInvMassCut));
  
  
      
  if(PairInOutCut)
        die->GetPairFilter().AddCuts( LMCL->GetPairCutsInOut(selectedPairInOutCut));
  
  

  /*

  if(PairCut){
    if (rejectionStep) {
      die->GetPairPreFilterLegs().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
      die->GetPairPreFilter().AddCuts( LMCL->GetPairPreFilterCuts(selectedPairCut));
      die->GetPairFilter().AddCuts( LMCL->GetPairCuts(selectedPairCut));
    }
    else {
      //      die->GetPairFilter().AddCuts( LMCL->GetPairCutsInvMass(selectedPairCut));
      die->GetPairFilter().AddCuts( LMCL->GetPairCuts(selectedPID));


      //  die->GetPairFilter().AddCuts( LMCL->GetPairCuts4(selectedPairMCut));
    }
  }
  */

  /*
  if (rejectionStep) {
    if (ESDanalysis) {
      die->GetTrackFilter().AddCuts( LMCL->GetESDTrackCutsAna(selectedPID) );
      die->GetPairPreFilterLegs().AddCuts( LMCL->GetESDTrackCutsAna(selectedPID) );
    }

    //die->GetTrackFilter().AddCuts(LMCL->GetPIDCutsPre(selectedPID) );
    die->GetTrackFilter().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
    die->GetPairPreFilterLegs().AddCuts(LMCL->GetPIDCutsAna(selectedPID) );
    die->GetPairPreFilter().AddCuts(LMCL->GetPairCuts(selectedPID) );

    //    if(PairCut){
    //      die->GetPairFilter().AddCuts( LMCL->GetPairCutsInvMass(selectedPairCut));
    //    }
  }
  else { //No Prefilter, no Pairfilter
    
    if (ESDanalysis) {
      die->GetTrackFilter().AddCuts( LMCL->GetESDTrackCutsAna(selectedPID) );
    }
    
    die->GetTrackFilter().AddCuts( LMCL->GetTrackCutsAna(selectedPID) );
    die->GetTrackFilter().AddCuts( LMCL->GetPIDCutsAna(selectedPID) );
    die->GetEventFilter().AddCuts(LMCL->GetCentralityCuts(selectedCentrality));
    
    //    if(PairCut){
    //      die->GetPairFilter().AddCuts( LMCL->GetPairCutsInvMass(selectedPairCut));
    //    }
    die->GetPairFilter().AddCuts(LMCL->GetPairCuts2(selectedPID,kFALSE));
    
  }
  //Introduce NULL-check for pp?
  die->GetEventFilter().AddCuts(LMCL->GetCentralityCuts(selectedCentrality));
  */



  AliDielectronTrackRotator *rot= 0x0;
  /*AliDielectronTrackRotator *rot= LMCL->GetTrackRotator(selectedPID);
  die->SetTrackRotator(rot);
   */
  AliDielectronMixingHandler *mix=LMCL->GetMixingHandler(selectedPID);
  die->SetMixingHandler(mix);

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);

  // the last definition uses no cuts and only the QA histograms should be filled!
//  InitCF(die,cutDefinition);

  return die;
}

//______________________________________________________________________________________

void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
  AliDielectronHistos *histos=
	new AliDielectronHistos(die->GetName(),
		die->GetTitle());
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Pre;RejTrack;RejPair");

  //Event class
//  if (cutDefinition==nDie-1) 
	  	histos->AddClass("Event");

  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
	histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }

  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }

  //ME and track rot
  if (die->GetMixingHandler()) {
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
  }
  if (die->GetTrackRotator()) {
	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
  }

  //PreFilter Classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
	histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
  }


  //Create Classes for Rejected Tracks/Pairs:
  for (Int_t i=0; i<2; ++i){
	histos->AddClass(Form("RejTrack_%s",AliDielectron::TrackClassName(i)));
  }
  for (Int_t i=0; i<3; ++i){
	histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
  }

  /*
  //track rotation

  histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
  */
	//add histograms to event class
	histos->UserHistogram("Event","nEvents","Number of processed events after cuts;Number events",
		1,0.,1.,AliDielectronVarManager::kNevents);
	histos->UserHistogram("Event","Centrality","Centrality;Centrality [%]","0,10,20,40,80,100,101",
		AliDielectronVarManager::kCentrality);


	histos->UserHistogram("Event","v0ACrpH2","VZERO-AC;v0ACrpH2",
			      100,-2.0,2.0,
			      AliDielectronVarManager::kv0ACrpH2);
	histos->UserHistogram("Event","v0ArpH2","VZERO-A;v0ArpH2",
			      100,-2.0,2.0,
			      AliDielectronVarManager::kv0ArpH2);
	histos->UserHistogram("Event","v0CrpH2","VZERO-C;v0CrpH2",
			      100,-2.0,2.0,
			      AliDielectronVarManager::kv0CrpH2);
	histos->UserHistogram("Event","RadomRP","RandomRP;RandomRP",
			      100,-2.0,2.0,
			      AliDielectronVarManager::kRandomRP);


	histos->UserHistogram("Event","ZDCArpH1","ZDC-ZN-A;ZDCrpH1",
                              100,-3.5,3.5,
                              AliDielectronVarManager::kZDCACrpH1);

        histos->UserProfile("Event","ZDCrpResH1Prof","ZDC;ZDCrpResH1",
			    AliDielectronVarManager::kZDCrpResH1, 
			    10, 0, 100,
                              AliDielectronVarManager::kCentrality);
    
        histos->UserProfile("Event","v0ZDCrpResProf","ZDC;v0ZDCrpRes",
			     AliDielectronVarManager::kv0ZDCrpRes,
			     10, 0, 100,
			    AliDielectronVarManager::kCentrality);
	
	
	histos->UserHistogram("Event","RefMult","RefMultiplicity;Multiplixity",
			      100,-3.5,3.5,
			      AliDielectronVarManager::kRefMult);
	histos->UserHistogram("Event","kXvPrim","VertexX;vertex_x",
			      100,0.03,0.1,
			      AliDielectronVarManager::kXvPrim);
	histos->UserHistogram("Event","kYvPrim","VartexY;vertex_y",
			      100,0.2,0.3,
			      AliDielectronVarManager::kYvPrim);


	


  //add histograms to Track classes

        histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);
        histos->UserHistogram("Track","Px","Px;Px [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPx);
        histos->UserHistogram("Track","Py","Py;Py [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPy);
        histos->UserHistogram("Track","Pz","Pz;Pz [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPz);


        histos->UserHistogram("Track","Eta","Eta; Eta;#tracks",
                              200,-2,2,AliDielectronVarManager::kEta);
        histos->UserHistogram("Track","Phi","Phi; Phi;#tracks",
                              200,0.,3.15,AliDielectronVarManager::kPhi);


  /*
  histos->UserHistogram("Track","NclsSFracTPC","NclsSFracTPC; NclsSFracTPC;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff","TPCclsDiff; TPCclsDiff;#tracks",200,0,10.,AliDielectronVarManager::kTPCclsDiff);

  histos->UserHistogram("Track","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
	  400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
	  400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
	  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
	  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
	  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);

  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons;#tracks",
	  400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,kTRUE);
  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions;#tracks",
	  400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,kTRUE);

  histos->UserHistogram("Track","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
	  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
  histos->UserHistogram("Track","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
	  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);

  histos->UserHistogram("Track","TOFbeta","TOF beta;P [GeV];TOF beta;#tracks",
	  400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);


  histos->UserHistogram("Track","Eta","Eta; Eta;#tracks",
	  200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","Phi; Phi;#tracks",
	  200,0.,3.15,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
	  200,-2,2,200,0,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  histos->UserHistogram("Track","dXY_dZ","dXY dZ Map; dXY; dZ;#tracks",
	  200,-2,2,200,-2,2.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);


  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParZ);

	  histos->UserHistogram("Track","TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",100,0.,1.,AliDielectronVarManager::kNFclsTPCfCross);
	  histos->UserHistogram("Track","TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",159,0.,159.,AliDielectronVarManager::kNFclsTPCr);
	  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);
	  histos->UserHistogram("Track","ITSnCls","Number of Clusters ITS;ITS number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsITS);

	  histos->UserHistogram("Track","TPCchi2","TPC Chi2 value;TPC chi2;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
	  histos->UserHistogram("Track","ITSchi2","ITS Chi2 value;ITS chi2;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);

	  histos->UserHistogram("Track","TPCnCls_kNFclsTPCr","nTPC vs nTPCr;nTPC vs nTPCr;#tracks",159,0.,159.,159,0.,159.,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);

	  histos->UserHistogram("Track","kNFclsTPCr_pT","nTPCr vs pt;nTPCr vs pt;#tracks",159,0.,159.,200,0.,20.,AliDielectronVarManager::kNFclsTPCr,AliDielectronVarManager::kPt);

  */
	  //add histograms to Pair classes
	  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
		  500,0.0,5.00,AliDielectronVarManager::kM);
	  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
		  100,-2.,2.,AliDielectronVarManager::kY);
	  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
		  100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
	  //2D Histo Plot
	  histos->UserHistogram("Pair","InvMassPairPt","Inv.Mass vs PairPt;Inv. Mass [GeV], pT [GeV];#pairs",
		  500,0.0,5.0,500,0.,50.,AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

	  histos->UserHistogram("Pair","InvMassOpeningAngle","Opening Angle vs Inv.Mass;Inv. Mass [GeV];#pairs",
		  500,0.0,5.0,200,0.,6.3,AliDielectronVarManager::kM,AliDielectronVarManager::kOpeningAngle);


	  histos->UserHistogram("Pair","Pt","Pt;Pt [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPt);
	  histos->UserHistogram("Pair","Px","Px;Px [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPx);
	  histos->UserHistogram("Pair","Py","Py;Py [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPy);
	  histos->UserHistogram("Pair","Pz","Pz;Pz [GeV];#tracks",300,0,30.,AliDielectronVarManager::kPz);
	  histos->UserHistogram("Pair","Phi","Phi;Phi[rad];#counts",100,-3.15,3.15,AliDielectronVarManager::kPhi );


	  histos->UserHistogram("Pair","DeltaPhiv0ArpH2","Phi;Phi[rad];#counts",
				100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0ArpH2);
	  histos->UserHistogram("Pair","DeltaPhiv0CrpH2","Phi;Phi[rad];#counts",
				100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0CrpH2);
	  histos->UserHistogram("Pair","DeltaPhiv0ACrpH2","Phi;Phi[rad];#counts",
				100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0ACrpH2);
	  histos->UserHistogram("Pair","DeltaPhiRandomRP","Phi;Phi[rad];#counts",
				100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiRandomRP);


	  histos->UserHistogram("Pair","PairPlaneAngle2C","Phi;Phi[rad];#counts",
				100,0,1.6,AliDielectronVarManager::kPairPlaneAngle2C);
	  histos->UserHistogram("Pair","PairPlaneAngle3C","Phi;Phi[rad];#counts",
				100,0,1.6,AliDielectronVarManager::kPairPlaneAngle3C);
	  histos->UserHistogram("Pair","PairPlaneAngle4C","Phi;Phi[rad];#counts",
				100,0,1.6,AliDielectronVarManager::kPairPlaneAngle4C);
	  histos->UserHistogram("Pair","PairPlaneAngleRan","Phi;Phi[rad];#counts",
				100,0,1.6,AliDielectronVarManager::kPairPlaneAngle3Ran);


	  //2D Histo Plot

	  histos->UserHistogram("Pair","InvMAllPP1C","Inv.Mass vs PairPlaneAngle;Inv. Mass [GeV];Phi [rad]",500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle1C);

	  histos->UserHistogram("Pair","InvMAllPP2C","Inv.Mass vs PairPlaneAngle;Inv. Mass [GeV];Phi [rad]",500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle2C);

	  histos->UserHistogram("Pair","InvMAllPP3C","Inv.Mass vs PairPlaneAngle;Inv. Mass [GeV];Phi [rad]",500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle3C);

	  histos->UserHistogram("Pair","InvMAllPP4C","Inv.Mass vs PairPlaneAngle;Inv. Mass [GeV];Phi [rad]",500,0.0,0.50,100,0.,3.15,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneAngle4C);


          histos->UserHistogram("Pair","PtAllPP1C","Pair Pt vs PairPlaneAngle;Pt [GeV];Phi [rad]",
				 500,0.,10.0,100,0.,3.15,AliDielectronVarManager::kPt, AliDielectronVarManager::kPairPlaneAngle1C);

	  histos->UserHistogram("Pair","PtAllPP2C","Pair Pt vs PairPlaneAngle;Pt [GeV];Phi [rad]",
				 500,0.,10.0,100,0.,3.15,AliDielectronVarManager::kPt, AliDielectronVarManager::kPairPlaneAngle2C);

	  histos->UserHistogram("Pair","PtAllPP3C","Pair Pt vs PairPlaneAngle;Pt [GeV];Phi [rad]",
				 500,0.,10.0,100,0.,3.15,AliDielectronVarManager::kPt, AliDielectronVarManager::kPairPlaneAngle3C);

	  histos->UserHistogram("Pair","PtAllPP4C","Pair Pt vs PairPlaneAngle;Pt [GeV];Phi [rad]",
				 500,0.,10.0,100,0.,3.15,AliDielectronVarManager::kPt, AliDielectronVarManager::kPairPlaneAngle4C);

	  /*
        histos->UserHistogram("Pair","InvMassAllPairplaneMagInPro","Inner Product of Mag and ee plane vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                              1000, 0.0,1.0,100,-2.0,2.0,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagInPro);
        histos->UserHistogram("Pair","InvMassLowPairplaneMagInPro","ee plane Mag component vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                              300, 0.0,0.03,100,-2.0,2.0,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagInPro);
        histos->UserHistogram("Pair","InvMassMiddlePairplaneMagInPro","ee plane Mag component vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                               180,0.12, 0.3,100,-2.0,2.0,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagInPro);
        histos->UserHistogram("Pair","InvMassHighPairplaneMagInPro","ee plane Mag component vs Inv.Mass;Inv. Mass [GeV];Phi [rad]",
                               200, 0.3, 0.5,100,-2.0,2.0,AliDielectronVarManager::kM,AliDielectronVarManager::kPairPlaneMagInPro);
	  */

	  histos->UserHistogram("Pair","DeltaPhiv0CrpH2","Phi;Phi[rad];#counts",
				100,-3.15,3.15,AliDielectronVarManager::kDeltaPhiv0CrpH2);

	  histos->UserHistogram("Pair","PtAllPairplaneMagInPro","ee plane Mag component vs Pt;Pt [GeV];Phi [rad]",
				500,0.0,10.0,100,-2.0,2.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);
	  histos->UserHistogram("Pair","PtLowPairplaneMagInPro","ee plane Mag component vs Pt;Pt [GeV];Phi [rad]",
				100,0.0,1.0,100,-2.0,2.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);
	  histos->UserHistogram("Pair","PtMiddlePairplaneMagInPro","ee plane Mag component vs Pt;Pt [GeV];Phi [rad]",
				100,1.0,2.0,100,-2.0,2.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);
	  histos->UserHistogram("Pair","PtHighPairplaneMagInPro","ee plane Mag component vs Pt;Pt [GeV];Phi [rad]",
				200,2.0,10.0,100,-2.0,2.0,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);
	  

	  /*
                        histos->UserHistogram("Pair","AllInvMassPtPairplaneMagInPro","ee plane Mag component;Inv.Mass[GeV];Pt[GeV];Phi[red]",
                              1000,0.0 ,0.5,500,0.0,10.0,100,-2.0,2.0,
                                              AliDielectronVarManager::kM,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);
	  */


	  //              histos->UserHistogram("Pair","AllInvMassPtPairplaneMagInPro","ee plane Mag component vs Pt;Inv.Mass[GeV];Pt [GeV];Phi [rad]",
	  //                                    1000,0.0,0.5,500,0.0,10.0,100,-2.0,2.0,
	  //                                    AliDielectronVarManager::kM,AliDielectronVarManager::kPt,AliDielectronVarManager::kPairPlaneMagInPro);


	  /*


	  //add histograms to Track classes
	  histos->UserHistogram("Pre","Pt","Pt;Pt [GeV];#tracks",200,0,20.,AliDielectronVarManager::kPt);

	  histos->UserHistogram("Pre","ITS_dEdx_P","ITS_dEdx;P [GeV];ITS signal (arb units);#tracks",
		  400,0.0,20.,1000,0.,1000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kITSsignal,kTRUE);

	  histos->UserHistogram("Pre","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
		  400,0.0,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);


	  histos->UserHistogram("Pre","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
		  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
	  histos->UserHistogram("Pre","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
		  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,kTRUE);
	  histos->UserHistogram("Pre","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
		  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,kTRUE);

	  histos->UserHistogram("Pre","TRDpidPobEle_P","TRD PID probability Electrons;P [GeV];TRD prob Electrons;#tracks",
		  400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle,kTRUE);
	  histos->UserHistogram("Pre","TRDpidPobPio_P","TRD PID probability Pions;P [GeV];TRD prob Pions;#tracks",
		  400,0.0,20.,100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio,kTRUE);

	  histos->UserHistogram("Pre","TOFnSigmaKao_P","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
		  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao,kTRUE);
	  histos->UserHistogram("Pre","TOFnSigmaPro_P","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
		  400,0.0,20.,100,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro,kTRUE);

	  histos->UserHistogram("Pre","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
		  200,-2,2,200,0,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

	  histos->UserHistogram("Pre","dXY","dXY;dXY [cm];#tracks",200,-2.,2.,AliDielectronVarManager::kImpactParXY);

  histos->UserHistogram("Pre","ZVertex ","ZVertex ;ZVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kZv);
  histos->UserHistogram("Pre","XVertex ","XVertex ;XVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kXv);
  histos->UserHistogram("Pre","YVertex ","YVertex ;YVertex[cm];#tracks",20,-20,20,AliDielectronVarManager::kYv);

  histos->UserHistogram("Pre","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",159,0.,159.,AliDielectronVarManager::kNclsTPC);
	  */
  //add histograms to Pair classes For Rejected Pairs:
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,200,0,20);
  cf->AddVariable(AliDielectronVarManager::kM,201,-0.01,4.01); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);

  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,10.0,30.0,40.0,60.,80.,100.");

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,200,0.,20.,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kITSsignal,1000,0.0.,1000.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,500,0.0.,500.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);

  //only in this case write MC truth info
  if (MCenabled) {
	cf->SetStepForMCtruth();
	cf->SetStepsForMCtruthOnly();
	cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
	cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}

//--------------------------------------
void EnableMC() {
  MCenabled=kTRUE;
}


//  LocalWords:  cutDefinition
