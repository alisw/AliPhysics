
AliGenerator* AddMCGenPythia8_ModSettings(Float_t e_cms= 13000.
					  , Bool_t kCR= kTRUE
					  , Int_t kF= 1
					  , Int_t kProcess= 0 
					  , Double_t ptHardMin= 0.
					  , Double_t ptHardMax= 1.
					  , Double_t ihalfMassForKT= 0.0
					  , Double_t ihalfScaleForKT= 0.0
					  , Double_t iprimordialKThard= 1.0
					  , Double_t iprimordialKTsoft= 0.0
					  ){
  // Add Pythia 8 generator: 
  // -- kProcess=0  MB generation
  //-- kProcess=1  Jet production, pthard generation
  //   -- Color reconnection = ON/OFF
  //  -- Set k factor, default = 1; range of possible values in xmldoc/CouplingsAndScales.xml

  gSystem->Load("liblhapdf");
    
  AliGenerator *myGen  = NULL;
  myGen                = CreatePythia8Gen(e_cms, kCR, kF, kProcess, ptHardMin, ptHardMax, ihalfMassForKT, ihalfScaleForKT, iprimordialKThard, iprimordialKTsoft);
    
  return myGen;
}

AliGenerator* CreatePythia8Gen( Float_t e_cms 
				, Bool_t kCR 
				, Int_t kF
				, Int_t kProcess
				, Double_t ptHardMin
				, Double_t ptHardMax
				, Double_t ihalfMassForKT
				, Double_t ihalfScaleForKT
				, Double_t iprimordialKThard
				, Double_t iprimordialKTsoft
				){
    
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libpythia8");
  gSystem->Load("libAliPythia8");
  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));


  AliGenPythiaPlus* myGener = new AliGenPythiaPlus(AliPythia8::Instance());

  if(kProcess==0)
    myGener->SetProcess(kPyMbDefault);

  if(kProcess==1){
    myGener->SetProcess(kPyJets);
    if(ptHardMin > 0.)
      myGener->SetPtHard(ptHardMin,ptHardMax);
  }

  //Center-of-mass energy 
  myGener->SetEnergyCMS(e_cms); // in GeV

  //Event list
  myGener->SetEventListRange(-1, -1);

  //Set tune
  // ... default (Monash 2013)

  //Set random seed (based on time)
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString("Random:seed = 0");

  //Set Color reconnection: ON / OFF
  if(kCR)
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
  else
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        
  //kFactor
  (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));


  //Additional settings

  (AliPythia8::Instance())->ReadString(Form("BeamRemnants:halfMassForKT = %f",ihalfMassForKT));
  (AliPythia8::Instance())->ReadString(Form("BeamRemnants:halfScaleForKT = %f",ihalfScaleForKT));
  (AliPythia8::Instance())->ReadString(Form("BeamRemnants:primordialKThard = %f",iprimordialKThard));
  (AliPythia8::Instance())->ReadString(Form("BeamRemnants:primordialKTsoft = %f",iprimordialKTsoft));    
    
/*    // for later use

    for( Long_t iCfs=0; iCfs < lNumConfStrings; iCfs++)
      myGener->ReadString(lConfigStrings[iCfs].Data());
*/
    
    return myGener;
}

