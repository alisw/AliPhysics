//Based on AddMCGen_Pythia8_MB_13TeV_Monash_NoCR.C
AliGenerator* CreatePythia8Gen(Float_t e_cms,       
			       Int_t tune    ,      
			       Bool_t kCR     ,     
			       Int_t kF        ,    
			       Int_t kProcess   ,   
			       Double_t ptHardMin,  
			       Double_t ptHardMax , 
			       Bool_t longlived ) {
  
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libpythia8");
  gSystem->Load("libAliPythia8");
  gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
  gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
  gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

  AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
  AliGenPythiaPlus* gener = new AliGenPythiaPlus(pythia);

  //Specify process
  if(kProcess==0){
    (AliPythia8::Instance())->ReadString("SoftQCD:inelastic = on");
  }else if(kProcess==1){
    //    (AliPythia8::Instance())->ReadString("PromptPhoton:all = on");
    (AliPythia8::Instance())->ReadString("PromptPhoton:qg2qgamma = on");//compton scattering
    (AliPythia8::Instance())->ReadString("PromptPhoton:qqbar2ggamma = on");//annihilation
  }else if(kProcess==2){//SoftQCD:elastic = on
    (AliPythia8::Instance())->ReadString("SoftQCD:inelastic = off");
    (AliPythia8::Instance())->ReadString("PromptPhoton:qg2qgamma = off");//compton scattering
    (AliPythia8::Instance())->ReadString("PromptPhoton:qqbar2ggamma = off");//annihilation
  }
  
  //Center of Mass energy
  //(AliPythia8::Instance())->ReadString("Beams:eCM = 13000.");//does not work
  gener->SetEnergyCMS(e_cms);//use this instead
  //Pythia8 Monash 2013 tune 
  (AliPythia8::Instance())->ReadString(Form("Tune:pp = %i",tune));

  //random seed based on time
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString("Random:seed = 0");

  //Color Reconnection
  if(kCR)             
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
  else
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        
  (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));

  return gener;
}

AliGenerator* AddMCGen_Pythia8_13TeV_Monash(Float_t e_cms       = 13000., 
					    Int_t tune          = 14, 
					    Bool_t kCR          = kFALSE, 
					    Int_t kF            = 1, 
					    Int_t kProcess      = 0, 
					    Double_t ptHardMin  = 0, 
					    Double_t ptHardMax  = 1.,
					    Bool_t longlived = kFALSE
) {

  gSystem->Load("liblhapdf");
  AliGenerator *genP = CreatePythia8Gen(e_cms, tune, kCR, kF, kProcess, ptHardMin, ptHardMax,longlived);  
  return genP;

}

  //Joshua's first shot
  /* pythia.readString("SoftQCD:inelastic = on");        ... kPyMbDefault */
  /* pythia.readString("Beams:eCM = 13000.");            ... set e_cms */
  /* pythia.readString("Tune:ee = 7");                   ... ??  */
  /* pythia.readString("Tune:pp = 14");                  ... kPyMbDefault */
  /* pythia.readString("ParticleDecays:limitTau0 = on"); ... not touched  */
  /* pythia.readString("ParticleDecays:tau0Max = 10");   ... ??  */
  /* pythia.readString("Random:setSeed = on");           ... ok */
  /* pythia.readString("Random:seed = 0");               ... ok  */
  /* pythia.readString("PromptPhoton:all = on");         ... --> need to set manually */

  //Should not be relevant 
  //  (AliPythia8::Instance())->ReadString("Tune:ee = 7");
  //Joshua's comment on this
  //"This  should not matter as this is the tune for e+e- collisions.
  //I think I included it as it is the default Monash 2013 tune as stated here: https://arxiv.org/pdf/1404.5630.pdf"

  //Should not be relevant 
  //  (AliPythia8::Instance())->ReadString("ParticleDecays:limitTau0 = on");
  //  (AliPythia8::Instance())->ReadString("ParticleDecays:tau0Max = 10");
  //Joshua's comment on this
  //"limit the particle decays in particular for strange decays.
  //It specifies that particles with a certain decay lenght (I now set it to 10mm) will not decay
  //(see https://pythia.org/latest-manual/ParticleDecays.html) So this should not be relevant for looking at direct photons."

