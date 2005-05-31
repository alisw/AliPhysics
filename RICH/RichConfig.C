class RichConfig : public TGMainFrame
{
RQ_OBJECT()
public:
 
          RichConfig(const char*sFileName);
         ~RichConfig()                    {Info("ctor","");Cleanup();}
//protected:
  enum EVersOpts  {kNo=-1,kVer0,kVer1,kTestBeam,kTestRadio,kOnTop,kDeclust=301,kSagita,kFeedback};
  enum EGenTypes  {kGunAlongZ=100,kBox1,kBox7,kGun7,kPythia7,kHijing,kHijingPara,kRichLib,kSignalHijing,kHijingPara2Proton};
  
  enum EDetectors {kPIPE=1,kITS,kTPC,kTRD,kTOF,kFRAME,kMAG,kCRT,kHALL,kPHOS,kSTART,kFMD,kABSO,kPMD,kDIPO,kEMCAL,kVZERO,kMUON,kZDC,kSHILD,kRICH};
  enum EProcesses {kDCAY=1,kPAIR,kCOMP,kPHOT,kPFIS,kDRAY,kANNI,kBREM,kMUNU,kCKOV,kHADR,kLOSS,kMULS,kRAYL,kALL};
  enum EBatchFlags{kNoRaw=0,kRawDdl,kRawDate,kRawRoot,kRecon,kFillEsd};
//SLOTS  
  void    AddDetector(Int_t id)          {if(id==kTRD || id==kTOF) fDetBG->SetButton(kFRAME);}
  void    RemoveDetector(Int_t id)       {if(id==kFRAME) {fDetBG->SetButton(kTRD,kFALSE);fDetBG->SetButton(kTOF,kFALSE);}}
  
  void    AddProcess(Int_t id)           {if(id==kALL) for(int i=1;i<=fProcBG->GetCount();i++) fProcBG->SetButton(i,kButtonDown);}
  void    RemoveProcess(Int_t id)        {if(id==kALL) for(int i=1;i<=fProcBG->GetCount();i++) fProcBG->SetButton(i,kButtonUp);}
  void    GeneratorSlot(Int_t type);
  void    VerSlot(Int_t ver)             {if(ver==kTestBeam){ fMagCB->SetState(kButtonUp); fGenTypeCO->Select(kGunAlongZ);}}        
  
  Float_t Eta2Theta(Float_t arg)    const{return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));}
  void    CreateConfig();
  void    CreateRichBatch();
  void    Exit();
     
  TGComboBox   *fRichCO;  TGButtonGroup *fRichBG;                                                                              //RICH
  TGButton     *fMagCB;                                                                                                     //MAG
  TGComboBox   *fGenTypeCO,*fGenPidCO,*fGenMinMomCO,*fGenMaxMomCO,*fGenChamberCO; TGGroupFrame *fGenGF,*fGenOptGF,*fGenPrimGF; TGNumberEntry *fGenPrimNE;             //GEN
  TGButtonGroup *fDetBG;                                                                                                       //DETECTORS
  TGButtonGroup *fProcBG;                                                                                                      //PROCESSES
  TGButtonGroup *fRawBG;                                                                                                       //RAW 
  TGButtonGroup *fRecoBG;                                                                                                      //RAW   
  char         *fFileName;
};//class RichConfig
//__________________________________________________________________________________________________	 
RichConfig::RichConfig(const char *sFileName):TGMainFrame(gClient->GetRoot(),700,400)
{// creates configuration file  
  fFileName=sFileName;

//  --------------------  
//  | |--------------------
//  | ||-------|
//  | ||       |
//  | ||-------|
//  | |  
//  | ||-------|   
  
//  | |--------------------
//  | ---------------------         
  
//  Connect("CloseWindow()","RichConfig",this,"CloseWindow()");   
  AddFrame(pMainHF=new TGHorizontalFrame(this,100,200)); //main horizontal frame
//Left vertical frame containing RICH and MagField  
  pMainHF->AddFrame(pLeftVF=new TGVerticalFrame(pMainHF,100,200));  
//RICH
  pLeftVF->AddFrame(pRichGF=new TGGroupFrame(pLeftVF,"RICH"));
  pRichGF->AddFrame(fRichCO=new TGComboBox(pRichGF,100));
  fRichCO->Connect("Selected(Int_t)" ,"RichConfig",this,"VerSlot(Int_t)");
    fRichCO->AddEntry("no RICH"      ,kNo);
    fRichCO->AddEntry("debug RICH 0" ,kVer0);
    fRichCO->AddEntry("normal RICH 1",kVer1);
    fRichCO->AddEntry("RICH on top"  ,kOnTop);
    fRichCO->AddEntry("test beam"    ,kTestBeam);
    fRichCO->AddEntry("radio source" ,kTestRadio);
    fRichCO->Select(kVer1); 
    fRichCO->Resize(150,20);
  pRichGF->AddFrame(fRichBG=new TGButtonGroup(pRichGF,"Options"));
    new TGCheckButton(fRichBG,"Decluster?"   ,kDeclust);      fRichBG->SetButton(kDeclust);
    new TGCheckButton(fRichBG,"Wire sagita?" ,kSagita);       fRichBG->SetButton(kSagita);
    new TGCheckButton(fRichBG,"Feedbacks?"   ,kFeedback);     fRichBG->SetButton(kFeedback); 
// Magnetic Field
  pLeftVF->AddFrame(pMagGF=new TGGroupFrame(pLeftVF,"Magnetic Field"));
  pMagGF->AddFrame(fMagCB=new TGCheckButton(pMagGF,"On/Off"));
  fMagCB->SetState(kButtonDown);
//Generator  
  pMainHF->AddFrame(fGenGF=new TGGroupFrame(pMainHF,"Generator"));
  fGenGF->AddFrame(fGenTypeCO=new TGComboBox(fGenGF,100));
  fGenTypeCO->Connect("Selected(Int_t)" ,"RichConfig",this,"GeneratorSlot(Int_t)");
  fGenTypeCO->AddEntry("gun along z"                  ,kGunAlongZ);  
  fGenTypeCO->AddEntry("box gun to single chamber"    ,kBox1);  
  fGenTypeCO->AddEntry("gun to all chambers"          ,kGun7);
  fGenTypeCO->AddEntry("box to all chambers"          ,kBox7);
  fGenTypeCO->AddEntry("7 guns+Pythia"                ,kPythia7);
  fGenTypeCO->AddEntry("HIJING"                       ,kHijing);
  fGenTypeCO->AddEntry("HIJING para"                  ,kHijingPara);
  fGenTypeCO->AddEntry("2 p+HIJING"                   ,kHijingPara2Proton);
  fGenTypeCO->AddEntry("Sr90 source"                  ,kSr90);
  fGenTypeCO->AddEntry("RICH lib"                     ,kRichLib);
  fGenTypeCO->AddEntry("RICH lib+HIJING"              ,kSignalHijing);
  fGenTypeCO->Select(kHijing);
  fGenTypeCO->Resize(160,20);
  
  
  fGenGF->AddFrame(fGenOptGF=new TGGroupFrame(fGenGF,"Options")); //PID for guns

  
  fGenOptGF->AddFrame(fGenPidCO=new TGComboBox(fGenOptGF,100)); //PID for guns
  fGenPidCO->AddEntry("Pion+"     ,kPiPlus);
  fGenPidCO->AddEntry("Pion-"     ,kPiMinus);
  fGenPidCO->AddEntry("Kaon+"     ,kKPlus);
  fGenPidCO->AddEntry("Kaon-"     ,kKMinus);
  fGenPidCO->AddEntry("K0s"       ,kK0Short);    
  fGenPidCO->AddEntry("Proton"    ,kProton);
  fGenPidCO->AddEntry("ProtonBar" ,kProtonBar);
  fGenPidCO->AddEntry("Lambda"    ,kLambda0);
  fGenPidCO->AddEntry("LambdaBar" ,kLambda0Bar);
  fGenPidCO->Select(kProton);  fGenPidCO->Resize(160,20);

  fGenOptGF->AddFrame(fGenMinMomCO=new TGComboBox(fGenOptGF,100)); //particle energy for guns  
  for(Int_t i=5;i<=95;i+=5)
    fGenMinMomCO->AddEntry(Form("%3.1f GeV",0.1*i), i);
  fGenMinMomCO->Select(15);  fGenMinMomCO->Resize(160,20);
  
  fGenOptGF->AddFrame(fGenMaxMomCO=new TGComboBox(fGenOptGF,100)); //particle energy for guns
  for(Int_t i=10;i<=95;i+=5)
    fGenMaxMomCO->AddEntry(Form("%3.1f GeV",0.1*i), i);
  fGenMaxMomCO->Select(40);  fGenMaxMomCO->Resize(160,20);
  
  fGenOptGF->AddFrame(fGenChamberCO=new TGComboBox(fGenOptGF,100)); //chamber number in case of gun1
  for(int i=1;i<=7;i++) fGenChamberCO->AddEntry(Form("Chamber %i",i),i);
  fGenChamberCO->Select(4); fGenChamberCO->Resize(160,20);
  
  fGenOptGF->AddFrame(fGenPrimGF=new TGGroupFrame(fGenOptGF,"Number of primaries"));
  fGenPrimGF->AddFrame(fGenPrimNE=new TGNumberEntry(fGenPrimGF,500));  
//Detectors
  pMainHF->AddFrame(fDetBG=new TGButtonGroup(pMainHF,"Detectors"));
  fDetBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"AddDetector(Int_t)");
  fDetBG->Connect("Released(Int_t)","RichConfig",this,"RemoveDetector(Int_t)");
    new TGCheckButton(fDetBG,"PIPE"  ,kPIPE));         
    new TGCheckButton(fDetBG,"ITS"   ,kITS));         
    new TGCheckButton(fDetBG,"TPC"   ,kTPC));
    new TGCheckButton(fDetBG,"TRD"   ,kTRD));
    new TGCheckButton(fDetBG,"TOF"   ,kTOF));         
    new TGCheckButton(fDetBG,"FRAME" ,kFRAME));         
    new TGCheckButton(fDetBG,"MAG"   ,kMAG));         
    new TGCheckButton(fDetBG,"CRT"   ,kCRT));         
    new TGCheckButton(fDetBG,"HALL"  ,kHALL));         
    new TGCheckButton(fDetBG,"PHOS"  ,kPHOS));         
    new TGCheckButton(fDetBG,"START" ,kSTART));         
    new TGCheckButton(fDetBG,"FMD"   ,kFMD));         
    new TGCheckButton(fDetBG,"ABSO"  ,kABSO));         
    new TGCheckButton(fDetBG,"PMD"   ,kPMD));         
    new TGCheckButton(fDetBG,"DIPO"  ,kDIPO));         
    new TGCheckButton(fDetBG,"EMCAL" ,kEMCAL));         
    new TGCheckButton(fDetBG,"VZERO" ,kVZERO));         
    new TGCheckButton(fDetBG,"MUON"  ,kMUON));         
    new TGCheckButton(fDetBG,"ZDC"   ,kZDC));         
    new TGCheckButton(fDetBG,"SHILD" ,kSHILD));         
//Processes  
  pMainHF->AddFrame(fProcBG=new TGButtonGroup(pMainHF,"Processes"));
  fProcBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"AddProcess(Int_t)");
  fProcBG->Connect("Released(Int_t)","RichConfig",this,"RemoveProcess(Int_t)");
    new TGCheckButton(fProcBG,"ALL  ON/OFF"                 ,kALL)) ;
    new TGCheckButton(fProcBG,"DCAY Decay"                  ,kDCAY));  fProcBG->SetButton(kDCAY);       
    new TGCheckButton(fProcBG,"PAIR Pair production"        ,kPAIR));  fProcBG->SetButton(kPAIR);       
    new TGCheckButton(fProcBG,"COMP Compton"                ,kCOMP));  fProcBG->SetButton(kCOMP);
    new TGCheckButton(fProcBG,"PHOT"                        ,kPHOT));  fProcBG->SetButton(kPHOT);
    new TGCheckButton(fProcBG,"PFIS Photofission"           ,kPFIS));  
    new TGCheckButton(fProcBG,"DRAY Delta electrons"        ,kDRAY));  
    new TGCheckButton(fProcBG,"ANNI Annihilation"           ,kANNI));  fProcBG->SetButton(kANNI);       
    new TGCheckButton(fProcBG,"BREM Bremstraslung"          ,kBREM));  fProcBG->SetButton(kBREM);       
    new TGCheckButton(fProcBG,"MUNU Muon-Nuclear"           ,kMUNU));  fProcBG->SetButton(kMUNU);       
    new TGCheckButton(fProcBG,"CKOV Cerenkovs"              ,kCKOV));  fProcBG->SetButton(kCKOV);       
    new TGCheckButton(fProcBG,"HADR Hadronic interactions " ,kHADR));  fProcBG->SetButton(kHADR);       
    new TGCheckButton(fProcBG,"LOSS Energy losses"          ,kLOSS));  fProcBG->SetButton(kLOSS);       
    new TGCheckButton(fProcBG,"MULS Multiple scattering"    ,kMULS));  fProcBG->SetButton(kMULS);       
    new TGCheckButton(fProcBG,"RAYL"                        ,kRAYL));  fProcBG->SetButton(kRAYL);       
//RichBatch part
  pMainHF->AddFrame(pBatchVF=new TGVerticalFrame(pMainHF));    
  pBatchVF->AddFrame(fRawBG=new TGButtonGroup(pBatchVF,"RichBatch"));
    new TGRadioButton(fRawBG,"No RAW files"      ,kNoRaw)) ;   fRawBG->SetButton(kNoRaw);
    new TGRadioButton(fRawBG,"RAW DDL"           ,kRawDdl)) ;
    new TGRadioButton(fRawBG,"RAW DATE"          ,kRawDate)) ;
    new TGRadioButton(fRawBG,"RAW ROOT"          ,kRawRoot)) ;
  fRawBG->AddFrame(fRecoBG=new TGButtonGroup(fRawBG,""));
    new TGCheckButton(fRecoBG,"Reconstruct"    ,kRecon)) ; fRecoBG->SetButton(kRecon);  
    new TGCheckButton(fRecoBG,"Fill ESD"       ,kFillEsd)) ; fRecoBG->SetButton(kFillEsd);
//File    
  AddFrame(pFileHF=new TGHorizontalFrame(this,100,200));
  pFileHF->AddFrame(pCrtTB=new TGTextButton(pFileHF,"Create"));
  pCrtTB->Connect("Clicked()","RichConfig",this,"Exit()");                                 
  pFileHF->AddFrame(new TGLabel(pFileHF,Form(" config file as %s",fFileName)));  
  
  MapSubwindows();   
  //Layout(); 
  Resize(GetDefaultSize());
  SetWindowName("Create AliROOT scripts");
  MapWindow();      
}//KirCondig::ctor
//__________________________________________________________________________________________________
void RichConfig::GeneratorSlot(Int_t type)
{//enable-disable generator widegets according to generator type
  if(type==kHijing) fGenGF->HideFrame(fGenOptGF); 
  
  if(type==kHijingPara){
    fGenGF->ShowFrame(fGenOptGF);
    fGenOptGF->ShowFrame(fGenPrimGF);        
    fGenOptGF->HideFrame(fGenPidCO);
    fGenOptGF->HideFrame(fGenMinMomCO);
    fGenOptGF->HideFrame(fGenMaxMomCO);
    fGenOptGF->HideFrame(fGenChamberCO);
  }
  
  if(type==kGunAlongZ) {
    fGenGF->ShowFrame(fGenOptGF);
    fGenOptGF->ShowFrame(fGenPidCO);
    fGenOptGF->ShowFrame(fGenMinMomCO);
    fGenOptGF->HideFrame(fGenMaxMomCO);
    fGenOptGF->HideFrame(fGenChamberCO);
    fGenOptGF->HideFrame(fGenPrimGF);    
  }
}//GeneratorSlot()        
//__________________________________________________________________________________________________          
void RichConfig::CreateConfig()
{   
  FILE *fp=fopen(fFileName,"w"); if(!fp){Info("CreateConfigFile","Cannot open output file:%sn",fFileName);return;}
  
  fprintf(fp,"void Config()\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  ::Info(\"RICH private config\",\"Start\");\n\n"); 
//Random
  fprintf(fp,"  gRandom->SetSeed(123456);//put 0 to use system time\n\n");    
//Geant  
  fprintf(fp,"  gSystem->Load(\"libgeant321\");\n");
  fprintf(fp,"  new TGeant3TGeo(\"C++ Interface to Geant3\");\n\n");
//File
  fprintf(fp,"  AliRunLoader *pAL=AliRunLoader::Open(\"galice.root\",AliConfig::GetDefaultEventFolderName(),\"recreate\");\n");    
  fprintf(fp,"  pAL->SetCompressionLevel(2);\n");
  fprintf(fp,"  pAL->SetNumberOfEventsPerFile(1000);\n");
  fprintf(fp,"  gAlice->SetRunLoader(pAL);\n\n");
//Decayer  
  fprintf(fp,"  TVirtualMCDecayer *pDecayer=new AliDecayerPythia();\n");
  fprintf(fp,"  pDecayer->SetForceDecay(kAll);\n"); 
  fprintf(fp,"  pDecayer->Init();\n"); 
  fprintf(fp,"  gMC->SetExternalDecayer(pDecayer);\n\n");
//Physics
  if(fProcBG->GetButton(kDCAY)->GetState()) fprintf(fp,"  gMC->SetProcess(\"DCAY\",1);");  else fprintf(fp,"  gMC->SetProcess(\"DCAY\",0);");
  if(fProcBG->GetButton(kPAIR)->GetState()) fprintf(fp,"  gMC->SetProcess(\"PAIR\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"PAIR\",0);\n");
  if(fProcBG->GetButton(kCOMP)->GetState()) fprintf(fp,"  gMC->SetProcess(\"COMP\",1);");  else fprintf(fp,"  gMC->SetProcess(\"COMP\",0);");
  if(fProcBG->GetButton(kPHOT)->GetState()) fprintf(fp,"  gMC->SetProcess(\"PHOT\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"PHOT\",0);\n");
  if(fProcBG->GetButton(kPFIS)->GetState()) fprintf(fp,"  gMC->SetProcess(\"PFIS\",1);");  else fprintf(fp,"  gMC->SetProcess(\"PFIS\",0);");
  if(fProcBG->GetButton(kDRAY)->GetState()) fprintf(fp,"  gMC->SetProcess(\"DRAY\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"DRAY\",0);\n");
  if(fProcBG->GetButton(kANNI)->GetState()) fprintf(fp,"  gMC->SetProcess(\"ANNI\",1);");  else fprintf(fp,"  gMC->SetProcess(\"ANNI\",0);");
  if(fProcBG->GetButton(kBREM)->GetState()) fprintf(fp,"  gMC->SetProcess(\"BREM\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"BREM\",0);\n");
  if(fProcBG->GetButton(kMUNU)->GetState()) fprintf(fp,"  gMC->SetProcess(\"MUNU\",1);");  else fprintf(fp,"  gMC->SetProcess(\"MUNU\",0);");
  if(fProcBG->GetButton(kCKOV)->GetState()) fprintf(fp,"  gMC->SetProcess(\"CKOV\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"CKOV\",0);\n");
  if(fProcBG->GetButton(kHADR)->GetState()) fprintf(fp,"  gMC->SetProcess(\"HADR\",1);");  else fprintf(fp,"  gMC->SetProcess(\"HADR\",0);");
  if(fProcBG->GetButton(kLOSS)->GetState()) fprintf(fp,"  gMC->SetProcess(\"LOSS\",2);\n");else fprintf(fp,"  gMC->SetProcess(\"LOSS\",0);\n");
  if(fProcBG->GetButton(kMULS)->GetState()) fprintf(fp,"  gMC->SetProcess(\"MULS\",1);");  else fprintf(fp,"  gMC->SetProcess(\"MULS\",0);");
  if(fProcBG->GetButton(kRAYL)->GetState()) fprintf(fp,"  gMC->SetProcess(\"RAYL\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"RAYL\",0);\n");
  
  fprintf(fp,"\n");  
  fprintf(fp,"  gMC->SetCut(\"CUTGAM\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"CUTELE\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTNEU\",0.001);\n"); fprintf(fp,"  gMC->SetCut(\"CUTHAD\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTMUO\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"BCUTE\" ,0.001);\n");
  fprintf(fp,"  gMC->SetCut(\"BCUTM\" ,0.001); ");  fprintf(fp,"  gMC->SetCut(\"DCUTE\" ,0.001); ");
  fprintf(fp,"  gMC->SetCut(\"DCUTM\" ,0.001);\n"); fprintf(fp,"  gMC->SetCut(\"PPCUTM\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"TOFMAX\",1e10);\n\n"); 
//Field
  if(fMagCB->GetState()==kButtonDown) fprintf(fp,"  gAlice->SetField();\n\n");else fprintf(fp,"  gAlice->SetField(0);\n\n");
  fprintf(fp,"  pAL->CdGAFile();\n\n");                                 //????       
//BODY-ALIC 
  fprintf(fp,"  new AliBODY(\"BODY\",\"Alice envelop\");\n\n");
//RICH
  if(!fRichBG->GetButton(kDeclust)->GetState()) fprintf(fp,"  AliRICHParam::SetDeclustering(kFALSE);\n");
  if(!fRichBG->GetButton(kSagita)->GetState())  fprintf(fp,"  AliRICHParam::SetWireSag(kFALSE);\n")     ;
  switch(fRichCO->GetSelected()){
    case kNo:                                                                                                  break;
    case kVer0:  fprintf(fp,"  AliRICH *pRICH=new AliRICHv0(\"RICH\",\"Rich with debug StepManager\");\n\n"); break;      
    case kVer1:  fprintf(fp,"  AliRICH *pRICH=new AliRICHv1(\"RICH\",\"Normal\");\n\n");                 break;      
  }//switch RICH
//Generator
  switch(fGenTypeCO->GetSelected()){
    case kHijingPara: 
      fprintf(fp,"  AliGenHIJINGpara *pGen=new AliGenHIJINGpara(%i);\n",(int)fGenPrimNE->GetNumber());
      fprintf(fp,"  pGen->SetMomentumRange(0,999); pGen->SetThetaRange(%f,%f); pGen->SetPhiRange(0,360);\n",Eta2Theta(8),Eta2Theta(-8));
      fprintf(fp,"  pGen->SetOrigin(0,0,0);  pGen->SetSigma(0,0,0);\n");
      fprintf(fp,"  pGen->Init();\n");
    break;
    case kBox1:     
      fprintf(fp,"  AliGenBox *pGen=new AliGenBox(1);\n");  
      fprintf(fp,"  pGen->SetPart(%i); pGen->SetOrigin(0,0,0);\n",fGenPidCO->GetSelected());  
      fprintf(fp,"  pGen->SetMomentumRange(%.1f,%.1f); \n",float(fGenMinMomCO->GetSelected())/10,   float(fGenMaxMomCO->GetSelected())/10);
      fprintf(fp,"  pGen->SetThetaRange(pRICH->C(%i)->ThetaD()-3,pRICH->C(%i)->ThetaD()-1); \n",fGenChamberCO->GetSelected(),fGenChamberCO->GetSelected());
      fprintf(fp,"  pGen->SetPhiRange(pRICH->C(%i)->PhiD()-1,pRICH->C(%i)->PhiD()+1); \n",fGenChamberCO->GetSelected(),fGenChamberCO->GetSelected());    
      fprintf(fp,"  pGen->Init();\n");
    break;    
    case kBox7:   
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenBox *pBox=new AliGenBox(1);\n");
      fprintf(fp,"    pBox->SetMomentumRange(%.1f,%.1f); \n",float(fGenMinMomCO->GetSelected())/10,   float(fGenMaxMomCO->GetSelected())/10);
      fprintf(fp,"    pBox->SetPart(%i); pBox->SetOrigin(0,0,0);\n",fGenPidCO->GetSelected());
      fprintf(fp,"    pBox->SetThetaRange(pRICH->C(i)->ThetaD()-3,pRICH->C(i)->ThetaD()-1); \n");
      fprintf(fp,"    pBox->SetPhiRange(pRICH->C(i)->PhiD()-1,pRICH->C(i)->PhiD()+1); \n");    
      fprintf(fp,"    pCocktail->AddGenerator(pBox,Form(\"Box %i\",i),1);\n  }\n");  
      fprintf(fp,"  pCocktail->Init();\n");
    break;
    case kGunAlongZ:   
      fprintf(fp,"  AliGenFixed *pGen=new AliGenFixed(1);\n");
      fprintf(fp,"  pGen->SetPart(%i); pGen->SetOrigin(0,0,0); pGen->SetMomentum(%.1f); pGen->SetTheta(0);\n",fGenPidCO->GetSelected(),float(fGenMinMomCO->GetSelected())/10);
      fprintf(fp,"  pGen->Init();\n");
    break;
    case kGun7:   
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenFixed *pFixed=new AliGenFixed(1);\n");
      fprintf(fp,"    pFixed->SetPart(%i); pFixed->SetMomentum(1.0+i*0.5); pFixed->SetOrigin(0,0,0);\n",fGenPidCO->GetSelected());
      fprintf(fp,"    pFixed->SetPhi(pRICH->C(i)->PhiD()); pFixed->SetTheta(pRICH->C(i)->ThetaD()-2);\n");                             
      fprintf(fp,"    pCocktail->AddGenerator(pFixed,Form(\"Fixed %i\",i),1);\n  }\n");  
      fprintf(fp,"  pCocktail->Init();\n");
    break;
    case kPythia7:      
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenFixed *pFixed=new AliGenFixed(1);\n");
      fprintf(fp,"    pFixed->SetPart(%i); pFixed->SetMomentum(2.5+i*0.4); pFixed->SetOrigin(0,0,0);\n",fGenPidCO->GetSelected());
      fprintf(fp,"    pFixed->SetPhiRange(pRICH->C(i)->PhiD()); pFixed->SetThetaRange(pRICH->C(i)->ThetaD()-2);\n");                             
      fprintf(fp,"    pCocktail->AddGenerator(pFixed,Form(\"Fixed %i\",i),1);\n  }\n");  
      fprintf(fp,"  AliGenPythia *pPythia = new AliGenPythia(-1);\n");
      fprintf(fp,"  pPythia->SetMomentumRange(0,999999); pPythia->SetPhiRange(20,80); pPythia->SetThetaRange(75,115);\n");
      fprintf(fp,"  pPythia->SetYRange(-12,12);  pPythia->SetPtRange(0,1000);  pPythia->SetStrucFunc(kCTEQ4L);\n");
      fprintf(fp,"  pPythia->SetProcess(kPyMb);  pPythia->SetEnergyCMS(14000);\n");      
      fprintf(fp,"  pCocktail->AddGenerator(pPythia,\"Pythia\",1);\n");  
      fprintf(fp,"  pCocktail->Init();\n");
    case kSr90:  
      fprintf(fp,"  AliGenRadioactive *pGen=new AliGenRadioactive(kSr90,%i);\n",(int)fGenNprimEntry->GetNumber());
      fprintf(fp,"  pGen->SetOrigin(0,0,0);  pGen->SetSigma(0.1,0.1,0.05);\n");      
      fprintf(fp,"  pGen->Init();\n");
    break;  
    case kHijing:  
      fprintf(fp,"  AliGenHijing *pGen=new AliGenHijing(-1); pGen->SetEnergyCMS(5500); pGen->SetReferenceFrame(\"CMS\");\n");
      fprintf(fp,"  pGen->SetProjectile(\"A\", 208, 82); pGen->SetTarget(\"A\", 208, 82);\n");      
      fprintf(fp,"  pGen->SetJetQuenching(0); pGen->SetShadowing(0);\n");
      fprintf(fp,"  pGen->KeepFullEvent(); pGen->SetSelectAll(0); \n");
      fprintf(fp,"  pGen->SetImpactParameterRange(0., 5.);\n");
      fprintf(fp,"  pGen->Init();\n");
    break;  
    case kHijingPara2Proton:  
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenFixed *pProton1=new AliGenFixed(1);\n");
      fprintf(fp,"    AliGenFixed *pProton2=new AliGenFixed(1);\n");
      fprintf(fp,"    pProton1->SetPart(kProton); pProton1->SetMomentum(2.0+i*0.5); pProton1->SetOrigin(0,0,0);\n");
      fprintf(fp,"    pProton2->SetPart(kProton); pProton2->SetMomentum(2.0+i*0.5); pProton2->SetOrigin(0,0,0);\n");
      fprintf(fp,"    pProton1->SetPhiRange(pRICH->C(i)->PhiD()); pProton1->SetThetaRange(pRICH->C(i)->ThetaD());\n");                             
      fprintf(fp,"    pProton2->SetPhiRange(pRICH->C(i)->PhiD()); pProton2->SetThetaRange(pRICH->C(i)->ThetaD()-4);\n");                             
      fprintf(fp,"    pCocktail->AddGenerator(pProton1,Form(\"proton1 %i\",i),1);\n");  
      fprintf(fp,"    pCocktail->AddGenerator(pProton2,Form(\"proton2 %i\",i),1);\n  }\n");  
      
      fprintf(fp,"  AliGenHIJINGpara *pHijingPara=new AliGenHIJINGpara(510);\n");
      fprintf(fp,"  pHijingPara->SetMomentumRange(0,99); pHijingPara->SetThetaRange(60,120); pHijingPara->SetPhiRange(0,60);\n");
      fprintf(fp,"  pHijingPara->SetOrigin(0,0,0);  pHijingPara->SetSigma(0,0,0);\n");
      fprintf(fp,"  pCocktail->AddGenerator(pHijingPara,\"HijingPara\",1);\n");  
      fprintf(fp,"  pCocktail->Init();\n");
    break;  
     case kRichLib:
      fprintf(fp,"  AliGenParam *pGen=new AliGenParam(2,new AliGenRICHlib,%i,\"EXP\"); \n",fGenPidCO->GetSelected());
      fprintf(fp,"  pGen->SetPtRange(%3.1f,%3.1f); \n",fGenMinMomCO->GetSelected()/10,fGenMaxMomCO->GetSelected()/10);
      fprintf(fp,"  pGen->SetYRange(-0.6,0.6); \n");
      fprintf(fp,"  pGen->SetPhiRange(0.,60.); \n");
      fprintf(fp,"  pGen->SetForceDecay(kAll); \n");
      fprintf(fp,"  pGen->Init();\n");
    break;
   case kSignalHijing:
      fprintf(fp,"  AliGenParam *pRichLib=new AliGenParam(2,new AliGenRICHlib,%i,\"FLAT\"); \n",fGenPidCO->GetSelected());
      fprintf(fp,"  pRichLib->SetPtRange(%3.1f,%3.1f); \n",fGenMinMomCO->GetSelected()/10,fGenMaxMomCO->GetSelected()/10);
      fprintf(fp,"  pRichLib->SetYRange(-0.6,0.6); \n");
      fprintf(fp,"  pRichLib->SetPhiRange(0.,60.); \n");
      fprintf(fp,"  pRichLib->SetForceDecay(kPhiKK); \n");

      fprintf(fp,"  AliGenHijing *pHijing=new AliGenHijing(-1); pHijing->SetEnergyCMS(5500); pHijing->SetReferenceFrame(\"CMS\");\n");
      fprintf(fp,"  pHijing->SetProjectile(\"A\", 208, 82); pHijing->SetTarget(\"A\", 208, 82);\n");
      fprintf(fp,"  pHijing->SetJetQuenching(0); pHijing->SetShadowing(0);\n");
      fprintf(fp,"  pHijing->KeepFullEvent(); pHijing->SetSelectAll(0); \n");
      fprintf(fp,"  pHijing->SetImpactParameterRange(0., 5.);\n");

      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  pCocktail->AddGenerator(pRichLib,\"RichLib\",1);\n");
      fprintf(fp,"  pCocktail->AddGenerator(pHijing,\"Hijing\",1);\n");
      fprintf(fp,"  pCocktail->Init();\n");
    break;
  }
//central before RICH detectors                  
  if(fDetBG->GetButton(kPIPE )->GetState()) fprintf(fp,"\n  new AliPIPEv0(\"PIPE\",\"Beam Pipe\");\n");
  if(fDetBG->GetButton(kSHILD)->GetState()) fprintf(fp,"\n  new AliSHILv2(\"SHIL\",\"Shielding Version 2\");\n");  
  if(fDetBG->GetButton(kITS  )->GetState()){
    fprintf(fp,"\n  AliITSvPPRasymmFMD *pIts =new AliITSvPPRasymmFMD(\"ITS\",\"ITS PPR detailed version\");\n");
    fprintf(fp,"  pIts->SetMinorVersion(2); pIts->SetReadDet(kTRUE);\n");
    fprintf(fp,"  pIts->SetThicknessDet1(200.); pIts->SetThicknessDet2(200.);\n");
    fprintf(fp,"  pIts->SetThicknessChip1(200.); pIts->SetThicknessChip2(200.);\n");
    fprintf(fp,"  pIts->SetRails(0); pIts->SetCoolingFluid(1);\n");
    fprintf(fp,"  pIts->SetEUCLID(0);\n");
  }  
  if(fDetBG->GetButton(kTPC  )->GetState()) {fprintf(fp,"\n  AliTPC *pTpc=new AliTPCv2(\"TPC\",\"Default\"); pTpc->SetSecAU(-1);pTpc->SetSecAL(-1);\n");}  
  if(fDetBG->GetButton(kFRAME)->GetState())  fprintf(fp,"\n  AliFRAMEv2 *pFrame=new AliFRAMEv2(\"FRAME\",\"Space Frame\"); pFrame->SetHoles(1);\n");
  if(fDetBG->GetButton(kTRD  )->GetState()) {
    fprintf(fp,"\n  AliTRD *pTrd=new AliTRDv1(\"TRD\",\"TRD slow simulator\");\n");
    fprintf(fp,"  pTrd->SetGasMix(1); pTrd->SetPHOShole(); pTrd->SetRICHhole();pTrd->CreateTR();\n");
  }  
  if(fDetBG->GetButton(kTOF  )->GetState()) fprintf(fp,"\n  new AliTOFv4T0(\"TOF\", \"normal TOF\");\n");
//central after RICH detectors  
  if(fDetBG->GetButton(kMAG  )->GetState()) fprintf(fp,"\n  new AliMAG(\"MAG\",\"Magnet\");\n");  
  if(fDetBG->GetButton(kHALL )->GetState()) fprintf(fp,"\n  new AliHALL(\"HALL\",\"Alice Hall\");\n");
//forward detectors  
  if(fDetBG->GetButton(kFMD  )->GetState()) fprintf(fp,"\n  new AliFMDv1(\"FMD\",\"normal FMD\");\n");
  if(fDetBG->GetButton(kABSO )->GetState()) fprintf(fp,"\n  new AliABSOv0(\"ABSO\",\"Muon absorber\");\n");
  if(fDetBG->GetButton(kDIPO )->GetState()) fprintf(fp,"\n  new AliDIPOv2(\"DIPO\",\"Dipole version 2\");\n");
  if(fDetBG->GetButton(kMUON )->GetState()) fprintf(fp,"\n  new AliMUONv1(\"MUON\",\"default\");\n");
  if(fDetBG->GetButton(kPMD  )->GetState()) fprintf(fp,"\n  new AliPMDv1(\"PMD\",\"normal PMD\");\n");
  if(fDetBG->GetButton(kSTART)->GetState()) fprintf(fp,"\n  new AliSTARTv1(\"START\",\"START Detector\");\n");
  if(fDetBG->GetButton(kVZERO)->GetState()) fprintf(fp,"\n  new AliVZEROv2(\"VZERO\",\"normal VZERO\");\n");
  if(fDetBG->GetButton(kZDC  )->GetState()) fprintf(fp,"\n  new AliZDCv2(\"ZDC\",\"normal ZDC\");\n");
//different phase space detectors  
  if(fDetBG->GetButton(kPHOS )->GetState()) fprintf(fp,"\n  new AliPHOSv1(\"PHOS\",\"IHEP\");\n");
  if(fDetBG->GetButton(kEMCAL)->GetState()) fprintf(fp,"\n  new AliEMCALv1(\"EMCAL\",\"G56_2_55_19_104_14\");\n");
  if(fDetBG->GetButton(kCRT  )->GetState()) fprintf(fp,"\n  new AliCRTv0(\"CRT\",\"normal ACORDE\");\n");

  fprintf(fp,"\n  ::Info(\"RICH private config\",\"Stop\");\n"); 
  fprintf(fp,"}\n");
out:  
  fclose(fp);  
//  CloseWindow();
}//CreateConfig
//__________________________________________________________________________________________________
void RichConfig::CreateRichBatch()
{//creates RichBatch.C file
  FILE *fp=fopen("RichBatch.C","w"); if(!fp){Info("CreateRichBatch","Cannot open output file: RichBatch.C");return;}
//header and debug   
  fprintf(fp,"void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  gSystem->Exec(\"rm -rf *.root hlt hough fort* raw* ZZZ*\");\n");
  fprintf(fp,"  if(isDebug)   AliLog::SetGlobalDebugLevel(AliLog::kDebug);\n");
  fprintf(fp,"  gBenchmark->Start(\"ALICE\");\n  TDatime time;\n\n");
//simulation section  
  fprintf(fp,"  AliSimulation     *pSim=new AliSimulation(sConfigFileName);\n");
  if(fGenTypeCO->GetSelected()==kBox1||fGenTypeCO->GetSelected()==kGun7) {
    fprintf(fp,"  pSim->SetRegionOfInterest(kTRUE);\n");
    fprintf(fp,"  pSim->SetMakeSDigits(\"TOF RICH\");\n");
    fprintf(fp,"  pSim->SetMakeDigitsFromHits(\"ITS TPC TRD\");\n");
  }
//RAW data generation  
  if     (fRawBG->GetButton(kRawDdl)->GetState())  fprintf(fp,"  pSim->SetWriteRawData(\"ALL\");\n");
  else if(fRawBG->GetButton(kRawDate)->GetState()) fprintf(fp,"  pSim->SetWriteRawData(\"ALL\",\".date\");\n");
  else if(fRawBG->GetButton(kRawRoot)->GetState()) fprintf(fp,"  pSim->SetWriteRawData(\"ALL\",\".root\");\n");
  
  fprintf(fp,"  pSim->Run(iNevents);\n  delete pSim;\n\n");
//reconstraction section  
  if(fRecoBG->GetButton(kRecon)->GetState()){
    fprintf(fp,"  AliReconstruction *pRec=new AliReconstruction;\n");
    fprintf(fp,"  pRec->SetRunLocalReconstruction(\"ITS TPC TRD TOF RICH\");\n");
    if(fRecoBG->GetButton(kFillEsd)->GetState())
      fprintf(fp,"  pRec->SetFillESD(\"ITS TPC TRD TOF RICH\");\n");
    fprintf(fp,"  pRec->Run();\n");         
    fprintf(fp,"  delete pRec;\n\n");
  }
//benchmarks  
  fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/RichBatch.C>: Start time: \";time.Print();\n");
  fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/RichBatch.C>: Stop  time: \";time.Set();  time.Print();\n");
  fprintf(fp,"  gBenchmark->Show(\"ALICE\");\n");
//indicate end of job  
  fprintf(fp,"  gSystem->Exec(\"touch ZZZ______finished_______ZZZ\");\n");
  fprintf(fp,"  gSystem->Exec(\"playwave ~/bin/end.wav\");\n");
  fprintf(fp,"}\n");
  
  fclose(fp);  
}//CreateRichBatch()
//__________________________________________________________________________________________________
void RichConfig::Exit()
{
//slot to be invoked by clik on the Create button  
  CreateConfig();
  CreateRichBatch();
  SendCloseMessage();
}
//__________________________________________________________________________________________________
RichConfig *rc;
void RichConfig()
{
   rc=new RichConfig("Config.C");
}   

