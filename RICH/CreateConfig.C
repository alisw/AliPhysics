class RichConfig
{
RQ_OBJECT()
public:
  
          RichConfig(const char*sFileName);
         ~RichConfig()                    {Info("ctor","");fMainFM->Cleanup(); delete fMainFM;}
protected:
  enum ERichVers  {kNoRich=-1,kNormalRich,kTestBeam,kTestRadio,kAerogel,kOnTop};
  enum EDetectors {kPIPE=1,kITS=2,kTPC=4,kTRD=8,kTOF=16,kFRAME=32,kMAG=64,kCRT=128,kHALL=256,kPHOS=512,kSTART=1024,kFMD=2048,kABSO=4096,
                   kPMD=8192,kDIPO=16384,kEMCAL=32768,kVZERO=65536,kMUON=131072,kZDC=262144,kSHILD=524288};
  enum EProcesses {kDCAY=1,kPAIR=2,kCOMP=4,kPHOT=8,kPFIS=16,kDRAY=32,kANNI=64,kBREM=128,kMUNU=256,kCKOV=512,kHADR=1024,kLOSS=2048,kMULS=4096,
                   kRAYL=8192};
  enum EGenTypes  {kGun1=100,kGun7,kPythia7,kHijing,kHijingPara};
  
  void    AddDetector(Int_t id)          {fDetectors+=id; if(id==kTRD || id==kTOF) fDetButGrp->SetButton(kFRAME);}
  void    RemoveDetector(Int_t id)       {fDetectors-=id; if(id==kFRAME) {fDetButGrp->SetButton(kTRD,kFALSE);fDetButGrp->SetButton(kTOF,kFALSE);}}
  void    AddProcess(Int_t id)           {fProcesses+=id;}
  void    RemoveProcess(Int_t id)        {fProcesses-=id;}
  Bool_t  IsDetectorOn(Int_t id)    const{return fDetectors&id;}
  Bool_t  IsProcessOn(Int_t id)     const{return fProcesses&id;}
  Float_t Eta2Theta(Float_t arg)    const{return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));}
  void    CreateConfig();
  void    CreateRichBatch();
  void    Exit();
  
  TGMainFrame  *fMainFM;//main window poiter
  TGComboBox   *fRichVerCombo;  TGButton *fRichDeclusterBC,*fRichSagBC,*fRichDebugBC;//RICH
  TGButton     *fMagFldBC;                               //MAG
  TGComboBox   *fGenTypeCombo,*fGenPartIdCombo,*fGenMinMomCombo,*fGenMaxMomCombo,*fGenChamberCombo; TGNumberEntry *fGenNprimEntry;//GEN
  Int_t         fDetectors; TGButtonGroup *fDetButGrp;       //DETECTORS
  Int_t         fProcesses;
  char         *fFileName;
};//class RichConfig
//__________________________________________________________________________________________________	 
RichConfig::RichConfig(const char *sFileName)
{// creates configuration file  
  fFileName=sFileName;
  fDetectors=0;
// Create main frame       
  fMainFM=new TGMainFrame(gClient->GetRoot(),500,400);
//  fMain->Connect("CloseWindow()","RichConfig",this,"CloseWindow()");   
  fMainFM->AddFrame(pHorFrame=new TGHorizontalFrame(fMainFM,100,200));
//RICH
  pHorFrame->AddFrame(pVerFrame=new TGVerticalFrame(pHorFrame,100,200));  
  pVerFrame->AddFrame(pRichFG=new TGGroupFrame(pHorFrame,"RICH"));
  pRichFG->AddFrame(fRichVerCombo=new TGComboBox(pRichFG,100));
  fRichVerCombo->AddEntry("no RICH"      ,kNoRich);
  fRichVerCombo->AddEntry("normal RICH"  ,kNormalRich);
  fRichVerCombo->AddEntry("RICH on top"  ,kOnTop);
  fRichVerCombo->AddEntry("test beam"    ,kTestBeam);
  fRichVerCombo->AddEntry("radio source" ,kTestRadio);
  fRichVerCombo->AddEntry("aerogel"      ,kAerogel);
  fRichVerCombo->Select(kNormalRich);  fRichVerCombo->Resize(150,20);
  pRichFG->AddFrame(fRichDeclusterBC=new TGCheckButton(pRichFG,"Declustering"));       fRichDeclusterBC->SetState(kButtonDown);
  pRichFG->AddFrame(fRichSagBC      =new TGCheckButton(pRichFG,"Wire sagita?"));       fRichSagBC      ->SetState(kButtonDown);
  pRichFG->AddFrame(fRichDebugBC    =new TGCheckButton(pRichFG,"Debug StepManager?"));
//Generator  
  pVerFrame->AddFrame(pGenGrpFrm=new TGGroupFrame(pHorFrame,"Generator"));
  pGenGrpFrm->AddFrame(fGenTypeCombo=new TGComboBox(pGenGrpFrm,100));
  fGenTypeCombo->AddEntry("gun to single chamber",kGun1);  
  fGenTypeCombo->AddEntry("gun to all chambers",kGun7);
  fGenTypeCombo->AddEntry("7 guns on top of Pythia",kPythia7);
  fGenTypeCombo->AddEntry("HIJING",kHijing);
  fGenTypeCombo->AddEntry("parametrized HIJING",kHijingPara);
  fGenTypeCombo->AddEntry("Sr90 source",kSr90);
  fGenTypeCombo->Select(kHijingPara);
  fGenTypeCombo->Resize(160,20);
  
  pGenGrpFrm->AddFrame(fGenPartIdCombo=new TGComboBox(pGenGrpFrm,100)); //PID for guns
  fGenPartIdCombo->AddEntry("Pion+",kPiPlus);
  fGenPartIdCombo->AddEntry("Pion-",kPiMinus);
  fGenPartIdCombo->AddEntry("Kaon+",kKPlus);
  fGenPartIdCombo->AddEntry("Kaon-",kKMinus);
  fGenPartIdCombo->AddEntry("Proton",kProton);
  fGenPartIdCombo->AddEntry("AntiProton",kProtonBar);
  fGenPartIdCombo->Select(kProton);  fGenPartIdCombo->Resize(160,20);

  pGenGrpFrm->AddFrame(fGenMinMomCombo=new TGComboBox(pGenGrpFrm,100)); //particle energy for guns
  fGenMinMomCombo->AddEntry("0.5 GeV", 5);
  fGenMinMomCombo->AddEntry("1.0 GeV",10);
  fGenMinMomCombo->AddEntry("1.5 GeV",15);
  fGenMinMomCombo->AddEntry("2.0 Gev",20);
  fGenMinMomCombo->AddEntry("2.5 GeV",25);
  fGenMinMomCombo->AddEntry("3.0 GeV",30);
  fGenMinMomCombo->AddEntry("3.5 GeV",35);
  fGenMinMomCombo->AddEntry("4.0 GeV",40);
  fGenMinMomCombo->AddEntry("4.5 GeV",45);
  fGenMinMomCombo->Select(10);  fGenMinMomCombo->Resize(160,20);
  
  pGenGrpFrm->AddFrame(fGenMaxMomCombo=new TGComboBox(pGenGrpFrm,100)); //particle energy for guns
  fGenMaxMomCombo->AddEntry("1.0 GeV",10);
  fGenMaxMomCombo->AddEntry("1.5 GeV",15);
  fGenMaxMomCombo->AddEntry("2.0 Gev",20);
  fGenMaxMomCombo->AddEntry("2.5 GeV",25);
  fGenMaxMomCombo->AddEntry("3.0 GeV",30);
  fGenMaxMomCombo->AddEntry("3.5 GeV",35);
  fGenMaxMomCombo->AddEntry("4.0 GeV",40);
  fGenMaxMomCombo->AddEntry("4.5 GeV",45);
  fGenMaxMomCombo->Select(40);  fGenMaxMomCombo->Resize(160,20);
  
  pGenGrpFrm->AddFrame(fGenChamberCombo=new TGComboBox(pGenGrpFrm,100)); //chamber number in case of gun1
  for(int i=1;i<=7;i++) fGenChamberCombo->AddEntry(Form("Chamber %i",i),i);
  fGenChamberCombo->Select(4); fGenChamberCombo->Resize(160,20);
  
  pGenGrpFrm->AddFrame(pGenNprimFrm=new TGGroupFrame(pGenGrpFrm,"Number of primaries"));//number of primiries in case of HIJING
  pGenNprimFrm->AddFrame(fGenNprimEntry=new TGNumberEntry(pGenNprimFrm,500));
  
// Magnetic Field
  pVerFrame->AddFrame(pFldGrpFrm=new TGGroupFrame(pHorFrame,"Magnetic Field"));
  pFldGrpFrm->AddFrame(fMagFldBC=new TGCheckButton(pFldGrpFrm,"On/Off"));
  fMagFldBC->SetState(kButtonDown);
//Detectors
  pHorFrame->AddFrame(fDetButGrp=new TGButtonGroup(pHorFrame,"Detectors"));
  fDetButGrp->Connect("Pressed(Int_t)","RichConfig",this,"AddDetector(Int_t)");
  fDetButGrp->Connect("Released(Int_t)","RichConfig",this,"RemoveDetector(Int_t)");
  new TGCheckButton(fDetButGrp,"PIPE"  ,kPIPE));         
  new TGCheckButton(fDetButGrp,"ITS"   ,kITS));         
  new TGCheckButton(fDetButGrp,"TPC"   ,kTPC));
  new TGCheckButton(fDetButGrp,"TRD"   ,kTRD));
  new TGCheckButton(fDetButGrp,"TOF"   ,kTOF));         
  new TGCheckButton(fDetButGrp,"FRAME" ,kFRAME));         
  new TGCheckButton(fDetButGrp,"MAG"   ,kMAG));         
  new TGCheckButton(fDetButGrp,"CRT"   ,kCRT));         
  new TGCheckButton(fDetButGrp,"HALL"  ,kHALL));         
  new TGCheckButton(fDetButGrp,"PHOS"  ,kPHOS));         
  new TGCheckButton(fDetButGrp,"START" ,kSTART));         
  new TGCheckButton(fDetButGrp,"FMD"   ,kFMD));         
  new TGCheckButton(fDetButGrp,"ABSO"  ,kABSO));         
  new TGCheckButton(fDetButGrp,"PMD"   ,kPMD));         
  new TGCheckButton(fDetButGrp,"DIPO"  ,kDIPO));         
  new TGCheckButton(fDetButGrp,"EMCAL" ,kEMCAL));         
  new TGCheckButton(fDetButGrp,"VZERO" ,kVZERO));         
  new TGCheckButton(fDetButGrp,"MUON"  ,kMUON));         
  new TGCheckButton(fDetButGrp,"ZDC"   ,kZDC));         
  new TGCheckButton(fDetButGrp,"SHILD" ,kSHILD));         
//Processes  
  pHorFrame->AddFrame(pProcBG=new TGButtonGroup(pHorFrame,"Processes"));
  pProcBG->Connect("Pressed(Int_t)","RichConfig",this,"AddProcess(Int_t)");
  pProcBG->Connect("Released(Int_t)","RichConfig",this,"RemoveProcess(Int_t)");
  new TGCheckButton(pProcBG,"DCAY Decay",kDCAY))                 ;  pProcBG->SetButton(kDCAY);       
  new TGCheckButton(pProcBG,"PAIR Pair production",kPAIR))       ;  pProcBG->SetButton(kPAIR);       
  new TGCheckButton(pProcBG,"COMP Compton",kCOMP))               ;  pProcBG->SetButton(kCOMP);
  new TGCheckButton(pProcBG,"PHOT",kPHOT))                       ;  pProcBG->SetButton(kPHOT);
  new TGCheckButton(pProcBG,"PFIS Photofission",kPFIS))          ;  
  new TGCheckButton(pProcBG,"DRAY Delta electrons",kDRAY))       ;  
  new TGCheckButton(pProcBG,"ANNI Annihilation",kANNI))          ;  pProcBG->SetButton(kANNI);       
  new TGCheckButton(pProcBG,"BREM Bremstraslung",kBREM))         ;  pProcBG->SetButton(kBREM);       
  new TGCheckButton(pProcBG,"MUNU",kMUNU))                       ;  pProcBG->SetButton(kMUNU);       
  new TGCheckButton(pProcBG,"CKOV Cerenkovs",kCKOV))             ;  pProcBG->SetButton(kCKOV);       
  new TGCheckButton(pProcBG,"HADR Hadronic interactions ",kHADR));  pProcBG->SetButton(kHADR);       
  new TGCheckButton(pProcBG,"LOSS Energy losses",kLOSS))         ;  pProcBG->SetButton(kLOSS);       
  new TGCheckButton(pProcBG,"MULS Multiple scattering",kMULS))   ;  pProcBG->SetButton(kMULS);       
  new TGCheckButton(pProcBG,"RAYL",kRAYL))                       ;  pProcBG->SetButton(kRAYL);       
//File    
  fMainFM->AddFrame(pFileHorFrm=new TGHorizontalFrame(fMainFM,100,200));
  pFileHorFrm->AddFrame(pCreateB=new TGTextButton(pFileHorFrm,"Create"));
  pCreateB->Connect("Clicked()","RichConfig",this,"Exit()");                                 
  pFileHorFrm->AddFrame(new TGLabel(pFileHorFrm,Form(" config file as %s",fFileName)));  
  
  fMainFM->MapSubwindows();   
  fMainFM->Layout();
  fMainFM->SetWindowName("Create AliROOT scripts");
  fMainFM->MapWindow();      
}//KirCondig::ctor
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
  fprintf(fp,"  new TGeant3(\"C++ Interface to Geant3\");\n\n");
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
  if(IsProcessOn(kDCAY)) fprintf(fp,"  gMC->SetProcess(\"DCAY\",1);");  else fprintf(fp,"  gMC->SetProcess(\"DCAY\",0);");
  if(IsProcessOn(kPAIR)) fprintf(fp,"  gMC->SetProcess(\"PAIR\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"PAIR\",0);\n");
  if(IsProcessOn(kCOMP)) fprintf(fp,"  gMC->SetProcess(\"COMP\",1);");  else fprintf(fp,"  gMC->SetProcess(\"COMP\",0);");
  if(IsProcessOn(kPHOT)) fprintf(fp,"  gMC->SetProcess(\"PHOT\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"PHOT\",0);\n");
  if(IsProcessOn(kPFIS)) fprintf(fp,"  gMC->SetProcess(\"PFIS\",1);");  else fprintf(fp,"  gMC->SetProcess(\"PFIS\",0);");
  if(IsProcessOn(kDRAY)) fprintf(fp,"  gMC->SetProcess(\"DRAY\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"DRAY\",0);\n");
  if(IsProcessOn(kANNI)) fprintf(fp,"  gMC->SetProcess(\"ANNI\",1);");  else fprintf(fp,"  gMC->SetProcess(\"ANNI\",0);");
  if(IsProcessOn(kBREM)) fprintf(fp,"  gMC->SetProcess(\"BREM\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"BREM\",0);\n");
  if(IsProcessOn(kMUNU)) fprintf(fp,"  gMC->SetProcess(\"MUNU\",1);");  else fprintf(fp,"  gMC->SetProcess(\"MUNU\",0);");
  if(IsProcessOn(kCKOV)) fprintf(fp,"  gMC->SetProcess(\"CKOV\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"CKOV\",0);\n");
  if(IsProcessOn(kHADR)) fprintf(fp,"  gMC->SetProcess(\"HADR\",1);");  else fprintf(fp,"  gMC->SetProcess(\"HADR\",0);");
  if(IsProcessOn(kLOSS)) fprintf(fp,"  gMC->SetProcess(\"LOSS\",2);\n");else fprintf(fp,"  gMC->SetProcess(\"LOSS\",0);\n");
  if(IsProcessOn(kMULS)) fprintf(fp,"  gMC->SetProcess(\"MULS\",1);");  else fprintf(fp,"  gMC->SetProcess(\"MULS\",0);");
  if(IsProcessOn(kRAYL)) fprintf(fp,"  gMC->SetProcess(\"RAYL\",1);\n");else fprintf(fp,"  gMC->SetProcess(\"RAYL\",0);\n");
  
  fprintf(fp,"\n");  
  fprintf(fp,"  gMC->SetCut(\"CUTGAM\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"CUTELE\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTNEU\",0.001);\n"); fprintf(fp,"  gMC->SetCut(\"CUTHAD\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTMUO\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"BCUTE\" ,0.001);\n");
  fprintf(fp,"  gMC->SetCut(\"BCUTM\" ,0.001); ");  fprintf(fp,"  gMC->SetCut(\"DCUTE\" ,0.001); ");
  fprintf(fp,"  gMC->SetCut(\"DCUTM\" ,0.001);\n");  fprintf(fp,"  gMC->SetCut(\"PPCUTM\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"TOFMAX\",1e10);\n\n"); 
//Field
  if(fMagFldBC->GetState()==kButtonDown) fprintf(fp,"  gAlice->SetField();\n\n");else fprintf(fp,"  gAlice->SetField(0);\n\n");
  fprintf(fp,"  pAL->CdGAFile();\n\n");                                 //????       
//BODY-ALIC 
  fprintf(fp,"  new AliBODY(\"BODY\",\"Alice envelop\");\n\n");
//RICH
  if(fRichVerCombo->GetSelected() != kNoRich){  
    TString richStr="RICH "; 
    if(fRichDeclusterBC->GetState()!=kButtonDown){richStr+="no declustering "  ; fprintf(fp,"  AliRICHParam::SetDeclustering(kFALSE);\n");}
    if(fRichSagBC->GetState()      !=kButtonDown){richStr+="no sagita "        ; fprintf(fp,"  AliRICHParam::SetWireSag(kFALSE);\n")     ;}
    switch(fRichVerCombo->GetSelected()){//RICH
      case kOnTop:      richStr+="top position";fprintf(fp,"  AliRICHParam::SetAngleRot(0);\n")      ;break;   
      case kAerogel:    richStr+="aerogel"     ;fprintf(fp,"  AliRICHParam::SetAerogel(kTRUE);\n") ;break;   
      case kTestRadio:  richStr+="radioactive" ;fprintf(fp,"  AliRICHParam::SetRadioSrc(kTRUE);\n");
                                                fprintf(fp,"  AliRICHParam::SetTestBeam(kTRUE);\n");break;   
      case kTestBeam:   richStr+="test beam"   ;fprintf(fp,"  AliRICHParam::SetTestBeam(kTRUE);\n");break;   
      case kNormalRich: richStr+="normal"      ;                                                    break;   
    }
    if(fRichDebugBC->GetState()    ==kButtonDown){
      richStr+="+debug StepManager"; 
      fprintf(fp,"  AliRICH *pRICH=new AliRICHv0(\"RICH\",\"%s\");\n\n",richStr.Data());
    }else{
      fprintf(fp,"  AliRICH *pRICH=new AliRICHv1(\"RICH\",\"%s\");\n\n",richStr.Data());
    }
  }//Rich not selected
//Generator
  switch(fGenTypeCombo->GetSelected()){
    case kHijingPara: 
      fprintf(fp,"  AliGenHIJINGpara *pGen=new AliGenHIJINGpara(%i);\n",(int)fGenNprimEntry->GetNumber());
      fprintf(fp,"  pGen->SetMomentumRange(0,999); pGen->SetThetaRange(%f,%f); pGen->SetPhiRange(0,360);\n",Eta2Theta(8),Eta2Theta(-8));
      fprintf(fp,"  pGen->SetOrigin(0,0,0);  pGen->SetSigma(0,0,0);\n");
      fprintf(fp,"  pGen->Init();\n");
    break;
    case kGun1:     
      fprintf(fp,"  AliGenBox *pGen=new AliGenBox(1);\n");  
      fprintf(fp,"  pGen->SetPart(%i); pGen->SetOrigin(0,0,0);\n",fGenPartIdCombo->GetSelected());  
      fprintf(fp,"  pGen->SetMomentumRange(%3.1f,%3.1f); \n",float(fGenMinMomCombo->GetSelected())/10,   float(fGenMaxMomCombo->GetSelected())/10);
      fprintf(fp,"  pGen->SetThetaRange(pRICH->C(%i)->ThetaD()-3,pRICH->C(%i)->ThetaD()-1); \n",fGenChamberCombo->GetSelected(),fGenChamberCombo->GetSelected());
      fprintf(fp,"  pGen->SetPhiRange(pRICH->C(%i)->PhiD()-1,pRICH->C(%i)->PhiD()+1); \n",fGenChamberCombo->GetSelected(),fGenChamberCombo->GetSelected());    
      fprintf(fp,"  pGen->Init();\n");
    break;    
    case kGun7:   
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenFixed *pFixed=new AliGenFixed(1);\n");
      fprintf(fp,"    pFixed->SetPart(%i); pFixed->SetMomentum(1.0+i*0.5); pFixed->SetOrigin(0,0,0);\n",fGenPartIdCombo->GetSelected());
      fprintf(fp,"    pFixed->SetPhi(pRICH->C(i)->PhiD()); pFixed->SetTheta(pRICH->C(i)->ThetaD()-2);\n");                             
      fprintf(fp,"    pCocktail->AddGenerator(pFixed,Form(\"Fixed %i\",i),1);\n  }\n");  
      fprintf(fp,"  pCocktail->Init();\n");
    break;
    case kPythia7:      
      fprintf(fp,"  AliGenCocktail *pCocktail=new AliGenCocktail();\n");
      fprintf(fp,"  for(int i=1;i<=7;i++){\n");
      fprintf(fp,"    AliGenFixed *pFixed=new AliGenFixed(1);\n");
      fprintf(fp,"    pFixed->SetPart(%i); pFixed->SetMomentum(2.5+i*0.4); pFixed->SetOrigin(0,0,0);\n",fGenPartIdCombo->GetSelected());
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
  }
//central before RICH detectors                  
  if(IsDetectorOn(kPIPE)) fprintf(fp,"\n  new AliPIPEv0(\"PIPE\",\"Beam Pipe\");\n");
  if(IsDetectorOn(kSHILD))fprintf(fp,"\n  new AliSHILv2(\"SHIL\",\"Shielding Version 2\");\n");  
  if(IsDetectorOn(kITS)){
    fprintf(fp,"\n  AliITSvPPRasymmFMD *pIts =new AliITSvPPRasymmFMD(\"ITS\",\"ITS PPR detailed version\");\n");
    fprintf(fp,"  pIts->SetMinorVersion(2); pIts->SetReadDet(kTRUE);\n");
    fprintf(fp,"  pIts->SetThicknessDet1(200.); pIts->SetThicknessDet2(200.);\n");
    fprintf(fp,"  pIts->SetThicknessChip1(200.); pIts->SetThicknessChip2(200.);\n");
    fprintf(fp,"  pIts->SetRails(0); pIts->SetCoolingFluid(1);\n");
    fprintf(fp,"  pIts->SetEUCLID(0);\n");
  }  
  if(IsDetectorOn(kTPC)) {fprintf(fp,"\n  AliTPC *pTpc=new AliTPCv2(\"TPC\",\"Default\"); pTpc->SetSecAU(-1);pTpc->SetSecAL(-1);\n");}  
  if(IsDetectorOn(kFRAME))fprintf(fp,"\n  AliFRAMEv2 *pFrame=new AliFRAMEv2(\"FRAME\",\"Space Frame\"); pFrame->SetHoles(1);\n");
  if(IsDetectorOn(kTRD)){
    fprintf(fp,"\n  AliTRD *pTrd=new AliTRDv1(\"TRD\",\"TRD slow simulator\");\n");
    fprintf(fp,"  pTrd->SetGasMix(1); pTrd->SetPHOShole(); pTrd->SetRICHhole();pTrd->CreateTR();\n");
  }  
  if(IsDetectorOn(kTOF))   fprintf(fp,"\n  new AliTOFv4T0(\"TOF\", \"normal TOF\");\n");
//central after RICH detectors  
  if(IsDetectorOn(kMAG))  fprintf(fp,"\n  new AliMAG(\"MAG\",\"Magnet\");\n");  
  if(IsDetectorOn(kHALL)) fprintf(fp,"\n  new AliHALL(\"HALL\",\"Alice Hall\");\n");
//forward detectors  
  if(IsDetectorOn(kFMD))   fprintf(fp,"\n  new AliFMDv1(\"FMD\",\"normal FMD\");\n");
  if(IsDetectorOn(kABSO))  fprintf(fp,"\n  new AliABSOv0(\"ABSO\",\"Muon absorber\");\n");
  if(IsDetectorOn(kDIPO))  fprintf(fp,"\n  new AliDIPOv2(\"DIPO\",\"Dipole version 2\");\n");
  if(IsDetectorOn(kMUON))  fprintf(fp,"\n  new AliMUONv1(\"MUON\",\"default\");\n");
  if(IsDetectorOn(kPMD))   fprintf(fp,"\n  new AliPMDv1(\"PMD\",\"normal PMD\");\n");
  if(IsDetectorOn(kSTART)) fprintf(fp,"\n  new AliSTARTv1(\"START\",\"START Detector\");\n");
  if(IsDetectorOn(kVZERO)) fprintf(fp,"\n  new AliVZEROv2(\"VZERO\",\"normal VZERO\");\n");
  if(IsDetectorOn(kZDC))   fprintf(fp,"\n  new AliZDCv2(\"ZDC\",\"normal ZDC\");\n");
//different phase space detectors  
  if(IsDetectorOn(kPHOS))  fprintf(fp,"\n  new AliPHOSv1(\"PHOS\",\"IHEP\");\n");
  if(IsDetectorOn(kEMCAL)) fprintf(fp,"\n  new AliEMCALv1(\"EMCAL\",\"G56_2_55_19_104_14\");\n");
  if(IsDetectorOn(kCRT))   fprintf(fp,"\n  new AliCRTv0(\"CRT\",\"normal ACORDE\");\n");

  fprintf(fp,"\n  ::Info(\"RICH private config\",\"Stop\");\n"); 
  fprintf(fp,"}\n");
out:  
  fclose(fp);  
//  fMain->CloseWindow();
}//CreateConfig
//__________________________________________________________________________________________________
void RichConfig::CreateRichBatch()
{//creates RichBatch.C file
  FILE *fp=fopen("RichBatch.C","w"); if(!fp){Info("CreateRichBatch","Cannot open output file: RichBatch.C");return;}
  
  fprintf(fp,"void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  gSystem->Exec(\"rm -rf *.root hlt hough fort* ZZZ*\");\n");
  fprintf(fp,"  if(isDebug)   AliLog::SetGlobalDebugLevel(AliLog::kDebug);\n");
  fprintf(fp,"  TStopwatch sw;TDatime time;\n\n");
  
  fprintf(fp,"  AliSimulation     *pSim=new AliSimulation(sConfigFileName);\n");
  if(fGenTypeCombo->GetSelected()==kGun1||fGenTypeCombo->GetSelected()==kGun7) {
    fprintf(fp,"  pSim->SetRegionOfInterest(kTRUE);\n");
    fprintf(fp,"  pSim->SetMakeSDigits(\"TOF RICH\");\n");
    fprintf(fp,"  pSim->SetMakeDigitsFromHits(\"ITS TPC TRD\");\n");
  }
  fprintf(fp,"  pSim->Run(iNevents); delete pSim;\n\n");
  if(fRichVerCombo->GetSelected()==kNormalRich){
    fprintf(fp,"  AliReconstruction *pRec=new AliReconstruction;\n");
    fprintf(fp,"  pRec->SetRunLocalReconstruction(\"ITS TPC TRD TOF RICH\");\n");
    fprintf(fp,"  pRec->SetFillESD(\"ITS TPC TRD TOF RICH\");\n");
    fprintf(fp,"  pRec->Run();         delete pRec;\n\n");
  }
  
  fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/RichBatch.C>: Start time: \";time.Print();\n");
  fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/RichBatch.C>: Stop  time: \";time.Set();  time.Print();\n");
  fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/RichBatch.C>: Time  used: \";sw.Print();\n\n");
  
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
  fMainFM->SendCloseMessage();
}
//__________________________________________________________________________________________________
RichConfig *p;
void CreateConfig()
{
   p=new RichConfig("Config.C");
}   

