class KirConfig
{
RQ_OBJECT()
public:
  enum EDetectors {kPIPE=1,kITS=2,kTPC=4,kTRD=8,kTOF=16,kFRAME=32,kMAG=64,kCRT=128,kHALL=256,kPHOS=512,kSTART=1024,kFMD=2048,kABSO=4096,
                   kPMD=8192,kDIPO=16384,kEMCAL=32768,kVZERO=65536,kMUON=131072,kZDC=262144,kSHILD=524288};
  enum EProcesses {kDCAY=1,kPAIR=2,kCOMP=4,kPHOT=8,kPFIS=16,kDRAY=32,kANNI=64,kBREM=128,kMUNU=256,kCKOV=512,kHADR=1024,kLOSS=2048,kMULS=4096,
                   kRAYL=8192};
  enum EGenTypes  {kGun1,kGun7,kHijing,kHijingPara};
      KirConfig();
     ~KirConfig()    {Info("ctor","");fMain->CleanUp(); delete fMain;}
  void   AddDetector(Int_t id)          {fDetectors+=id;}
  void   RemoveDetector(Int_t id)       {fDetectors-=id;}
  void   AddProcess(Int_t id)           {fProcesses+=id;}
  void   RemoveProcess(Int_t id)        {fProcesses-=id;}
  Bool_t IsDetectorOn(Int_t id)    const{return fDetectors&id;}
  Bool_t IsProcessOn(Int_t id)     const{return fProcesses&id;}
  void   CreateConfigFile();
protected:
  TGMainFrame  *fMain;//main window poiter
  TGComboBox   *fRichVersionCombo;
  TGButton     *fRichTopChkBtn, *fMagFldChkBtn;
  TGComboBox   *fGenTypeCombo,*fPartIdCombo,*fMomCombo;
  Int_t         fDetectors;
  Int_t         fProcesses;
};//class KirConfig
	 

KirConfig::KirConfig()
{   
// Create main frame       
  fMain=new TGMainFrame(gClient->GetRoot(),500,400);
//  fMain->Connect("CloseWindow()","KirConfig",this,"CloseWindow()");   
  fMain->AddFrame(pHorFrame=new TGHorizontalFrame(fMain,100,200));
//RICH
  pHorFrame->AddFrame(pVerFrame=new TGVerticalFrame(pHorFrame,100,200));  
  pVerFrame->AddFrame(pRichGrpFrm=new TGGroupFrame(pHorFrame,"RICH"));
  pRichGrpFrm->AddFrame(fRichVersionCombo=new TGComboBox(pRichGrpFrm,100));
  fRichVersionCombo->AddEntry("version 0",0);
  fRichVersionCombo->AddEntry("version 1",1);
  fRichVersionCombo->AddEntry("version 3",3);
  fRichVersionCombo->Select(3);
  fRichVersionCombo->Resize(140,20);
  pRichGrpFrm->AddFrame(fRichTopChkBtn=new TGCheckButton(pRichGrpFrm,"Rotate to Top?"));
//Generator  
  pVerFrame->AddFrame(pGenGrpFrm=new TGGroupFrame(pHorFrame,"Generator"));
  pGenGrpFrm->AddFrame(fGenTypeCombo=new TGComboBox(pGenGrpFrm,100));
  fGenTypeCombo->AddEntry("gun to central chamber",kGun1);
  fGenTypeCombo->AddEntry("gun to all chambers",kGun7);
  fGenTypeCombo->AddEntry("HIJING",kHijing);
  fGenTypeCombo->AddEntry("parametrized HIJING",kHijingPara);
  fGenTypeCombo->Select(kGun1);
  fGenTypeCombo->Resize(160,20);
  
  pGenGrpFrm->AddFrame(fPartIdCombo=new TGComboBox(pGenGrpFrm,100));
  fPartIdCombo->AddEntry("Pion+",kPiPlus);
  fPartIdCombo->AddEntry("Pion-",kPiMinus);
  fPartIdCombo->AddEntry("Kaon+",kKPlus);
  fPartIdCombo->AddEntry("Kaon-",kKMinus);
  fPartIdCombo->AddEntry("Proton",kProton);
  fPartIdCombo->AddEntry("AntiProton",kProtonBar);
  fPartIdCombo->Select(kPiPlus);
  fPartIdCombo->Resize(160,20);

  pGenGrpFrm->AddFrame(fMomCombo=new TGComboBox(pGenGrpFrm,100));
  fMomCombo->AddEntry("1.5 GeV",15);
  fMomCombo->AddEntry("2.0 Gev",20);
  fMomCombo->AddEntry("2.5 GeV",25);
  fMomCombo->AddEntry("3.0 GeV",30);
  fMomCombo->AddEntry("3.5 GeV",35);
  fMomCombo->AddEntry("4.0 GeV",40);
  fMomCombo->AddEntry("4.5 GeV",45);
  fMomCombo->Select(40);
  fMomCombo->Resize(160,20);
// Magnetic Field
  pVerFrame->AddFrame(pFldGrpFrm=new TGGroupFrame(pHorFrame,"Magnetic Field"));
  pFldGrpFrm->AddFrame(fMagFldChkBtn=new TGCheckButton(pFldGrpFrm,"On/Off?"));
  fMagFldChkBtn->SetState(kButtonDown);
//Detectors
  pHorFrame->AddFrame(pDetButGrp=new TGButtonGroup(pHorFrame,"Detectors"));
  pDetButGrp->Connect("Pressed(Int_t)","KirConfig",this,"AddDetector(Int_t)");
  pDetButGrp->Connect("Released(Int_t)","KirConfig",this,"RemoveDetector(Int_t)");
  new TGCheckButton(pDetButGrp,"PIPE"  ,kPIPE));         
  new TGCheckButton(pDetButGrp,"ITS"   ,kITS));         
  new TGCheckButton(pDetButGrp,"TPC"   ,kTPC));
  new TGCheckButton(pDetButGrp,"TRD"   ,kTRD));
  new TGCheckButton(pDetButGrp,"TOF"   ,kTOF));         
  new TGCheckButton(pDetButGrp,"FRAME" ,kFRAME));         
  new TGCheckButton(pDetButGrp,"MAG"   ,kMAG));         
  new TGCheckButton(pDetButGrp,"CRT"   ,kCRT));         
  new TGCheckButton(pDetButGrp,"HALL"  ,kHALL));         
  new TGCheckButton(pDetButGrp,"PHOS"  ,kPHOS));         
  new TGCheckButton(pDetButGrp,"START" ,kSTART));         
  new TGCheckButton(pDetButGrp,"FMD"   ,kFMD));         
  new TGCheckButton(pDetButGrp,"ABSO"  ,kABSO));         
  new TGCheckButton(pDetButGrp,"PMD"   ,kPMD));         
  new TGCheckButton(pDetButGrp,"DIPO"  ,kDIPO));         
  new TGCheckButton(pDetButGrp,"EMCAL" ,kEMCAL));         
  new TGCheckButton(pDetButGrp,"VZERO" ,kVZERO));         
  new TGCheckButton(pDetButGrp,"MUON"  ,kMUON));         
  new TGCheckButton(pDetButGrp,"ZDC"   ,kZDC));         
  new TGCheckButton(pDetButGrp,"SHILD" ,kSHILD));         
//Processes  
  pHorFrame->AddFrame(pProcButGrp=new TGButtonGroup(pHorFrame,"Processes"));
  pProcButGrp->Connect("Pressed(Int_t)","KirConfig",this,"AddProcess(Int_t)");
  pProcButGrp->Connect("Released(Int_t)","KirConfig",this,"RemoveProcess(Int_t)");
  new TGCheckButton(pProcButGrp,"DCAY",kDCAY));  pProcButGrp->SetButton(kDCAY);       
  new TGCheckButton(pProcButGrp,"PAIR",kPAIR));  pProcButGrp->SetButton(kPAIR);       
  new TGCheckButton(pProcButGrp,"COMP",kCOMP));  pProcButGrp->SetButton(kCOMP);
  new TGCheckButton(pProcButGrp,"PHOT",kPHOT));  pProcButGrp->SetButton(kPHOT);
  new TGCheckButton(pProcButGrp,"PFIS",kPFIS));  
  new TGCheckButton(pProcButGrp,"DRAY",kDRAY));  
  new TGCheckButton(pProcButGrp,"ANNI",kANNI));  pProcButGrp->SetButton(kANNI);       
  new TGCheckButton(pProcButGrp,"BREM",kBREM));  pProcButGrp->SetButton(kBREM);       
  new TGCheckButton(pProcButGrp,"MUNU",kMUNU));  pProcButGrp->SetButton(kMUNU);       
  new TGCheckButton(pProcButGrp,"CKOV",kCKOV));  pProcButGrp->SetButton(kCKOV);       
  new TGCheckButton(pProcButGrp,"HADR",kHADR));  pProcButGrp->SetButton(kHADR);       
  new TGCheckButton(pProcButGrp,"LOSS",kLOSS));  pProcButGrp->SetButton(kLOSS);       
  new TGCheckButton(pProcButGrp,"MULS",kMULS));  pProcButGrp->SetButton(kMULS);       
  new TGCheckButton(pProcButGrp,"RAYL",kRAYL));  pProcButGrp->SetButton(kRAYL);       
    
  fMain->AddFrame(pFileHorFrm=new TGHorizontalFrame(fMain,100,200));
  pFileHorFrm->AddFrame(pCreateBtn=new TGTextButton(pFileHorFrm,"Create"));
  pCreateBtn->Connect("Clicked()","KirConfig",this,"CreateConfigFile()");                                 
  
  fMain->MapSubwindows();   
  fMain->Layout();
  fMain->SetWindowName("Create AliROOT Config script");
  fMain->MapWindow();      
}//KirCondig::ctor
          
void KirConfig::CreateConfigFile()
{   
  char *sName="ConfigRich.C";
  FILE *fp=fopen(sName,"w"); if(!fp){Info("CreateConfigFile","Cannot open output file:%sn",sName);return;}
  
  fprintf(fp,"void Config()\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  ::Info(\"RICH private config\",\"Start\");\n\n"); 
//Random
  fprintf(fp,"  gRandom->SetSeed(123456);//put 0 to use system time\n\n");    
//Geant  
  fprintf(fp,"  gSystem->Load(\"libgeant321\");\n");
  fprintf(fp,"  new TGeant3(\"C++ Interface to Geant3\");\n\n");
//File
  fprintf(fp,"  AliRunLoader *pAL=AliRunLoader::Open(\"galice.root\",AliConfig::fgkDefaultEventFolderName,\"recreate\");\n");    
  fprintf(fp,"  pAL->SetCompressionLevel(2);\n");
  fprintf(fp,"  pAL->SetNumberOfEventsPerFile(50);\n");
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
  
  fprintf(fp,"  gMC->SetCut(\"CUTGAM\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"CUTELE\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTNEU\",0.001);\n"); fprintf(fp,"  gMC->SetCut(\"CUTHAD\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"CUTMUO\",0.001); ");  fprintf(fp,"  gMC->SetCut(\"BCUTE\" ,0.001);\n");
  fprintf(fp,"  gMC->SetCut(\"BCUTM\" ,0.001); ");  fprintf(fp,"  gMC->SetCut(\"DCUTE\" ,0.001); ");
  fprintf(fp,"  gMC->SetCut(\"DCUTM\" ,0.001);\n");  fprintf(fp,"  gMC->SetCut(\"PPCUTM\",0.001); ");
  fprintf(fp,"  gMC->SetCut(\"TOFMAX\",1e10);\n\n"); 
//Field
  if(fMagFldChkBtn->GetState()==kButtonDown) fprintf(fp,"  gAlice->SetField();\n\n");else fprintf(fp,"  gAlice->SetField(0);\n\n");
  fprintf(fp,"  pAL->CdGAFile();\n\n");                                 //????       
//BODY-ALIC 
  fprintf(fp,"  new AliBODY(\"BODY\",\"Alice envelop\");\n\n");
//RICH  
  if(fRichTopChkBtn->GetState()==kButtonDown) fprintf(fp,"  AliRICHParam::AngleRot(0);\n");
  switch(fRichVersionCombo->GetSelected()){//RICH
    case 0:   fprintf(fp,"  pRICH=new AliRICHv0(\"RICH\",\"RICH version 0\");\n\n"); break;   
    case 1:   fprintf(fp,"  pRICH=new AliRICHv1(\"RICH\",\"RICH version 1\");\n\n"); break;   
    case 3:   fprintf(fp,"  pRICH=new AliRICHv3(\"RICH\",\"RICH version 3\");\n\n"); break;   
  }
//Generator
  switch(fGenTypeCombo->GetSelected()){
    case kHijingPara:   fprintf(fp,"  pGenerator->AddGenerator(new AliGenHIJINGpara(85700),\"HIJING para\",1);\n");  break;
    case kGun7:   break;
    case kGun1:     
      fprintf(fp,"  AliGenFixed *pGen=new AliGenFixed(1);\n");  
      fprintf(fp,"  pGen->SetPart(%i); pGen->SetMomentum(%3.1f); pGen->SetOrigin(0,0,0);\n",fPartIdCombo->GetSelected(),float(fMomCombo->GetSelected())/10);  
      fprintf(fp,"  pGen->SetPhiRange(pRICH->C(4)->PhiD()); pGen->SetThetaRange(pRICH->C(4)->ThetaD()-2);\n");            
      fprintf(fp,"  pGen->Init();\n");
      break;    
  }
//Other detectors                  
  if(IsDetectorOn(kMAG))  fprintf(fp,"\n  new AliMAG(\"MAG\",\"Magnet\");\n");
  if(IsDetectorOn(kABSO)) fprintf(fp,"\n  new AliABSOv0(\"ABSO\",\"Muon absorber\");\n");
  if(IsDetectorOn(kDIPO)) fprintf(fp,"\n  new AliDIPOv2(\"DIPO\",\"Dipole version 2\");\n");
  if(IsDetectorOn(kHALL)) fprintf(fp,"\n  new AliHALL(\"HALL\",\"Alice Hall\");\n");
  if(IsDetectorOn(kFRAME))fprintf(fp,"\n  AliFRAMEv2 *pFrame=new AliFRAMEv2(\"FRAME\",\"Space Frame\"); pFrame->SetHoles(1);\n");
  if(IsDetectorOn(kSHILD))fprintf(fp,"\n  new AliSHILv2(\"SHIL\",\"Shielding Version 2\");\n");
  if(IsDetectorOn(kPIPE)) fprintf(fp,"\n  new AliPIPEv0(\"PIPE\",\"Beam Pipe\");\n");
  
  if(IsDetectorOn(kITS)){
    fprintf(fp,"\n  AliITSvPPRasymmFMD *pIts =new AliITSvPPRasymmFMD(\"ITS\",\"ITS PPR detailed version\");\n");
    fprintf(fp,"  pIts->SetMinorVersion(2); pIts->SetReadDet(kTRUE);\n");
    fprintf(fp,"  pIts->SetThicknessDet1(200.); pIts->SetThicknessDet2(200.);\n");
    fprintf(fp,"  pIts->SetThicknessChip1(200.); pIts->SetThicknessChip2(200.);\n");
    fprintf(fp,"  pIts->SetRails(0); pIts->SetCoolingFluid(1);\n");
    fprintf(fp,"  pIts->SetEUCLID(0);\n");
  }
  
  if(IsDetectorOn(kTPC)){fprintf(fp,"\n  AliTPC *pTpc=new AliTPCv2(\"TPC\",\"Default\"); pTpc->SetSecAU(-1);pTpc->SetSecAL(-1);\n");}
  if(IsDetectorOn(kZDC))fprintf(fp,"\n  new AliZDCv2(\"ZDC\",\"normal ZDC\");\n");
  
  if(IsDetectorOn(kTRD)){
    fprintf(fp,"\n  AliTRD *pTrd=new AliTRDv1(\"TRD\",\"TRD slow simulator\");\n");
    fprintf(fp,"  pTrd->SetGasMix(1); pTrd->SetPHOShole(); pTrd->SetRICHhole();pTrd->CreateTR();\n");
  }
  
  if(IsDetectorOn(kFMD))   fprintf(fp,"\n  new AliFMDv1(\"FMD\",\"normal FMD\");\n");
  if(IsDetectorOn(kMUON))  fprintf(fp,"\n  new AliMUONv1(\"MUON\",\"default\");\n");
  if(IsDetectorOn(kPHOS))  fprintf(fp,"\n  new AliPHOSv1(\"PHOS\",\"IHEP\");\n");
  if(IsDetectorOn(kPMD))   fprintf(fp,"\n  new AliPMDv1(\"PMD\",\"normal PMD\");\n");
  if(IsDetectorOn(kSTART)) fprintf(fp,"\n  new AliSTARTv1(\"START\",\"START Detector\");\n");
  if(IsDetectorOn(kEMCAL)) fprintf(fp,"\n  new AliEMCALv1(\"EMCAL\",\"G56_2_55_19_104_14\");\n");
  if(IsDetectorOn(kCRT))   fprintf(fp,"\n  new AliCRTv0(\"CRT\",\"normal ACORDE\");\n");
  if(IsDetectorOn(kVZERO)) fprintf(fp,"\n  new AliVZEROv2(\"VZERO\",\"normal VZERO\");\n");

  fprintf(fp,"\n  ::Info(\"RICH private config\",\"Stop\");\n"); 
  fprintf(fp,"}\n");
  fclose(fp);
  gApplication->Terminate(0);
}//CreateConfigFile

KirConfig *p;
void CreateConfig()
{
   p=new KirConfig;
}   

