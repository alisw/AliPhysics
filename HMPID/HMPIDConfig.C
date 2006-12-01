#include <TPDGCode.h>

class RichConfig : public TGMainFrame
{
RQ_OBJECT()
public:
          RichConfig(const char*sFileName);
         ~RichConfig()                    {Info("ctor","");Cleanup();}
         
  enum EVersOpts  {kNo=101,kVer0,kVer1,kVer2,kTest, kDeclust=301,kSagita,kFeedback,kSecRad,kQe0=400,kQeNorm,kOptics};
  enum EGenTypes  {kGunZ=1,kGun1,kGun7,kBox,kHijing,kHijingPara,kPythia,kRichLib,kNotUsed=999};
  
  enum EDetectors {kPIPE=1,kITS,kTPC,kTRD,kTOF,kFRAME,kMAG,kACORDE,kHALL,kPHOS,kT0,kFMD,kABSO,kPMD,kDIPO,kEMCAL,kVZERO,kMUON,kZDC,kSHILD};
  enum EProcesses {kDCAY=1,kPAIR,kCOMP,kPHOT,kPFIS,kDRAY,kANNI,kBREM,kMUNU,kCKOV,kHADR,kLOSS,kMULS,kRAYL,kALL};
  enum EBatchFlags{kPrim=555,kTransport,kAll,kOnly,kRawDdl,kRawDate,kRawRoot,kVertex,kTrack,kHlt,kEsd,kAlign};
  enum EMagField  {kFld0,kFld2,kFld4,kFld5,kFld_2,kFld_4,kFld_5};
  
  Float_t Eta2Theta      (Float_t arg)  const{return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));}
  
  void GuiRich (TGHorizontalFrame *pMainF);  void WriteRich (FILE *pF);  void RichSlot   (Int_t ver);
  void GuiPhys (TGHorizontalFrame *pMainF);  void WritePhys (FILE *pF);  void PhysAddSlot(Int_t id); void PhysRemSlot(Int_t id);
  void GuiGen  (TGCompositeFrame  *pMainF);  void WriteGen  (FILE *pF);  void GenAddSlot (Int_t id); void GenRemSlot (Int_t id);
  void GuiDet  (TGHorizontalFrame *pMainF);  void WriteDet  (FILE *pF);  void DetAddSlot (Int_t id); void DetRemSlot (Int_t id);
  void GuiBatch(TGHorizontalFrame *pMainF);  void WriteBatch(        );  void SlotBatch  (Int_t id);
  
  void    WriteConfig();
  void    ExitSlot();
     
  TGButtonGroup *fVerBG,*fOptBG,*fQeBG; //Rich widgets  
  TGButtonGroup *fMagBG;                //mag field widgets   
  TGButton      *fDecayerB;             //external decayer widgets
  
  TGComboBox        *fGenPidCO,*fGenPminCO,*fGenPmaxCO,*fGenChamCO,*fGenNprimCO; //generator widgets combos
  TGCompositeFrame  *fGenF;                                                     //generator widgets frames     
  TGButtonGroup     *fGenBG;                                                    //generator widgets
  
  TGButtonGroup *fDetBG,*fProcBG;                                      //detectors and processes widgets
  TGButtonGroup *fSimBG,*fSDigBG,*fDigBG,*fRawBG,*fClusBG,*fRecoBG;    //batch control widgets

   char         *fFileName;
};//class RichConfig
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RichConfig::RichConfig(const char *sFileName):TGMainFrame(gClient->GetRoot(),700,400)
{
  fFileName=sFileName;

  AddFrame(pMainF=new TGHorizontalFrame(this,100,200)); //main horizontal frame
  
  GuiRich (pMainF);
  GuiGen  (pMainF);   
  GuiDet  (pMainF);
  GuiPhys (pMainF);
  GuiBatch(pMainF);  
  
//File    
  AddFrame(pFileHF=new TGHorizontalFrame(this,100,200));
  pFileHF->AddFrame(pCrtTB=new TGTextButton(pFileHF,"Create"));
  pCrtTB->Connect("Clicked()","RichConfig",this,"ExitSlot()");                                 
  pFileHF->AddFrame(new TGLabel(pFileHF,Form(" config file as %s",fFileName)));  
  
  MapSubwindows();   
  Resize(GetDefaultSize());
  SetWindowName("Create AliROOT scripts");
  MapWindow(); 
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GuiRich(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(pLeftVF=new TGVerticalFrame(pMainHF,100,200));  
  pLeftVF->AddFrame(pRichGF=new TGGroupFrame(pLeftVF,"HMPID"));
  pRichGF->AddFrame(fVerBG=new TGButtonGroup(pRichGF,""));  fVerBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"RichVerSlot(Int_t)");
    new TGRadioButton(fVerBG,   "No"         ,kNo       ));   
    new TGRadioButton(fVerBG,   "ver0"       ,kVer0     ));
    new TGRadioButton(fVerBG,   "ver1"       ,kVer1     ));  fVerBG->SetButton(kVer1);
    new TGRadioButton(fVerBG,   "ver2"       ,kVer2     ));
  pRichGF->AddFrame(fOptBG=new TGButtonGroup(pRichGF,""));  fOptBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"RichVerSlot(Int_t)");
    new TGCheckButton(fOptBG,"Test beam position?"  ,kTest);       
    new TGCheckButton(fOptBG,"Second radiator?"     ,kSecRad);       
    new TGCheckButton(fOptBG,"Decluster?"           ,kDeclust);      fOptBG->SetButton(kDeclust);
    new TGCheckButton(fOptBG,"Wire sagita?"         ,kSagita);       fOptBG->SetButton(kSagita);
    new TGCheckButton(fOptBG,"Feedbacks?"           ,kFeedback);     fOptBG->SetButton(kFeedback); 
    new TGCheckButton(fOptBG,"Plot optics?"         ,kOptics);     
  pRichGF->AddFrame(fQeBG=new TGButtonGroup(pRichGF,""));
    new TGRadioButton(fQeBG,"QE=0"                 ,kQe0);       
    new TGRadioButton(fQeBG,"QE normal"            ,kQeNorm);       fQeBG->SetButton(kQeNorm);
  pLeftVF->AddFrame(fDecayerB=new TGCheckButton(pLeftVF,"Puthia decayer On/Off"));  fDecayerB->SetState(kButtonDown);
  pLeftVF->AddFrame(fMagBG=new TGButtonGroup(pLeftVF,"Mag field")); 
    new TGRadioButton(fMagBG,   "0.5 T"      ,kFld5   ));
    new TGRadioButton(fMagBG,   "0.4 T"      ,kFld4   ));
    new TGRadioButton(fMagBG,   "0.2 T"      ,kFld2   ));  fMagBG->SetButton(kFld2);
    new TGRadioButton(fMagBG,   "0 T"        ,kFld0   ));   
    new TGRadioButton(fMagBG,   "-0.2 T"     ,kFld_2  ));
    new TGRadioButton(fMagBG,   "-0.4 T"     ,kFld_4  ));
    new TGRadioButton(fMagBG,   "-0.5 T"     ,kFld_5  ));
}//Rich()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::RichVerSlot(Int_t ver)  
{
  if(ver==kTest)    {fMagBG->SetButton(kFld0); fGenBG->SetButton(kGunZ); }
} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WriteRich(FILE *pF)
{
  if(!fVerBG->GetButton(kNo)->GetState()){
    TString title;
    if( fOptBG->GetButton(kSecRad)->GetState())             title+=" Radiator2 ";
    if(!fOptBG->GetButton(kSagita)->GetState())             fprintf(pF,"  AliHMPIDParam::fgIsWireSagita=kFALSE;\n");
    if( fOptBG->GetButton(kTest)  ->GetState())             title+=" TestBeam ";
    if( fOptBG->GetButton(kOptics)->GetState())             title+=" ShowOptics ";
    if(title.Length()==0) title="Default";
    
    if     (fVerBG->GetButton(kVer0)->GetState())           fprintf(pF,"  new AliHMPIDv0(\"Gel %s\");\n\n",title.Data());    
    else if(fVerBG->GetButton(kVer1)->GetState())           fprintf(pF,"  new AliHMPIDv1(\"HMPID\",\"%s\");\n\n",title.Data());   
    else if(fVerBG->GetButton(kVer2)->GetState())           fprintf(pF,"  new AliHMPIDv2(\"Tic %s\");\n\n",title.Data());   
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GuiPhys(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(fProcBG=new TGButtonGroup(pMainHF,"Procs"));
  fProcBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"PhysAddSlot(Int_t)");
  fProcBG->Connect("Released(Int_t)","RichConfig",this,"PhysRemSlot(Int_t)");
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
    new TGCheckButton(fProcBG,"RAYL Rayleigh scattering"    ,kRAYL));  fProcBG->SetButton(kRAYL);       
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::PhysAddSlot(Int_t id)
{//slot is invoked when any of physics check button is pressed
  if(id==kALL)  //if the pressed check button is kALL then switch on all the buttons 
    for(int i=1;i<=fProcBG->GetCount();i++) 
      fProcBG->SetButton(i,kButtonDown);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  RichConfig::PhysRemSlot(Int_t id)
{//slot is invoked when any of physics check button is released
  if(id==kALL) //if the released check button is kALL then switch off all the buttons 
    for(int i=1;i<=fProcBG->GetCount();i++) 
      fProcBG->SetButton(i,kButtonUp);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WritePhys(FILE *pF)
{
  if(fProcBG->GetButton(kDCAY)->GetState()) fprintf(pF,"  gMC->SetProcess(\"DCAY\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"DCAY\",0);  ");
  if(fProcBG->GetButton(kPAIR)->GetState()) fprintf(pF,"  gMC->SetProcess(\"PAIR\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"PAIR\",0);  ");
  if(fProcBG->GetButton(kCOMP)->GetState()) fprintf(pF,"  gMC->SetProcess(\"COMP\",1);\n");else fprintf(pF,"  gMC->SetProcess(\"COMP\",0);\n");
  if(fProcBG->GetButton(kPHOT)->GetState()) fprintf(pF,"  gMC->SetProcess(\"PHOT\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"PHOT\",0);  ");
  if(fProcBG->GetButton(kPFIS)->GetState()) fprintf(pF,"  gMC->SetProcess(\"PFIS\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"PFIS\",0);  ");
  if(fProcBG->GetButton(kDRAY)->GetState()) fprintf(pF,"  gMC->SetProcess(\"DRAY\",1);\n");else fprintf(pF,"  gMC->SetProcess(\"DRAY\",0);\n");
  if(fProcBG->GetButton(kANNI)->GetState()) fprintf(pF,"  gMC->SetProcess(\"ANNI\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"ANNI\",0);  ");
  if(fProcBG->GetButton(kBREM)->GetState()) fprintf(pF,"  gMC->SetProcess(\"BREM\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"BREM\",0);  ");
  if(fProcBG->GetButton(kMUNU)->GetState()) fprintf(pF,"  gMC->SetProcess(\"MUNU\",1);\n");else fprintf(pF,"  gMC->SetProcess(\"MUNU\",0);\n");
  if(fProcBG->GetButton(kCKOV)->GetState()) fprintf(pF,"  gMC->SetProcess(\"CKOV\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"CKOV\",0);  ");
  if(fProcBG->GetButton(kHADR)->GetState()) fprintf(pF,"  gMC->SetProcess(\"HADR\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"HADR\",0);  ");
  if(fProcBG->GetButton(kLOSS)->GetState()) fprintf(pF,"  gMC->SetProcess(\"LOSS\",2);\n");else fprintf(pF,"  gMC->SetProcess(\"LOSS\",0);\n");
  if(fProcBG->GetButton(kMULS)->GetState()) fprintf(pF,"  gMC->SetProcess(\"MULS\",1);  ");else fprintf(pF,"  gMC->SetProcess(\"MULS\",0);  ");
  if(fProcBG->GetButton(kRAYL)->GetState()) fprintf(pF,"  gMC->SetProcess(\"RAYL\",1);\n");else fprintf(pF,"  gMC->SetProcess(\"RAYL\",0);\n");
  
  fprintf(pF,"\n");  
  fprintf(pF,"  gMC->SetCut(\"CUTGAM\",0.001);  "); fprintf(pF,"  gMC->SetCut(\"CUTELE\",0.001);  "); fprintf(pF,"  gMC->SetCut(\"CUTNEU\",0.001);\n"); 
  fprintf(pF,"  gMC->SetCut(\"CUTHAD\",0.001);  "); fprintf(pF,"  gMC->SetCut(\"CUTMUO\",0.001);  "); fprintf(pF,"  gMC->SetCut(\"BCUTE\" ,0.001);\n");
  fprintf(pF,"  gMC->SetCut(\"BCUTM\" ,0.001);  "); fprintf(pF,"  gMC->SetCut(\"DCUTE\" ,0.001);  "); fprintf(pF,"  gMC->SetCut(\"DCUTM\" ,0.001);\n"); 
  fprintf(pF,"  gMC->SetCut(\"PPCUTM\",0.001);  "); fprintf(pF,"  gMC->SetCut(\"TOFMAX\",1e10);\n\n"); 
}//WritePhys()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GuiGen(TGCompositeFrame *pMainF)
{//generator configurator implemented as group of radio buttons
  pMainF->AddFrame(fGenF=new TGVerticalFrame(pMainF));//add generator vertical frame to horizontal main frame
  fGenF->AddFrame(fGenBG=new TGButtonGroup(fGenF,"Gene"));//add type button group to vertical generator frame
  fGenBG->Connect("Pressed(Int_t)","RichConfig",this,"GenAddSlot(Int_t)"); fGenBG->Connect("Released(Int_t)","RichConfig",this,"GenRemSlot(Int_t)");
    new TGCheckButton(fGenBG,"gun along z"           ,kGunZ);
    new TGCheckButton(fGenBG,"gun to 1 chamber"      ,kGun1);
    new TGCheckButton(fGenBG,"gun to 7 chambers"     ,kGun7);
    new TGCheckButton(fGenBG,"box HMPID phase space"  ,kBox );  
    new TGCheckButton(fGenBG,"HIJING"                ,kHijing);
    new TGCheckButton(fGenBG,"HIJING para"           ,kHijingPara);
    new TGCheckButton(fGenBG,"Pythia"                ,kPythia);
    new TGCheckButton(fGenBG,"HMPID lib"              ,kRichLib);
//N prims for box and Hijing para    
  fGenF->AddFrame(fGenNprimCO=new TGComboBox(fGenF,100)); //add N prims combo to generator vertical frame
  fGenNprimCO->AddEntry("not used"    ,kNotUsed);
  fGenNprimCO->AddEntry("N prim=1"    ,1);
  fGenNprimCO->AddEntry("N prim=2"    ,2);
  fGenNprimCO->AddEntry("N prim=5"    ,5);
  fGenNprimCO->AddEntry("N prim=100"  ,100);
  fGenNprimCO->AddEntry("N prim=500"  ,500);
  fGenNprimCO->AddEntry("N prim=1000" ,1000);
  fGenNprimCO->AddEntry("N prim=10000",10000);
  fGenNprimCO->AddEntry("N prim=80000",80000);  fGenNprimCO->Resize(160,20);   fGenNprimCO->Select(kNotUsed);
//PID    
  fGenF->AddFrame(fGenPidCO=new TGComboBox(fGenF,100)); //add pid combo to generator vertical frame
  fGenPidCO->AddEntry("not used"   ,kNotUsed);
  fGenPidCO->AddEntry("electron"   ,kElectron);
  fGenPidCO->AddEntry("positron"   ,kPositron);
  fGenPidCO->AddEntry("pion +"     ,kPiPlus);
  fGenPidCO->AddEntry("pion -"     ,kPiMinus);
  fGenPidCO->AddEntry("kaon +"     ,kKPlus);
  fGenPidCO->AddEntry("kaon -"     ,kKMinus);
  fGenPidCO->AddEntry("kaon 0"     ,kK0Short);    
  fGenPidCO->AddEntry("proton"     ,kProton);
  fGenPidCO->AddEntry("antiproton" ,kProtonBar);
  fGenPidCO->AddEntry("lambda"     ,kLambda0);
  fGenPidCO->AddEntry("antilambda" ,kLambda0Bar);  fGenPidCO->Resize(160,20); fGenPidCO->Select(kNotUsed);
//Pmin  
  fGenF->AddFrame(fGenPminCO=new TGComboBox(fGenF,100)); //add Pmin combo to generator vertical frame  
  fGenPminCO->AddEntry("not used",kNotUsed);
  for(Int_t i=5;i<=295;i+=5) fGenPminCO->AddEntry(Form("Pmin=%3.1f GeV",0.1*i), i); fGenPminCO->Resize(160,20); fGenPminCO->Select(kNotUsed);
//Pmax  
  fGenF->AddFrame(fGenPmaxCO=new TGComboBox(fGenF,100)); //add Pmax combo to generator vertical frame
  fGenPmaxCO->AddEntry("not used",kNotUsed);
  for(Int_t i=10;i<=295;i+=5) fGenPmaxCO->AddEntry(Form("Pmax=%3.1f GeV",0.1*i), i);  fGenPmaxCO->Resize(160,20); fGenPmaxCO->Select(kNotUsed);
//Chamber number  
  fGenF->AddFrame(fGenChamCO=new TGComboBox(fGenF,100)); //add chamber number combo to generator vertical frame
  fGenChamCO->AddEntry("not used",kNotUsed);
  for(int i=1;i<=7;i++) fGenChamCO->AddEntry(Form("Chamber %i",i),i);  fGenChamCO->Resize(160,20); fGenChamCO->Select(kNotUsed);
}//GuiGen()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GenAddSlot(Int_t id)
{//is invoked when any of generator check button is pressed
  if(id==kGunZ){
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO); fGenF->HideFrame(fGenPmaxCO);
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);                                
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  fGenPidCO->Select(kProton);    fGenPminCO->Select(25);
  } 
  if(id==kGun1){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO); fGenF->HideFrame(fGenPmaxCO); 
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);  fGenPidCO->Select(kProton);    fGenPminCO->Select(25); fGenChamCO->Select(4);
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  
  }
  if(id==kGun7){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO);
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);                                fGenPidCO->Select(kProton); 
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  fGenF->HideFrame(fGenPmaxCO); fGenPminCO->Select(25);
                                                   fGenF->HideFrame(fGenChamCO); 
  } 
  if(id==kBox){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);                                fGenNprimCO->Select(500);
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);                                fGenPidCO  ->Select(kProton);
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);                                fGenPminCO ->Select(15);       fGenPmaxCO->Select(15); 
  }
  if(id==kHijing){
    fGenBG->GetButton(kHijing)->ChangeBackground(0xff0000);
    fGenBG->GetButton(kHijingPara)->SetEnabled(kFALSE);
    fGenBG->GetButton(kPythia)->SetEnabled(kFALSE);
  }
  if(id==kHijingPara){
    fGenBG->GetButton(kHijing)->SetEnabled(kFALSE);   fGenNprimCO->Select(500);
    fGenBG->GetButton(kPythia)->SetEnabled(kFALSE);
  }
  if(id==kPythia){
    fGenBG->GetButton(kHijing)->SetEnabled(kFALSE);
    fGenBG->GetButton(kHijingPara)->SetEnabled(kFALSE);
  }
}//GenAddSlot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GenRemSlot(Int_t id)
{//is invoked when any of generator check button is released
  if(id==kGunZ){
    fGenBG->GetButton(kGun1)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);                                                                  
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed);
  } 
  if(id==kGun1){
    fGenBG->GetButton(kGunZ)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed);  
  }
  if(id==kGun7){
    fGenBG->GetButton(kGunZ)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);
    fGenBG->GetButton(kGun1)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed); 
                                             fGenF->ShowFrame(fGenChamCO); 
  } 
  if(id==kBox){
    fGenBG->GetButton(kGunZ)->SetEnabled();                                     fGenNprimCO->Select(kNotUsed);
    fGenBG->GetButton(kGun1)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPminCO ->Select(kNotUsed);   fGenPmaxCO->Select(kNotUsed);
                                                                                fGenChamCO ->Select(kNotUsed);
  }
  if(id==kHijing){
    fGenBG->GetButton(kHijing)->ChangeBackground(0xbebebe);
    fGenBG->GetButton(kHijingPara)->SetEnabled();
    fGenBG->GetButton(kPythia)->SetEnabled();
  }
  if(id==kHijingPara){
    fGenBG->GetButton(kHijing)->SetEnabled();
    fGenBG->GetButton(kPythia)->SetEnabled();
  }
  if(id==kPythia){
    fGenBG->GetButton(kHijing)->SetEnabled();
    fGenBG->GetButton(kHijingPara)->SetEnabled();
  }
}//GenRemSlot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WriteGen(FILE *pF)
{
  Int_t pid=fGenPidCO->GetSelected();
  Float_t pmin=0.1*fGenPminCO->GetSelected(); //particle momentum, GeV
  Float_t pmax=0.1*fGenPmaxCO->GetSelected(); //particle momentum, GeV
  
  fprintf(pF,"  AliGenCocktail *pG=new AliGenCocktail();\n\n");
  
  if(fGenBG->GetButton(kGunZ)->GetState()==kButtonDown)//1 particle along Z axis 
    fprintf(pF,"  AliGenFixed *pGz=new AliGenFixed(1); pGz->SetPart(%i); pGz->SetMomentum(%.1f); pGz->SetOrigin(0,0,-200); pG->AddGenerator(pGz,\"Gz\",1);\n",pid,p);
  
  if(fGenBG->GetButton(kGun1)->GetState()==kButtonDown){//1 gun towards 1 HMPID chamber
    switch(fGenChamCO->GetSelected()){
     case 1: fprintf(pF,"  AliGenFixed *pG1=new AliGenFixed(1); pG1->SetPart(%i); pG1->SetMomentum(%.1f);\n",pid,pmin); 
             fprintf(pF,"               pG1->SetTheta(109.5);   pG1->SetPhi(10);  pG->AddGenerator(pG1,\"g1\",1);\n"); break;
     case 2: fprintf(pF,"  AliGenFixed *pG2=new AliGenFixed(1); pG2->SetPart(%i); pG2->SetMomentum(%.1f);\n",pid,pmin);
	     fprintf(pF,"               pG2->SetTheta( 90.0);   pG2->SetPhi(10);  pG->AddGenerator(pG2,\"g2\",1);\n"); break;
     case 3: fprintf(pF,"  AliGenFixed *pG3=new AliGenFixed(1); pG3->SetPart(%i); pG3->SetMomentum(%.1f);\n",pid,pmin);
	     fprintf(pF,"               pG3->SetTheta(109.5);   pG3->SetPhi(30);  pG->AddGenerator(pG3,\"g3\",1);\n"); break;
     case 4: fprintf(pF,"  AliGenFixed *pG4=new AliGenFixed(1); pG4->SetPart(%i); pG4->SetMomentum(%.1f);\n",pid,pmin);
             fprintf(pF,"               pG4->SetTheta( 90.0);   pG4->SetPhi(30);  pG->AddGenerator(pG4,\"g4\",1);\n"); break;
     case 5: fprintf(pF,"  AliGenFixed *pG5=new AliGenFixed(1); pG5->SetPart(%i); pG5->SetMomentum(%.1f);\n",pid,pmin);
             fprintf(pF,"               pG5->SetTheta( 70.5);   pG5->SetPhi(30);  pG->AddGenerator(pG5,\"g5\",1);\n"); break;
     case 6: fprintf(pF,"  AliGenFixed *pG6=new AliGenFixed(1); pG6->SetPart(%i); pG6->SetMomentum(%.1f);\n",pid,pmin);
             fprintf(pF,"               pG6->SetTheta( 90.0);   pG6->SetPhi(50);  pG->AddGenerator(pG6,\"g6\",1);\n"); break;
     case 7: fprintf(pF,"  AliGenFixed *pG7=new AliGenFixed(1); pG7->SetPart(%i); pG7->SetMomentum(%.1f);\n",pid,pmin);
             fprintf(pF,"               pG7->SetTheta( 70.5);   pG7->SetPhi(50);  pG->AddGenerator(pG7,\"g7\",1);\n"); break;
    }    
  }
  
  if(fGenBG->GetButton(kGun7)->GetState()==kButtonDown){//7 guns towards 7 HMPID chambers
    fprintf(pF,"  AliGenFixed *pG1=new AliGenFixed(1); pG1->SetPart(%i); pG1->SetMomentum(%.1f);pG1->SetTheta(109.5-3); pG1->SetPhi(10);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG1,\"g1\",1);\n");
    fprintf(pF,"  AliGenFixed *pG2=new AliGenFixed(1); pG2->SetPart(%i); pG2->SetMomentum(%.1f);pG2->SetTheta( 90.0-3); pG2->SetPhi(10);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG2,\"g2\",1);\n");
    fprintf(pF,"  AliGenFixed *pG3=new AliGenFixed(1); pG3->SetPart(%i); pG3->SetMomentum(%.1f);pG3->SetTheta(109.5-3); pG3->SetPhi(30);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG3,\"g3\",1);\n");
    fprintf(pF,"  AliGenFixed *pG4=new AliGenFixed(1); pG4->SetPart(%i); pG4->SetMomentum(%.1f);pG4->SetTheta( 90.0-3); pG4->SetPhi(30);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG4,\"g4\",1);\n");
    fprintf(pF,"  AliGenFixed *pG5=new AliGenFixed(1); pG5->SetPart(%i); pG5->SetMomentum(%.1f);pG5->SetTheta( 70.0-3); pG5->SetPhi(30);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG5,\"g5\",1);\n");
    fprintf(pF,"  AliGenFixed *pG6=new AliGenFixed(1); pG6->SetPart(%i); pG6->SetMomentum(%.1f);pG6->SetTheta( 90.0-3); pG6->SetPhi(50);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG6,\"g6\",1);\n");
    fprintf(pF,"  AliGenFixed *pG7=new AliGenFixed(1); pG7->SetPart(%i); pG7->SetMomentum(%.1f);pG7->SetTheta( 70.0-3); pG7->SetPhi(50);\n",pid,pmin); 
    fprintf(pF,"               pG->AddGenerator(pG7,\"g7\",1);\n");
  }  
    
  if(fGenBG->GetButton(kBox)->GetState()==kButtonDown){// box towards HMPID phase space
    fprintf(pF,"  AliGenBox *pB=new AliGenBox(%i);      pB->SetPart(%i);       pB->SetMomentumRange(%.1f,%.1f);\n",(int)fGenNprimCO->GetSelected(),pid,pmin,pmax); 
    fprintf(pF,"             pB->SetThetaRange(65,115); pB->SetPhiRange(5,55); pG->AddGenerator(pB,\"b\",1);\n");
  }     
  
  if(fGenBG->GetButton(kHijing)->GetState()==kButtonDown){//normal HIJING
    fprintf(pF,"  AliGenHijing *pH=new AliGenHijing(-1);           pH->SetEnergyCMS(5500);        pH->SetReferenceFrame(\"CMS\");\n");
    fprintf(pF,"                pH->SetProjectile(\"A\", 208, 82); pH->SetTarget(\"A\", 208, 82); pH->SetJetQuenching(0);\n");      
    fprintf(pF,"                pH->SetShadowing(0);               pH->KeepFullEvent();           pH->SetSelectAll(0);\n");
    fprintf(pF,"                pH->SetImpactParameterRange(0, 5); //fermi\n");
    fprintf(pF,"  pG->AddGenerator(pH,\"h\",1);\n\n");
  }
  
  if(fGenBG->GetButton(kHijingPara)->GetState()==kButtonDown){//parametrized HIJING 
    fprintf(pF,"  AliGenHIJINGpara *pHP=new AliGenHIJINGpara(%i);\n",(int)fGenNprimCO->GetSelected());
    fprintf(pF,"  pHP->SetMomentumRange(0,999); pHP->SetThetaRange(%f,%f); pHP->SetPhiRange(0,360);\n",Eta2Theta(8),Eta2Theta(-8));
    fprintf(pF,"  pG->AddGenerator(pHP,\"hp\",1);\n\n");
  }
      
  if(fGenBG->GetButton(kPythia)->GetState()==kButtonDown){//Pythia
    fprintf(pF,"  AliGenPythia *pP=new AliGenPythia(-1);\n");
    fprintf(pF,"  pP->SetMomentumRange(0,999); pP->SetPhiRange(20,80); pP->SetThetaRange(75,115);\n");
    fprintf(pF,"  pP->SetYRange(-12,12);  pP->SetPtRange(0,1000);      pP->SetStrucFunc(kCTEQ4L);\n");
    fprintf(pF,"  pP->SetProcess(kPyMb);  pP->SetEnergyCMS(14000);\n");      
    fprintf(pF,"  pG->AddGenerator(pP,\"p\",1);\n\n");  
  }
  fprintf(pF,"  pG->Init();\n\n");
}//WriteGenerator()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GuiDet(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(fDetBG=new TGButtonGroup(pMainHF,"Dets"));
  fDetBG->Connect("Pressed(Int_t)" ,"RichConfig",this,"DetAddSlot(Int_t)");
  fDetBG->Connect("Released(Int_t)","RichConfig",this,"DetRemSlot(Int_t)");
    new TGCheckButton(fDetBG,"PIPE"  ,kPIPE));     new TGCheckButton(fDetBG,"ITS"   ,kITS));   new TGCheckButton(fDetBG,"TPC"   ,kTPC));
    new TGCheckButton(fDetBG,"TRD"   ,kTRD));      new TGCheckButton(fDetBG,"TOF"   ,kTOF));   new TGCheckButton(fDetBG,"FRAME" ,kFRAME));         
    new TGCheckButton(fDetBG,"MAG"   ,kMAG));      new TGCheckButton(fDetBG,"ACORDE"   ,kACORDE));   new TGCheckButton(fDetBG,"HALL"  ,kHALL));         
    new TGCheckButton(fDetBG,"PHOS"  ,kPHOS));     new TGCheckButton(fDetBG,"T0" ,kT0)); new TGCheckButton(fDetBG,"FMD"   ,kFMD));         
    new TGCheckButton(fDetBG,"ABSO"  ,kABSO));     new TGCheckButton(fDetBG,"PMD"   ,kPMD));   new TGCheckButton(fDetBG,"DIPO"  ,kDIPO));         
    new TGCheckButton(fDetBG,"EMCAL" ,kEMCAL));    new TGCheckButton(fDetBG,"VZERO" ,kVZERO)); new TGCheckButton(fDetBG,"MUON"  ,kMUON));         
    new TGCheckButton(fDetBG,"ZDC"   ,kZDC));      new TGCheckButton(fDetBG,"SHILD" ,kSHILD));         
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::DetAddSlot(Int_t id)
{//slot is invoked when any detector check button is pressed
  if(id==kTRD || id==kTOF)     //TRD and TOF geometries depend on FRAME 
    fDetBG->SetButton(kFRAME);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::DetRemSlot(Int_t id)
{//slot is invoked when any detector check button is released
  if(id==kFRAME){
    fDetBG->SetButton(kTRD,kFALSE); //switch off TRD and TOF when FRAME is switched off by user hence TRD&TOF depend on FRAME
    fDetBG->SetButton(kTOF,kFALSE);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WriteDet(FILE *pF)
{
//CORE detectors                  
  if(fDetBG->GetButton(kPIPE )->GetState()) fprintf(pF,"\n  new AliPIPEv0(\"PIPE\",\"Beam Pipe\");\n");
  if(fDetBG->GetButton(kSHILD)->GetState()) fprintf(pF,"\n  new AliSHILv2(\"SHIL\",\"Shielding Version 2\");\n");  
  if(fDetBG->GetButton(kITS  )->GetState()){
    fprintf(pF,"\n  AliITSvPPRasymmFMD *pIts =new AliITSvPPRasymmFMD(\"ITS\",\"ITS PPR detailed version\");\n");
    fprintf(pF,"  pIts->SetMinorVersion(2); pIts->SetReadDet(kFALSE);\n");
    fprintf(pF,"  pIts->SetThicknessDet1(200.); pIts->SetThicknessDet2(200.);\n");
    fprintf(pF,"  pIts->SetThicknessChip1(150.); pIts->SetThicknessChip2(150.);\n");
    fprintf(pF,"  pIts->SetRails(0); pIts->SetCoolingFluid(1);\n");
    fprintf(pF,"  pIts->SetEUCLID(0);\n");
  }  
  if(fDetBG->GetButton(kTPC  )->GetState())  fprintf(pF,"\n  new AliTPCv2(\"TPC\",\"Default\");\n");
  if(fDetBG->GetButton(kFRAME)->GetState())  fprintf(pF,"\n  AliFRAMEv2 *pFrame=new AliFRAMEv2(\"FRAME\",\"Space Frame\"); pFrame->SetHoles(1);\n");
  if(fDetBG->GetButton(kTRD  )->GetState())  fprintf(pF,"\n  AliTRD *pTrd=new AliTRDv1(\"TRD\",\"TRD slow simulator\");\n");  
  if(fDetBG->GetButton(kTOF  )->GetState()) fprintf(pF,"\n  new AliTOFv4T0(\"TOF\", \"normal TOF\");\n");
//central detectors  behind HMPID
  if(fDetBG->GetButton(kMAG  )->GetState()) fprintf(pF,"\n  new AliMAG(\"MAG\",\"Magnet\");\n");  
  if(fDetBG->GetButton(kHALL )->GetState()) fprintf(pF,"\n  new AliHALL(\"HALL\",\"Alice Hall\");\n");
//forward detectors  
  if(fDetBG->GetButton(kFMD  )->GetState()) fprintf(pF,"\n  new AliFMDv1(\"FMD\",\"normal FMD\");\n");
  if(fDetBG->GetButton(kABSO )->GetState()) fprintf(pF,"\n  new AliABSOv0(\"ABSO\",\"Muon absorber\");\n");
  if(fDetBG->GetButton(kDIPO )->GetState()) fprintf(pF,"\n  new AliDIPOv2(\"DIPO\",\"Dipole version 2\");\n");
  if(fDetBG->GetButton(kMUON )->GetState()) fprintf(pF,"\n  new AliMUONv1(\"MUON\",\"default\");\n");
  if(fDetBG->GetButton(kPMD  )->GetState()) fprintf(pF,"\n  new AliPMDv1(\"PMD\",\"normal PMD\");\n");
  if(fDetBG->GetButton(kT0)->GetState()) fprintf(pF,"\n  new AliT0v1(\"T0\",\"T0 Detector\");\n");
  if(fDetBG->GetButton(kVZERO)->GetState()) fprintf(pF,"\n  new AliVZEROv2(\"VZERO\",\"normal VZERO\");\n");
  if(fDetBG->GetButton(kZDC  )->GetState()) fprintf(pF,"\n  new AliZDCv2(\"ZDC\",\"normal ZDC\");\n");
//different phase space detectors  
  if(fDetBG->GetButton(kPHOS )->GetState()) fprintf(pF,"\n  new AliPHOSv1(\"PHOS\",\"IHEP\");\n");
  if(fDetBG->GetButton(kEMCAL)->GetState()) fprintf(pF,"\n  new AliEMCALv1(\"EMCAL\",\"G56_2_55_19_104_14\");\n");
  if(fDetBG->GetButton(kACORDE  )->GetState()) fprintf(pF,"\n  new AliACORDEv0(\"ACORDE\",\"normal ACORDE\");\n");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::GuiBatch(TGHorizontalFrame *pMainF)
{
  TGGroupFrame *pSimF=new TGGroupFrame(pMainF,"Simu"); pMainF->AddFrame(pSimF);    
  pSimF->AddFrame(fSimBG=new TGButtonGroup(pSimF,""    ));
    new TGCheckButton(fSimBG,   "Primary"         ,kPrim     ));   fSimBG->SetButton(kPrim);
    new TGCheckButton(fSimBG,   "Transport"       ,kTransport));   fSimBG->SetButton(kTransport);
  pSimF->AddFrame (fSDigBG=new TGButtonGroup(pSimF,""  ));
    new TGRadioButton(fSDigBG,  "No SDigits"      ,kNo       ));   
    new TGRadioButton(fSDigBG,  "SDigits CORE"    ,kAll      ));  
    new TGRadioButton(fSDigBG,  "SDigits HMPID"    ,kOnly     ));   fSDigBG->SetButton(kOnly);
  pSimF->AddFrame (fDigBG=new TGButtonGroup(pSimF,""   ));
    new TGRadioButton(fDigBG,   "No Digits"       ,kNo       ));   
    new TGRadioButton(fDigBG,   "Digits CORE"     ,kAll      ));  
    new TGRadioButton(fDigBG,   "Digits HMPID"     ,kOnly     ));   fDigBG->SetButton(kOnly); 
  pSimF->AddFrame(fRawBG=new TGButtonGroup(pSimF,""    ));
    new TGRadioButton(fRawBG,   "No RAW files"    ,kNo       ));   fRawBG->SetButton(kNo);
    new TGRadioButton(fRawBG,   "RAW DDL"         ,kRawDdl   ));
    new TGRadioButton(fRawBG,   "RAW DATE"        ,kRawDate  ));
    new TGRadioButton(fRawBG,   "RAW ROOT"        ,kRawRoot  ));
    
  TGCompositeFrame *pRecF=new TGGroupFrame(pMainF,"Reco"); pMainF->AddFrame(pRecF);    
  pRecF->AddFrame(fClusBG=new TGButtonGroup(pRecF,""   ));
    new TGRadioButton(fClusBG,  "No Clusters"     ,kNo       ));  
    new TGRadioButton(fClusBG,  "Clusters CORE"   ,kAll      ));  
    new TGRadioButton(fClusBG,  "Clusters HMPID"   ,kOnly     ));   fClusBG->SetButton(kOnly);
    
  pRecF->AddFrame(fRecoBG=new TGButtonGroup(pRecF,""   ));                       fRecoBG->Connect("Pressed(Int_t)","RichConfig",this,"SlotBatch(Int_t)");
    new TGCheckButton(fRecoBG,  "Load Align data" ,kAlign    ));  
    new TGCheckButton(fRecoBG,  "Find ESD tracks" ,kTrack    ));  
    new TGCheckButton(fRecoBG,  "Find vertex"     ,kVertex   ));  
    new TGCheckButton(fRecoBG,  "Fill ESD"        ,kEsd      ));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::SlotBatch(Int_t id)
{//slot is invoked when any button in fRecoBG is pressed
  if(id==kTrack){
    fSDigBG->SetButton(kAll); 
    fDigBG->SetButton(kAll);
    fClusBG->SetButton(kAll);
    fDetBG->SetButton(kITS);
    fDetBG->SetButton(kTPC);
    fDetBG->SetButton(kTRD);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WriteBatch()
{//creates Batch.C file
  char *sBatchName="HMPID Batch";
  FILE *fp=fopen(Form("%s.C",sBatchName),"w"); if(!fp){Info("CreateBatch","Cannot open output file: %s.C",sBatchName);return;}
  
                                                       fprintf(fp,"void %s(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)\n{\n",sBatchName);
                                                       
  if(fSimBG->GetButton(kPrim)->GetState())             fprintf(fp,"  gSystem->Exec(\"rm -rf *.root hlt hough fort* raw* ZZZ*\");\n");
  else                                                 fprintf(fp,"  gSystem->Exec(\"rm -rf  hlt hough raw* ZZZ*\");\n");
  
                                                       fprintf(fp,"  if(isDebug)   AliLog::SetGlobalDebugLevel(AliLog::kDebug);\n");
                                                       fprintf(fp,"  gBenchmark->Start(\"ALICE\");\n  TDatime time;\n\n");
//simulation section  
  if(fSimBG->GetButton(kPrim)->GetState()){
                                                       fprintf(fp,"  AliSimulation     *pSim=new AliSimulation(sConfigFileName);\n");
                                                       
  if(fSimBG->GetButton(kPrim)->GetState())             fprintf(fp,"  pSim->SetRunGeneration(kTRUE);           //initial kinematics generated\n");
  else                                                 fprintf(fp,"  pSim->SetRunGeneration(kFALSE);          //no initial kinematics generated\n");
//transport and hit creations  
  if(fSimBG->GetButton(kTransport)->GetState()) 
                                                       fprintf(fp,"  pSim->SetRunSimulation(kTRUE);                  //transport and hits creation\n");
  else
                                                       fprintf(fp,"  pSim->SetRunSimulation(kFALSE);                 //no transport and hits creation\n");
//sdigits  
  if     (fSDigBG->GetButton(kNo)    ->GetState())     fprintf(fp,"  pSim->SetMakeSDigits(\"\");                     //no sdigits created\n");
  else if(fSDigBG->GetButton(kAll)   ->GetState())     fprintf(fp,"  pSim->SetMakeSDigits(\"ITS TPC TRD TOF HMPID\"); //sdigits created for these detectors\n");
  else if(fSDigBG->GetButton(kOnly)  ->GetState())     fprintf(fp,"  pSim->SetMakeSDigits(\"HMPID\");                 //sdigits created for HMPID only\n");
//digits  
  if     (fDigBG->GetButton(kNo)    ->GetState())      fprintf(fp,"  pSim->SetMakeDigits(\"\");                      //no digits created\n");
  else if(fDigBG->GetButton(kAll)   ->GetState())      fprintf(fp,"  pSim->SetMakeDigits(\"ITS TPC TRD TOF HMPID\");  //digits created for these detectors\n");
  else if(fDigBG->GetButton(kOnly)  ->GetState())      fprintf(fp,"  pSim->SetMakeDigits(\"HMPID\");                  //digits created for HMPID only\n");
//raw data generation  
  if     (fRawBG->GetButton(kRawDdl)->GetState())      fprintf(fp,"  pSim->SetWriteRawData(\"ITS TPC TRD TOF HMPID\");//raw data in DDL  format for these detectors\n");
  else if(fRawBG->GetButton(kRawDate)->GetState())     fprintf(fp,"  pSim->SetWriteRawData(\"ITS TPC TRD TOF HMPID\",\".date\");//raw data in DATE format for these detectors\n");
  else if(fRawBG->GetButton(kRawRoot)->GetState())     fprintf(fp,"  pSim->SetWriteRawData(\"ITS TPC TRD TOF HMPID\",\".root\");//raw data in ROOT format for these detectors\n");
  
                                                       fprintf(fp,"  pSim->Run(iNevents);\n  delete pSim;\n\n");
  }                                                     


//reconstraction section - cluster finder
                                                       fprintf(fp,"  AliReconstruction *pRec=new AliReconstruction;\n");
  if(fRecoBG->GetButton(kAlign)->GetState())           fprintf(fp,"  pRec->SetLoadAlignData(\"ITS TPC TRD TOF HMPID\");          //Misalign data for these detectors\n");     
  else                                                 fprintf(fp,"  pRec->SetLoadAlignData(\"\");                              //Misalign data not loaded\n");     
  
    if     (fClusBG->GetButton(kAll)   ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"ITS TPC TRD TOF HMPID\"); //clusters created for these detectors\n");
    else if(fClusBG->GetButton(kOnly)  ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"HMPID\");                 //clusters created for HMPID only\n");
    else if(fClusBG->GetButton(kNo)    ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"\");                     //clusters are not created\n");
//reconstraction section - vertex finder
    if     (fRecoBG->GetButton(kVertex)->GetState())   fprintf(fp,"  pRec->SetRunVertexFinder(kTRUE);                           //run vertex finder\n");
    else                                               fprintf(fp,"  pRec->SetRunVertexFinder(kFALSE);                          //do not run vertex finder\n");    
//reconstraction section - tracks finder    
    if     (fRecoBG->GetButton(kTrack) ->GetState())   fprintf(fp,"  pRec->SetRunTracking(\"ITS TPC TRD TOF HMPID\");            //run tracking for these detectors\n");
    else                                               fprintf(fp,"  pRec->SetRunTracking(\"\");                                //do not run tracking\n");    
//reconstraction section - fill ESD (prob vector creator)    
    if     (fRecoBG->GetButton(kEsd)   ->GetState())   fprintf(fp,"  pRec->SetFillESD(\"ITS TPC TRD TOF HMPID\");                //run fill ESD (prob vect)\n");      
    else                                               fprintf(fp,"  pRec->SetFillESD(\"\");                                    //do not run fill ESD (prob vect)\n");
                                                       fprintf(fp,"  pRec->Run();delete pRec;\n\n");         
//benchmarks  
                                                       fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/Batch.C>: Start time: \";time.Print();\n");
                                                       fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <my/Batch.C>: Stop  time: \";time.Set();  time.Print();\n");
                                                       fprintf(fp,"  gBenchmark->Show(\"ALICE\");\n");
  
                                                       fprintf(fp,"  gSystem->Exec(\"play -c 2 ~/my/end.wav\");\n");
                                                       fprintf(fp,"  gSystem->Exec(\"touch ZZZ______finished_______ZZZ\");\n}\n");
  fclose(fp);  
}//WriteBatch()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::ExitSlot()
{
//slot to be invoked by clik on the Create button  
  WriteConfig();
  WriteBatch();
  SendCloseMessage();
  gApplication->Terminate(0);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void RichConfig::WriteConfig()
{   
  FILE *pF=fopen(fFileName,"w"); if(!pF){Info("CreateConfigFile","Cannot open output file:%sn",fFileName);return;}
  
  fprintf(pF,"void Config()\n");
  fprintf(pF,"{\n");
  fprintf(pF,"\n  ::Info(\"\\n\\n\\n----------> HMPID private config\",\"Start\");\n"); 
//Random
  fprintf(pF,"  gRandom->SetSeed(123456);//put 0 to use system time\n\n");    
//Geant  
  fprintf(pF,"  gSystem->Load(\"libgeant321\");\n");
  fprintf(pF,"  new TGeant3TGeo(\"C++ Interface to Geant3\");\n\n");
//File
  fprintf(pF,"  AliRunLoader *pAL=AliRunLoader::Open(\"galice.root\",AliConfig::GetDefaultEventFolderName(),\"recreate\");\n");    
  fprintf(pF,"  pAL->SetCompressionLevel(2);\n");
  fprintf(pF,"  pAL->SetNumberOfEventsPerFile(1000);\n");
  fprintf(pF,"  gAlice->SetRunLoader(pAL);\n\n");
//Decayer  
  if(fDecayerB->GetState()==kButtonDown){   
    fprintf(pF,"  TVirtualMCDecayer *pDecayer=new AliDecayerPythia();\n");
    fprintf(pF,"  pDecayer->SetForceDecay(kAll);\n"); 
    fprintf(pF,"  pDecayer->Init();\n"); 
    fprintf(pF,"  gMC->SetExternalDecayer(pDecayer);\n\n");
  }
  WritePhys(pF); //physics processes
  
//Field
       if(fMagBG->GetButton(kFld0)->GetState())     fprintf(pF,"  gAlice->SetField(0);        //no field\n\n");
  else if(fMagBG->GetButton(kFld2)->GetState())     fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,1,10,0));//0.2 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld4)->GetState())     fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,1,10,1));//0.4 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld5)->GetState())     fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,1,10,2));//0.5 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_2)->GetState())    fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,-1,10,0));//-0.2 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_4)->GetState())    fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,-1,10,1));//-0.4 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_5)->GetState())    fprintf(pF,"  gAlice->SetField(new AliMagFMaps(\"Maps\",\"Maps\",2,-1,10,2));//-0.5 Tesla field\n\n");
  
  fprintf(pF,"  pAL->CdGAFile();\n\n");                                 //????       
//Generator 
  WriteGen(pF);//generator
//BODY-ALIC 
  fprintf(pF,"  new AliBODY(\"BODY\",\"Alice envelop\");\n\n");
//HMPID
  WriteRich(pF);  //private HMPID part
  WriteDet(pF);  //other detectors
//end of Config.C file:  
  fprintf(pF,"\n  ::Info(\"----------> HMPID private config\",\"Stop\\n\\n\\n\");\n"); 
  fprintf(pF,"}\n");
  fclose(pF);  
}//WriteConfig()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HMPIDConfig()
{
   new RichConfig("Config.C");
}   
