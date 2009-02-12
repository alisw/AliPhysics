//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Change log:  21st June 2007 by Levente Molnar
//             Detector description upgrade: ITS,ToF,Absorber, Dipole, V0, Emcal
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


#include <TPDGCode.h>

Bool_t gGun1      =kFALSE;
Bool_t gGun7      =kFALSE;
Bool_t gGunZ      =kFALSE;
Bool_t gBox       =kFALSE;
Bool_t gHijing    =kFALSE;
Bool_t gHijingPara=kFALSE;
Bool_t gPythia    =kFALSE;

class HmpConfig : public TGMainFrame
{
RQ_OBJECT()
public:
          HmpConfig(const char*sFileName);
         ~HmpConfig()                    {Info("ctor","");Cleanup();}
         
  enum EVersOpts  {kNo=101,kVer0,kVer1,kVer2,kVer3,kTest, kDeclust=301,kSagita,kFeedback,kElNoise,kQe0=400,kQeNorm,kFlatIdx,kOptics};
  enum EGenTypes  {kGunZ=1,kGun1,kGun7,kBox,kHijing,kHijingPara,kPythia,kHmpLib,kNotUsed=999};
  
  enum EDetectors {kPIPE=1,kITS,kTPC,kTRD,kTOF,kFRAME,kMAG,kACORDE,kHALL,kPHOS,kT0,kFMD,kABSO,kPMD,kDIPO,kEMCAL,kVZERO,kMUON,kZDC,kSHILD};
  enum EProcesses {kDCAY=1,kPAIR,kCOMP,kPHOT,kPFIS,kDRAY,kANNI,kBREM,kMUNU,kCKOV,kHADR,kLOSS,kMULS,kRAYL,kALL};
  enum EBatchFlags{kAll=55,kHmp,kDdl,kDat,kRoo,kVtx,kTrk,kHlt,kPid,kAln,kRecoPar};
  enum EMagField  {kFld0,kFld2,kFld4,kFld5,kFld_2,kFld_4,kFld_5};
  
  Float_t Eta2Theta      (Float_t arg)  const{return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));}
  
  void GuiHmp (TGHorizontalFrame *pMainF);  void WriteHmp (FILE *pF);  void HmpSlot   (Int_t ver);
  void GuiPhys (TGHorizontalFrame *pMainF);  void WritePhys (FILE *pF);  void PhysAddSlot(Int_t); void PhysRemSlot(Int_t);
  void GuiGen  (TGCompositeFrame  *pMainF);  void WriteGen  (FILE *pF);  void GenAddSlot (Int_t); void GenRemSlot (Int_t);
  void GuiDet  (TGHorizontalFrame *pMainF);  void WriteDet  (FILE *pF);  void DetAddSlot (Int_t); void DetRemSlot (Int_t);
  void GuiBatch(TGHorizontalFrame *pMainF);  void WriteBatch(        );  void SlotBatch  (Int_t); void SlotRec    (Bool_t); 
  
  void    WriteConfig();
  void    ExitSlot();
     
  TGButtonGroup *fVerBG,*fOptBG,*fQeBG; //Hmp widgets  
  TGButtonGroup *fMagBG;                //mag field widgets   
  TGButton      *fDecayerB;             //external decayer widgets
  
  TGComboBox        *fGenPidCO,*fGenPminCO,*fGenPmaxCO,*fGenChamCO,*fGenNprimCO; //generator widgets combos
  TGCompositeFrame  *fGenF;                                                     //generator widgets frames     
  TGButtonGroup     *fGenBG;                                                    //generator widgets
  
  TGButtonGroup *fDetBG,*fProcBG;                                               //detectors and processes sections
  TGCheckButton *fSimB,*fTraB;   TGButtonGroup *fSdiBG,*fDigBG,*fRawBG;         //sim section
  TGCheckButton *fRecB;          TGButtonGroup *fInpBG,*fCluBG,*fTrkBG;         //rec section
      //batch control widgets

   char         *fFileName;
};//class HmpConfig
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
HmpConfig::HmpConfig(const char *sFileName):TGMainFrame(gClient->GetRoot(),700,400)
{
  fFileName=sFileName;

  AddFrame(pMainF=new TGHorizontalFrame(this,100,200)); //main horizontal frame
  
  GuiHmp (pMainF);
  GuiGen  (pMainF);   
  GuiDet  (pMainF);
  GuiPhys (pMainF);
  GuiBatch(pMainF);  
  
//File    
  AddFrame(pFileHF=new TGHorizontalFrame(this,100,200));
  pFileHF->AddFrame(pCrtTB=new TGTextButton(pFileHF,"Create"));
  pCrtTB->Connect("Clicked()","HmpConfig",this,"ExitSlot()");                                 
  pFileHF->AddFrame(new TGLabel(pFileHF,Form(" config file as %s",fFileName)));  
  
  MapSubwindows();   
  Resize(GetDefaultSize());
  SetWindowName("Create AliROOT scripts");
  MapWindow(); 
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::GuiHmp(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(pLeftVF=new TGVerticalFrame(pMainHF,100,200));  
  pLeftVF->AddFrame(pHmpGF=new TGGroupFrame(pLeftVF,"HMPID"));
  pHmpGF->AddFrame(fVerBG=new TGButtonGroup(pHmpGF,""));  fVerBG->Connect("Pressed(Int_t)" ,"HmpConfig",this,"HmpVerSlot(Int_t)");
    new TGRadioButton(fVerBG,   "No"         ,kNo       );   
    new TGRadioButton(fVerBG,   "ver0"       ,kVer0     );
    new TGRadioButton(fVerBG,   "ver1"       ,kVer1     );  
    new TGRadioButton(fVerBG,   "ver2"       ,kVer2     ); 
    new TGRadioButton(fVerBG,   "ver3"       ,kVer3     ); fVerBG->SetButton(kVer3);
  pHmpGF->AddFrame(fOptBG=new TGButtonGroup(pHmpGF,""));  fOptBG->Connect("Pressed(Int_t)" ,"HmpConfig",this,"HmpVerSlot(Int_t)");
    new TGCheckButton(fOptBG,"Test run position"   ,kTest);       
    new TGCheckButton(fOptBG,"Unfold cluster    "  ,kDeclust);      fOptBG->SetButton(kDeclust);
    new TGCheckButton(fOptBG,"Wire sagitta      "  ,kSagita);       fOptBG->SetButton(kSagita);
    new TGCheckButton(fOptBG,"Photon feedback   "  ,kFeedback);     fOptBG->SetButton(kFeedback); 
    new TGCheckButton(fOptBG,"Electronic noise  "  ,kElNoise);     // fOptBG->SetButton(kElNoise); 
    new TGCheckButton(fOptBG,"C6F14 N=1.292     "  ,kFlatIdx);     
    new TGCheckButton(fOptBG,"Plot optics       "  ,kOptics);     
  pHmpGF->AddFrame(fQeBG=new TGButtonGroup(pHmpGF,""));
    new TGRadioButton(fQeBG,"QE=0"                 ,kQe0);       
    new TGRadioButton(fQeBG,"QE normal"            ,kQeNorm);       fQeBG->SetButton(kQeNorm);
  pLeftVF->AddFrame(fDecayerB=new TGCheckButton(pLeftVF,"Puthia decayer On/Off"));  fDecayerB->SetState(kButtonDown);
  pLeftVF->AddFrame(fMagBG=new TGButtonGroup(pLeftVF,"Mag field")); 
    new TGRadioButton(fMagBG,   "0.5 T"      ,kFld5   );
    new TGRadioButton(fMagBG,   "0.4 T"      ,kFld4   );
    new TGRadioButton(fMagBG,   "0.2 T"      ,kFld2   );  fMagBG->SetButton(kFld2);
    new TGRadioButton(fMagBG,   "0 T"        ,kFld0   );   
    new TGRadioButton(fMagBG,   "-0.2 T"     ,kFld_2  );
    new TGRadioButton(fMagBG,   "-0.4 T"     ,kFld_4  );
    new TGRadioButton(fMagBG,   "-0.5 T"     ,kFld_5  );
}//Hmp()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::HmpVerSlot(Int_t ver)  
{
  if(ver==kTest)    {fMagBG->SetButton(kFld0); fGenBG->SetButton(kGunZ); }
} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WriteHmp(FILE *pF)
{
  if(!fVerBG->GetButton(kNo)->GetState()){
    TString title;
    if(!fOptBG->GetButton(kSagita)  ->GetState())             fprintf(pF,"  AliHMPIDParam::fgIsWireSagita=kFALSE;\n");
    if(!fOptBG->GetButton(kFeedback)->GetState())             title+=" NoFeedBack ";
    if( fOptBG->GetButton(kElNoise) ->GetState())             fprintf(pF,"  AliHMPIDDigitizer::DoNoise(kTRUE);\n");
    if( fOptBG->GetButton(kTest)    ->GetState())             title+=" TestBeam ";
    if( fOptBG->GetButton(kOptics)  ->GetState())             title+=" ShowOptics ";
    if( fOptBG->GetButton(kFlatIdx) ->GetState())             title+=" FlatIdx ";
    if(title.Length()==0) title="Default";
    
    if     (fVerBG->GetButton(kVer0)->GetState())           fprintf(pF,"  new AliHMPIDv0(\"Gel %s\");\n\n",title.Data());    
    else if(fVerBG->GetButton(kVer1)->GetState())           fprintf(pF,"  new AliHMPIDv1(\"HMPID\",\"%s\");\n\n",title.Data());   
    else if(fVerBG->GetButton(kVer2)->GetState())           fprintf(pF,"  new AliHMPIDv2(\"HMPID\",\"%s\");\n\n",title.Data());   
    else if(fVerBG->GetButton(kVer3)->GetState())           fprintf(pF,"  new AliHMPIDv3(\"HMPID\",\"%s\");\n\n",title.Data());   
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::GuiPhys(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(fProcBG=new TGButtonGroup(pMainHF,"Procs"));
  fProcBG->Connect("Pressed(Int_t)" ,"HmpConfig",this,"PhysAddSlot(Int_t)");
  fProcBG->Connect("Released(Int_t)","HmpConfig",this,"PhysRemSlot(Int_t)");
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
void HmpConfig::PhysAddSlot(Int_t id)
{//slot is invoked when any of physics check button is pressed
  if(id==kALL)  //if the pressed check button is kALL then switch on all the buttons 
    for(int i=1;i<=fProcBG->GetCount();i++) 
      fProcBG->SetButton(i,kButtonDown);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void  HmpConfig::PhysRemSlot(Int_t id)
{//slot is invoked when any of physics check button is released
  if(id==kALL) //if the released check button is kALL then switch off all the buttons 
    for(int i=1;i<=fProcBG->GetCount();i++) 
      fProcBG->SetButton(i,kButtonUp);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WritePhys(FILE *pF)
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
void HmpConfig::GuiGen(TGCompositeFrame *pMainF)
{//generator configurator implemented as group of radio buttons
  pMainF->AddFrame(fGenF=new TGVerticalFrame(pMainF));//add generator vertical frame to horizontal main frame
  fGenF->AddFrame(fGenBG=new TGButtonGroup(fGenF,"Gene"));//add type button group to vertical generator frame
  fGenBG->Connect("Pressed(Int_t)","HmpConfig",this,"GenAddSlot(Int_t)"); fGenBG->Connect("Released(Int_t)","HmpConfig",this,"GenRemSlot(Int_t)");
    new TGCheckButton(fGenBG,"gun along z"           ,kGunZ);
    new TGCheckButton(fGenBG,"gun to 1 chamber"      ,kGun1);
    new TGCheckButton(fGenBG,"gun to 7 chambers"     ,kGun7);
    new TGCheckButton(fGenBG,"box HMPID phase space"  ,kBox );  
    new TGCheckButton(fGenBG,"HIJING"                ,kHijing);
    new TGCheckButton(fGenBG,"HIJING para"           ,kHijingPara);
    new TGCheckButton(fGenBG,"Pythia"                ,kPythia);
    new TGCheckButton(fGenBG,"HMPID lib"              ,kHmpLib);
//N prims for box and Hijing para    
  fGenF->AddFrame(fGenNprimCO=new TGComboBox(fGenF,100)); //add N prims combo to generator vertical frame
  fGenNprimCO->AddEntry("not used"    ,kNotUsed);
  fGenNprimCO->AddEntry("N prim=1"    ,1);
  fGenNprimCO->AddEntry("N prim=2"    ,2);
  fGenNprimCO->AddEntry("N prim=5"    ,5);
  fGenNprimCO->AddEntry("N prim=10"   ,10);
  fGenNprimCO->AddEntry("N prim=20"   ,20);
  fGenNprimCO->AddEntry("N prim=50"   ,50);
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
void HmpConfig::GenAddSlot(Int_t id)
{//is invoked when any of generator check button is pressed
  if(id==kGunZ){
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO); fGenF->HideFrame(fGenPmaxCO);
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);                                
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  fGenPidCO->Select(kProton);    fGenPminCO->Select(25);
                                                   fGenF->HideFrame(fGenChamCO);
    gGunZ=!gGunZ;
  } 
  if(id==kGun1){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO); fGenF->HideFrame(fGenPmaxCO); 
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);  fGenPidCO->Select(kProton);    fGenPminCO->Select(25); fGenChamCO->Select(4);
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  
    gGun1=!gGun1;
  }
  if(id==kGun7){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);  fGenF->HideFrame(fGenNprimCO);
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);                                fGenPidCO->Select(kProton); 
    fGenBG->GetButton(kBox )->SetEnabled(kFALSE);  fGenF->HideFrame(fGenPmaxCO); fGenPminCO->Select(25);
                                                   fGenF->HideFrame(fGenChamCO); 
    gGun7=!gGun7;
  } 
  if(id==kBox){
    fGenBG->GetButton(kGunZ)->SetEnabled(kFALSE);                                fGenNprimCO->Select(500);
    fGenBG->GetButton(kGun1)->SetEnabled(kFALSE);                                fGenPidCO  ->Select(kProton);
    fGenBG->GetButton(kGun7)->SetEnabled(kFALSE);                                fGenPminCO ->Select(15);       fGenPmaxCO->Select(15); 
                                                   fGenF->HideFrame(fGenChamCO);
    gBox=!gBox;
  }
  if(id==kHijing){
//    fGenBG->GetButton(kHijing)->ChangeBackground(0xff0000);
    fGenBG->GetButton(kHijingPara)->SetEnabled(kFALSE);
    fGenBG->GetButton(kPythia)->SetEnabled(kFALSE);
    gHijing=!gHijing;
  }
  if(id==kHijingPara){
    fGenBG->GetButton(kHijing)->SetEnabled(kFALSE);   fGenNprimCO->Select(500);
    fGenBG->GetButton(kPythia)->SetEnabled(kFALSE);
    gHijingPara=!gHijingPara;
  }
  if(id==kPythia){
    fGenBG->GetButton(kHijing)->SetEnabled(kFALSE);
    fGenBG->GetButton(kHijingPara)->SetEnabled(kFALSE);
    gPythia=!gPythia;
  }
}//GenAddSlot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::GenRemSlot(Int_t id)
{//is invoked when any of generator check button is released
  if(id==kGunZ&&!gGunZ){
    fGenBG->GetButton(kGun1)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);                                                                  
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed);
  } 
  if(id==kGun1&&!gGun1){
    fGenBG->GetButton(kGunZ)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed);  
                                                                                fGenChamCO ->Select(kNotUsed); 
  }
  if(id==kGun7&&!gGun7){
    fGenBG->GetButton(kGunZ)->SetEnabled();  fGenF->ShowFrame(fGenNprimCO);
    fGenBG->GetButton(kGun1)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kBox )->SetEnabled();  fGenF->ShowFrame(fGenPmaxCO);      fGenPminCO ->Select(kNotUsed); 
                                             fGenF->ShowFrame(fGenChamCO); 
  } 
  if(id==kBox&&!gBox){
    fGenBG->GetButton(kGunZ)->SetEnabled();                                     fGenNprimCO->Select(kNotUsed);
    fGenBG->GetButton(kGun1)->SetEnabled();                                     fGenPidCO  ->Select(kNotUsed);
    fGenBG->GetButton(kGun7)->SetEnabled();                                     fGenPminCO ->Select(kNotUsed);   fGenPmaxCO->Select(kNotUsed);
                                                                                fGenChamCO ->Select(kNotUsed);
  }
  if(id==kHijing&&!gHijing){
//    fGenBG->GetButton(kHijing)->ChangeBackground(0xbebebe);
    fGenBG->GetButton(kHijingPara)->SetEnabled();
    fGenBG->GetButton(kPythia)->SetEnabled();
  }
  if(id==kHijingPara&&!gHijingPara){
    fGenBG->GetButton(kHijing)->SetEnabled();                                   fGenNprimCO->Select(kNotUsed);
    fGenBG->GetButton(kPythia)->SetEnabled();
  }
  if(id==kPythia&&!gPythia){
    fGenBG->GetButton(kHijing)->SetEnabled();
    fGenBG->GetButton(kHijingPara)->SetEnabled();
  }
}//GenRemSlot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WriteGen(FILE *pF)
{
  Int_t pid=fGenPidCO->GetSelected();
  Float_t pmin=0.1*fGenPminCO->GetSelected(); //particle momentum, GeV
  Float_t pmax=0.1*fGenPmaxCO->GetSelected(); //particle momentum, GeV
  
  if(fGenBG->GetButton(kHijing)->GetState()==kButtonDown     ||
     fGenBG->GetButton(kHijingPara)->GetState()==kButtonDown ||
     fGenBG->GetButton(kPythia)->GetState()==kButtonDown) {;}
  else {fprintf(pF,"  AliGenCocktail *pG=new AliGenCocktail();\n\n");}
  
  if(fGenBG->GetButton(kGunZ)->GetState()==kButtonDown)//1 particle along Z axis 
    fprintf(pF,"  AliGenFixed *pGz=new AliGenFixed(1); pGz->SetPart(%i); pGz->SetMomentum(%.1f); pGz->SetOrigin(0,0,-200); pG->AddGenerator(pGz,\"Gz\",1);\n",pid,pmin);
  
  if(fGenBG->GetButton(kGun1)->GetState()==kButtonDown){//1 gun towards 1 HMPID chamber
    switch(fGenChamCO->GetSelected()){
     case 1: fprintf(pF,"  AliGenFixed *pG1=new AliGenFixed(1); pG1->SetPart(%i); pG1->SetMomentum(%.1f);\n",pid,pmin); 
             fprintf(pF,"               pG1->SetTheta(109.5);   pG1->SetPhi(10);  pG->AddGenerator(pG1,\"g1\",1);\n"); break;
     case 2: fprintf(pF,"  AliGenFixed *pG2=new AliGenFixed(1); pG2->SetPart(%i); pG2->SetMomentum(%.1f);\n",pid,pmin);
	     fprintf(pF,"               pG2->SetTheta( 90.0);   pG2->SetPhi(10);  pG->AddGenerator(pG2,\"g2\",1);\n"); break;
     case 3: fprintf(pF,"  AliGenFixed *pG3=new AliGenFixed(1); pG3->SetPart(%i); pG3->SetMomentum(%.1f);\n",pid,pmin);
	     fprintf(pF,"               pG3->SetTheta(109.5);   pG3->SetPhi(30);  pG->AddGenerator(pG3,\"g3\",1);\n"); break;
     case 4: fprintf(pF,"  AliGenFixed *pG4=new AliGenFixed(1); pG4->SetPart(%i); pG4->SetMomentum(%.1f);\n",pid,pmin);
             fprintf(pF,"               pG4->SetTheta( 87.0);   pG4->SetPhi(30);  pG->AddGenerator(pG4,\"g4\",1);\n"); break;
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
    fprintf(pF,"                pH->SetProjectile(\"A\", 208, 82);   pH->SetTarget(\"A\", 208, 82);   pH->SetJetQuenching(0);\n");      
    fprintf(pF,"                pH->SetShadowing(0);               pH->KeepFullEvent();           pH->SetSelectAll(0);\n");
    fprintf(pF,"                pH->SetImpactParameterRange(0, 5); //fermi\n");
    fprintf(pF,"                AliGenerator *gener = (AliGenerator*)pH; \n");
    fprintf(pF,"                gener->SetOrigin(0, 0, 0);    // vertex position \n");
    fprintf(pF,"                gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position \n");
    fprintf(pF,"                gener->SetCutVertexZ(1.);     // Truncate at 1 sigma \n");
    fprintf(pF,"                gener->SetVertexSmear(kPerEvent);  \n");
    fprintf(pF,"                gener->SetTrackingFlag(1); \n");
    fprintf(pF,"                gener->Init(); \n\n");
  }
  
  if(fGenBG->GetButton(kHijingPara)->GetState()==kButtonDown){//parametrized HIJING 
    fprintf(pF,"  AliGenHIJINGpara *pHP=new AliGenHIJINGpara(%i);\n",(int)fGenNprimCO->GetSelected());
    fprintf(pF,"  pHP->SetMomentumRange(0,999); pHP->SetThetaRange(%f,%f); pHP->SetPhiRange(0,360);\n",Eta2Theta(8),Eta2Theta(-8));
    fprintf(pF,"                AliGenerator *gener = (AliGenerator*)pHP; \n");
    fprintf(pF,"                gener->SetOrigin(0, 0, 0);    // vertex position \n");
    fprintf(pF,"                gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position \n");
    fprintf(pF,"                gener->SetCutVertexZ(1.);     // Truncate at 1 sigma \n");
    fprintf(pF,"                gener->SetVertexSmear(kPerEvent);  \n");
    fprintf(pF,"                gener->SetTrackingFlag(1); \n");
    fprintf(pF,"                gener->Init(); \n\n");
  }
      
  if(fGenBG->GetButton(kPythia)->GetState()==kButtonDown){//Pythia
    fprintf(pF,"  AliGenPythia *pP=new AliGenPythia(-1);\n");
    fprintf(pF,"  pP->SetMomentumRange(0,999);\n");
    fprintf(pF,"  pP->SetYRange(-12,12);  pP->SetPtRange(0,1000);      pP->SetStrucFunc(kCTEQ4L);\n");
    fprintf(pF,"  pP->SetProcess(kPyMb);  pP->SetEnergyCMS(14000);\n");      
    fprintf(pF,"                AliGenerator *gener = (AliGenerator*)pP; \n");
    fprintf(pF,"                gener->SetOrigin(0, 0, 0);    // vertex position \n");
    fprintf(pF,"                gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position \n");
    fprintf(pF,"                gener->SetCutVertexZ(1.);     // Truncate at 1 sigma \n");
    fprintf(pF,"                gener->SetVertexSmear(kPerEvent);  \n");
    fprintf(pF,"                gener->SetTrackingFlag(1); \n");
    fprintf(pF,"                gener->Init(); \n\n");
  }
  if(fGenBG->GetButton(kHijing)->GetState()==kButtonDown     ||
     fGenBG->GetButton(kHijingPara)->GetState()==kButtonDown ||
     fGenBG->GetButton(kPythia)->GetState()==kButtonDown) continue;
  else fprintf(pF,"  pG->Init();\n\n");
}//WriteGenerator()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::GuiDet(TGHorizontalFrame *pMainHF)
{
  pMainHF->AddFrame(fDetBG=new TGButtonGroup(pMainHF,"Dets"));
  fDetBG->Connect("Pressed(Int_t)" ,"HmpConfig",this,"DetAddSlot(Int_t)");
  fDetBG->Connect("Released(Int_t)","HmpConfig",this,"DetRemSlot(Int_t)");
    new TGCheckButton(fDetBG,"PIPE"  ,kPIPE));     new TGCheckButton(fDetBG,"ITS"   ,kITS));   new TGCheckButton(fDetBG,"TPC"   ,kTPC));
    new TGCheckButton(fDetBG,"TRD"   ,kTRD));      new TGCheckButton(fDetBG,"TOF"   ,kTOF));   new TGCheckButton(fDetBG,"FRAME" ,kFRAME));         
    new TGCheckButton(fDetBG,"MAG"   ,kMAG));      new TGCheckButton(fDetBG,"ACORDE",kACORDE));new TGCheckButton(fDetBG,"HALL"  ,kHALL));         
    new TGCheckButton(fDetBG,"PHOS"  ,kPHOS));     new TGCheckButton(fDetBG,"T0"    ,kT0));    new TGCheckButton(fDetBG,"FMD"   ,kFMD));         
    new TGCheckButton(fDetBG,"ABSO"  ,kABSO));     new TGCheckButton(fDetBG,"PMD"   ,kPMD));   new TGCheckButton(fDetBG,"DIPO"  ,kDIPO));         
    new TGCheckButton(fDetBG,"EMCAL" ,kEMCAL));    new TGCheckButton(fDetBG,"VZERO" ,kVZERO)); new TGCheckButton(fDetBG,"MUON"  ,kMUON));         
    new TGCheckButton(fDetBG,"ZDC"   ,kZDC));      new TGCheckButton(fDetBG,"SHILD" ,kSHILD));         
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::DetAddSlot(Int_t id)
{//slot is invoked when any detector check button is pressed
  if(id==kTRD || id==kTOF)     //TRD and TOF geometries depend on FRAME 
    fDetBG->SetButton(kFRAME);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::DetRemSlot(Int_t id)
{//slot is invoked when any detector check button is released
  if(id==kFRAME){
    fDetBG->SetButton(kTRD,kFALSE); //switch off TRD and TOF when FRAME is switched off by user hence TRD&TOF depend on FRAME
    fDetBG->SetButton(kTOF,kFALSE);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WriteDet(FILE *pF)
{
//CORE detectors                  
  if(fDetBG->GetButton(kPIPE )->GetState()) fprintf(pF,"\n  new AliPIPEv3(\"PIPE\",\"Beam Pipe\");\n");
  if(fDetBG->GetButton(kSHILD)->GetState()) fprintf(pF,"\n  new AliSHILv3(\"SHIL\",\"Shielding Version 2\");\n");  
  if(fDetBG->GetButton(kITS  )->GetState()){
    fprintf(pF,"\n AliITSvPPRasymmFMD *ITS  = new AliITSvPPRasymmFMD(\"ITS\",\"ITS PPR detailed version with asymmetric services\");\n");
    fprintf(pF,"\n  ITS->SetMinorVersion(2);\n");  
    fprintf(pF,"\n  ITS->SetReadDet(kFALSE);\n");	
    fprintf(pF,"\n  ITS->SetThicknessDet1(200.);\n"); 
    fprintf(pF,"\n  ITS->SetThicknessDet2(200.); \n");
    fprintf(pF,"\n  ITS->SetThicknessChip1(150.);\n");
    fprintf(pF,"\n  ITS->SetThicknessChip2(150.);\n");
    fprintf(pF,"\n  ITS->SetRails(0);\n");	
    fprintf(pF,"\n  ITS->SetCoolingFluid(1);\n");
    fprintf(pF,"\n  ITS->SetEUCLID(0);\n");
  }  
  
  if(fDetBG->GetButton(kTPC  )->GetState())  fprintf(pF,"\n  new AliTPCv2(\"TPC\",\"Default\");\n");
  if(fDetBG->GetButton(kFRAME)->GetState())  fprintf(pF,"\n  AliFRAMEv2 *pFrame=new AliFRAMEv2(\"FRAME\",\"Space Frame\"); pFrame->SetHoles(1);\n");
  if(fDetBG->GetButton(kTRD  )->GetState())  fprintf(pF,"\n  AliTRD *pTrd=new AliTRDv1(\"TRD\",\"TRD slow simulator\");\n");  
  if(fDetBG->GetButton(kTOF  )->GetState())  {
    fprintf(pF,"\n  AliTOF *TOF = new AliTOFv6T0(\"TOF\", \"normal TOF\");\n");
    fprintf(pF,"\n  Int_t TOFSectors[18]={-1,0,0,-1,-1,-1,0,0,-1,0,0,0,0,-1,-1,0,0,0};\n");
    fprintf(pF,"\n  TOF->SetTOFSectors(TOFSectors);\n");
  }
//central detectors  behind HMPID
  if(fDetBG->GetButton(kMAG  )->GetState()) fprintf(pF,"\n  new AliMAG(\"MAG\",\"Magnet\");\n");  
  if(fDetBG->GetButton(kHALL )->GetState()) fprintf(pF,"\n  new AliHALL(\"HALL\",\"Alice Hall\");\n");
//forward detectors  
  if(fDetBG->GetButton(kFMD  )->GetState()) fprintf(pF,"\n  AliFMD *FMD = new AliFMDv1(\"FMD\",\"normal FMD\");\n");
  if(fDetBG->GetButton(kABSO )->GetState()) fprintf(pF,"\n  AliABSO *ABSO = new AliABSOv3(\"ABSO\",\"Muon absorber\");\n");
  if(fDetBG->GetButton(kDIPO )->GetState()) fprintf(pF,"\n  AliDIPO *DIPO = new AliDIPOv3(\"DIPO\",\"Dipole version 3\");\n");
  if(fDetBG->GetButton(kMUON )->GetState()) fprintf(pF,"\n  AliMUON *MUON = new AliMUONv1(\"MUON\",\"default\");\n");
  if(fDetBG->GetButton(kPMD  )->GetState()) fprintf(pF,"\n  AliPMD *PMD = new AliPMDv1(\"PMD\",\"normal PMD\");\n");
  if(fDetBG->GetButton(kT0)->GetState())    fprintf(pF,"\n  AliT0 *T0 = new AliT0v1(\"T0\",\"T0 Detector\");\n");
  if(fDetBG->GetButton(kVZERO)->GetState()) fprintf(pF,"\n  AliVZERO *VZERO = new AliVZEROv7(\"VZERO\",\"normal VZERO\");\n");
  if(fDetBG->GetButton(kZDC  )->GetState()) fprintf(pF,"\n  AliZDC *ZDC = new AliZDCv2(\"ZDC\",\"normal ZDC\");\n");
//different phase space detectors  
  if(fDetBG->GetButton(kPHOS )->GetState()) fprintf(pF,"\n  AliPHOS *PHOS = new AliPHOSv1(\"PHOS\",\"IHEP\");\n");
  if(fDetBG->GetButton(kEMCAL)->GetState()) fprintf(pF,"\n  AliEMCAL *EMCAL = new AliEMCALv2(\"EMCAL\",\"SHISH_77_TRD1_2X2_FINAL_110DEG\");\n");
  if(fDetBG->GetButton(kACORDE  )->GetState()) fprintf(pF,"\n  AliACORDE *ACORDE = new AliACORDEv1(\"ACORDE\",\"normal ACORDE\");\n");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::GuiBatch(TGHorizontalFrame *pMainF)
{
  TGGroupFrame *pSimF=new TGGroupFrame(pMainF,"Simu"); pMainF->AddFrame(pSimF);    
  pSimF->AddFrame(fSimB =new TGCheckButton(pSimF,"Enable"));      fSimB->Connect("Toggled(Bool_t)","HmpConfig",this,"SlotSim(Bool_t)"); 
  pSimF->AddFrame(fTraB =new TGCheckButton(pSimF,"Transport")); 
  pSimF->AddFrame(fSdiBG=new TGButtonGroup(pSimF,""));
  pSimF->AddFrame(fDigBG=new TGButtonGroup(pSimF,""));
  pSimF->AddFrame(fRawBG=new TGButtonGroup(pSimF,""));
    new TGRadioButton(fSdiBG,  "No SDigits"      ,kNo   );   
    new TGRadioButton(fSdiBG,  "SDigits CORE"    ,kAll  );  
    new TGRadioButton(fSdiBG,  "SDigits HMPID"   ,kHmp  );   
    new TGRadioButton(fDigBG,  "No Digits"       ,kNo   );   
    new TGRadioButton(fDigBG,  "Digits CORE"     ,kAll  );  
    new TGRadioButton(fDigBG,  "Digits HMPID"    ,kHmp  );   
    new TGRadioButton(fRawBG,  "No RAW files"    ,kNo   );   
    new TGRadioButton(fRawBG,  "RAW DDL"         ,kDdl  );
    new TGRadioButton(fRawBG,  "RAW DATE"        ,kDat  );
    new TGRadioButton(fRawBG,  "RAW ROOT"        ,kRoo  );
    
    
  TGCompositeFrame *pRecF=new TGGroupFrame(pMainF,"Reco"); pMainF->AddFrame(pRecF);    
  pRecF->AddFrame(fRecB =new TGCheckButton(pRecF,"Enable"));  fRecB->Connect("Toggled(Bool_t)","HmpConfig",this,"SlotRec(Bool_t)"); 
  pRecF->AddFrame(fInpBG=new TGButtonGroup(pRecF,""   ));       
  pRecF->AddFrame(fCluBG=new TGButtonGroup(pRecF,""   ));       
  pRecF->AddFrame(fTrkBG=new TGButtonGroup(pRecF,""   ));     fTrkBG->Connect("Pressed(Int_t)","HmpConfig",this,"SlotBatch(Int_t)");  
    new TGRadioButton(fInpBG,  "From sim"        ,kNo     );  
    new TGRadioButton(fInpBG,  "From DDL"        ,kDdl    );  
    new TGRadioButton(fInpBG,  "From DATE"       ,kDat    );  
    new TGRadioButton(fInpBG,  "From ROOT"       ,kRoo    );  
    new TGRadioButton(fCluBG,  "No Clusters"     ,kNo     );  
    new TGRadioButton(fCluBG,  "Clusters CORE"   ,kAll    );  
    new TGRadioButton(fCluBG,  "Clusters HMPID"  ,kHmp    );   

    new TGCheckButton(fTrkBG,  "Apply RecoParam" ,kRecoPar);  

    new TGCheckButton(fTrkBG,  "Load Align data" ,kAln    );  
    new TGCheckButton(fTrkBG,  "Prim vertex"     ,kVtx    );  
    new TGCheckButton(fTrkBG,  "ESD tracks"      ,kTrk    );  
    new TGCheckButton(fTrkBG,  "Assign PID"      ,kPid    );
    
  fSimB->SetState(kButtonDown); fTraB->SetState(kButtonDown); fSdiBG->SetButton(kHmp); fDigBG->SetButton(kHmp); fRawBG->SetButton(kNo);
  fRecB->SetState(kButtonDown);                               fInpBG->SetButton(kNo);  fCluBG->SetButton(kHmp); 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::SlotSim(Bool_t isChk)
{
  if(isChk){
    fTraB                  ->SetState(kButtonDown);
    
    fSdiBG->GetButton(kNo )->SetState(kButtonEngaged);
    fSdiBG->GetButton(kAll)->SetState(kButtonEngaged);
    fSdiBG->GetButton(kHmp)->SetState(kButtonDown);
    
    fDigBG->GetButton(kNo )->SetState(kButtonEngaged);
    fDigBG->GetButton(kAll)->SetState(kButtonEngaged);
    fDigBG->GetButton(kHmp)->SetState(kButtonDown);
    
    fRawBG->GetButton(kNo )->SetState(kButtonDown);
    fRawBG->GetButton(kDdl)->SetState(kButtonEngaged);
    fRawBG->GetButton(kDat)->SetState(kButtonEngaged);
    fRawBG->GetButton(kRoo)->SetState(kButtonEngaged);    
  }else{
    fTraB                  ->SetState(kButtonDisabled);
    
    fSdiBG->GetButton(kNo )->SetState(kButtonDisabled);
    fSdiBG->GetButton(kAll)->SetState(kButtonDisabled);
    fSdiBG->GetButton(kHmp)->SetState(kButtonDisabled);
    
    fDigBG->GetButton(kNo )->SetState(kButtonDisabled);
    fDigBG->GetButton(kAll)->SetState(kButtonDisabled);
    fDigBG->GetButton(kHmp)->SetState(kButtonDisabled);
    
    fRawBG->GetButton(kNo )->SetState(kButtonDisabled);
    fRawBG->GetButton(kDdl)->SetState(kButtonDisabled);
    fRawBG->GetButton(kDat)->SetState(kButtonDisabled);
    fRawBG->GetButton(kRoo)->SetState(kButtonDisabled);    
  }  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::SlotRec(Bool_t isChk)
{
  if(isChk){
    fInpBG->GetButton(kNo) ->SetState(kButtonDown);
    fInpBG->GetButton(kDdl)->SetState(kButtonEngaged);
    fInpBG->GetButton(kDat)->SetState(kButtonEngaged);
    fInpBG->GetButton(kRoo)->SetState(kButtonEngaged);
    
    fCluBG->GetButton(kNo) ->SetState(kButtonEngaged);
    fCluBG->GetButton(kAll)->SetState(kButtonEngaged);
    fCluBG->GetButton(kHmp)->SetState(kButtonDown);
    
    fTrkBG->GetButton(kRecoPar)->SetState(kButtonEngaged);
    fTrkBG->GetButton(kAln)->SetState(kButtonEngaged);
    fTrkBG->GetButton(kVtx)->SetState(kButtonEngaged);
    fTrkBG->GetButton(kTrk)->SetState(kButtonEngaged);
    fTrkBG->GetButton(kPid)->SetState(kButtonEngaged);
  }else{
    fInpBG->GetButton(kNo) ->SetState(kButtonDisabled);
    fInpBG->GetButton(kDdl)->SetState(kButtonDisabled);
    fInpBG->GetButton(kDat)->SetState(kButtonDisabled);
    fInpBG->GetButton(kRoo)->SetState(kButtonDisabled);
    
    fCluBG->GetButton(kNo) ->SetState(kButtonDisabled);
    fCluBG->GetButton(kAll)->SetState(kButtonDisabled);
    fCluBG->GetButton(kHmp)->SetState(kButtonDisabled);
    
    fTrkBG->GetButton(kRecoPar)->SetState(kButtonDisabled);
    fTrkBG->GetButton(kAln)->SetState(kButtonDisabled);
    fTrkBG->GetButton(kVtx)->SetState(kButtonDisabled);
    fTrkBG->GetButton(kTrk)->SetState(kButtonDisabled);
    fTrkBG->GetButton(kPid)->SetState(kButtonDisabled);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::SlotBatch(Int_t id)
{//slot is invoked when any button in fTrkBG is pressed
  if(id==kTrk){
    fSdiBG->SetButton(kAll); 
    fDigBG->SetButton(kAll);
    fCluBG->SetButton(kAll);
    fDetBG->SetButton(kITS);
    fDetBG->SetButton(kTPC);
    fDetBG->SetButton(kTRD);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WriteBatch()
{//creates Batch.C file
  TString det;
  if(fDetBG->GetButton(kITS  )->GetState())  det+="ITS ";
  if(fDetBG->GetButton(kTPC  )->GetState())  det+="TPC ";
  if(fDetBG->GetButton(kTRD  )->GetState())  det+="TRD ";
  if(fDetBG->GetButton(kTOF  )->GetState())  det+="TOF ";
  if(!fVerBG->GetButton(kNo)->GetState())    det+="HMPID ";
  char *sBatchName="sim";
  FILE *fp=fopen(Form("%s.C",sBatchName),"w"); if(!fp){Info("CreateSim","Cannot open output file: %s.C",sBatchName);return;}
  
                                                    fprintf(fp,"void %s(Int_t iNevt=1,Bool_t isDbg=kFALSE,char *sCfg=\"Config.C\")\n{\n",sBatchName);
                                                    fprintf(fp,"  gSystem->Exec(\"rm -rf hlt hough gphysi* fort* ZZZ* raw*\");  //remove garbage\n"); 
                                                    fprintf(fp,"  gBenchmark->Start(\"ALICE\"); TDatime time;      //start benchmarking\n\n");
                                                       
                                                    fprintf(fp,"  if(isDbg) AliLog::SetGlobalDebugLevel(AliLog::kDebug);\n");
//simulation section  
  if(fSimB->GetState()){                            fprintf(fp,"  gSystem->Exec(\"rm -rf *.root \");               //remove previous simulation\n"); 
                                                    fprintf(fp,"  gSystem->Exec(\"rm -rf *.date \");               //remove previous simulation\n"); 
                                                    fprintf(fp,"  AliSimulation *pSim=new AliSimulation(sCfg);   //init simulation\n");
                                                       
    if(fTraB->GetState())                           fprintf(fp,"  pSim->SetRunSimulation(kTRUE);                 //transport and hits creation\n");
    else                                            fprintf(fp,"  pSim->SetRunSimulation(kFALSE);                //no transport and hits creation\n");
    
    if     (fSdiBG->GetButton(kNo )->GetState())    fprintf(fp,"  pSim->SetMakeSDigits(\"\");                      //no sdigits\n");
    else if(fSdiBG->GetButton(kAll)->GetState())    fprintf(fp,"  pSim->SetMakeSDigits(\"%s\");                   //sdigits for all\n",det.Data());
    else if(fSdiBG->GetButton(kHmp)->GetState())    fprintf(fp,"  pSim->SetMakeSDigits(\"HMPID\");                 //sdigits for HMPID\n");
    
    if     (fDigBG->GetButton(kNo )->GetState())    fprintf(fp,"  pSim->SetMakeDigits(\"\");                       //no digits\n");
    else if(fDigBG->GetButton(kAll)->GetState())    fprintf(fp,"  pSim->SetMakeDigits(\"%s\");                    //digits for all\n",det.Data());
    else if(fDigBG->GetButton(kHmp)->GetState())    fprintf(fp,"  pSim->SetMakeDigits(\"HMPID\");                  //digits for HMPID\n");
    
    if     (fRawBG->GetButton(kDdl)->GetState())    fprintf(fp,"  pSim->SetWriteRawData(\"%s\");                  //raw data as DDL\n",det.Data());
    else if(fRawBG->GetButton(kDat)->GetState())    fprintf(fp,"  pSim->SetWriteRawData(\"%s\",\"raw.date\");     //raw data as DATE\n",det.Data());
    else if(fRawBG->GetButton(kRoo)->GetState())    fprintf(fp,"  pSim->SetWriteRawData(\"%s\",\"raw.root\");     //raw data as ROOT\n",det.Data());

                                                    fprintf(fp,"  pSim->SetRunHLT(\"\");                           //no HLT stuff\n");   
                                                    fprintf(fp,"  pSim->SetRunQA(\":\");                           //no QA\n");
                                                    fprintf(fp,"  pSim->Run(iNevt);                                //run iNevt events\n  delete pSim;\n\n");
  }//sim section
                                                    fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <sim.C>: Start time: \";time.Print();\n");
                                                    fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <sim.C>: Stop  time: \";time.Set();  time.Print();\n");
                                                    fprintf(fp,"  gBenchmark->Show(\"ALICE\");\n");
  
                                                    fprintf(fp,"  gSystem->Exec(\"touch ZZZ______finished_______SSS\");\n");
                                                    fprintf(fp,"  gSystem->Exec(\"aliroot -q rec.C &\");\n}\n");
  fclose(fp);  
// rec section
  char *sBatchName="rec";
  FILE *fp=fopen(Form("%s.C",sBatchName),"w"); if(!fp){Info("CreateRec","Cannot open output file: %s.C",sBatchName);return;}
  
                                                    fprintf(fp,"void %s()\n{\n",sBatchName);
                                                    fprintf(fp,"  gSystem->Exec(\"rm -rf *RRR \");  //remove garbage\n"); 
  if(fRecB->GetState()){
/*
       if(fMagBG->GetButton(kFld0)->GetState())     fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 0,1,10.,2));//no field\n");
  else if(fMagBG->GetButton(kFld2)->GetState())     fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k2kG));//0.2 Tesla field\n");
  else if(fMagBG->GetButton(kFld4)->GetState())     fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k4kG));//0.4 Tesla field\n");
  else if(fMagBG->GetButton(kFld5)->GetState())     fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k5kG));//0.5 Tesla field\n");
  else if(fMagBG->GetButton(kFld_2)->GetState())    fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k2kG));//-0.2 Tesla field\n");
  else if(fMagBG->GetButton(kFld_4)->GetState())    fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k4kG));//-0.4 Tesla field\n");
  else if(fMagBG->GetButton(kFld_5)->GetState())    fprintf(fp,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k5kG));//-0.5 Tesla field\n");
*/  
                                                    fprintf(fp,"  AliReconstruction *pRec=new AliReconstruction;\n");
                                                    fprintf(fp,"  gBenchmark->Start(\"ALICE\"); TDatime time;      //start benchmarking\n\n");
                                                    
    //---------------------------------------------
    if     (fTrkBG->GetButton(kRecoPar)->GetState())
      { 
         fprintf(fp,"  AliHMPIDRecoParam * hmpidRecoParam = AliHMPIDRecoParam::GetUserModeParam(); //Get the HMPID reco param\n"); 
         fprintf(fp,"  hmpidRecoParam->SetUserCutMode(kTRUE);                                      //Switch to RecoParam\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(0,4);                                            //eg cut for UserCutSigma (Values: ch0)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(1,4);                                            //eg cut for UserCutSigma (Values: ch1)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(2,4);                                            //eg cut for UserCutSigma (Values: ch2)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(3,4);                                            //eg cut for UserCutSigma (Values: ch3)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(4,4);                                            //eg cut for UserCutSigma (Values: ch4)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(5,4);                                            //eg cut for UserCutSigma (Values: ch5)\n");
         fprintf(fp,"  hmpidRecoParam->SetUserCut(6,4);                                            //eg cut for UserCutSigma (Values: ch6)\n");
         fprintf(fp,"  AliHMPIDReconstructor::SetRecoParam(hmpidRecoParam);                        //Pass the RecoPar to the Reconstructor\n");
      }
    //---------------------------------------------                                                    
                                                    
    if     (fInpBG->GetButton(kNo )->GetState())    fprintf(fp,"  pRec->SetInput(\"\");                       //from digits\n");   
    else if(fInpBG->GetButton(kDdl)->GetState())    fprintf(fp,"  pRec->SetInput(\"./\");                     //from raw data in DDL format\n");                                            
    else if(fInpBG->GetButton(kDat)->GetState())    fprintf(fp,"  pRec->SetInput(\"raw.date\");          //from raw data in DATE format\n");                                            
    else if(fInpBG->GetButton(kRoo)->GetState())    fprintf(fp,"  pRec->SetInput(\"raw.root\");               //from raw data in ROOT format\n");                                            
    
    if     (fCluBG->GetButton(kAll) ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"%s\");   //clusters for all detectors\n",det.Data());
    else if(fCluBG->GetButton(kHmp) ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"HMPID\"); //clusters for HMPID only\n");
    else if(fCluBG->GetButton(kNo)  ->GetState())   fprintf(fp,"  pRec->SetRunLocalReconstruction(\"\");      //no clusters\n");
    
    
    
    if     (fTrkBG->GetButton(kAln)->GetState())    fprintf(fp,"  pRec->SetLoadAlignData(\"%s\");            //with misalignment\n",det.Data());     
    else                                            fprintf(fp,"  pRec->SetLoadAlignData(\"\");               //no misalignment\n");     
  
    if     (fTrkBG->GetButton(kVtx)->GetState())    fprintf(fp,"  pRec->SetRunVertexFinder(kTRUE);          //primary vertex\n");
    else                                            fprintf(fp,"  pRec->SetRunVertexFinder(kFALSE);         //no primary vertex\n");    
    
    if     (fTrkBG->GetButton(kTrk)   ->GetState()) fprintf(fp,"  pRec->SetRunTracking(\"%s\");              //tracking\n",det.Data());
    else                                            fprintf(fp,"  pRec->SetRunTracking(\"\");                 //no tracking\n");    
    
    if     (fTrkBG->GetButton(kPid)   ->GetState()) fprintf(fp,"  pRec->SetFillESD(\"%s\");                  //prob vect\n",det.Data());      
    else                                            fprintf(fp,"  pRec->SetFillESD(\"\");                     //no prob vect\n");
    
                                                    fprintf(fp,"  pRec->Run();delete pRec;\n\n");         
  }//rec part                                                       
//benchmarks  
                                                    fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <rec.C>: Start time: \";time.Print();\n");
                                                    fprintf(fp,"  cout<<\"!!!!!!!!!!!!Info in <rec.C>: Stop  time: \";time.Set();  time.Print();\n");
                                                    fprintf(fp,"  gBenchmark->Show(\"ALICE\");\n");
  
                                                    fprintf(fp,"  gSystem->Exec(\"touch ZZZ______finished_______RRR\");\n}\n");
  fclose(fp);  
}//WriteBatch()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::ExitSlot()
{
//slot to be invoked by clik on the Create button  
  if(fSimB->GetState()) WriteConfig();
  WriteBatch();
  SendCloseMessage();
  gApplication->Terminate(0);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpConfig::WriteConfig()
{   
  FILE *pF=fopen(fFileName,"w"); if(!pF){Info("CreateConfigFile","Cannot open output file:%sn",fFileName);return;}
  
  fprintf(pF,"void Config()\n");
  fprintf(pF,"{\n");
  fprintf(pF,"\n  ::Info(\"\\n\\n\\n----------> HMPID private config\",\"Start\");\n"); 
//Random
  fprintf(pF,"  gRandom->SetSeed(123456);//put 0 to use system time\n\n");    
//File
  fprintf(pF,"  AliRunLoader *pAL=AliRunLoader::Open(\"galice.root\",AliConfig::GetDefaultEventFolderName(),\"recreate\");\n");    
  fprintf(pF,"  pAL->SetCompressionLevel(2);\n");
  fprintf(pF,"  pAL->SetNumberOfEventsPerFile(1000);\n");
  fprintf(pF,"  gAlice->SetRunLoader(pAL);\n\n");
//Decayer  
  if(fDecayerB->GetState()==kButtonDown){
    fprintf(pF,"  gSystem->Load(\"liblhapdf.so\"); // Parton density functions \n");
    fprintf(pF,"  gSystem->Load(\"libEGPythia6.so\");   // TGenerator interface \n");
    fprintf(pF,"  gSystem->Load(\"libpythia6.so\");   // Pythia \n");
    fprintf(pF,"  gSystem->Load(\"libAliPythia6.so\");   // ALICE specifics implementations \n\n");
//Geant  
    fprintf(pF,"  gSystem->Load(\"libgeant321\");\n");
    fprintf(pF,"  new TGeant3TGeo(\"C++ Interface to Geant3\");\n\n");
  
    fprintf(pF,"  AliDecayer *pDecayer=new AliDecayerPythia();\n");
    fprintf(pF,"  pDecayer->SetForceDecay(kAll);\n"); 
    fprintf(pF,"  pDecayer->Init();\n"); 
    fprintf(pF,"  gMC->SetExternalDecayer(pDecayer);\n\n");
  }
  WritePhys(pF); //physics processes
  
//Field
       if(fMagBG->GetButton(kFld0)->GetState())     fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 0,1,10,AliMagF::k2kG));// NO field\n");
  else if(fMagBG->GetButton(kFld2)->GetState())     fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k2kG));//0.2 Tesla field\n");
  else if(fMagBG->GetButton(kFld4)->GetState())     fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k4kG));//0.4 Tesla field\n");
  else if(fMagBG->GetButton(kFld5)->GetState())     fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2, 1,1,10,AliMagF::k5kG));//0.5 Tesla field\n");
  else if(fMagBG->GetButton(kFld_2)->GetState())    fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k2kG));//-0.2 Tesla field\n");
  else if(fMagBG->GetButton(kFld_4)->GetState())    fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k4kG));//-0.4 Tesla field\n");
  else if(fMagBG->GetButton(kFld_5)->GetState())    fprintf(pF,"  TGeoGlobalMagField::Instance()->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,AliMagF::k5kG));//-0.5 Tesla field\n");
  
  fprintf(pF,"  pAL->CdGAFile();\n\n");                                 //????       
//Generator 
  WriteGen(pF);//generator
//BODY-ALIC 
  fprintf(pF,"  new AliBODY(\"BODY\",\"Alice envelop\");\n\n");
//HMPID
  WriteHmp(pF);  //private HMPID part
  WriteDet(pF);  //other detectors
//end of Config.C file:  
  fprintf(pF,"\n  ::Info(\"----------> HMPID private config\",\"Stop\\n\\n\\n\");\n");
//
  fprintf(pF,"}\n");
  fclose(pF);  
}//WriteConfig()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hconfig()
{
   new HmpConfig("Config.C");
}   
