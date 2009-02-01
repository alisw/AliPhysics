#include <TGFrame.h>
#include <TGToolBar.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGStatusBar.h>
#include <TApplication.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGComboBox.h>
#include <TPDGCode.h>
#include <TMath.h>
#include <Riostream.h>
#include <AliESD.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TPolyMarker.h>
#include "AliHMPIDTracker.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDRecon.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDCluster.h"
#include <TLegend.h>
#include <TLatex.h>
                           
class HmpGui: public TGMainFrame 
{
public:
                  HmpGui();
  virtual        ~HmpGui()                                {Cleanup();}
          void    ToolBar      (                      );
          void    TstTab       (TGCompositeFrame *pTab);
          void    RunTab       (TGCompositeFrame *pTab);          
  virtual void    CloseWindow  (                      )     {gApplication->Terminate(0);}          //executed when user click cross button of main window  
          void    DetSlot1     (Int_t id              )     {if(id==kAllDet) for(Int_t i=1;i<=fDetBG->GetCount();i++) fDetBG->SetButton(i,kButtonDown);}
          void    DetSlot2     (Int_t id              )     {if(id==kAllDet) for(Int_t i=1;i<=fDetBG->GetCount();i++) fDetBG->SetButton(i,kButtonUp);}
          void    ProSlot1     (Int_t id              )     {if(id==kAllPro) for(Int_t i=1;i<=fProBG->GetCount();i++) fProBG->SetButton(i,kButtonDown);}
          void    ProSlot2     (Int_t id              )     {if(id==kAllPro) for(Int_t i=1;i<=fProBG->GetCount();i++) fProBG->SetButton(i,kButtonUp);}
          void    CreateConfigC(); 
          Float_t Eta2Theta    (Float_t arg           )const{return (180./TMath::Pi())*2.*TMath::ATan(TMath::Exp(-arg));}
          
          void    SimEvt       (                      );    
  static  void    SimEsd       (AliESD *pEsd          );
          void    Render       (                      );       
          void    DoTstGo      (                                   ); 
          void    DoZoom       (Int_t evt,Int_t x,Int_t y,TObject *);
          void    MakeStorage  (                                   ); //creates all containers and render elements
          void    DrawChamber(Int_t iCh);
          void    DrawLegend();
protected:
  enum EVersOpts  {kVerNo=101,kVer0,kVer1,kVer2,kTest, kDeclust=301,kSagita,kFeedback,kElNoise,kQe0=400,kQeNorm,kFlatIdx,kOptics};
  enum EGenTypes  {kGunZ=1,kGun1,kGun7,kBox,kHijing,kHijingPara,kPythia,kHmpLib,kNotUsed=999};
  enum EDetectors {kAllDet=1,kPIPE,kITS,kTPC,kTRD,kTOF,kFRAME,kMAG,kACRD,kHALL,kPHOS,kT0,kFMD,kABSO,kPMD,kDIPO,kEMCAL,kVZERO,kMUON,kZDC,kSHILD};
  enum EProcesses {kAllPro=1,kDCAY,kPAIR,kCOMP,kPHOT,kPFIS,kDRAY,kANNI,kBREM,kMUNU,kCKOV,kHADR,kLOSS,kMULS,kRAYL};
  enum EBatchFlags{kNo,kAll,kHmp,kDdl,kDat,kRoo,kVtx,kTrk,kHlt,kPid,kAln};
  enum EMagField  {kFld0,kFld2,kFld4,kFld5,kFld_2,kFld_4,kFld_5};
  
  AliESD          *fEsd; TPolyMarker *fRenTxC[7]; TPolyMarker *fRenRin[7]; 
  TClonesArray *fHitLst; TPolyMarker *fRenMip[7]; TPolyMarker *fRenCko[7]; TPolyMarker *fRenFee[7];
  TClonesArray *fSdiLst;
  TObjArray    *fDigLst; TPolyMarker *fRenDig[7];
  TObjArray    *fCluLst; TPolyMarker *fRenClu[7];
  
  TGStatusBar         *fStatBar;                      // Status bar
  TCanvas             *fCanvas;                       // Hists canvas
  TGButton            *fTstGoB;                      //buttons of Tst Tab
  TGButton            *fSimB,*fRecB,*fTraB;           
  TGButtonGroup       *fVerBG,*fOptBG,*fQeBG ,*fMagBG,*fGenBG,*fDetBG,*fProBG,*fSdiBG,*fDigBG,*fRawBG,*fInpBG,*fCluBG,*fTrkBG; //button groups for Conf Tab
  TGComboBox          *fPidCO,*fNprCO,*fPmiCO,*fPmaCO,*fChaCO;                                                                 //combo bo 
  
  ClassDef(HmpGui,0)
};

ClassImp(HmpGui)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
HmpGui::HmpGui():  TGMainFrame(gClient->GetRoot(), 1000, 700)
{//default ctor  
  ToolBar();
  
  TGTab *pTab; AddFrame(pTab=new TGTab(this,580,360),new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));  
  TstTab(pTab->AddTab("Test"));
  RunTab(pTab->AddTab("Run"));           
  
  AddFrame(fStatBar=new TGStatusBar(this, 50, 10, kHorizontalFrame),new TGLayoutHints(kLHintsBottom| kLHintsExpandX, 0, 0, 0, 0));  
  Int_t parts[] = {45, 45, 10};  fStatBar->SetParts(parts, 3);
  fStatBar->SetText("Waiting for commands...",0);  fStatBar->SetText("No file",1);

  SetCleanup(kDeepCleanup);
  SetWindowName("HMPID Control panel"); 
  MapSubwindows();
  Resize(GetDefaultSize()); // this is used here to launch layout algorithm
  MapWindow();      
  
  MakeStorage();
}//ctor()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::ToolBar()
{//creates toolbar on main frame
  TGToolBar *pTlb;  
  AddFrame(pTlb=new TGToolBar(this, 60, 20, kHorizontalFrame | kRaisedFrame),new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,0,0));
  
  ToolBarData_t b; b.fStayDown=kFALSE;
  b.fPixmap="x_pic.xpm";b.fTipText="Run" ;b.fId=1;b.fButton=NULL;pTlb->AddButton(this,&b, 0);b.fButton->Connect("Pressed()","HmpGui",this,"DoRun()");  
  b.fPixmap="y_pic.xpm";b.fTipText="Info";b.fId++;b.fButton=NULL;pTlb->AddButton(this,&b,40);b.fButton->Connect("Pressed()","HmpGui",this,"DoInfo()");  
  b.fPixmap="z_pic.xpm";b.fTipText="Help";b.fId++;b.fButton=NULL;pTlb->AddButton(this,&b, 0);b.fButton->Connect("Pressed()","HmpGui",this,"DoHelp()");  
}//ToolBar()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::TstTab(TGCompositeFrame *pTab)
{// Create embedded canvas widget
  TRootEmbeddedCanvas *pEmbCan; pTab->AddFrame(pEmbCan=new TRootEmbeddedCanvas("pEmbCnv",pTab,980,660),new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY,2,2,2,2));
  TGVerticalFrame *pV1;         pTab->AddFrame(pV1=new TGVerticalFrame(pTab));    
                                pV1 ->AddFrame(fTstGoB=new TGTextButton(pV1,"Go!"));   fTstGoB ->Connect("Clicked()","HmpGui",this,"DoTstGo()");
  fCanvas = pEmbCan->GetCanvas();
  fCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","hmpGui",this,"ZoomSlot(Int_t,Int_t,Int_t,TObject*)");
  fCanvas->Divide(3,3,0,0);
}//TstTab()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::RunTab(TGCompositeFrame *pTab)
{
  TGCompositeFrame *pH1; pTab->AddFrame(pH1   =new TGHorizontalFrame(pTab,100,200));  
  TGCompositeFrame *pV1; pH1 ->AddFrame(pV1   =new TGVerticalFrame  (pTab,100,200));  //version and options
  TGCompositeFrame *pV2; pH1 ->AddFrame(pV2   =new TGVerticalFrame  (pTab,100,200));  //generator 
                         pH1 ->AddFrame(fDetBG=new TGButtonGroup    (pH1 ,"Dets"));        
                         pH1 ->AddFrame(fProBG=new TGButtonGroup    (pH1 ,"Procs"));       
  TGCompositeFrame *pV3; pH1 ->AddFrame(pV3   =new TGVerticalFrame  (pH1,100,200));   
  TGGroupFrame     *pS;  pV3 ->AddFrame(pS    =new TGGroupFrame     (pV3,"Simu"));   //simu
  TGGroupFrame     *pR;  pV3 ->AddFrame(pR    =new TGGroupFrame     (pV3,"Reco"));   //reco
  
  pV1 ->AddFrame(fVerBG=new TGButtonGroup(pV1 ,"version"));   pV2->AddFrame(fGenBG=new TGButtonGroup(pV2 ,"Type")); pS->AddFrame(fSimB =new TGCheckButton(pS,"Enable"));
  pV1 ->AddFrame(fOptBG=new TGButtonGroup(pV1 ,"Options"));   pV2->AddFrame(fNprCO=new TGComboBox   (pV2,100));     pS->AddFrame(fTraB =new TGCheckButton(pS,"Trans")); 
  pV1 ->AddFrame(fQeBG =new TGButtonGroup(pV1 ,"QE"));        pV2->AddFrame(fPidCO=new TGComboBox   (pV2,100));     pS->AddFrame(fSdiBG=new TGButtonGroup(pS,""));            
  pV1 ->AddFrame(fMagBG=new TGButtonGroup(pV1 ,"Mag field")); pV2->AddFrame(fPmiCO=new TGComboBox   (pV2,100));     pS->AddFrame(fDigBG=new TGButtonGroup(pS,""));
                                                              pV2->AddFrame(fPmaCO=new TGComboBox   (pV2,100));     pS->AddFrame(fRawBG=new TGButtonGroup(pS,""));  
                                                              pV2->AddFrame(fChaCO=new TGComboBox   (pV2,100));     pR->AddFrame(fRecB =new TGCheckButton(pR,"Enable")); 
                                                                                                                    pR->AddFrame(fInpBG=new TGButtonGroup(pR,""   ));       
                                                                                                                    pR->AddFrame(fCluBG=new TGButtonGroup(pR,""   ));       
                                                                                                                    pR->AddFrame(fTrkBG=new TGButtonGroup(pR,""   ));     
 
                          
                         
                         
                         
                         
  fDetBG->Connect("Pressed(Int_t)" ,"HmpGui",this,"DetSlot1(Int_t)");  fDetBG->Connect("Released(Int_t)","HmpGui",this,"DetSlot2(Int_t)");
  fProBG->Connect("Pressed(Int_t)" ,"HmpGui",this,"ProSlot1(Int_t)");  fProBG->Connect("Released(Int_t)","HmpGui",this,"ProSlot2(Int_t)");
  fGenBG->Connect("Pressed(Int_t)" ,"HmpGui",this,"GenSlot1(Int_t)");  fGenBG->Connect("Released(Int_t)","HmpGui",this,"GenSlot2(Int_t)");
  fSimB ->Connect("Toggled(Bool_t)","HmpGui",this,"SimSlot(Bool_t)");
                                 
  new TGRadioButton(fVerBG,"No",kVerNo); new TGRadioButton(fVerBG,"ver0", kVer0);  new TGRadioButton(fVerBG,"ver1",kVer1); new TGRadioButton(fVerBG,"ver2", kVer2);
  
  new TGCheckButton(fOptBG,"Test pos"  ,kTest);     new TGCheckButton(fOptBG,"Unfold clus",kDeclust);  new TGCheckButton(fOptBG,"Sagitta"      ,kSagita);       
  new TGCheckButton(fOptBG,"Phot feed" ,kFeedback); new TGCheckButton(fOptBG,"Elec noise" ,kElNoise);  new TGCheckButton(fOptBG,"N=1.292",kFlatIdx);     
  new TGCheckButton(fOptBG,"Plot optics",kOptics);     
  
  new TGRadioButton(fQeBG,"QE=0",kQe0);         new TGRadioButton(fQeBG,"QE normal",kQeNorm);      
  
  new TGRadioButton(fMagBG,"0.5 T" ,kFld5   );  new TGRadioButton(fMagBG,"0.4 T" ,kFld4   );  new TGRadioButton(fMagBG,"0.2 T" ,kFld2); 
  new TGRadioButton(fMagBG,"0 T"   ,kFld0   );  new TGRadioButton(fMagBG,"-0.2 T",kFld_2  );  new TGRadioButton(fMagBG,"-0.4 T",kFld_4);
  new TGRadioButton(fMagBG,"-0.5 T",kFld_5  );
  
  new TGCheckButton(fGenBG,"Z gun",kGunZ);         new TGCheckButton(fGenBG,"gun 1" ,kGun1);   new TGCheckButton(fGenBG,"gun 7",kGun7);
  new TGCheckButton(fGenBG,"HMPID space",kBox );   new TGCheckButton(fGenBG,"HIJING",kHijing); new TGCheckButton(fGenBG,"HIJING para",kHijingPara);
  new TGCheckButton(fGenBG,"Pythia"     ,kPythia); new TGCheckButton(fGenBG,"HMPID lib",kHmpLib);
  
  fNprCO->AddEntry("not used"    ,kNotUsed);  fNprCO->AddEntry("N prim=1"    ,1);    fNprCO->AddEntry("N prim=2"    ,2);     fNprCO->AddEntry("N prim=5"    ,5);
  fNprCO->AddEntry("N prim=100"  ,100);       fNprCO->AddEntry("N prim=500"  ,500);  fNprCO->AddEntry("N prim=1000" ,1000);  fNprCO->AddEntry("N prim=10000",10000);
  fNprCO->AddEntry("N prim=80000",80000);  
  fNprCO->Resize(160,20);   
  
  fPidCO->AddEntry("not used"   ,kNotUsed);
  fPidCO->AddEntry("e-",kElectron); fPidCO->AddEntry("pi+",kPiPlus);  fPidCO->AddEntry("K+"     ,kKPlus);    fPidCO->AddEntry("p+"     ,kProton);
  fPidCO->AddEntry("e+",kPositron); fPidCO->AddEntry("pi-",kPiMinus); fPidCO->AddEntry("K-"     ,kKMinus);   fPidCO->AddEntry("p-" ,kProtonBar);
  fPidCO->AddEntry("K0",kK0Short);  fPidCO->AddEntry("lambda"     ,kLambda0); fPidCO->AddEntry("antilambda" ,kLambda0Bar); 
  fPidCO->Resize(160,20); fPidCO->Select(kNotUsed);
  
  fPmiCO->AddEntry("not used",kNotUsed); for(Int_t i= 5;i<=295;i+=5) fPmiCO->AddEntry(Form("Pmin=%3.1f GeV",0.1*i), i); fPmiCO->Resize(160,20); 
  fPmaCO->AddEntry("not used",kNotUsed); for(Int_t i=10;i<=295;i+=5) fPmaCO->AddEntry(Form("Pmax=%3.1f GeV",0.1*i), i); fPmaCO->Resize(160,20);   
  fChaCO->AddEntry("not used",kNotUsed); for(Int_t i= 1;i<=  7;i++ ) fChaCO->AddEntry(Form("Chamber %i",i),i);          fChaCO->Resize(160,20); 
  
  new TGCheckButton(fProBG,"ALL  ON/OFF"                 ,kAllPro );
  new TGCheckButton(fProBG,"DCAY Decay"                  ,kDCAY);  fProBG->SetButton(kDCAY);       
  new TGCheckButton(fProBG,"PAIR Pair production"        ,kPAIR);  fProBG->SetButton(kPAIR);       
  new TGCheckButton(fProBG,"COMP Compton"                ,kCOMP);  fProBG->SetButton(kCOMP);
  new TGCheckButton(fProBG,"PHOT Photoelectric"          ,kPHOT);  fProBG->SetButton(kPHOT);
  new TGCheckButton(fProBG,"PFIS Photofission"           ,kPFIS);  
  new TGCheckButton(fProBG,"DRAY Delta electrons"        ,kDRAY);  
  new TGCheckButton(fProBG,"ANNI Annihilation"           ,kANNI);  fProBG->SetButton(kANNI);       
  new TGCheckButton(fProBG,"BREM Bremstraslung"          ,kBREM);  fProBG->SetButton(kBREM);       
  new TGCheckButton(fProBG,"MUNU Muon-Nuclear"           ,kMUNU);  fProBG->SetButton(kMUNU);       
  new TGCheckButton(fProBG,"CKOV Cerenkovs"              ,kCKOV);  fProBG->SetButton(kCKOV);       
  new TGCheckButton(fProBG,"HADR Hadronic interactions " ,kHADR);  fProBG->SetButton(kHADR);       
  new TGCheckButton(fProBG,"LOSS Energy losses"          ,kLOSS);  fProBG->SetButton(kLOSS);       
  new TGCheckButton(fProBG,"MULS Multiple scattering"    ,kMULS);  fProBG->SetButton(kMULS);       
  new TGCheckButton(fProBG,"RAYL Rayleigh scattering"    ,kRAYL);  fProBG->SetButton(kRAYL);       

  new TGCheckButton(fDetBG,"ALL"   ,kAllDet ); new TGCheckButton(fDetBG,"PIPE" ,kPIPE);  new TGCheckButton(fDetBG,"ITS"  ,kITS);   new TGCheckButton(fDetBG,"TPC" ,kTPC);
  new TGCheckButton(fDetBG,"TRD"   ,kTRD ); new TGCheckButton(fDetBG,"TOF"  ,kTOF);   new TGCheckButton(fDetBG,"FRAME",kFRAME); new TGCheckButton(fDetBG,"MAG" ,kMAG);  
  new TGCheckButton(fDetBG,"ACORDE",kACRD); new TGCheckButton(fDetBG,"HALL" ,kHALL);  new TGCheckButton(fDetBG,"PHOS" ,kPHOS);  new TGCheckButton(fDetBG,"T0"  ,kT0); 
  new TGCheckButton(fDetBG,"FMD"   ,kFMD);  new TGCheckButton(fDetBG,"ABSO",kABSO);   new TGCheckButton(fDetBG,"PMD"   ,kPMD);  new TGCheckButton(fDetBG,"DIPO",kDIPO);
  new TGCheckButton(fDetBG,"EMCAL" ,kEMCAL); new TGCheckButton(fDetBG,"VZERO",kVZERO); new TGCheckButton(fDetBG,"MUON" ,kMUON);  new TGCheckButton(fDetBG,"ZDC" ,kZDC);
  new TGCheckButton(fDetBG,"SHILD" ,kSHILD);         
  
  new TGRadioButton(fSdiBG,  "No Sdis"  ,kNo); new TGRadioButton(fSdiBG,"Sdis CORE",kAll); new TGRadioButton(fSdiBG,"Sdis HMPID" ,kHmp);   
  new TGRadioButton(fDigBG,  "No Digs"  ,kNo); new TGRadioButton(fDigBG,"Digs CORE",kAll); new TGRadioButton(fDigBG,"Digs HMPID" ,kHmp);   
  new TGRadioButton(fRawBG,  "No RAW"   ,kNo); new TGRadioButton(fRawBG,"RAW DDL"  ,kDdl); new TGRadioButton(fRawBG,"RAW DATE"   ,kDat);  new TGRadioButton(fRawBG,"RAW ROOT",kRoo);
    
  new TGRadioButton(fInpBG,"From sim"  ,kNo ); new TGRadioButton(fInpBG,"From DDL"   ,kDdl); new TGRadioButton(fInpBG,"From DATE" ,kDat    ); new TGRadioButton(fInpBG,"From ROOT" ,kRoo);  
  new TGRadioButton(fCluBG,"No Clus"   ,kNo ); new TGRadioButton(fCluBG,"Clus CORE"  ,kAll); new TGRadioButton(fCluBG,"Clus HMPID",kHmp    );   
  new TGCheckButton(fTrkBG,"Load Align",kAln); new TGCheckButton(fTrkBG,"Prim vertex",kVtx); new TGCheckButton(fTrkBG,"ESD tracks",kTrk    ); new TGCheckButton(fTrkBG,"Assign PID",kPid);
    
  fSimB->SetState(kButtonDown); fTraB->SetState(kButtonDown); fSdiBG->SetButton(kHmp); fDigBG->SetButton(kHmp); fRawBG->SetButton(kNo);
  fRecB->SetState(kButtonDown);                               fInpBG->SetButton(kNo);  fCluBG->SetButton(kHmp); 
  
  fVerBG->SetButton(kVer1); fOptBG->SetButton(kDeclust); fOptBG->SetButton(kSagita); fOptBG->SetButton(kFeedback); fQeBG->SetButton(kQeNorm); fMagBG->SetButton(kFld2);
  fNprCO->Select(kNotUsed); fPmiCO->Select(kNotUsed);  fPmaCO->Select(kNotUsed);  fChaCO->Select(kNotUsed);
}//RunTab()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::CreateConfigC()
{
  ofstream f; f.open("Config.C");
  
  f<<Form("void Config()\n"); f<<Form("{\n");  f<<Form("\n  ::Info(\"\\n\\n\\n----------> HMPID private config\",\"Start\");\n"); //start & info
  f<<Form("  gRandom->SetSeed(123456);//put 0 to use system time\n\n");                                                           //random seed  
  f<<Form("  gSystem->Load(\"libgeant321\");\n");  f<<Form("  new TGeant3TGeo(\"C++ Interface to Geant3\");\n\n");                //geant lib 
  f<<Form("  AliRunLoader *pAL=AliRunLoader::Open(\"galice.root\",AliConfig::GetDefaultEventFolderName(),\"recreate\");\n");      //run loader
  f<<Form("  pAL->SetCompressionLevel(2);\n");  f<<Form("  pAL->SetNumberOfEventsPerFile(1000);\n");  f<<Form("  gAlice->SetRunLoader(pAL);\n\n");
  f<<Form("  TVirtualMCDecayer *pDecayer=new AliDecayerPythia();\n");
  f<<Form("  pDecayer->SetForceDecay(kAll);\n"); 
  f<<Form("  pDecayer->Init();\n"); 
  f<<Form("  gMC->SetExternalDecayer(pDecayer);\n\n");
  
//Field
       if(fMagBG->GetButton(kFld0)->GetState())     f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",0,1,1,10,0));       //no field\n\n");
  else if(fMagBG->GetButton(kFld2)->GetState())     f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,1,1,10,0));//0.2 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld4)->GetState())     f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,1,1,10,1));//0.4 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld5)->GetState())     f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,1,1,10,2));//0.5 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_2)->GetState())    f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,0));//-0.2 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_4)->GetState())    f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,1));//-0.4 Tesla field\n\n");
  else if(fMagBG->GetButton(kFld_5)->GetState())    f<<Form("  gAlice->SetField(new AliMagF(\"Maps\",\"Maps\",2,-1,1,10,2));//-0.5 Tesla field\n\n");
  
  if(fProBG->GetButton(kDCAY)->GetState()) f<<Form("  gMC->SetProcess(\"DCAY\",1);  ");else f<<Form("  gMC->SetProcess(\"DCAY\",0);  ");
  if(fProBG->GetButton(kPAIR)->GetState()) f<<Form("  gMC->SetProcess(\"PAIR\",1);  ");else f<<Form("  gMC->SetProcess(\"PAIR\",0);  ");
  if(fProBG->GetButton(kCOMP)->GetState()) f<<Form("  gMC->SetProcess(\"COMP\",1);\n");else f<<Form("  gMC->SetProcess(\"COMP\",0);\n");
  if(fProBG->GetButton(kPHOT)->GetState()) f<<Form("  gMC->SetProcess(\"PHOT\",1);  ");else f<<Form("  gMC->SetProcess(\"PHOT\",0);  ");
  if(fProBG->GetButton(kPFIS)->GetState()) f<<Form("  gMC->SetProcess(\"PFIS\",1);  ");else f<<Form("  gMC->SetProcess(\"PFIS\",0);  ");
  if(fProBG->GetButton(kDRAY)->GetState()) f<<Form("  gMC->SetProcess(\"DRAY\",1);\n");else f<<Form("  gMC->SetProcess(\"DRAY\",0);\n");
  if(fProBG->GetButton(kANNI)->GetState()) f<<Form("  gMC->SetProcess(\"ANNI\",1);  ");else f<<Form("  gMC->SetProcess(\"ANNI\",0);  ");
  if(fProBG->GetButton(kBREM)->GetState()) f<<Form("  gMC->SetProcess(\"BREM\",1);  ");else f<<Form("  gMC->SetProcess(\"BREM\",0);  ");
  if(fProBG->GetButton(kMUNU)->GetState()) f<<Form("  gMC->SetProcess(\"MUNU\",1);\n");else f<<Form("  gMC->SetProcess(\"MUNU\",0);\n");
  if(fProBG->GetButton(kCKOV)->GetState()) f<<Form("  gMC->SetProcess(\"CKOV\",1);  ");else f<<Form("  gMC->SetProcess(\"CKOV\",0);  ");
  if(fProBG->GetButton(kHADR)->GetState()) f<<Form("  gMC->SetProcess(\"HADR\",1);  ");else f<<Form("  gMC->SetProcess(\"HADR\",0);  ");
  if(fProBG->GetButton(kLOSS)->GetState()) f<<Form("  gMC->SetProcess(\"LOSS\",2);\n");else f<<Form("  gMC->SetProcess(\"LOSS\",0);\n");
  if(fProBG->GetButton(kMULS)->GetState()) f<<Form("  gMC->SetProcess(\"MULS\",1);  ");else f<<Form("  gMC->SetProcess(\"MULS\",0);  ");
  if(fProBG->GetButton(kRAYL)->GetState()) f<<Form("  gMC->SetProcess(\"RAYL\",1);\n");else f<<Form("  gMC->SetProcess(\"RAYL\",0);\n");
  
  f<<Form("\n");  
  f<<Form("  gMC->SetCut(\"CUTGAM\",0.001);  "); f<<Form("  gMC->SetCut(\"CUTELE\",0.001);  "); f<<Form("  gMC->SetCut(\"CUTNEU\",0.001);\n"); 
  f<<Form("  gMC->SetCut(\"CUTHAD\",0.001);  "); f<<Form("  gMC->SetCut(\"CUTMUO\",0.001);  "); f<<Form("  gMC->SetCut(\"BCUTE\" ,0.001);\n");
  f<<Form("  gMC->SetCut(\"BCUTM\" ,0.001);  "); f<<Form("  gMC->SetCut(\"DCUTE\" ,0.001);  "); f<<Form("  gMC->SetCut(\"DCUTM\" ,0.001);\n"); 
  f<<Form("  gMC->SetCut(\"PPCUTM\",0.001);  "); f<<Form("  gMC->SetCut(\"TOFMAX\",1e10);\n\n");
  f<<Form("  pAL->CdGAFile();\n\n");                                 //????       
  
  f<<Form("  AliGenCocktail *pG=new AliGenCocktail();\n\n"); //Generator                                 
  Int_t pid=fPidCO->GetSelected();  Float_t pmin=0.1*fPmiCO->GetSelected();   Float_t pmax=0.1*fPmaCO->GetSelected();
  
  if(fGenBG->GetButton(kGunZ)->GetState()==kButtonDown)//1 particle along Z axis 
    f<<Form("  AliGenFixed *pGz=new AliGenFixed(1); pGz->SetPart(%i); pGz->SetMomentum(%.1f); pGz->SetOrigin(0,0,-200); pG->AddGenerator(pGz,\"Gz\",1);\n",pid,pmin);
  
  if(fGenBG->GetButton(kGun1)->GetState()==kButtonDown){//1 gun towards 1 HMPID chamber
    switch(fChaCO->GetSelected()){
     case 1: f<<Form("  AliGenFixed *pG1=new AliGenFixed(1); pG1->SetPart(%i); pG1->SetMomentum(%.1f);\n",pid,pmin); 
             f<<Form("               pG1->SetTheta(109.5);   pG1->SetPhi(10);  pG->AddGenerator(pG1,\"g1\",1);\n"); break;
     case 2: f<<Form("  AliGenFixed *pG2=new AliGenFixed(1); pG2->SetPart(%i); pG2->SetMomentum(%.1f);\n",pid,pmin);
	     f<<Form("               pG2->SetTheta( 90.0);   pG2->SetPhi(10);  pG->AddGenerator(pG2,\"g2\",1);\n"); break;
     case 3: f<<Form("  AliGenFixed *pG3=new AliGenFixed(1); pG3->SetPart(%i); pG3->SetMomentum(%.1f);\n",pid,pmin);
	     f<<Form("               pG3->SetTheta(109.5);   pG3->SetPhi(30);  pG->AddGenerator(pG3,\"g3\",1);\n"); break;
     case 4: f<<Form("  AliGenFixed *pG4=new AliGenFixed(1); pG4->SetPart(%i); pG4->SetMomentum(%.1f);\n",pid,pmin);
             f<<Form("               pG4->SetTheta( 87.0);   pG4->SetPhi(30);  pG->AddGenerator(pG4,\"g4\",1);\n"); break;
     case 5: f<<Form("  AliGenFixed *pG5=new AliGenFixed(1); pG5->SetPart(%i); pG5->SetMomentum(%.1f);\n",pid,pmin);
             f<<Form("               pG5->SetTheta( 70.5);   pG5->SetPhi(30);  pG->AddGenerator(pG5,\"g5\",1);\n"); break;
     case 6: f<<Form("  AliGenFixed *pG6=new AliGenFixed(1); pG6->SetPart(%i); pG6->SetMomentum(%.1f);\n",pid,pmin);
             f<<Form("               pG6->SetTheta( 90.0);   pG6->SetPhi(50);  pG->AddGenerator(pG6,\"g6\",1);\n"); break;
     case 7: f<<Form("  AliGenFixed *pG7=new AliGenFixed(1); pG7->SetPart(%i); pG7->SetMomentum(%.1f);\n",pid,pmin);
             f<<Form("               pG7->SetTheta( 70.5);   pG7->SetPhi(50);  pG->AddGenerator(pG7,\"g7\",1);\n"); break;
    }    
  }
  
  if(fGenBG->GetButton(kGun7)->GetState()==kButtonDown){//7 guns towards 7 HMPID chambers
    f<<Form("  AliGenFixed *pG1=new AliGenFixed(1); pG1->SetPart(%i); pG1->SetMomentum(%.1f);pG1->SetTheta(109.5-3); pG1->SetPhi(10);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG1,\"g1\",1);\n");
    f<<Form("  AliGenFixed *pG2=new AliGenFixed(1); pG2->SetPart(%i); pG2->SetMomentum(%.1f);pG2->SetTheta( 90.0-3); pG2->SetPhi(10);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG2,\"g2\",1);\n");
    f<<Form("  AliGenFixed *pG3=new AliGenFixed(1); pG3->SetPart(%i); pG3->SetMomentum(%.1f);pG3->SetTheta(109.5-3); pG3->SetPhi(30);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG3,\"g3\",1);\n");
    f<<Form("  AliGenFixed *pG4=new AliGenFixed(1); pG4->SetPart(%i); pG4->SetMomentum(%.1f);pG4->SetTheta( 90.0-3); pG4->SetPhi(30);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG4,\"g4\",1);\n");
    f<<Form("  AliGenFixed *pG5=new AliGenFixed(1); pG5->SetPart(%i); pG5->SetMomentum(%.1f);pG5->SetTheta( 70.0-3); pG5->SetPhi(30);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG5,\"g5\",1);\n");
    f<<Form("  AliGenFixed *pG6=new AliGenFixed(1); pG6->SetPart(%i); pG6->SetMomentum(%.1f);pG6->SetTheta( 90.0-3); pG6->SetPhi(50);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG6,\"g6\",1);\n");
    f<<Form("  AliGenFixed *pG7=new AliGenFixed(1); pG7->SetPart(%i); pG7->SetMomentum(%.1f);pG7->SetTheta( 70.0-3); pG7->SetPhi(50);\n",pid,pmin); 
    f<<Form("               pG->AddGenerator(pG7,\"g7\",1);\n");
  }  
    
  if(fGenBG->GetButton(kBox)->GetState()==kButtonDown){// box towards HMPID phase space
    f<<Form("  AliGenBox *pB=new AliGenBox(%i);      pB->SetPart(%i);       pB->SetMomentumRange(%.1f,%.1f);\n",(int)fNprCO->GetSelected(),pid,pmin,pmax); 
    f<<Form("             pB->SetThetaRange(65,115); pB->SetPhiRange(5,55); pG->AddGenerator(pB,\"b\",1);\n");
  }     
  
  if(fGenBG->GetButton(kHijing)->GetState()==kButtonDown){//normal HIJING
    f<<Form("  AliGenHijing *pH=new AliGenHijing(-1);           pH->SetEnergyCMS(14000);        pH->SetReferenceFrame(\"CMS\");\n");
    f<<Form("                pH->SetProjectile(\"P\", 1, 1); pH->SetTarget(\"P\", 1, 1  ); pH->SetJetQuenching(0);\n");      
    f<<Form("                pH->SetShadowing(0);               pH->KeepFullEvent();           pH->SetSelectAll(0);\n");
    f<<Form("                pH->SetImpactParameterRange(0, 5); //fermi\n");
    f<<Form("  pG->AddGenerator(pH,\"h\",1);\n\n");
  }
  
  if(fGenBG->GetButton(kHijingPara)->GetState()==kButtonDown){//parametrized HIJING 
    f<<Form("  AliGenHIJINGpara *pHP=new AliGenHIJINGpara(%i);\n",(int)fNprCO->GetSelected());
    f<<Form("  pHP->SetMomentumRange(0,999); pHP->SetThetaRange(%f,%f); pHP->SetPhiRange(0,360);\n",Eta2Theta(8),Eta2Theta(-8));
    f<<Form("  pG->AddGenerator(pHP,\"hp\",1);\n\n");
  }
      
  if(fGenBG->GetButton(kPythia)->GetState()==kButtonDown){//Pythia
    f<<Form("  AliGenPythia *pP=new AliGenPythia(-1);\n");
    f<<Form("  pP->SetMomentumRange(0,999); pP->SetPhiRange(20,80); pP->SetThetaRange(75,115);\n");
    f<<Form("  pP->SetYRange(-12,12);  pP->SetPtRange(0,1000);      pP->SetStrucFunc(kCTEQ4L);\n");
    f<<Form("  pP->SetProcess(kPyMb);  pP->SetEnergyCMS(14000);\n");      
    f<<Form("  pG->AddGenerator(pP,\"p\",1);\n\n");  
  }
  f<<Form("  pG->Init();\n\n");


                                            f<<Form("\n  new AliBODY           (\"BODY\"  ,\"Alice envelop\");");             //BODY-ALIC 
  if(fDetBG->GetButton(kPIPE )->GetState()) f<<Form("\n  new AliPIPEv3         (\"PIPE\"  ,\"Beam Pipe\");");
  if(fDetBG->GetButton(kSHILD)->GetState()) f<<Form("\n  new AliSHILv3         (\"SHIL\"  ,\"Shielding Version 2\");");  
  if(fDetBG->GetButton(kITS  )->GetState()) f<<Form("\n  new AliITSvPPRasymmFMD(\"ITS\"   ,\"ITS PPR\");");
  if(fDetBG->GetButton(kTPC  )->GetState()) f<<Form("\n  new AliTPCv2          (\"TPC\"   ,\"Default\");");
  if(fDetBG->GetButton(kFRAME)->GetState()) f<<Form("\n  new AliFRAMEv2        (\"FRAME\" ,\"Space Frame\");");
  if(fDetBG->GetButton(kTRD  )->GetState()) f<<Form("\n  new AliTRDv1          (\"TRD\"   ,\"TRD slow simulator\");");  
  if(fDetBG->GetButton(kTOF  )->GetState()) f<<Form("\n  new AliTOFv6T0        (\"TOF\"   , \"normal TOF\");");
  if(fDetBG->GetButton(kMAG  )->GetState()) f<<Form("\n  new AliMAG            (\"MAG\"   ,\"Magnet\");");                                
  if(fDetBG->GetButton(kHALL )->GetState()) f<<Form("\n  new AliHALL           (\"HALL\"  ,\"Alice Hall\");");
  if(fDetBG->GetButton(kFMD  )->GetState()) f<<Form("\n  new AliFMDv1          (\"FMD\"   ,\"normal FMD\");");            
  if(fDetBG->GetButton(kABSO )->GetState()) f<<Form("\n  new AliABSOv3         (\"ABSO\"  ,\"Muon absorber\");");
  if(fDetBG->GetButton(kDIPO )->GetState()) f<<Form("\n  new AliDIPOv3         (\"DIPO\"  ,\"Dipole version 3\");");
  if(fDetBG->GetButton(kMUON )->GetState()) f<<Form("\n  new AliMUONv1         (\"MUON\"  ,\"default\");");
  if(fDetBG->GetButton(kPMD  )->GetState()) f<<Form("\n  new AliPMDv1          (\"PMD\"   ,\"normal PMD\");");
  if(fDetBG->GetButton(kT0   )->GetState()) f<<Form("\n  new AliT0v1           (\"T0\"    ,\"T0 Detector\");");
  if(fDetBG->GetButton(kVZERO)->GetState()) f<<Form("\n  new AliVZEROv7        (\"VZERO\" ,\"normal VZERO\");");
  if(fDetBG->GetButton(kZDC  )->GetState()) f<<Form("\n  new AliZDCv2          (\"ZDC\"   ,\"normal ZDC\");");
  if(fDetBG->GetButton(kPHOS )->GetState()) f<<Form("\n  new AliPHOSv1         (\"PHOS\"  ,\"IHEP\");");              
  if(fDetBG->GetButton(kEMCAL)->GetState()) f<<Form("\n  new AliEMCALv2        (\"EMCAL\" ,\"SHISH_77_TRD1_2X2_FINAL_110DEG\");");
  if(fDetBG->GetButton(kACRD )->GetState()) f<<Form("\n  new AliACORDEv1       (\"ACORDE\",\"normal ACORDE\");");
  
  f<<Form("\n  ::Info(\"----------> HMPID private config\",\"Stop\\n\\n\\n\");\n");   f<<Form("}\n");  //end of Config.C file:  
  f.close(); 
}//CreateConfigC()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::SimEvt()
{//simulates event
  fEsd->Clear();  fHitLst->Clear();  fSdiLst->Clear();  for(Int_t ch=0;ch<=6;ch++) {fDigLst->At(ch)->Clear(); fCluLst->At(ch)->Clear();}
  SimEsd(fEsd);
//  SimHits(fEsd,fHitLst);
//                 AliHMPIDv1::Hit2Sdi(fHitLst,fSdiLst);                               
//          AliHMPIDDigitizer::Sdi2Dig(fSdiLst,fDigLst);     
//      AliHMPIDReconstructor::Dig2Clu(fDigLst,fCluLst);
//            AliHMPIDTracker::Recon(fEsd,fCluLst,(TObjArray*)pNmeanEnt->GetObject());
//  AliCDBManager* pCDB = AliCDBManager::Instance();  pCDB->SetDefaultStorage("local://$HOME"); pCDB->SetRun(0);
//  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean");
            
}//SimEvt()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::SimEsd(AliESD *pEsd)
{
  TParticle part; TLorentzVector mom;
  for(Int_t iTrk=0;iTrk<100;iTrk++){//tracks loop
    part.SetPdgCode(kProton);    part.SetProductionVertex(0,0,0,0);
    Double_t eta= -0.2+gRandom->Rndm()*0.4;                //rapidity is random [-0.2,+0.2]
    Double_t phi= gRandom->Rndm()*60.*TMath::DegToRad();   //phi is random      [ 0  , 60 ] degrees    
    mom.SetPtEtaPhiM(5,eta,phi,part.GetMass());   part.SetMomentum(mom);
    pEsd->AddTrack(new AliESDtrack(&part));
  }//tracks loop  
}//SimEsd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimHits(AliESD *pEsd, TClonesArray *pHits)
{//used by SimulateEvent to simulate hits out from provided ESD
  const Int_t kCerenkov=50000050;  const Int_t kFeedback=50000051;
  
  AliHMPIDRecon rec;
  Float_t eMip=200e-9,ePho=7.5e-9; 
  Int_t hc=0; 
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//tracks loop
    AliESDtrack *pTrk=pEsd->GetTrack(iTrk);
    Float_t xRa,yRa;
    Int_t ch=AliHMPIDTracker::IntTrkCha(pTrk,xRa,yRa);
    if(ch<0) continue; //this track does not hit HMPID
    Float_t beta = pTrk->GetP()/(TMath::Sqrt(pTrk->GetP()*pTrk->GetP()+0.938*0.938));
    Float_t ckov=TMath::ACos(1./(beta*1.292));

    Float_t theta,phi,xPc,yPc,; pTrk->GetHMPIDtrk(xPc,yPc,theta,phi); rec.SetTrack(xRa,yRa,theta,phi); 
    
    if(!AliHMPIDParam::IsInDead(xPc,yPc)) new((*pHits)[hc++]) AliHMPIDHit(ch,eMip,kProton  ,iTrk,xPc,yPc);                 //mip hit
    Int_t nPhots = (Int_t)(20.*TMath::Power(TMath::Sin(ckov),2)/TMath::Power(TMath::Sin(TMath::ACos(1./1.292)),2));
    for(int i=0;i<nPhots;i++){
      TVector2 pos;
      pos=rec.TracePhot(ckov,gRandom->Rndm()*TMath::TwoPi());
      if(!AliHMPIDParam::IsInDead(pos.X(),pos.Y())) new((*pHits)[hc++]) AliHMPIDHit(ch,ePho,kCerenkov,iTrk,pos.X(),pos.Y());
    }                      //photon hits  
    for(int i=0;i<3;i++){//feedback photons
      Float_t x=gRandom->Rndm()*160; Float_t y=gRandom->Rndm()*150;
      if(!AliHMPIDParam::IsInDead(x,y)) new((*pHits)[hc++]) AliHMPIDHit(ch,ePho,kFeedback,iTrk,x,y);                 //feedback hits  
    }//photon hits loop                      
  }//tracks loop    
}//SimHits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::DoTstGo()
{

  SimEvt();
  Render();
  
  for(Int_t iCh=0;iCh<=6;iCh++){//chambers loop    
    switch(iCh){
      case 6: fCanvas->cd(1); break; case 5: fCanvas->cd(2); break;
      case 4: fCanvas->cd(4); break; case 3: fCanvas->cd(5); break; case 2: fCanvas->cd(6); break;
                                     case 1: fCanvas->cd(8); break; case 0: fCanvas->cd(9); break;
    }
    gPad->SetEditable(kTRUE); gPad->Clear(); DrawChamber(iCh);
    fRenTxC[iCh]->Draw();        
    fRenMip[iCh]->Draw();        
    fRenFee[iCh]->Draw();        
    fRenCko[iCh]->Draw();       
    fRenRin[iCh]->Draw("CLP");   
    fRenDig[iCh]->Draw();        
    fRenClu[iCh]->Draw();        
    gPad->SetEditable(kFALSE);
  }//chambers loop  
  fCanvas->Modified();fCanvas->Update();
  
}//TstGoSlot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::Render()
{//render all containers to polymarker structures, one per chamber
  for(Int_t iCh=0;iCh<=6;iCh++){                                                                                     //chambers loop   
    fRenTxC[iCh]->SetPolyMarker(0);     fRenRin[iCh]->SetPolyMarker(0);                                                //clear all renders alements for this chamber 
    fRenMip[iCh]->SetPolyMarker(0);     fRenCko[iCh]->SetPolyMarker(0);     fRenFee[iCh]->SetPolyMarker(0);  
    fRenDig[iCh]->SetPolyMarker(0);     fRenClu[iCh]->SetPolyMarker(0); 
    TClonesArray *pDigCh=(TClonesArray*)fDigLst->At(iCh); TClonesArray *pCluCh=(TClonesArray*)fCluLst->At(iCh);      
    for(Int_t iDig=0;iDig<pDigCh->GetEntries();iDig++){AliHMPIDDigit   *pDig = (AliHMPIDDigit*)  pDigCh->At(iDig); fRenDig[iCh]->SetNextPoint(pDig->LorsX(),pDig->LorsY());}
    for(Int_t iClu=0;iClu<pCluCh->GetEntries();iClu++){AliHMPIDCluster *pClu = (AliHMPIDCluster*)pCluCh->At(iClu); fRenClu[iCh]->SetNextPoint(pClu->X()    ,pClu->Y());    }
  }//chambers loop
  for(Int_t iHit=0;iHit<fHitLst->GetEntries();iHit++){       //hits loop
    AliHMPIDHit *pHit = (AliHMPIDHit*)fHitLst->At(iHit); Int_t ch=pHit->Ch(); Float_t x=pHit->LorsX(); Float_t y=pHit->LorsY();    //get current hit        
    switch(pHit->Pid()){
      case 50000050: fRenCko[ch]->SetNextPoint(x,y);break; 
      case 50000051: fRenFee[ch]->SetNextPoint(x,y);break;
      default:       fRenMip[ch]->SetNextPoint(x,y);break;
    }//switch hit PID      
  }//hits loop
  AliHMPIDRecon rec;
  for(Int_t iTrk=0;iTrk<fEsd->GetNumberOfTracks();iTrk++){//tracks loop to collect cerenkov rings and intersection points
    AliESDtrack *pTrk=fEsd->GetTrack(iTrk);    Int_t ch=pTrk->GetHMPIDcluIdx(); //get track and chamber intersected by it
    if(ch<0) continue;                                                          //this track does not intersect any chamber
    Float_t thRa,phRa,xRa,yRa; pTrk->GetHMPIDtrk(xRa,yRa,thRa,phRa);            //get info on current track
    ch/=1000000;                            
    Float_t xPc=0,yPc=0; AliHMPIDTracker::IntTrkCha(pTrk,xPc,yPc);              //find again intersection of track with PC--> it is not stored in ESD!
    fRenTxC[ch]->SetNextPoint(xPc,yPc);                                         //add this intersection point
    Float_t ckov=pTrk->GetHMPIDsignal();                                        //get ckov angle stored for this track  
    if(ckov>0){
      rec.SetTrack(xRa,yRa,thRa,phRa);
      for(Int_t j=0;j<100;j++){ 
        TVector2 pos; pos=rec.TracePhot(ckov,j*0.0628);
       if(!AliHMPIDParam::IsInDead(pos.X(),pos.Y())) fRenRin[ch]->SetNextPoint(pos.X(),pos.Y());
      }      
    }//if ckov is valid
  }//tracks loop  
}//Render()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::DoZoom(Int_t evt, Int_t px, Int_t py, TObject *)
{
  if(evt!=5 && evt!=6) return; //5- zoom in 6-zoom out
  const Int_t minZoom=64;
  const Int_t maxZoom=2;
  static Int_t zoom=minZoom; //zoom level
  if(evt==5&&zoom==maxZoom) return; 
  if(evt==6&&zoom==minZoom) return; 
  
 // if(!obj->IsA()->InheritsFrom("TPad")) return;  //current object is not pad
  TVirtualPad *pPad=gPad->GetSelectedPad();
  if(pPad->GetNumber()==3 || pPad->GetNumber()==7) return; //current pad is wrong

 // Printf("evt=%i (%i,%i) %s",evt,px,py,obj->GetName());
    
  Float_t x=pPad->AbsPixeltoX(px); Float_t y=pPad->AbsPixeltoY(py); 
 
  if(evt==5){ zoom=zoom/2;     pPad->Range(x-zoom*2,y-zoom*2,x+zoom*2,y+zoom*2);} //zoom in
  else      { zoom=zoom*2;     pPad->Range(x-zoom*2,y-zoom*2,x+zoom*2,y+zoom*2);} //zoom out 
  if(zoom==minZoom) pPad->Range(-10,-10,AliHMPIDParam::SizeAllX()+5,AliHMPIDParam::SizeAllY()+5);
  ((TCanvas *)gTQSender)->SetTitle(Form("zoom x%i",minZoom/zoom));
  pPad->Modified();
  pPad->Update();                                              
}//DoZoom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::MakeStorage()
{//creates all render elements
  fHitLst=new TClonesArray("AliHMPIDHit");
  fSdiLst=new TClonesArray("AliHMPIDDigit");
  fDigLst=new TObjArray(7); 
  fCluLst=new TObjArray(7); 
  fEsd   =new AliESD;
  for(Int_t ch=0;ch<7;ch++){
    fDigLst->AddAt(new TClonesArray("AliHMPIDDigit"),ch);       fDigLst->SetOwner(kTRUE);
    fCluLst->AddAt(new TClonesArray("AliHMPIDCluster"),ch);     fCluLst->SetOwner(kTRUE); 
    fRenMip[ch]=new TPolyMarker; fRenMip[ch]->SetMarkerStyle(kOpenTriangleUp);  fRenMip[ch]->SetMarkerColor(kRed);
    fRenCko[ch]=new TPolyMarker; fRenCko[ch]->SetMarkerStyle(kOpenCircle);      fRenCko[ch]->SetMarkerColor(kRed);
    fRenFee[ch]=new TPolyMarker; fRenFee[ch]->SetMarkerStyle(kOpenDiamond);     fRenFee[ch]->SetMarkerColor(kRed);
    fRenDig[ch]=new TPolyMarker; fRenDig[ch]->SetMarkerStyle(kOpenSquare);      fRenDig[ch]->SetMarkerColor(kGreen);
    fRenClu[ch]=new TPolyMarker; fRenClu[ch]->SetMarkerStyle(kStar);            fRenClu[ch]->SetMarkerColor(kBlue);
    fRenTxC[ch]=new TPolyMarker; fRenTxC[ch]->SetMarkerStyle(kPlus);            fRenTxC[ch]->SetMarkerColor(kRed);      fRenTxC[ch]->SetMarkerSize(3);
    fRenRin[ch]=new TPolyMarker; fRenRin[ch]->SetMarkerStyle(kFullDotSmall);    fRenRin[ch]->SetMarkerColor(kMagenta);
  }
}//MakeStorage()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::DrawChamber(Int_t iCh) 
{//used by Draw() to draw chamber structure
  gPad->Range(-10,-10,AliHMPIDParam::SizeAllX()+5,AliHMPIDParam::SizeAllY()+5); 
  if(iCh>=0){TLatex txt; txt.SetTextSize(0.1); txt.DrawLatex(-5,-5,Form("%i",iCh));}
  
  for(Int_t iPc=AliHMPIDParam::kMinPc;iPc<=AliHMPIDParam::kMaxPc;iPc++){
    TBox *pBox=new TBox(AliHMPIDParam::MinPcX(iPc),AliHMPIDParam::MinPcY(iPc),
                        AliHMPIDParam::MaxPcX(iPc),AliHMPIDParam::MaxPcY(iPc));
    pBox->SetFillStyle(0);  pBox->Draw();
  }//PC loop      
}//DrawChamber()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpGui::DrawLegend()
{//used by Draw() to draw legend
  Int_t nTxC=0,nMip=0,nCko=0,nFee=0,nDig=0,nClu=0;
  for(Int_t ch=0;ch<7;ch++){
    nTxC+=fRenTxC[ch]->GetN();
    nMip+=fRenMip[ch]->GetN();
    nCko+=fRenCko[ch]->GetN();
    nFee+=fRenFee[ch]->GetN();
    nDig+=fRenDig[ch]->GetN();
    nClu+=fRenClu[ch]->GetN();
  }
  TLegend *pLeg=new TLegend(0.2,0.2,0.8,0.8);
//  pLeg->SetHeader(Form("Event %i Total %i",fEvt,fNevt));
  pLeg->AddEntry(fRenTxC[0],Form("TRKxPC %i"     ,nTxC),"p");
  pLeg->AddEntry(fRenMip[0],Form("Mip hits %i"   ,nMip),"p");    
  pLeg->AddEntry(fRenCko[0],Form("Ckov hits %i"  ,nCko),"p");    
  pLeg->AddEntry(fRenFee[0],Form("Feed hits %i"  ,nFee),"p");    
  pLeg->AddEntry(fRenDig[0],Form("Digs %i"       ,nDig),"p");    
  pLeg->AddEntry(fRenClu[0],Form("Clus %i"       ,nClu),"p");    
  pLeg->Draw();
}//DrawLegend()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void Hnew(){   new HmpGui;    }   
