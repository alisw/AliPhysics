const static Int_t   iRICH  =  3;//0-1-3
const static Bool_t  IsRichUp=kFALSE;
const static Int_t   kEventsPerFile=50;

enum  EGenTypes {kGun0,kGun1,kGun7,kPP7};
const static EGenTypes kGen=kGun1;

Int_t   iPIPE  =  0;//central before RICH
Int_t   iITS   =  0;
Int_t   iTPC   =  0;
Int_t   iTRD   =  0;
Int_t   iTOF   =  0;
Int_t   iFRAME =  0;

Int_t   iMAG   =  0;//central after RICH
Int_t   iCRT   =  0;
Int_t   iHALL  =  0;

Int_t   iPHOS  =  0;

Int_t   iSTART =  0;//forward
Int_t   iFMD   =  0;
Int_t   iSHIL  =  0;
Int_t   iABSO  =  1;
Int_t   iPMD   =  0;
Int_t   iDIPO  =  0;
Int_t   iEMCAL =  0;
Int_t   iVZERO =  0;
Int_t   iMUON  =  0;
Int_t   iZDC   =  0;

AliRICH *gRICH;

void Config()
{   
  ::Info("kir","Start RICH VERSION %i debug %i",iRICH,gAlice->GetDebug());
  
  int n,pid,chamber;
  double p;
  Geant3(); File(); Decayer(); Field(); Other();
  

  if(IsRichUp) AliRICHParam::AngleRot(0); 
  AliRICH *pRICH;
  switch(iRICH){
    case 0:
      gRICH=new AliRICHv0("RICH","coarse RICH for material budget");
      break;
    case 1:
      gRICH=new AliRICHv1("RICH","detailed RICH for full simulation");
      break;
    case 3:
      gRICH=new AliRICHv3("RICH","old parametrised RICH with rotation");
      break;
  }   
  switch(kGen){
	case kPP7:          Pythia7(pid=kPiPlus,p=4);                break;
        case kGun7:         Gun7(pid=kPiPlus,p=4);                   break;
        case kGun1:         Gun1(n=1,pid=kPiPlus,p=4,chamber=4);     break;
        case kGun0:         Gun1(n=1,pid=kNeutron,p=4,chamber=4);    break;
        default:            Fatal("Config","No generator");          break;
  }  
  ::Info("kir","Stop.");
}//void Config()

void Gun1(Int_t iNprim,Int_t iPID,Double_t p,Int_t iChamber)
{
  Double_t theta=gRICH->C(iChamber)->ThetaD();
  Double_t phi  =gRICH->C(iChamber)->PhiD();
  theta-=2;
  ::Info("kir-Gun1","%i primarie(s) of %i PID  with p=%f GeV at (%f,%f)",iNprim,iPID,p,theta,phi);
   
  AliGenFixed *pGen=new AliGenFixed(iNprim);
  pGen->SetMomentum(p);
  pGen->SetPhiRange(phi);
  pGen->SetThetaRange(theta);
  pGen->SetOrigin(0,0,0);                 
  pGen->SetPart(iPID);    
  pGen->Init();
}//Gun()     

void Para()
{
  Int_t iNprim=85700;
  ::Info("kir-Para","%i primaries",iNprim);
   
  AliGenHIJINGpara *pGen=new AliGenHIJINGpara(iNprim);
  pGen->SetMomentumRange(0,999);                      //GeV
  pGen->SetPhiRange(0,360);                           //degree
  pGen->SetThetaRange(Eta2Theta(8),Eta2Theta(-8));    //degree  
  pGen->SetOrigin(0,0,0);                             //IP, cm 
  pGen->SetSigma(0,0,0);                              //IP sigma, cm
  pGen->Init();
}


void Gun7(Int_t iPID,Double_t p)
{
  ::Info("kir-Gun7","7 primaries of %i PID  with p=%f GeV",iPID,p);
  AliGenCocktail *pCocktail=new AliGenCocktail();
  for(int i=1;i<=7;i++){
    AliGenFixed *pFixed=new AliGenFixed(1);
    pFixed->SetMomentum(p);
    pFixed->SetPhiRange(gRICH->C(i)->PhiD());
    pFixed->SetThetaRange(gRICH->C(i)->ThetaD()+2);
    pFixed->SetOrigin(0,0,0);                 
    pFixed->SetPart(iPID);    
    pCocktail->AddGenerator(pFixed,Form("Fixed %i",i),1);
  }
  pCocktail->Init();
}

void Pythia7(Int_t iPID,Double_t p)
{
  ::Info("kir-Pythia","7 primaries of %i PID  with p=%f GeV plus Pythia",iPID,p);
  AliGenCocktail *pCocktail=new AliGenCocktail();
  for(int i=1;i<=7;i++){
    AliGenFixed *pFixed=new AliGenFixed(1);
    pFixed->SetMomentum(p);
    pFixed->SetPhiRange(gRICH->C(i)->PhiD());
    pFixed->SetThetaRange(gRICH->C(i)->ThetaD()+2);
    pFixed->SetOrigin(0,0,0);                 
    pFixed->SetPart(iPID);    
    pCocktail->AddGenerator(pFixed,Form("Fixed %i",i),1);
  }
  
  AliGenPythia *pPythia = new AliGenPythia(-1);
  pPythia->SetMomentumRange(0,999999);
  pPythia->SetPhiRange(0,360);
  pPythia->SetThetaRange(0., 180.);
  pPythia->SetYRange(-12,12);
  pPythia->SetPtRange(0,1000);
  pPythia->SetStrucFunc(kCTEQ4L);
  pPythia->SetProcess(kPyMb);
  pPythia->SetEnergyCMS(14000.);  
    
  pCocktail->AddGenerator(pPythia,"Pythia",1);  
  pCocktail->Init();
}



void Scan()
{
}


void File()
{
  ::Info("kir-File","Create galice.root, %i events per file",kEventsPerFile);

  AliRunLoader *pRL = AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,"recreate");
  if(!pRL)
    Fatal("my/AliceConfig.C::File","Can not instatiate the Run Loader");
  
  pRL->SetCompressionLevel(2);
  pRL->SetNumberOfEventsPerFile(kEventsPerFile);
  gAlice->SetRunLoader(pRL);
}

void Geant3()
{
  ::Info("kir-Geant3","Initialize the actual MC code.");

  gSystem->Load("libgeant321.so");   
  new TGeant3("C++ Interface to Geant");
}//void Geant3()

void Decayer()
{
  ::Info("kir-Decayer","Initialise external decayer.");
  TVirtualMCDecayer *pDecayer = new AliDecayerPythia();
  pDecayer->SetForceDecay(kAll);
  pDecayer->Init();
  gMC->SetExternalDecayer(pDecayer);

  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",0);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1);
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  gMC->SetProcess("RAYL",1);

  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut);
  gMC->SetCut("BCUTM",  cut);
  gMC->SetCut("DCUTE",  cut);
  gMC->SetCut("DCUTM",  cut);
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax);
}

void Field()
{
  ::Info("kir-Field","Set default magnetic field L3 0.4T.");
  //AliMagFMaps* pField = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  gAlice->SetField();
}

void Other()
{
  ::Info("kir-Other","Init all other detectors");
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");

  if(iMAG){
    new AliMAG("MAG", "Magnet");
  }

  if(iABSO){
    new AliABSOv0("ABSO", "Muon Absorber");
  }

  if(iDIPO){
    new AliDIPOv2("DIPO", "Dipole version 2");
  }

  if(iHALL){
    new AliHALL("HALL", "Alice Hall");
  }

  if(iFRAME){
    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    FRAME->SetHoles(1);
  }

  if(iSHIL){
    new AliSHILv2("SHIL", "Shielding Version 2");
  }

  if(iPIPE){
    AliPIPEv0("PIPE", "Beam Pipe");
  }
 
  if(iITS){
    AliITSvPPRasymm *ITS  = new AliITSvPPRasymm("ITS","New ITS PPR detailed version with asymmetric services");
    ITS->SetMinorVersion(2);       // don't touch this parameter if you're not an ITS developer
    ITS->SetReadDet(kFALSE);       // don't touch this parameter if you're not an ITS developer
    ITS->SetThicknessDet1(200.);   // detector thickness on layer 1 must be in the range [100,300]
    ITS->SetThicknessDet2(200.);   // detector thickness on layer 2 must be in the range [100,300]
    ITS->SetThicknessChip1(200.);  // chip thickness on layer 1 must be in the range [150,300]
    ITS->SetThicknessChip2(200.);  // chip thickness on layer 2 must be in the range [150,300]
    ITS->SetRails(0);              // 1 --> rails in ; 0 --> rails out
    ITS->SetCoolingFluid(1);       // 1 --> water ; 0 --> freon
    ITS->SetEUCLID(0);  
  }

  if(iTPC){
    AliTPC *TPC = new AliTPCv2("TPC", "Default");
    TPC->SetSecAU(-1);
    TPC->SetSecAL(-1);
  }
  
  if(iTOF){
    new AliTOFv2FHoles("TOF", "TOF with Holes");
  }
  
  if(iZDC){
    new AliZDCv2("ZDC", "normal ZDC");
  }

  if(iTRD){
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    TRD->SetGasMix(1);
    TRD->SetPHOShole();
    TRD->SetRICHhole();
    AliTRDsim *TRDsim = TRD->CreateTR();
  }

  if(iFMD){
    AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    FMD->SetRingsSi1(256);
    FMD->SetRingsSi2(128);
    FMD->SetSectorsSi1(20);
    FMD->SetSectorsSi2(40);      
  }

  if(iMUON){
    new AliMUONv1("MUON", "default");
  }
  if(iPHOS){
    new AliPHOSv1("PHOS", "IHEP");
  }

  if(iPMD){
    new AliPMDv1("PMD", "normal PMD");
  }

  if(iSTART){
    new AliSTARTv1("START", "START Detector");
  }

  if(iEMCAL){
    new AliEMCALv1("EMCAL", "EMCALArch1a");
  }

  if(iCRT){
    new AliCRTv0("CRT", "normal ACORDE");
  }

  if(iVZERO){
    new AliVZEROv2("VZERO", "normal VZERO");
  }
  
  ::Info("kir-Other","Stop.");
}//Other()

Float_t Eta2Theta(Float_t arg)
{
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

