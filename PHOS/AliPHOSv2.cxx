//-*-C++-*-
//_________________________________________________________________________
// Manager and hit classes for PHOS
//*-- Author : Maxim Volkov, RRC KI
// AliPHOSv2 derives directly from AliDetector, because too much functionality
// has been put in AliPHOS for my liking.
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TBRIK.h"
#include "TNode.h"

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// --- galice header files ---
#include "AliPHOSv2.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMC.h" 

///////////////////////////////////////////////////////////////////////////////

ClassImp(AliPHOSv2)

///////////////////////////////////////////////////////////////////////////////
  
AliPHOSv2::~AliPHOSv2(void)
{

  //  fNV       = 0;
  //  fNH       = 0;
  fIshunt   = 0;
  
  delete fHits;

}

///////////////////////////////////////////////////////////////////////////////

AliPHOSv2::AliPHOSv2()
{

  //  fNV       = 0;
  //  fNH       = 0;
  fIshunt   = 0;
  fHits     = 0;
  fDigits   = 0;

}

///////////////////////////////////////////////////////////////////////////////

AliPHOSv2::AliPHOSv2(const char *name, const char *title):
  AliDetector(name,title)
{
  
  //Begin_Html
  /*
    <img src="picts/aliphos.gif">
  */
  //End_Html
  
  fHits   = new TClonesArray("AliPHOShitv2",  405);
  
  //  fNV         =  4;
  //  fNH         =  4;
  fIshunt     =  1; // All hits are associated with primary particles

  DefPars(); // Set geometry parameters
  
  //  SetMarkerColor(kGreen);
  //  SetMarkerStyle(2);
  //  SetMarkerSize(0.4);
  
}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::DefPars(void)
{

  // Initialization of GEANT3 geometry parameters
  fXtlSize[0]=2.2;
  fXtlSize[1]=18.0;
  fXtlSize[2]=2.2;

  fWrapThickness=0.01;

  fPINSize[0]=1.0;
  fPINSize[1]=0.1;
  fPINSize[2]=1.0;

  fCPVThickness=1.0;

  fPHOSFoam[0]=214.6;
  fPHOSFoam[1]=80.0;
  fPHOSFoam[2]=260.0;

  fPHOStxwall[0]=209.0;
  fPHOStxwall[1]=71.0;
  fPHOStxwall[2]=250.0;

  fPHOSAir[0]=206.0;
  fPHOSAir[1]=66.0;
  fPHOSAir[2]=244.0;

  fRadius[0]=447.0;
  fRadius[1]=460.0;

  fPHOSextra[0]=0.005; // Titanium cover thickness
  fPHOSextra[1]=6.95; // Crystal support height
  fPHOSextra[2]=4.0; // Thermo Insulating outer cover Upper plate thickness
  fPHOSextra[3]=5.0; // Upper Polystyrene Foam plate thickness
  fPHOSextra[4]=2.0; // Thermo insulating Crystal Block wall thickness
  fPHOSextra[5]=0.06; // Upper Cooling Plate thickness
  fPHOSextra[6]=10.0; // Al Support Plate thickness
  fPHOSextra[7]=3.0; // Lower Thermo Insulating Plate thickness
  fPHOSextra[8]=1.0; // Lower Textolit Plate thickness
  fPHOSextra[9]=0.03; // 1/2 total gap between adjacent crystals

  fNphi=88;
  fNz=104;

  fNModules=4;

  fPHOSAngle[0]=0.0; // Module position angles are set in CreateGeometry()
  fPHOSAngle[1]=0.0;
  fPHOSAngle[2]=0.0;
  fPHOSAngle[3]=0.0;

}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::Init(void)
{

  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" PHOS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");

  // Here the PHOS initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");

}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::BuildGeometry()
{

  // Stolen completely from A. Zvyagine

  TNode *Node, *Top;

  const int kColorPHOS = kRed;

  Top=gAlice->GetGeometry()->GetNode("alice");

  // PHOS
  Float_t pphi=12.9399462;
  new TRotMatrix("rot988","rot988",90,-3*pphi,90,90-3*pphi,0,0);
  new TRotMatrix("rot989","rot989",90,-  pphi,90,90-  pphi,0,0);
  new TRotMatrix("rot990","rot990",90,   pphi,90,90+  pphi,0,0);
  new TRotMatrix("rot991","rot991",90, 3*pphi,90,90+3*pphi,0,0);
  new TBRIK("S_PHOS","PHOS box","void",107.3,40,130);
  Top->cd();
  Node=new TNode("PHOS1","PHOS1","S_PHOS",-317.824921,-395.014343,0,"rot988");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);
  Top->cd();
  Node=new TNode("PHOS2","PHOS2","S_PHOS",-113.532333,-494.124908,0,"rot989");
  fNodes->Add(Node);
  Node->SetLineColor(kColorPHOS);
  Top->cd();
  Node=new TNode("PHOS3","PHOS3","S_PHOS", 113.532333,-494.124908,0,"rot990");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);
  Top->cd();
  Node=new TNode("PHOS4","PHOS4","S_PHOS", 317.824921,-395.014343,0,"rot991");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);

}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::CreateMaterials()
{

  // DEFINITION OF AVAILABLE PHOS MATERIALS
  
  AliMC* pMC=AliMC::GetMC();
  
  Int_t   ISXFLD=gAlice->Field()->Integ();
  Float_t SXMGMX=gAlice->Field()->Max();

  // --- The PbWO4 crystals ---
  Float_t AX[3]={207.19, 183.85, 16.0};
  Float_t ZX[3]={82.0, 74.0, 8.0};
  Float_t WX[3]={1.0, 1.0, 4.0};
  Float_t DX=8.28;

  // --- Titanium ---
  Float_t ATIT[3]={47.88, 26.98, 54.94};
  Float_t ZTIT[3]={22.0, 13.0, 25.0};
  Float_t WTIT[3]={69.0, 6.0, 1.0};
  Float_t DTIT=4.5;

  // --- The polysterene scintillator (CH) ---
  Float_t AP[2]={12.011, 1.00794};
  Float_t ZP[2]={6.0, 1.0};
  Float_t WP[2]={1.0, 1.0};
  Float_t DP=1.032;

  // --- Tyvek (CnH2n) ---
  Float_t AT[2]={12.011, 1.00794};
  Float_t ZT[2]={6.0, 1.0};
  Float_t WT[2]={1.0, 2.0};
  Float_t DT=0.331;

  // --- Polystyrene foam ---
  Float_t AF[2]={12.011, 1.00794};
  Float_t ZF[2]={6.0, 1.0};
  Float_t WF[2]={1.0, 1.0};
  Float_t DF=0.12;

  // --- Foam thermo insulation ---
  Float_t ATI[2]={12.011, 1.00794};
  Float_t ZTI[2]={6.0, 1.0};
  Float_t WTI[2]={1.0, 1.0};
  Float_t DTI=0.1;

  // --- Textolit ---
  Float_t ATX[4]={16.0, 28.09, 12.011, 1.00794};
  Float_t ZTX[4]={8.0, 14.0, 6.0, 1.0};
  Float_t WTX[4]={292.0, 68.0, 462.0, 736.0};
  Float_t DTX=1.75;

  Int_t *idtmed = gAlice->Idtmed();

  AliMixture(0, "PbWO4$", AX, ZX, DX, -3, WX);
  AliMixture(1, "Polystyrene$", AP, ZP, DP, -2, WP);
  AliMaterial(2, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0);
  // ---                     Absorption length^ is ignored ---
  AliMixture(3, "Tyvek$", AT, ZT, DT, -2, WT);
  AliMixture(4, "Foam$", AF, ZF, DF, -2, WF);
  AliMixture(5, "Titanium$", ATIT, ZTIT, DTIT, -3, WTIT);
  AliMaterial(6, "Si$", 28.09, 14., 2.33, 9.36, 42.3, 0, 0);
  AliMixture(7, "Thermo Insul.$", ATI, ZTI, DTI, -2, WTI);
  AliMixture(8, "Textolit$", ATX, ZTX, DTX, -4, WTX);
  AliMaterial(99, "Air$", 14.61, 7.3, 0.001205, 30420., 67500., 0, 0);
  
  AliMedium(700, "PHOS Xtal    $", 0, 1,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0);
  AliMedium(701, "CPV scint.   $", 1, 1,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0);
  AliMedium(702, "Al parts     $", 2, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0);
  AliMedium(703, "Tyvek wrapper$", 3, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0);
  AliMedium(704, "Polyst. foam $", 4, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0);
  AliMedium(705, "Titan. cover $", 5, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.0001, 0.0001, 0, 0);
  AliMedium(706, "Si PIN       $", 6, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.01, 0.01, 0, 0);
  AliMedium(707, "Thermo Insul.$", 7, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0);
  AliMedium(708, "Textolit     $", 8, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0);
  AliMedium(799, "Air          $", 99, 0,
	    ISXFLD, SXMGMX, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0);

  // --- Set decent energy thresholds for gamma and electron tracking
  pMC->Gstpar(idtmed[700],"CUTGAM",0.5E-4);
  pMC->Gstpar(idtmed[700],"CUTELE",1.0E-4);
  // --- Generate explicitly delta rays in the titan cover ---
  pMC->Gstpar(idtmed[705],"LOSS",3.);
  pMC->Gstpar(idtmed[705],"DRAY",1.);
  // --- and in aluminium parts ---
  pMC->Gstpar(idtmed[702],"LOSS",3.);
  pMC->Gstpar(idtmed[702],"DRAY",1.);

}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::CreateGeometry()
{

  AliMC* pMC = AliMC::GetMC();

  AliPHOSv2 *PHOS_tmp = (AliPHOSv2*)gAlice->GetModule("PHOS");
  if(PHOS_tmp==NULL){
    
    fprintf(stderr,"PHOS detector not found!\n");
    return;
    
  }

  // --- Dimensions of volumes ---
  Float_t DPHOS[3], DPTXW[3], DPAIR[3];
  Float_t DPUFP[3], DPUCP[3], DPASP[3], DPTIP[3], DPTXP[3];
  Float_t DPTCB[3], DPCBL[3], DPSTC[3], DPPAP[3], DPXTL[3], DPSUP[3],
    DPPIN[3];
  Float_t DPCPV[3], DPCPA[3];

  Float_t R, YO, XP1, YP1, PPHI, angle;
  Int_t IDROTM[99];
  Int_t i;
  
  Double_t const RADDEG=180.0/kPI;
  //  Double_t const DEGRAD=kPI/180.0;

  // --- Dimensions of PbWO4 crystal ---
  //      PARAMETER(XTL_X=2.2,XTL_Y=18.,XTL_Z=2.2)
  Float_t XTL_X=GetCrystalSize(0);
  Float_t XTL_Y=GetCrystalSize(1);
  Float_t XTL_Z=GetCrystalSize(2);

  // --- Tyvek wrapper thickness
  //      PARAMETER(PAP_THICK=0.01)
  Float_t PAP_THICK=GetWrapThickness();

  // --- Steel (titanium) cover thickness ---
  //      PARAMETER(STE_THICK=0.005)
  Float_t STE_THICK=GetPHOSextra(0);

  // --- Crystal support height ---
  //      PARAMETER(SUP_Y=6.95)
  Float_t SUP_Y=GetPHOSextra(1);

  // --- PIN-diode dimensions ---
  //      PARAMETER(PIN_X=1.4,PIN_Y=0.4,PIN_Z=1.4)
  Float_t PIN_X=GetPINSize(0);
  Float_t PIN_Y=GetPINSize(1);
  Float_t PIN_Z=GetPINSize(2);

  // --- CPV thickness ---
  //      PARAMETER(CPV_Y=0.5)
  Float_t CPV_Y=GetCPVThickness();

  // --- Foam Thermo Insulating outer cover dimensions ---
  //      PARAMETER(FTI_X=214.6,FTI_Y=80.,FTI_Z=260.)
  Float_t FTI_X=GetPHOSFoam(0);
  Float_t FTI_Y=GetPHOSFoam(1);
  Float_t FTI_Z=GetPHOSFoam(2);

  // --- Thermo Insulating outer cover Upper plate thickness ---
  //      PARAMETER(FTIU_THICK=4.)
  Float_t FTIU_THICK=GetPHOSextra(2);

  // --- Textolit Wall box dimentions ---
  //      PARAMETER(TXW_X=209.,TXW_Y=71.,TXW_Z=250.)
  Float_t TXW_X=GetPHOStxwall(0);
  Float_t TXW_Y=GetPHOStxwall(1);
  Float_t TXW_Z=GetPHOStxwall(2);

  // --- Inner AIR volume dimensions ---
  //      PARAMETER(AIR_X=206.,AIR_Y=66.,AIR_Z=244.)
  Float_t AIR_X=GetPHOSAir(0);
  Float_t AIR_Y=GetPHOSAir(1);
  Float_t AIR_Z=GetPHOSAir(2);

  // --- Upper Polystyrene Foam plate thickness ---
  //      PARAMETER(UFP_Y=5.)
  Float_t UFP_Y=GetPHOSextra(3);

  // --- Thermo insulating Crystal Block wall thickness ---
  //      PARAMETER(TCB_THICK=2.)
  Float_t TCB_THICK=GetPHOSextra(4);

  // --- Upper Cooling Plate thickness ---
  //      PARAMETER(UCP_Y=0.06)
  Float_t UCP_Y=GetPHOSextra(5);

  // --- Al Support Plate thickness ---
  //      PARAMETER(ASP_Y=10.)
  Float_t ASP_Y=GetPHOSextra(6);

  // --- Lower Thermo Insulating Plate thickness ---
  //      PARAMETER(TIP_Y=3.)
  Float_t TIP_Y=GetPHOSextra(7);

  // --- Lower Textolit Plate thickness ---
  //      PARAMETER(TXP_Y=1.)
  Float_t TXP_Y=GetPHOSextra(8);

  // --- 1/2 total gap between adjacent crystals
  Float_t TOTAL_GAP=GetPHOSextra(9);

  // --- Distance from IP to Foam Thermo Insulating top plate ---
  //      PARAMETER(FTI_R=467.)
  Float_t FTI_R=GetRadius(0);

  // --- Distance from IP to Crystal Block top Surface (needs to be 460.) ---
  //      PARAMETER(CBS_R=480.)
  Float_t CBS_R=GetRadius(1);

  // Get pointer to the array containing media indeces
  Int_t *IDTMED = gAlice->Idtmed();

  // --- Define PHOS box volume, fill with thermo insulating foam ---
  DPHOS[0]=FTI_X/2.0;
  DPHOS[1]=FTI_Y/2.0;
  DPHOS[2]=FTI_Z/2.0;
  pMC->Gsvolu("PHOS", "BOX ", IDTMED[706], DPHOS, 3);

  // --- Define Textolit Wall box, position inside PHOS ---
  DPTXW[0]=TXW_X/2.0;
  DPTXW[1]=TXW_Y/2.0;
  DPTXW[2]=TXW_Z/2.0;
  pMC->Gsvolu("PTXW", "BOX ", IDTMED[707], DPTXW, 3);
  YO=(FTI_Y-TXW_Y)/2.0-FTIU_THICK;
  pMC->Gspos("PTXW", 1, "PHOS", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Upper Polystyrene Foam Plate, place inside PTXW ---
  // --- immediately below Foam Thermo Insulation Upper plate ---
  DPUFP[0]=TXW_X/2.0;
  DPUFP[1]=UFP_Y/2.0;
  DPUFP[2]=TXW_Z/2.0;
  pMC->Gsvolu("PUFP", "BOX ", IDTMED[703], DPUFP, 3);
  YO=(TXW_Y-UFP_Y)/2.0;
  pMC->Gspos("PUFP", 1, "PTXW", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define air-filled box, place inside PTXW ---
  DPAIR[0]=AIR_X/2.0;
  DPAIR[1]=AIR_Y/2.0;
  DPAIR[2]=AIR_Z/2.0;
  pMC->Gsvolu("PAIR", "BOX ", IDTMED[798], DPAIR, 3);
  YO=(TXW_Y-AIR_Y)/2.0-UFP_Y;
  pMC->Gspos("PAIR", 1, "PTXW", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Thermo insulating Crystal Box, position inside PAIR ---
  DPTCB[0]=GetNphi()*(XTL_X+2*TOTAL_GAP)/2.0+TCB_THICK;
  DPTCB[1]=(XTL_Y+SUP_Y+PAP_THICK+STE_THICK)/2.0+TCB_THICK/2.0;
  DPTCB[2]=GetNz()*(XTL_Z+2*TOTAL_GAP)/2.0+TCB_THICK;
  pMC->Gsvolu("PTCB", "BOX ", IDTMED[706], DPTCB, 3);
  YO=AIR_Y/2.0-DPTCB[1]-
    (CBS_R-FTI_R-TCB_THICK-FTIU_THICK-UFP_Y);
  pMC->Gspos("PTCB", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Crystal BLock filled with air, position it inside PTCB ---
  DPCBL[0]=GetNphi()*(XTL_X+2*TOTAL_GAP)/2.0;
  DPCBL[1]=(XTL_Y+SUP_Y+PAP_THICK+STE_THICK)/2.0;
  DPCBL[2]=GetNz()*(XTL_Z+2*TOTAL_GAP)/2.0;
  pMC->Gsvolu("PCBL", "BOX ", IDTMED[798], DPCBL, 3);
  
  // --- Divide PCBL in X (phi) and Z directions --
  pMC->Gsdvn("PROW", "PCBL", Int_t (GetNphi()), 1);
  pMC->Gsdvn("PCEL", "PROW", Int_t (GetNz()), 3);
  YO=-TCB_THICK/2.0;
  pMC->Gspos("PCBL", 1, "PTCB", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define STeel (actually, it's titanium) Cover volume, place inside PCEL
  DPSTC[0]=(XTL_X+2*PAP_THICK)/2.0;
  DPSTC[1]=(XTL_Y+SUP_Y+PAP_THICK+STE_THICK)/2.0;
  DPSTC[2]=(XTL_Z+2*PAP_THICK+2*STE_THICK)/2.0;
  pMC->Gsvolu("PSTC", "BOX ", IDTMED[704], DPSTC, 3);
  pMC->Gspos("PSTC", 1, "PCEL", 0.0, 0.0, 0.0, 0, "ONLY");

  // --- Define Tyvek volume, place inside PSTC ---
  DPPAP[0]=XTL_X/2.0+PAP_THICK;
  DPPAP[1]=(XTL_Y+SUP_Y+PAP_THICK)/2.0;
  DPPAP[2]=XTL_Z/2.0+PAP_THICK;
  pMC->Gsvolu("PPAP", "BOX ", IDTMED[702], DPPAP, 3);
  YO=(XTL_Y+SUP_Y+PAP_THICK)/2.0-(XTL_Y+SUP_Y+PAP_THICK+STE_THICK)/2.0;
  pMC->Gspos("PPAP", 1, "PSTC", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define PbWO4 crystal volume, place inside PPAP ---
  DPXTL[0]=XTL_X/2.0;
  DPXTL[1]=XTL_Y/2.0;
  DPXTL[2]=XTL_Z/2.0;
  pMC->Gsvolu("PXTL", "BOX ", IDTMED[699], DPXTL, 3);
  YO=(XTL_Y+SUP_Y+PAP_THICK)/2.0-XTL_Y/2.0-PAP_THICK;
  pMC->Gspos("PXTL", 1, "PPAP", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define crystal support volume, place inside PPAP ---
  DPSUP[0]=XTL_X/2.0+PAP_THICK;
  DPSUP[1]=SUP_Y/2.0;
  DPSUP[2]=XTL_Z/2.0+PAP_THICK;
  pMC->Gsvolu("PSUP", "BOX ", IDTMED[798], DPSUP, 3);
  YO=SUP_Y/2.0-(XTL_Y+SUP_Y+PAP_THICK)/2.0;
  pMC->Gspos("PSUP", 1, "PPAP", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define PIN-diode volume and position it inside crystal support ---
  // --- right behind PbWO4 crystal
  DPPIN[0]=PIN_X/2.0;
  DPPIN[1]=PIN_Y/2.0;
  DPPIN[2]=PIN_Z/2.0;
  pMC->Gsvolu("PPIN", "BOX ", IDTMED[705], DPPIN, 3);
  YO=SUP_Y/2.0-PIN_Y/2.0;
  pMC->Gspos("PPIN", 1, "PSUP", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Upper Cooling Panel, place it on top of PTCB ---
  DPUCP[0]=DPTCB[0];
  DPUCP[1]=UCP_Y/2.0;
  DPUCP[2]=DPTCB[2];
  pMC->Gsvolu("PUCP", "BOX ", IDTMED[701], DPUCP,3);
  YO=(AIR_Y-UCP_Y)/2.0-(CBS_R-FTI_R-TCB_THICK-FTIU_THICK-UFP_Y-UCP_Y);
  pMC->Gspos("PUCP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Al Support Plate, position it inside PAIR ---
  // --- right beneath PTCB ---
  DPASP[0]=AIR_X/2.0;
  DPASP[1]=ASP_Y/2.0;
  DPASP[2]=AIR_Z/2.0;
  pMC->Gsvolu("PASP", "BOX ", IDTMED[701], DPASP, 3);
  YO=(AIR_Y-ASP_Y)/2.0-(CBS_R-FTI_R-FTIU_THICK-UFP_Y+DPCBL[1]*2);
  pMC->Gspos("PASP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Thermo Insulating Plate, position it inside PAIR ---
  // --- right beneath PASP ---
  DPTIP[0]=AIR_X/2.0;
  DPTIP[1]=TIP_Y/2.0;
  DPTIP[2]=AIR_Z/2.0;
  pMC->Gsvolu("PTIP", "BOX ", IDTMED[706], DPTIP, 3);
  YO=(AIR_Y-TIP_Y)/2.0-(CBS_R-FTI_R-FTIU_THICK-UFP_Y+DPCBL[1]*2+ASP_Y);
  pMC->Gspos("PTIP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define Textolit Plate, position it inside PAIR ---
  // --- right beneath PTIP ---
  DPTXP[0]=AIR_X/2.0;
  DPTXP[1]=TXP_Y/2.0;
  DPTXP[2]=AIR_Z/2.0;
  pMC->Gsvolu("PTXP", "BOX ", IDTMED[707], DPTXP, 3);
  YO=(AIR_Y-TXP_Y)/2.0-
    (CBS_R-FTI_R-FTIU_THICK-UFP_Y+DPCBL[1]*2+ASP_Y+TIP_Y);
  pMC->Gspos("PTXP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY");

  // --- Define CPV volume, DON'T PLACE IT YET ---
  // --- Divide in X and Z direction (same way as PCBL) ---
  DPCPV[0]=DPCBL[0];
  DPCPV[1]=CPV_Y/2.0;
  DPCPV[2]=DPCBL[2];
  //  pMC->Gsvolu("PCPV", "BOX ", IDTMED[700], DPCPV, 3);
  pMC->Gsvolu("PCPV", "BOX ", IDTMED[798], DPCPV, 3);
  pMC->Gsdvn("PCRO", "PCPV", Int_t (GetNphi()), 1);
  pMC->Gsdvn("PCCE", "PCRO", Int_t (GetNz()), 3);

  // Define CPV sensitive pad. It has the same size as PCCE.
  DPCPA[0]=DPCBL[0]/GetNphi();
  DPCPA[1]=CPV_Y/2.0;
  DPCPA[2]=DPCBL[2]/GetNz();
  pMC->Gsvolu("PCPA", "BOX ", IDTMED[700], DPCPA, 3);
  pMC->Gspos("PCPA", 1, "PCCE", 0.0, 0.0, 0.0, 0, "ONLY");

  // --- Position various PHOS units in ALICE setup ---
  // --- PHOS itself first ---
  PPHI=TMath::ATan(FTI_X/(2.0*FTI_R));
  PPHI*=RADDEG;

  for(i=1; i<=GetNModules(); i++){

    angle=PPHI*2*(i-GetNModules()/2.0-0.5);
    AliMatrix(IDROTM[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0);

    // --- Position various PHOS units in ALICE setup ---
    // --- PHOS itself first ---
    R=FTI_R+FTI_Y/2.0;
    XP1=R*TMath::Sin(angle/RADDEG);
    YP1=-R*TMath::Cos(angle/RADDEG);
    pMC->Gspos("PHOS", i, "ALIC", XP1, YP1, 0.0, IDROTM[i-1], "ONLY");

    // --- Now position PCPV so that its plates are right on top of ---
    // --- corresponding PHOS modules (previously called cradles) ---
    R=FTI_R-CPV_Y/2.0;
    XP1=R*TMath::Sin(angle/RADDEG);
    YP1=-R*TMath::Cos(angle/RADDEG);
    pMC->Gspos("PCPV", i, "ALIC", XP1, YP1, 0.0, IDROTM[i-1], "ONLY");
    GetModuleAngle(i-1)=angle-90.0;

  }

  // --- Set volumes seen without their descendants for drawing ---
  pMC->Gsatt("PCEL", "SEEN", -2);
  pMC->Gsatt("PCCE", "SEEN", -2);

}

//////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::StepManager(void)
{

  AliMC *pMC=AliMC::GetMC();
  Int_t blrc[4]; // (box, layer, row, column) indices
  Float_t xyze[4]; // position wrt MRS and energy deposited

  Int_t *IDTMED=gAlice->Idtmed();

  if(pMC->GetMedium()==IDTMED[700]){ // We are inside a CPV sensitive pad

    pMC->TrackPosition(xyze);
    xyze[3]=pMC->Edep();
    
    pMC->CurrentVolOff(3, (Text_t*)NULL, blrc[0]);
    blrc[1]=1; // CPV corresponds to layer 1
    pMC->CurrentVolOff(2, (Text_t*)NULL, blrc[2]);
    pMC->CurrentVolOff(1, (Text_t*)NULL, blrc[3]);
    
    AddHit(gAlice->CurrentTrack(), blrc, xyze);

  }

  if(pMC->GetMedium()==IDTMED[699]){ // We are inside a PWO crystal

    pMC->TrackPosition(xyze);
    xyze[3]=pMC->Edep();

    pMC->CurrentVolOff(9, (Text_t*)NULL, blrc[0]);
    blrc[1]=2; // PWO crystals correspond to layer 2
    pMC->CurrentVolOff(4, (Text_t*)NULL, blrc[2]);
    pMC->CurrentVolOff(3, (Text_t*)NULL, blrc[3]);

    AddHit(gAlice->CurrentTrack(), blrc, xyze);

  }
}

////////////////////////////////////////////////////////////////////////////

void AliPHOSv2::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{

  Int_t hitCounter;
  TClonesArray &lhits = *fHits;
  AliPHOShitv2 *newHit,*curHit;
  
  newHit=new AliPHOShitv2(fIshunt, track, vol, hits);
  
  for(hitCounter=0;hitCounter<fNhits;hitCounter++){
    curHit=(AliPHOShitv2*)lhits[hitCounter];
    if(*curHit==*newHit){
      *curHit=*curHit+*newHit;
      delete newHit;
      return;
    }
  }
  
  new(lhits[fNhits++]) AliPHOShitv2(*newHit);
  delete newHit;

}

///////////////////////////////////////////////////////////////////////////////

ClassImp(AliPHOShitv2)

//////////////////////////////////////////////////////////////////////////////

AliPHOShitv2::AliPHOShitv2(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt, track)
{

   Int_t i;
   for (i=0;i<4;i++)fVolume[i]=vol[i];

   fX       = hits[0];
   fY       = hits[1];
   fZ       = hits[2];
   fELOS    = hits[3];

}

//////////////////////////////////////////////////////////////////////////////

Bool_t AliPHOShitv2::operator==(AliPHOShitv2 const &rValue) const
{

  Int_t volCounter;

  //  if(fDet!=rValue.GetDet()) return kFALSE;
  if(fTrack!=rValue.GetTrack()) return kFALSE;

  for(volCounter=0;volCounter<4;volCounter++)
    if(fVolume[volCounter]!=rValue.GetVolume(volCounter))return kFALSE;

  return kTRUE;

}

/////////////////////////////////////////////////////////////////////////////

AliPHOShitv2 const AliPHOShitv2::operator+(AliPHOShitv2 const &rValue) const
{

  AliPHOShitv2 added(*this);

  added.fELOS+=rValue.GetEnergy();
  return added;

}
 
//////////////////////////////////////////////////////////////////////////////
