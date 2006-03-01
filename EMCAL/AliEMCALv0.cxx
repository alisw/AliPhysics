/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

//_________________________________________________________________________
// Implementation version v0 of EMCAL Manager class 
// An object of this class does not produce hits nor digits
// It is the one to use if you do not want to produce outputs in TREEH or TREED
// This class places a Geometry of the EMCAL in the ALICE Detector as defined in AliEMCALGeometry.cxx                 
//*-- Author: Yves Schutz (SUBATECH)
//*-- and   : Sahal Yacoob (LBL / UCT)
//          : Aleksei Pavlinov (WSU)     SHASHLYK

// --- ROOT system ---

#include "TNode.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTRAP.h"
#include "TPGON.h"
#include "TTUBS.h"
#include "TGeometry.h"
#include "TVirtualMC.h"
#include "TArrayI.h"
#include "TROOT.h"
#include "TArrayF.h"
#include "TList.h"
#include "TVector2.h"

#include "AliEMCALShishKebabModule.h"
#include "AliEMCALShishKebabTrd1Module.h"

// --- Standard library ---

//#include <stdio.h>

// --- AliRoot header files ---

#include "AliEMCALv0.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliLog.h"

ClassImp(AliEMCALv0)

TArrayF ENVELOP1; // now global for BuildGeometry() - 30-aug-2004
int idAIR=1599, idPB = 1600, idSC = 1601, idSTEEL = 1603; // global
Int_t *idtmed=0, idrotm=0;
double sampleWidth=0.;
double parEMOD[5], smodPar0=0., smodPar1=0., smodPar2=0.;

//______________________________________________________________________
AliEMCALv0::AliEMCALv0(const char *name, const char *title):
  AliEMCAL(name,title)
{
  // ctor : title is used to identify the layout
  AliEMCALGeometry *geom = GetGeometry() ; 
  //geom->CreateListOfTrd1Modules(); 
  fShishKebabModules = geom->GetShishKebabTrd1Modules(); 
}

//______________________________________________________________________
void AliEMCALv0::BuildGeometry()
{
    // Display Geometry for display.C

    const Int_t kColorArm1   = kBlue ;

    AliEMCALGeometry * geom = GetGeometry();

    TString gn(geom->GetName());
    gn.ToUpper(); 

    // Define the shape of the Calorimeter 
    TNode * top = gAlice->GetGeometry()->GetNode("alice") ; // See AliceGeom/Nodes
    TNode * envelopNode = 0;
    char *envn = "Envelop1";
    if(!gn.Contains("SHISH") || gn.Contains("TRD2")){
      new TTUBS(envn, "Tubs that contains arm 1", "void", 
	      geom->GetEnvelop(0) -10, // rmin 
	      geom->GetEnvelop(1) +40 ,// rmax
	      geom->GetEnvelop(2)/2.0, // half length in Z
	      geom->GetArm1PhiMin(),   // minimum phi angle
	      geom->GetArm1PhiMax()    // maximum phi angle
	);
      top->cd();
      envelopNode = new TNode(envn, "Arm1 Envelop", "Envelop1", 0., 0., 0., "") ;
    } else {
      if(gn.Contains("WSUC")) {
        envelopNode = BuildGeometryOfWSUC();
      } else { // Shish-kebab now for compact, twist and TRD1 cases (ALIC)
        envn="Envelop2";
        TPGON *pgon = new TPGON(envn, "PGON that contains arm 1", "void", 
        geom->GetArm1PhiMin(),geom->GetArm1PhiMax()-geom->GetArm1PhiMin(),geom->GetNPhiSuperModule(), 2);
      // define section
        pgon->DefineSection(0, ENVELOP1[4],  ENVELOP1[5], ENVELOP1[6]);
        pgon->DefineSection(1, ENVELOP1[7],  ENVELOP1[5], ENVELOP1[6]);
        top->cd();
        envelopNode = new TNode(envn, "Arm1 Envelop2", envn, 0., 0., 0., "") ;
      }
    }
				
    envelopNode->SetLineColor(kColorArm1) ;
    fNodes->Add(envelopNode);
}

TNode *AliEMCALv0::BuildGeometryOfWSUC()
{ // June 8, 2005; see directory geant3/TGeant3/G3toRoot.cxx
  // enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan } - see $ROOTSYS/include/Gtypes.h
   AliEMCALGeometry * g = GetGeometry(); 
   TNode * top = gAlice->GetGeometry()->GetNode("alice") ; // See AliceGeom/Nodes
   top->cd();

   TNode *envelopNode = 0;
   char *name = "";
   /*
    name = "WSUC";
   new TBRIK(name, "WSUC(XEN1 in Geant)","void",ENVELOP1[0],ENVELOP1[1],ENVELOP1[2]);
    envelopNode = new TNode(name, "envelope for WSUC", name, 0., 0., 0., "");
   envelopNode->SetVisibility(0);
   */

   TNode *emod=0, *scmx=0;
   name = "SMOD"; // super module
   new TBRIK(name, "SMOD(SMOD in Geant)","void", smodPar0,smodPar1,smodPar2);
   if(envelopNode) envelopNode->cd();
   TNode *smod = new TNode(name, "SMOD", name, 0., 0., 0., "");   
   smod->SetLineColor(kBlue) ;
   if(envelopNode==0) envelopNode = smod;

   name = "EMOD"; // see CreateEMOD 
   TTRD1 *EMOD = new TTRD1(name, "EMOD(EMOD in Geant)","void", float(parEMOD[0]),
   float(parEMOD[1]),float(parEMOD[2]),float(parEMOD[3]));

   // SCMX = EMOD/4 for simplicity of drawing
   name = "SCMX";
   Float_t dz=0.,theta=0.,phi=0.,h1=0.,bl1=0.,tl1=0.,alpha1=0.,h2=0.,bl2=0.,tl2=0.,alpha2=0.;
   h1     = EMOD->GetDy()/2.;
   bl1    = EMOD->GetDx()/2.;
   tl1    = bl1;
   alpha1 = 0.;
   h2     = EMOD->GetDy()/2.;
   bl2    = EMOD->GetDx2()/2.;
   tl2    = bl2;
   alpha2 = 0.;

   dz       = EMOD->GetDz();
   double dr = TMath::Sqrt((h2-h1)*(h2-h1)+(bl2-bl1)*(bl2-bl1));
   theta    = TMath::ATan2(dr,2.*dz) * TMath::RadToDeg(); 
   phi      = 180.;

   TTRAP *SCMX = new TTRAP(name, "SCMX(SCMX as in Geant)","void",
   	     dz,theta,phi, h1,bl1,tl1,alpha1, h2,bl2,tl2,alpha2);
   //   SCMX->Dump();
   Float_t xShiftSCMX = (EMOD->GetDx() + EMOD->GetDx2())/4.;
   Float_t yShiftSCMX = EMOD->GetDy()/2.;
   printf(" xShiftSCMX %7.4f yShiftSCMX %7.4f \n",xShiftSCMX,yShiftSCMX);

   name = "EMOD"; // see CreateEMOD 
   smod->cd();
   
   AliEMCALShishKebabTrd1Module *mod=0;
   Double_t  angle=90., xpos=0.,ypos=0.,zpos=0., xposSCMX=0.,yposSCMX=0.,zposSCMX=0.;
   char rtmn[100], rtmt[100];
   TRotMatrix *rtm=0, *rtmSCMX=0;
   int numEmod=0;
   for(int iz=0; iz<g->GetNZ(); iz++) {
     mod   = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(iz);
     zpos = mod->GetPosZ()      - smodPar2;
     ypos = mod->GetPosXfromR() - smodPar1;

     angle = mod->GetThetaInDegree();
     sprintf(rtmn,"rmEmod%5.1f",angle);
     sprintf(rtmt,"rotation matrix for EMOD, iz=%i, angle = %6.3f",iz, angle);
     if(iz==0) rtm = new TRotMatrix(rtmn, rtmt,0.,0., 90.,0., 90.,90.); // z'(x); y'(y); x'(z)
     else      rtm = new TRotMatrix(rtmn, rtmt,90.-angle,270., 90.0,0.0, angle,90.);

     TGeometry *tg = gAlice->GetGeometry();
     for(int ix=0; ix<g->GetNPhi(); ix++) { // flat in phi
       xpos = g->GetPhiModuleSize()*(2*ix+1 - g->GetNPhi())/2.;
       sprintf(rtmt,"EMOD, iz %i, ix %i, angle %6.3f",iz,ix, angle);
       TString namNode=name;
       namNode += numEmod++;
       smod->cd();
       emod = new TNode(namNode.Data(), rtmt, (TShape*)EMOD, xpos,ypos,zpos,rtm);   
       //       emod->SetLineColor(kGreen) ;
       emod->SetVisibility(0); // SCMX will bi visible 
       if(SCMX) { // 4(2x2) sensetive volume inside EMOD
         emod->cd();
         zposSCMX = 0.;
         for(int jy=0; jy<2; jy++){ // division on y
           yposSCMX = yShiftSCMX *(2*jy - 1); 
           for(int jx=0; jx<2; jx++){ // division on x
             Double_t theta1=90.,phi1=0., theta2=90.,phi2=90., theta3=0.,phi3=0 ; 
             xposSCMX = xShiftSCMX *(2*jx - 1);              
             namNode = "SCMX";
             namNode += jy;
             namNode += jx;
             sprintf(rtmn,"rm%s",namNode.Data());
             sprintf(rtmt,"rotation matrix for %s inside EMOD",namNode.Data());
             rtmSCMX = tg->GetRotMatrix(rtmn);
             if(jx == 1) {
               phi1 = 180.; // x' = -x
               phi2 = 270.; // y' = -y
             }
             if(rtmSCMX == 0) rtmSCMX = new TRotMatrix(rtmn,rtmt, theta1,phi1, theta2,phi2, theta3,phi3);
             sprintf(rtmt,"%s inside %s", namNode.Data(), emod->GetName());
             scmx = new TNode(namNode.Data(), rtmt, (TShape*)SCMX, xposSCMX,yposSCMX,zposSCMX,rtmSCMX);   
             scmx->SetLineColor(kMagenta);
           }
         }
       }
     }
   }
   //   emod->Draw(); // for testing

   return envelopNode;
}

//______________________________________________________________________
void AliEMCALv0::CreateGeometry()
{
  // Create the EMCAL geometry for Geant
  // Geometry of a tower
  //|-----------------------------------------------------| XEN1
  //| |                                                 | |
  //| |    Al thickness = GetAlFrontThickness()         | |
  //| |                                                 | |
  //| |                                                 | |
  //| |                                                 | |
  //|  -------------------------------------------------  |
  //| |    Air Gap = GetGap2Active()                    | |
  //| |                                                 | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPST (ECAL e = GetECScintThick() )     | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPBX (ECAL e = GetECPbRadThick() )     | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPST (ECAL e = GetECScintThick()       | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPBX (ECAL e = GetECPbRadThick() )     | |
  //|  -------------------------------------------------  |
  //|    etc ..... GetNECLayers() - 1 times               |
  //|  -------------------------------------------------  |
  //| |    XUNLayer : XPST (ECAL e = GetECScintThick()  | |
  //|  -------------------------------------------------  |

    AliEMCALGeometry * geom = GetGeometry() ; 
    TString gn(geom->GetName());
    gn.ToUpper(); 

    if(!(geom->IsInitialized())){
	Error("CreateGeometry","EMCAL Geometry class has not been set up.");
    } // end if

    // Get pointer to the array containing media indices
    idtmed = fIdtmed->GetArray() - 1599 ;

    idrotm = 1;
    //  gMC->Matrix(nmat, theta1, phi1, theta2, phi2, theta3, phi3) - see AliModule
    AliMatrix(idrotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ; 

    // Create the EMCAL Mother Volume (a polygone) within which to place the Detector and named XEN1 

    Float_t envelopA[10];
    if(gn.Contains("TRD2")) { // TUBS
       envelopA[0] = geom->GetEnvelop(0) - 10.; // rmin 
       envelopA[1] = geom->GetEnvelop(1) + 12.; // rmax
       //       envelopA[2] = geom->ZFromEtaR(geom->GetEnvelop(1), geom->GetArm1EtaMin());
       envelopA[2] = 390.; // 6-feb-05
       envelopA[3] = geom->GetArm1PhiMin();
       envelopA[4] = geom->GetArm1PhiMax();
       gMC->Gsvolu("XEN1", "TUBS", idtmed[1599], envelopA, 5) ;   // Tubs filled with air 
       ENVELOP1.Set(5, envelopA);
    // Position the EMCAL Mother Volume (XEN1) in Alice (ALIC)  
       gMC->Gspos("XEN1", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
    } else if(gn.Contains("TRD1") && gn.Contains("WSUC") ) { // TRD1 for WSUC facility
      // 17-may-05 - just BOX
      envelopA[0] = 26;
      envelopA[1] = 15;
      envelopA[2] = 30;
      gMC->Gsvolu("XEN1", "BOX", idtmed[idSC], envelopA, 3) ;
      ENVELOP1.Set(3);
      for(int i=0; i<3; i++) ENVELOP1[i] = envelopA[i]; // 23-may-05  
      // Position the EMCAL Mother Volume (XEN1) in WSUC  
      gMC->Gspos("XEN1", 1, "WSUC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
    } else { 
      envelopA[0] = geom->GetArm1PhiMin();                         // minimum phi angle
      envelopA[1] = geom->GetArm1PhiMax() - geom->GetArm1PhiMin(); // angular range in phi
      envelopA[2] = geom->GetNPhi();                               // number of sections in phi
      envelopA[3] = 2;                                             // 2 z coordinates
      envelopA[4] = geom->ZFromEtaR(geom->GetEnvelop(1),
				   geom->GetArm1EtaMin());       // z coordinate 1
    //add some padding for mother volume
      envelopA[5] = geom->GetEnvelop(0) ;                          // rmin at z1
      envelopA[6] = geom->GetEnvelop(1) ;                          // rmax at z1
      envelopA[7] = geom->ZFromEtaR(geom->GetEnvelop(1),
				  geom->GetArm1EtaMax());        // z coordinate 2
      envelopA[8] = envelopA[5] ;                                  // radii are the same.
      envelopA[9] = envelopA[6] ;                                  // radii are the same.

      if(gn.Contains("SHISH")) envelopA[2] = geom->GetNPhiSuperModule();

      gMC->Gsvolu("XEN1", "PGON ", idtmed[1599], envelopA, 10) ;   // Polygone filled with air 
      ENVELOP1.Set(10, envelopA);
      if (gDebug==2) {
        printf("CreateGeometry: XEN1 = %f, %f\n", envelopA[5], envelopA[6]); 
        printf("CreateGeometry: XU0 = %f, %f\n", envelopA[5], envelopA[6]); 
      }
    // Position the EMCAL Mother Volume (XEN1) in Alice (ALIC)  
      gMC->Gspos(geom->GetNameOfEMCALEnvelope(), 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
    }

    if(gn.Contains("SHISH")){
      // COMPACT, TWIST, TRD2 or TRD1
      CreateShishKebabGeometry();
    }
}

//______________________________________________________________________
void AliEMCALv0::Init(void)
{
    // Just prints an information message
  
  if(AliLog::GetGlobalDebugLevel()>0) { 
    TString message("\n") ; 
    message += "*****************************************\n" ;
    
    // Here the EMCAL initialisation code (if any!)
    
    AliEMCALGeometry * geom = GetGeometry() ; 
    
    if (geom!=0) {   
      message += "AliEMCAL " ; 
      message += Version() ; 
      message += "EMCAL geometry initialized for " ; 
      message += geom->GetName()  ;
    }
    else {
      message += "AliEMCAL " ; 
      message += Version() ;  
      message += "EMCAL geometry initialization failed !" ; 
    }
    message += "\n*****************************************" ;
    printf(message.Data() ) ; 
  }
}

// 24-aug-04 by PAI
void AliEMCALv0::CreateShishKebabGeometry()
{  // TWIST, TRD1 and TRD2 
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper(); 
  // see AliModule::fIdtmed
  //  idtmed = fIdtmed->GetArray() - 1599 ; // see AliEMCAL::::CreateMaterials()
  //  int idAIR=1599, idPB = 1600, idSC = 1601, idSTEEL = 1603;
  //  idAL = 1602;
  Double_t par[10], xpos=0., ypos=0., zpos=0.;

  CreateSmod(g->GetNameOfEMCALEnvelope());

  CreateEmod("SMOD","EMOD"); // 18-may-05

  if(gn.Contains("110DEG")) CreateEmod("SM10","EMOD"); // 12-oct-05

  // Sensitive SC  (2x2 tiles)
  double parSCM0[5], *dummy = 0, parTRAP[11];
  Double_t trd1Angle = g->GetTrd1Angle()*TMath::DegToRad(), tanTmp = TMath::Tan(trd1Angle/2.);
   if(!gn.Contains("TRD")) { // standard module
    par[0] = (g->GetECPbRadThick()+g->GetECScintThick())*g->GetNECLayers()/2.;
    par[1] = par[2] = g->GetPhiTileSize();   // Symetric case
    gMC->Gsvolu("SCM0", "BOX", idtmed[idSC], par, 3); // 2x2 tiles
    gMC->Gspos("SCM0", 1, "EMOD", 0., 0., 0., 0, "ONLY") ;
    // Division to tile size
    gMC->Gsdvn("SCM1","SCM0", g->GetNPHIdiv(), 2); // y-axis
    gMC->Gsdvn("SCM2","SCM1", g->GetNETAdiv(), 3); // z-axis
  // put LED to the SCM2 
    par[0] = g->GetECPbRadThick()/2.;
    par[1] = par[2] = g->GetPhiTileSize()/2.; // Symetric case
    gMC->Gsvolu("PBTI", "BOX", idtmed[idPB], par, 3);

    printf(" Pb tiles \n");
    int nr=0;
    ypos = zpos = 0.0;
    xpos = -sampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;
    for(int ix=0; ix<g->GetNECLayers(); ix++){
      gMC->Gspos("PBTI", ++nr, "SCM2", xpos, ypos, zpos, 0, "ONLY") ;
    //    printf(" %i xpos %f \n", ix+1, xpos);
      xpos += sampleWidth;
    } 
    printf(" Number of Pb tiles in SCM2 %i \n", nr);
  } else if(gn.Contains("TRD1")) { // TRD1 - 30-sep-04
    if(gn.Contains("MAY05")){
      Double_t dzTmp = g->GetFrontSteelStrip()+g->GetPassiveScintThick();
      parSCM0[0] = parEMOD[0] + tanTmp*dzTmp; // dx1
      parSCM0[1] = parEMOD[1];                // dx2
      parSCM0[2] = parEMOD[2];                // dy
      for(int i=0; i<3; i++) parSCM0[i] -= g->GetLateralSteelStrip();
      parSCM0[3] = parEMOD[3] - dzTmp/2.; // dz

      gMC->Gsvolu("SCM0", "TRD1", idtmed[idAIR], parSCM0, 4);
      gMC->Gspos("SCM0", 1, "EMOD", 0., 0., dzTmp/2., 0, "ONLY") ;
    } else { // before MAY 2005
      double wallThickness = g->GetPhiModuleSize()/2. -  g->GetPhiTileSize(); // Need check
      for(int i=0; i<3; i++) parSCM0[i] = parEMOD[i] - wallThickness;
      parSCM0[3] = parEMOD[3];
      gMC->Gsvolu("SCM0", "TRD1", idtmed[idAIR], parSCM0, 4);
      gMC->Gspos("SCM0", 1, "EMOD", 0., 0., 0., 0, "ONLY") ;
    }

    if(g->GetNPHIdiv()==2 && g->GetNETAdiv()==2) {
    // Division to tile size - 1-oct-04
      printf("<I> Divide SCM0 on y-axis %i\n", g->GetNETAdiv());
      gMC->Gsdvn("SCMY","SCM0", g->GetNETAdiv(), 2); // y-axis
    // Trapesoid 2x2
      parTRAP[0] = parSCM0[3];    // dz
      parTRAP[1] = TMath::ATan2((parSCM0[1]-parSCM0[0])/2.,2.*parSCM0[3])*180./TMath::Pi(); // theta
      parTRAP[2] = 0.;           // phi
      // bottom
      parTRAP[3] = parSCM0[2]/2.; // H1
      parTRAP[4] = parSCM0[0]/2.; // BL1
      parTRAP[5] = parTRAP[4];    // TL1
      parTRAP[6] = 0.0;           // ALP1
      // top
      parTRAP[7] = parSCM0[2]/2.; // H2
      parTRAP[8] = parSCM0[1]/2.; // BL2
      parTRAP[9] = parTRAP[8];    // TL2
      parTRAP[10]= 0.0;           // ALP2
      printf(" ** TRAP ** \n");
      for(int i=0; i<11; i++) printf(" par[%2.2i] %9.4f\n", i, parTRAP[i]);

      gMC->Gsvolu("SCMX", "TRAP", idtmed[idSC], parTRAP, 11);
      xpos = +(parSCM0[1]+parSCM0[0])/4.;
      gMC->Gspos("SCMX", 1, "SCMY", xpos, 0.0, 0.0, 0, "ONLY") ;

      // Using rotation because SCMX should be the same due to Pb tiles
      xpos = -xpos; 
      AliMatrix(idrotm, 90.0,180., 90.0, 270.0, 0.0,0.0) ;
      gMC->Gspos("SCMX", 2, "SCMY", xpos, 0.0, 0.0, idrotm, "ONLY");
    // put LED to the SCM0 
      AliEMCALShishKebabTrd1Module *mod = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(0);
      gMC->Gsvolu("PBTI", "BOX", idtmed[idPB], dummy, 0);

      par[1] = parSCM0[2]/2;            // y 
      par[2] = g->GetECPbRadThick()/2.; // z

      int nr=0;
      ypos = 0.0;
      zpos = -sampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;
      double xCenterSCMX =  (parTRAP[4] +  parTRAP[8])/2.;
      if(!gn.Contains("NOPB")) { // for testing - 11-jul-05
        printf(" Pb tiles \n");
        for(int iz=0; iz<g->GetNECLayers(); iz++){
          par[0] = (parSCM0[0] + mod->GetTanBetta()*sampleWidth*iz)/2.;
          xpos   = par[0] - xCenterSCMX;
          gMC->Gsposp("PBTI", ++nr, "SCMX", xpos, ypos, zpos, 0, "ONLY", par, 3) ;
          printf(" %i xpos %f zpos %f par[0] %f \n", iz+1, xpos, zpos, par[0]);
          zpos += sampleWidth;
        } 
        printf(" Number of Pb tiles in SCMX %i \n", nr);
      }
    } else if(g->GetNPHIdiv()==3 && g->GetNETAdiv()==3) {
      Trd1Tower3X3(parSCM0);
    } else if(g->GetNPHIdiv()==4 && g->GetNETAdiv()==4) {
      Trd1Tower4X4();
    }
  } else if(gn.Contains("TRD2")) {    // TRD2 - 14-jan-05
    //    Scm0InTrd2(g, parEMOD, parSCM0); // First dessin 
    PbmoInTrd2(g, parEMOD, parSCM0); // Second dessin 
  }
}

void AliEMCALv0::CreateSmod(const char* mother)
{ // 18-may-05; mother="XEN1"; 
  // child="SMOD" from first to 10th, "SM10" (11th and 12th) (TRD1 case)
  // child="SMON" and "SMOP"("TRD2" case)
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper();

  Double_t par[3], parTubs[5], xpos=0., ypos=0., zpos=0., rpos=0., dphi=0., phi=0.0, phiRad=0.;
  Double_t par1C = 0.;
  //  ===== define Super Module from air - 14x30 module ==== ;
  sampleWidth = double(g->GetECPbRadThick()+g->GetECScintThick());
  printf("\n ## Super Module | sampleWidth %5.3f ## %s \n", sampleWidth, gn.Data());
  par[0] = g->GetShellThickness()/2.;
  par[1] = g->GetPhiModuleSize()*g->GetNPhi()/2.; 
  par[2] = g->GetEtaModuleSize()*15.; 
  idrotm=0;
  int nphism = g->GetNumberOfSuperModules()/2; // 20-may-05
  if(nphism>0) {
    dphi = (g->GetArm1PhiMax() - g->GetArm1PhiMin())/nphism;
    //    if(gn.Contains("110DEG")) dphi = (g->GetArm1PhiMax() - g->GetArm1PhiMin())/(nphism-1);
    rpos = (g->GetEnvelop(0) + g->GetEnvelop(1))/2.;
    printf(" rpos %8.2f : dphi %6.1f degree \n", rpos, dphi);
  }

  if (gn.Contains("TRD2")) { // tubs - 27-jan-05
    parTubs[0] = g->GetTubsR();                       // rmin
    parTubs[1] = parTubs[0] + g->GetShellThickness(); // rmax ?? 
    parTubs[2] = 380./2.;                             // DZ half length in z; 11-oct-04 - for 26 division
    parTubs[3] = -dphi/2.;                            // PHI1 starting angle of the segment;
    parTubs[4] = +dphi/2.;                            // PHI2 ending angle of the segment;

    gMC->Gsvolu("SMOP", "TUBS", idtmed[idAIR], parTubs, 5); // pozitive Z
    gMC->Gsvolu("SMON", "TUBS", idtmed[idAIR], parTubs, 5); // negative Z

    printf(" SMOP,N ** TUBS **\n"); 
    printf("tmed %i | Rmin %7.2f Rmax %7.2f dz %7.2f phi1,2 (%7.2f,%7.2f)\n", 
    idtmed[idAIR], parTubs[0],parTubs[1],parTubs[2], parTubs[3],parTubs[4]);
    // have to add 1 cm plastic before EMOD - time solution 
  } else if(gn.Contains("WSUC")) {
    par[0] = g->GetPhiModuleSize()*g->GetNPhi()/2.; 
    par[1] = g->GetShellThickness()/2.;
    par[2] = g->GetEtaModuleSize()*g->GetNZ()/2. + 5; 

    gMC->Gsvolu("SMOD", "BOX", idtmed[idAIR], par, 3);

    printf("SMOD in WSUC : tmed %i | dx %7.2f dy %7.2f dz %7.2f (SMOD, BOX)\n", 
    idtmed[idAIR], par[0],par[1],par[2]);
    smodPar0 = par[0]; 
    smodPar1 = par[1];
    smodPar2 = par[2];
    nphism   =  g->GetNumberOfSuperModules();
  } else {
    if     (gn.Contains("TWIST")) {
      par[2] += 0.4;      // for 27 division
    } else if(gn.Contains("TRD")) {
      par[2]  = 350./2.; // 11-oct-04 - for 26 division
      if(gn.Contains("TRD1")) {
        printf(" par[0] %7.2f (old) \n",  par[0]);  
        Float_t *parSM = g->GetSuperModulesPars(); 
        for(int i=0; i<3; i++) par[i] = parSM[i];
      }
    }
    gMC->Gsvolu("SMOD", "BOX", idtmed[idAIR], par, 3);
    printf("tmed %i | dx %7.2f dy %7.2f dz %7.2f (SMOD, BOX)\n", 
    idtmed[idAIR], par[0],par[1],par[2]);
    smodPar0 = par[0]; 
    smodPar2 = par[2];
    if(gn.Contains("110DEG")) { // 12-oct-05
      par1C = par[1];
      par[1] /= 2.;
      gMC->Gsvolu("SM10", "BOX", idtmed[idAIR], par, 3);
      printf(" Super module with name \"SM10\" was created too par[1] = %f\n", par[1]);
      par[1] = par1C;
    }
  // Steel plate
    if(g->GetSteelFrontThickness() > 0.0) { // 28-mar-05
      par[0] = g->GetSteelFrontThickness()/2.;
      gMC->Gsvolu("STPL", "BOX", idtmed[idSTEEL], par, 3);
      printf("tmed %i | dx %7.2f dy %7.2f dz %7.2f (STPL) \n", idtmed[idSTEEL], par[0],par[1],par[2]);
      xpos = -(g->GetShellThickness() - g->GetSteelFrontThickness())/2.;
      gMC->Gspos("STPL", 1, "SMOD", xpos, 0.0, 0.0, 0, "ONLY") ;
    }
  }

  int nr=0, nrsmod=0, i0=0;
  if(gn.Contains("TEST")) {nphism = 1;} // just only 2 super modules;

  // Turn whole super module
  int TurnSupMod = 1; // should be ONE; for testing =0
  for(int i=i0; i<nphism; i++) {
    if (gn.Contains("TRD2")) {      // tubs - 27-jan-05
      if(i==i0) {
        printf("** TRD2 ** ");
        if(TurnSupMod==1) printf(" No 3 degree rotation !!! ");
        printf("\n");
      }
      Double_t phic=0., phicRad=0.; // phi angle of arc center
      phic    = g->GetArm1PhiMin() + dphi*(2*i+1)/2.; //
      phicRad = phic*TMath::DegToRad();
      phi     = phic - g->GetTubsTurnAngle();
      phiRad  = phi*TMath::DegToRad();
      if(TurnSupMod==1) {
        TVector2  vc;     // position of super module center
        vc.SetMagPhi(parTubs[0], phicRad);
        TVector2  vcTurn; // position of super module center with turn
        vcTurn.SetMagPhi(parTubs[0], phiRad);
        TVector2 vcShift = vc - vcTurn;
        phic = phi;

        xpos = vcShift.X();
        ypos = vcShift.Y();
      } else { // 1-mar-05 ; just for testing - no turn od SMOD; looks good
        xpos = ypos = 0.0;
      }
      zpos = parTubs[2];
      AliMatrix(idrotm, 90.0, phic, 90.0, 90.0+phic, 0.0, 0.0);

      gMC->Gspos("SMOP", ++nr, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
      printf("SMOP %2i | %2i idrotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f \n", 
      i, nr, idrotm, phic, phicRad, xpos, ypos, zpos);
      printf(" phiy(90+phic)  %6.1f \n", 90. + phic);

      if(!gn.Contains("TEST1") && g->GetNumberOfSuperModules() > 1){
	//        double  phiy = 90. + phic + 180.;
	//        if(phiy>=360.) phiy -= 360.;
	//        printf(" phiy  %6.1f \n", phiy);
	//        AliMatrix(idrotm, 90.0, phic, 90.0, phiy, 180.0, 0.0);
        gMC->Gspos("SMON", nr, mother, xpos, ypos, -zpos, idrotm, "ONLY") ;
        printf("SMON %2i | %2i idrotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f \n", 
        i, nr, idrotm, phic, phicRad, xpos, ypos, -zpos);
      }
    } else if(gn.Contains("WSUC")) {
      xpos = ypos = zpos = 0.0;
      idrotm = 0;
      gMC->Gspos("SMOD", 1, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
      printf(" idrotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f \n", 
      idrotm, phi, phiRad, xpos, ypos, zpos);
      nr++;
    } else { // TRD1 
      TString smName("SMOD"); // 12-oct-05
      if(i==5 && gn.Contains("110DEG")) {
        smName = "SM10";
        nrsmod = nr;
        nr     = 0;
      }
      phi    = g->GetArm1PhiMin() + dphi*(2*i+1)/2.; // phi= 70, 90, 110, 130, 150, 170
      phiRad = phi*TMath::Pi()/180.;

      AliMatrix(idrotm, 90.0, phi, 90.0, 90.0+phi, 0.0, 0.0);

      xpos = rpos * TMath::Cos(phiRad);
      ypos = rpos * TMath::Sin(phiRad);
      zpos = smodPar2; // 21-sep-04
      if(i==5 && gn.Contains("110DEG")) {
        xpos += (par1C/2. * TMath::Sin(phiRad)); 
        ypos -= (par1C/2. * TMath::Cos(phiRad)); 
      }
      
      // 1th module in z-direction;
      gMC->Gspos(smName.Data(), ++nr, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
      printf(" %s : %2i idrotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f : i %i \n", 
      smName.Data(), nr, idrotm, phi, phiRad, xpos, ypos, zpos, i);
      // 2th module in z-direction;
      if(gn.Contains("TWIST") || gn.Contains("TRD")) {
      // turn arround X axis; 0<phi<360
        double phiy = 90. + phi + 180.;
        if(phiy>=360.) phiy -= 360.;
 
        AliMatrix(idrotm, 90.0, phi, 90.0, phiy, 180.0, 0.0);
        gMC->Gspos(smName.Data(), ++nr, mother, xpos, ypos, -zpos, idrotm, "ONLY");
        printf(" %s : %2i idrotm %3i phiy %6.1f  xpos %7.2f ypos %7.2f zpos %7.2f \n", 
        smName.Data(), nr, idrotm, phiy, xpos, ypos, -zpos);
      } else {
        gMC->Gspos("SMOD", ++nr, mother, xpos, ypos, -zpos, idrotm, "ONLY");
      }
    }
  }
  printf(" Number of Super Modules %i \n", nr+nrsmod);
}

void AliEMCALv0::CreateEmod(const char* mother, const char* child)
{ // 17-may-05; mother="SMOD"; child="EMOD"
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper(); 
  // Module definition
  Double_t par[10], parTubs[5], xpos=0., ypos=0., zpos=0., rpos=0.;
  Double_t parSCPA[5], zposSCPA=0.; // passive SC - 13-MAY-05, TRD1 case
  Double_t trd1Angle = g->GetTrd1Angle()*TMath::DegToRad(), tanTrd1 = TMath::Tan(trd1Angle/2.);
  Double_t tanTrd2y  = TMath::Tan(g->GetTrd2AngleY()*TMath::DegToRad()/2.);
  int nr=0;
  idrotm=0;
  if(!gn.Contains("TRD")) { // standard module
    par[0] = (sampleWidth*g->GetNECLayers())/2.; 
    par[1] = par[2] = g->GetPhiModuleSize()/2.;
    gMC->Gsvolu(child, "BOX", idtmed[idSTEEL], par, 3);

  } else if (gn.Contains("TRD1")){ // TRD1 system coordinate iz differnet
    if(strcmp(mother,"SMOD")==0) {
      parEMOD[0] = g->GetEtaModuleSize()/2.;   // dx1
      parEMOD[1] = g->Get2Trd1Dx2()/2.;        // dx2
      parEMOD[2] = g->GetPhiModuleSize()/2.;;  // dy
      parEMOD[3] = g->GetLongModuleSize()/2.;  // dz
      gMC->Gsvolu(child, "TRD1", idtmed[idSTEEL], parEMOD, 4);
      if(gn.Contains("WSUC") || gn.Contains("MAY05")){
        parSCPA[0] = g->GetEtaModuleSize()/2. + tanTrd1*g->GetFrontSteelStrip();   // dx1
        parSCPA[1] = parSCPA[0]               + tanTrd1*g->GetPassiveScintThick(); // dx2
        parSCPA[2] = g->GetPhiModuleSize()/2.;     // dy
        parSCPA[3] = g->GetPassiveScintThick()/2.; // dz
        gMC->Gsvolu("SCPA", "TRD1", idtmed[idSC], parSCPA, 4);
        zposSCPA   = -parEMOD[3] + g->GetFrontSteelStrip() + g->GetPassiveScintThick()/2.;
        gMC->Gspos ("SCPA", ++nr, child, 0.0, 0.0, zposSCPA, 0, "ONLY");
      }
    }
  } else if (gn.Contains("TRD2")){ // TRD2 as for TRD1 - 27-jan-05
    parEMOD[0] = g->GetEtaModuleSize()/2.;   // dx1
    parEMOD[1] = g->Get2Trd1Dx2()/2.;        // dx2
    parEMOD[2] = g->GetPhiModuleSize()/2.;   // dy1
    parEMOD[3] = parEMOD[2] + tanTrd2y*g->GetLongModuleSize();// dy2
    parEMOD[4] = g->GetLongModuleSize()/2.;  // dz
    gMC->Gsvolu(child, "TRD2", idtmed[idSTEEL], parEMOD, 5);
  }

  nr   = 0;
  if(gn.Contains("TWIST")) { // 13-sep-04
    fShishKebabModules = new TList;
    AliEMCALShishKebabModule *mod=0, *mTmp; // current module
    for(int iz=0; iz<g->GetNZ(); iz++) {
      //for(int iz=0; iz<4; iz++) {
      if(iz==0) {
        mod    = new AliEMCALShishKebabModule();
        idrotm = 0;
      } else {
        mTmp = new AliEMCALShishKebabModule(*mod);
        mod  = mTmp;
        Double_t  angle = mod->GetThetaInDegree();
        AliMatrix(idrotm, angle,0., 90.0,90.0, 90.-angle, 180.);
      }
      fShishKebabModules->Add(mod);

      xpos = mod->GetPosXfromR() + g->GetSteelFrontThickness() - smodPar0;
      zpos = mod->GetPosZ() - smodPar2;
      for(int iy=0; iy<g->GetNPhi(); iy++) {
        ypos = g->GetPhiModuleSize()*(2*iy+1 - g->GetNPhi())/2.;
        gMC->Gspos(child, ++nr, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
	//        printf(" %3i(%2i,2i) xpos %7.2f ypos %7.2f zpos %7.2f \n", nr,iy,iz, xpos, ypos, zpos);
      }
    }    
  } else if(gn.Contains("TRD")) { // 30-sep-04; 27-jan-05 - as for TRD1 as for TRD2
    // X->Z(0, 0); Y->Y(90, 90); Z->X(90, 0)
    AliEMCALShishKebabTrd1Module *mod=0, *mTmp; // current module

    for(int iz=0; iz<g->GetNZ(); iz++) {
      Double_t  angle=90., phiOK=0;
      if(gn.Contains("TRD1")) {
        mod = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(iz);
        angle = mod->GetThetaInDegree();
        if(!gn.Contains("WSUC")) { // ALICE 
          if(iz==0) AliMatrix(idrotm, 0.,0., 90.,90., 90.,0.); // z'(x); y'(y); x'(z)
          else      AliMatrix(idrotm, 90.-angle,180., 90.0,90.0, angle, 0.);
          phiOK = mod->GetCenterOfModule().Phi()*180./TMath::Pi(); 
	  //          printf(" %2i | angle | %6.3f - %6.3f = %6.3f(eta %5.3f)\n", 
          //iz+1, angle, phiOK, angle-phiOK, mod->GetEtaOfCenterOfModule());
          xpos = mod->GetPosXfromR() + g->GetSteelFrontThickness() - smodPar0;
          zpos = mod->GetPosZ() - smodPar2;

          int iyMax = g->GetNPhi();
          if(strcmp(mother,"SMOD") && gn.Contains("110DEG")) {
            iyMax /= 2;
          }
          for(int iy=0; iy<iyMax; iy++) { // flat in phi
            ypos = g->GetPhiModuleSize()*(2*iy+1 - iyMax)/2.;
            gMC->Gspos(child, ++nr, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
        //printf(" %2i xpos %7.2f ypos %7.2f zpos %7.2f idrotm %i\n", nr, xpos, ypos, zpos, idrotm);
            printf("%3.3i(%2.2i,%2.2i) ", nr,iy+1,iz+1);
          }
          printf("\n");
	} else {
          if(iz==0) AliMatrix(idrotm, 0.,0., 90.,0., 90.,90.); // (x')z; y'(x); z'(y)
          else      AliMatrix(idrotm, 90-angle,270., 90.0,0.0, angle,90.);
          phiOK = mod->GetCenterOfModule().Phi()*180./TMath::Pi(); 
          printf(" %2i | angle -phiOK | %6.3f - %6.3f = %6.3f(eta %5.3f)\n", 
          iz+1, angle, phiOK, angle-phiOK, mod->GetEtaOfCenterOfModule());
          zpos = mod->GetPosZ()      - smodPar2;
          ypos = mod->GetPosXfromR() - smodPar1;
          printf(" zpos %7.2f ypos %7.2f idrotm %i\n xpos ", zpos, xpos, idrotm);
          for(int ix=0; ix<g->GetNPhi(); ix++) { // flat in phi
            xpos = g->GetPhiModuleSize()*(2*ix+1 - g->GetNPhi())/2.;
            gMC->Gspos(child, ++nr, mother, xpos, ypos, zpos, idrotm, "ONLY") ;
            printf(" %7.2f ", xpos);
          }
          printf("\n");
        }
      } else if(gn.Contains("TRD2")){ // 1-feb-05 - TRD2;  curve in phi
        double angEtaRow = 0.;
	double theta1=0.,phi1=0., theta2=0.,phi2=0., theta3=0.,phi3=0.;
        angle=90.;
        if(iz==0) {
          mod   = new AliEMCALShishKebabTrd1Module();
        } else {
          mTmp  = new AliEMCALShishKebabTrd1Module(*mod);
          mod   = mTmp;
          angle = mod->GetThetaInDegree();
        }

        fShishKebabModules->Add(mod);
        phiOK = mod->GetCenterOfModule().Phi()*180./TMath::Pi(); 
        printf(" %i | theta | %6.3f - %6.3f = %6.3f\n", iz+1, angle, phiOK, angle-phiOK);

        zpos = mod->GetPosZ() - parTubs[2];
        rpos = parTubs[0] + mod->GetPosXfromR();

        angle     = mod->GetThetaInDegree();
        Double_t stepAngle = (parTubs[4] -  parTubs[3])/g->GetNPhi(); // 11-mar-04
	for(int iy=0; iy<g->GetNPhi(); iy++) {
          angEtaRow = parTubs[3] + stepAngle*(0.5+double(iy));
	  //          angEtaRow = 0;
          theta1  = 90. +  angle; phi1 = angEtaRow;      // x' ~-z;
          theta2  = 90.;          phi2 = 90. + angEtaRow;// y' ~ y;
          theta3  = angle;        phi3 = angEtaRow;      // z' ~ x;
          if(phi3 < 0.0) phi3 += 360.; 
          AliMatrix(idrotm, theta1,phi1, theta2,phi2, theta3,phi3);

          xpos = rpos * TMath::Cos(angEtaRow*TMath::DegToRad());
          ypos = rpos * TMath::Sin(angEtaRow*TMath::DegToRad());
          gMC->Gspos(child, ++nr, "SMOP", xpos, ypos, zpos, idrotm, "ONLY") ;
	  // SMON; 
	  phi1    = 180 + angEtaRow;
	  theta3  = 180.-theta3;  phi3 = angEtaRow;
          AliMatrix(idrotm, theta1,phi1, theta2,phi2, theta3,phi3);
          gMC->Gspos(child,  nr, "SMON", xpos, ypos, -zpos, idrotm, "ONLY") ;
          if(0) {
	    printf(" angEtaRow(phi) %7.2f |  angle(eta) %7.2f \n",  angEtaRow, angle);
	    printf("iy=%2i xpos %7.2f ypos %7.2f zpos %7.2f idrotm %i\n", iy, xpos, ypos, zpos, idrotm);
          }
        } // for(int iy=0; iy<g->GetNPhi(); iy++)
      }
    } 
  } else {
    xpos = g->GetSteelFrontThickness()/2.;
    for(int iz=0; iz<g->GetNZ(); iz++) {
      zpos = -smodPar2 + g->GetEtaModuleSize()*(2*iz+1)/2.;
      for(int iy=0; iy<g->GetNPhi(); iy++) {
        ypos = g->GetPhiModuleSize()*(2*iy+1 - g->GetNPhi())/2.;
        gMC->Gspos(child, ++nr, mother, xpos, ypos, zpos, 0, "ONLY") ;
      //printf(" %2i xpos %7.2f ypos %7.2f zpos %7.2f \n", nr, xpos, ypos, zpos);
      }
    }
  }
  printf(" Number of modules in Super Module %i \n", nr);
}

// 8-dec-04 by PAI
void AliEMCALv0::Trd1Tower3X3(const double parSCM0[4])
{// PB should be for whole SCM0 - ?
  double parTRAP[11], *dummy=0;
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()), scmx; 
  gn.ToUpper(); 
 // Division to tile size 
  Info("Trd1Tower3X3()","<I> Divide SCM0 on y-axis %i", g->GetNETAdiv());
  gMC->Gsdvn("SCMY","SCM0", g->GetNETAdiv(), 2); // y-axis
  double dx1=parSCM0[0], dx2=parSCM0[1], dy=parSCM0[2], dz=parSCM0[3];
  double ndiv=3., xpos=0.0;
  // should be defined once
  gMC->Gsvolu("PBTI", "BOX", idtmed[idPB], dummy, 0);
  if(gn.Contains("TEST")==0) { // one name for all trapesoid
    scmx = "SCMX"; 
    gMC->Gsvolu(scmx.Data(), "TRAP", idtmed[idSC], dummy, 0);
  }

  
  for(int ix=1; ix<=3; ix++) { // 3X3
    // ix=1
    parTRAP[0] = dz;
    parTRAP[1] = TMath::ATan2((dx2-dx1)/2.,2.*dz)*TMath::RadToDeg(); // theta
    parTRAP[2] = 0.;           // phi
    // bottom
    parTRAP[3] = dy/ndiv;      // H1
    parTRAP[4] = dx1/ndiv;     // BL1
    parTRAP[5] = parTRAP[4];   // TL1
    parTRAP[6] = 0.0;          // ALP1
    // top
    parTRAP[7] = dy/ndiv;      // H2
    parTRAP[8] = dx2/ndiv;     // BL2
    parTRAP[9] = parTRAP[8];   // TL2
    parTRAP[10]= 0.0;          // ALP2
    xpos = +(dx1+dx2)/3.;      // 6 or 3

    if      (ix==3) {
      parTRAP[1] = -parTRAP[1];
      xpos = -xpos;
    } else if(ix==2) { // central part is box but we treat as trapesoid due to numbering
      parTRAP[1] = 0.;
      parTRAP[8] = dx1/ndiv;     // BL2
      parTRAP[9] = parTRAP[8];   // TL2
      xpos = 0.0;
    }
    printf(" ** TRAP ** xpos %9.3f\n", xpos);
    for(int i=0; i<11; i++) printf(" par[%2.2i] %9.4f\n", i, parTRAP[i]);
    if(gn.Contains("TEST")){
      scmx = "SCX"; scmx += ix;
      gMC->Gsvolu(scmx.Data(), "TRAP", idtmed[idSC], parTRAP, 11);
      gMC->Gspos(scmx.Data(), 1, "SCMY", xpos, 0.0, 0.0, 0, "ONLY") ;
    } else {
      gMC->Gsposp(scmx.Data(), ix, "SCMY", xpos, 0.0, 0.0, 0, "ONLY", parTRAP, 11) ;
    }
    PbInTrap(parTRAP, scmx);
  }

  Info("Trd1Tower3X3()", "Ver. 1.0 : was tested.");
}

// 8-dec-04 by PAI
void AliEMCALv0::PbInTrap(const double parTRAP[11], TString n)
{// see for example CreateShishKebabGeometry(); just for case TRD1
  static int nr=0;
  printf(" Pb tiles : nrstart %i\n", nr);
  AliEMCALGeometry * g = GetGeometry(); 

  double par[3];
  double sampleWidth = double(g->GetECPbRadThick()+g->GetECScintThick());
  double xpos = 0.0, ypos = 0.0;
  double zpos = -sampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;

  double coef = (parTRAP[8] -  parTRAP[4]) / (2.*parTRAP[0]);
  double xCenterSCMX =  (parTRAP[4] +  parTRAP[8])/2.; // ??
  //  double tan = TMath::Tan(parTRAP[1]*TMath::DegToRad());

  par[1] = parTRAP[3];              // y 
  par[2] = g->GetECPbRadThick()/2.; // z
  for(int iz=0; iz<g->GetNECLayers(); iz++){
    par[0] = parTRAP[4] + coef*sampleWidth*iz;
    xpos   = par[0] - xCenterSCMX;
    if(parTRAP[1] < 0.) xpos = -xpos;
    gMC->Gsposp("PBTI", ++nr, n.Data(), xpos, ypos, zpos, 0, "ONLY", par, 3) ;
    printf(" %i xpos %9.3f zpos %9.3f par[0] %9.3f |", iz+1, xpos, zpos, par[0]);
    zpos += sampleWidth;
    if(iz%2>0) printf("\n");
  } 
  printf(" Number of Pb tiles in SCMX %i coef %9.7f \n", nr, coef);
  printf(" par[1] %9.3f  par[2] %9.3f ypos %9.3f \n", par[1], par[2], ypos); 
  Info("PbInTrap", "Ver. 1.0 : was tested.");
}

// 8-dec-04 by PAI
void AliEMCALv0::Trd1Tower4X4()
{// see for example CreateShishKebabGeometry()
}
// 3-feb-05
void AliEMCALv0::Scm0InTrd2(const AliEMCALGeometry * g, const Double_t parEMOD[5], Double_t parSCM0[5])
{
  // Passive material inside the detector
  double wallThickness = g->GetPhiModuleSize()/2. -  g->GetPhiTileSize(); //Need check
  printf(" wall thickness %7.5f \n", wallThickness);
  for(int i=0; i<4; i++) { // on pictures sometimes I can not see 0 -> be carefull!!
    parSCM0[i] = parEMOD[i] - wallThickness;
    printf(" %i parSCMO %7.3f parEMOD %7.3f : dif %7.3f \n", i, parSCM0[i],parEMOD[i], parSCM0[i]-parEMOD[i]);
  }
  parSCM0[4] = parEMOD[4];
  gMC->Gsvolu("SCM0", "TRD2", idtmed[idSC], parSCM0, 5); // idAIR -> idSC
  gMC->Gspos("SCM0", 1, "EMOD", 0., 0., 0., 0, "ONLY") ;
  // Division 
  if(g->GetNPHIdiv()==2 && g->GetNETAdiv()==2) {
    Division2X2InScm0(g, parSCM0);
  } else {
    Info("Scm0InTrd2"," no division SCM0 in this geometry |%s|\n", g->GetName());
    assert(0);
  }
}

void AliEMCALv0::Division2X2InScm0(const AliEMCALGeometry * g, const Double_t parSCM0[5])
{
  Double_t parTRAP[11], xpos=0.,ypos=0., dx1=0.0,dx2=0.,dy1=0.0,dy2=0.,dz=0.
  ,dr1=0.0,dr2=0.;
  idrotm=0;

  Info("Division2X2InScm0","Divide SCM0 on y-axis %i\n", g->GetNETAdiv());
  TString n("SCMX"), overLapFlagSCMY("ONLY"), overLapFlagSCMX("ONLY");
  n = "SCM0"; // for testing - 14-mar-05
  if(n=="SCM0"){
    PbInTrapForTrd2(parSCM0, n);
    // overLapFlagSCMY=overLapFlagSCMX="MANY"; // do not work
    return;
  }

  dy1 = parSCM0[2] , dy2 = parSCM0[3], dz = parSCM0[4];

  parTRAP[0] = parSCM0[4];    // dz
  parTRAP[1] = TMath::ATan2((dy2-dy1)/2.,2.*dz)*TMath::RadToDeg();
  parTRAP[2] = 90.;           // phi
  // bottom
  parTRAP[3] = parSCM0[2]/2.; // H1
  parTRAP[4] = parSCM0[0];    // BL1
  parTRAP[5] = parTRAP[4];    // TL1
  parTRAP[6] = 0.0;           // ALP1
  // top
  parTRAP[7] = parSCM0[3]/2.; // H2
  parTRAP[8] = parSCM0[1];    // BL2
  parTRAP[9] = parTRAP[8];    // TL2
  parTRAP[10]= 0.0;           // ALP2
  printf(" ** SCMY ** \n");
  for(int i=0; i<11; i++) printf(" par[%2.2i] %9.4f\n", i, parTRAP[i]);

  idrotm=0;
  gMC->Gsvolu("SCMY", "TRAP", idtmed[idSC], parTRAP, 11); // idAIR -> idSC
  ypos = +(parTRAP[3]+parTRAP[7])/2.; //
  printf(" Y shift SCMY inside SCM0 : %7.3f : opt %s\n", ypos, overLapFlagSCMY.Data()); 
  gMC->Gspos("SCMY", 1, "SCM0", 0.0, ypos, 0.0, idrotm,  overLapFlagSCMY.Data()) ;
  // Rotation SCMY around z-axis on 180 degree; x'=-x; y'=-y and z=z
  AliMatrix(idrotm, 90.0,180., 90.0, 270.0, 0.0,0.0) ;
  // We may have problem with numeration due to rotation - 4-feb-05
  gMC->Gspos("SCMY", 2, "SCM0", 0.0, -ypos, 0.0, idrotm,  overLapFlagSCMY.Data()); 

  Info("Division2X2InScm0","Divide SCMY on x-axis %i\n", g->GetNPHIdiv());
  dx1 = parSCM0[0]; 
  dx2 = parSCM0[1]; 
  dr1=TMath::Sqrt(dx1*dx1+dy1*dy1);
  dr2=TMath::Sqrt(dx2*dx2+dy2*dy2);

  parTRAP[0] = parSCM0[4];    // dz
  parTRAP[1] = TMath::ATan2((dr2-dr1)/2.,2.*dz)*TMath::RadToDeg(); // 
  parTRAP[2] = 45.;           // phi
  // bottom
  parTRAP[3] = parSCM0[2]/2.; // H1
  parTRAP[4] = parSCM0[0]/2.; // BL1
  parTRAP[5] = parTRAP[4];    // TL1
  parTRAP[6] = 0.0;           // ALP1
  // top
  parTRAP[7] = parSCM0[3]/2.; // H2
  parTRAP[8] = parSCM0[1]/2;  // BL2
  parTRAP[9] = parTRAP[8];    // TL2
  parTRAP[10]= 0.0;           // ALP2
  printf(" ** SCMX ** \n");
  for(int i=0; i<11; i++) printf(" par[%2.2i] %9.4f\n", i, parTRAP[i]);

  idrotm=0;
  gMC->Gsvolu("SCMX", "TRAP", idtmed[idSC], parTRAP, 11);
  xpos = (parTRAP[4]+parTRAP[8])/2.;
  printf(" X shift SCMX inside SCMX : %7.3f : opt %s\n", xpos, overLapFlagSCMX.Data()); 
  gMC->Gspos("SCMX", 1, "SCMY", xpos, 0.0, 0.0, idrotm,  overLapFlagSCMX.Data()) ;
  //  AliMatrix(idrotm, 90.0,270., 90.0, 0.0, 0.0,0.0); // x'=-y; y'=x; z'=z
  AliMatrix(idrotm, 90.0,90., 90.0, -180.0, 0.0,0.0);     // x'=y;  y'=-x; z'=z
  gMC->Gspos("SCMX", 2, "SCMY", -xpos, 0.0, 0.0, idrotm,  overLapFlagSCMX.Data()) ;
  // PB:
  if(n=="SCMX" && overLapFlagSCMY == "ONLY") {
    PbInTrapForTrd2(parTRAP, n);
  }
} 

// 4-feb-05 by PAI
void AliEMCALv0::PbInTrapForTrd2(const double *parTRAP, TString name)
{// TRD2 cases
  Double_t *dummy=0;
  TString pbShape("BOX"), pbtiChonly("ONLY");
  if(name=="SCM0") {
    pbShape    = "TRD2";
    //    pbtiChonly = "MANY";
  }
  gMC->Gsvolu("PBTI", pbShape.Data(), idtmed[idPB], dummy, 0);

  int nr=0;
  Info("PbInTrapForTrd2"," Pb tiles inside %s: shape %s :pbtiChonly %s\n nrstart %i\n", 
  name.Data(), pbShape.Data(), pbtiChonly.Data(), nr);
  AliEMCALGeometry * g = GetGeometry(); 

  double par[5], parPB[5], parSC[5];
  double sampleWidth = double(g->GetECPbRadThick()+g->GetECScintThick());
  double xpos = 0.0, ypos = 0.0;
  double zpos = -sampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;
  if(name == "SCMX") { // common trapezoid - 11 parameters
    double coef = (parTRAP[8] -  parTRAP[4]) / (2.*parTRAP[0]);
    double xCenterSCMX =  (parTRAP[4] +  parTRAP[8])/2.; // the same for y
    printf(" xCenterSCMX %8.5f : coef %8.7f \n", xCenterSCMX, coef);

    par[2] = g->GetECPbRadThick()/2.; // z
    for(int iz=0; iz<g->GetNECLayers(); iz++){
      par[0] = parTRAP[4] + coef*sampleWidth*iz;
      par[1] = par[0];
      xpos   = ypos = par[0] - xCenterSCMX;
    //if(parTRAP[1] < 0.) xpos = -xpos;
      gMC->Gsposp("PBTI", ++nr, name.Data(), xpos, ypos, zpos, 0, "ONLY", par, 3) ;
      printf(" %2.2i xpos %8.5f zpos %6.3f par[0,1] %6.3f |", iz+1, xpos, zpos, par[0]);
      if(iz%2>0) printf("\n");
      zpos += sampleWidth;
    } 
    printf(" Number of Pb tiles in SCMX %i coef %9.7f \n", nr, coef);
    printf(" par[1] %9.5f  par[2] %9.5f ypos %9.5f \n", par[1], par[2], ypos); 
  } else if(name == "SCM0") { // 1-mar-05 ; TRD2 - 5 parameters
    printf(" SCM0 par = ");
    for(int i=0; i<5; i++) printf(" %9.5f ", parTRAP[i]);
    printf("\n zpos %f \n",zpos);

    double tanx = (parTRAP[1] -  parTRAP[0]) / (2.*parTRAP[4]); //  tanx =  tany now
    double tany = (parTRAP[3] -  parTRAP[2]) / (2.*parTRAP[4]), ztmp=0.;
    parPB[4] = g->GetECPbRadThick()/2.;
    parSC[2] = g->GetECScintThick()/2.;
    for(int iz=0; iz<g->GetNECLayers(); iz++){
      ztmp     = sampleWidth*double(iz);
      parPB[0] = parTRAP[0] + tanx*ztmp;
      parPB[1] = parPB[0]   + tanx*g->GetECPbRadThick();
      parPB[2] = parTRAP[2] + tany*ztmp;
      parPB[3] = parPB[2]   + tany*g->GetECPbRadThick();
      gMC->Gsposp("PBTI", ++nr, name.Data(), xpos, ypos, zpos, 0, pbtiChonly.Data(), parPB, 5) ;
      printf("\n PBTI %2i | zpos %6.3f | par = ", nr, zpos);
      /*
      for(int i=0; i<5; i++) printf(" %9.5f ", parPB[i]);
      // individual SC tile
      parSC[0] = parPB[0];
      parSC[1] = parPB[1];
      gMC->Gsposp("SCTI", nr, name.Data(), xpos, ypos, zpos+g->GetECScintThick(), 
      0, pbtiChonly.Data(), parSC, 3) ;
      printf("\n SCTI     zpos %6.3f | par = ", zpos+g->GetECScintThick());
      for(int i=0; i<3; i++) printf(" %9.5f ", parPB[i]);
      */
      zpos  += sampleWidth;
    }
    printf("\n");
  }
  Info("PbInTrapForTrd2", "Ver. 0.03 : was tested.");
}

// 15-mar-05
void AliEMCALv0::PbmoInTrd2(const AliEMCALGeometry * g, const Double_t parEMOD[5], Double_t parPBMO[5])
{
  Info("PbmoInTrd2"," started : geometry %s ", g->GetName());
  double wallThickness = g->GetPhiModuleSize()/2. -  g->GetPhiTileSize();
  printf(" wall thickness %7.5f \n", wallThickness);
  for(int i=0; i<4; i++) {
    parPBMO[i] = parEMOD[i] - wallThickness;
    printf(" %i parPBMO %7.3f parEMOD %7.3f : dif %7.3f \n", 
    i, parPBMO[i],parEMOD[i], parPBMO[i]-parEMOD[i]);
  }
  parPBMO[4] = parEMOD[4];
  gMC->Gsvolu("PBMO", "TRD2", idtmed[idPB], parPBMO, 5);
  gMC->Gspos("PBMO", 1, "EMOD", 0., 0., 0., 0, "ONLY") ;
  // Division 
  if(g->GetNPHIdiv()==2 && g->GetNETAdiv()==2) {
    Division2X2InPbmo(g, parPBMO);
    printf(" PBMO, division 2X2 | geometry |%s|\n", g->GetName());
  } else {
    printf(" no division PBMO in this geometry |%s|\n", g->GetName());
    assert(0);
  }
}

void AliEMCALv0::Division2X2InPbmo(const AliEMCALGeometry * g, const Double_t parPBMO[5])
{
  Info("Division2X2InPbmo"," started : geometry %s ", g->GetName());
  //Double_t *dummy=0;
  //  gMC->Gsvolu("SCTI", "BOX", idtmed[idSC], dummy, 0);

  double parSC[3];
  double sampleWidth = double(g->GetECPbRadThick()+g->GetECScintThick());
  double xpos = 0.0, ypos = 0.0, zpos = 0.0, ztmp=0;;
  double tanx = (parPBMO[1] -  parPBMO[0]) / (2.*parPBMO[4]); //  tanx =  tany now
  double tany = (parPBMO[3] -  parPBMO[2]) / (2.*parPBMO[4]);
  char name[10], named[10], named2[10];

  printf(" PBMO par = ");
  for(int i=0; i<5; i++) printf(" %9.5f ", parPBMO[i]);
  printf("\n");

  parSC[2] = g->GetECScintThick()/2.;
  zpos = -sampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick() + g->GetECScintThick()/2.;
  printf(" parSC[2] %9.5f \n", parSC[2]);
  for(int iz=0; iz<g->GetNECLayers(); iz++){
    ztmp     = g->GetECPbRadThick() + sampleWidth*double(iz); // Z for previous PB
    parSC[0] =  parPBMO[0] + tanx*ztmp;
    parSC[1] =  parPBMO[2] + tany*ztmp;

    sprintf(name,"SC%2.2i", iz+1);
    gMC->Gsvolu(name, "BOX", idtmed[idSC], parSC, 3);
    gMC->Gspos(name, 1, "PBMO", xpos, ypos, zpos, 0, "ONLY") ;
    printf("%s | zpos %6.3f | parSC[0,1]=(%7.5f,%7.5f) -> ", 
    name, zpos, parSC[0], parSC[1]);
    
    sprintf(named,"SY%2.2i", iz+1);
    printf(" %s -> ", named);
    gMC->Gsdvn(named,name, 2, 2);

    sprintf(named2,"SX%2.2i", iz+1);
    printf(" %s \n", named2);
    gMC->Gsdvn(named2,named, 2, 1);

    zpos    += sampleWidth;
  }
}

AliEMCALShishKebabTrd1Module* AliEMCALv0::GetShishKebabModule(const Int_t neta)
{ // 28-oct-05
  AliEMCALShishKebabTrd1Module* trd1=0;
  if(fShishKebabModules && neta>=0 && neta<fShishKebabModules->GetSize()) {
    trd1 = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(neta);
  }
  return trd1;
}
