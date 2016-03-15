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
//          : Alexei Pavlinov (WSU)     SHASHLYK
//          : Adapted for DCAL by M.L. Wang CCNU & Subatech Oct-19-2012
// --- ROOT system ---
#include <cassert>

#include <TGeometry.h>
#include <TGeoPhysicalNode.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TVirtualMC.h>
#include <TArrayI.h>
#include <TROOT.h>
#include <TList.h>
#include <TVector2.h>
#include <cassert>

// --- Standard library ---

//#include <stdio.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliLog.h"
#include "AliGeomManager.h"

//--- EMCAL system---
#include "AliEMCALShishKebabTrd1Module.h"
#include "AliEMCALv0.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALSpaceFrame.h"



ClassImp(AliEMCALv0)

// EMCAL material: look to the AliEMCAL.cxx
enum
 {
  kIdAIR   = 1599, 
  kIdPB    = 1600, 
  kIdSC    = 1601, 
  kIdAL    = 1602, 
  kIdSTEEL = 1603,
  kIdPAPER = 1604
 };


//______________________________________________________________________
AliEMCALv0::AliEMCALv0()
  : AliEMCAL(),
    fShishKebabModules(),fEnvelop1(0),fIdRotm(0),fIdTmedArr(0),
    fSampleWidth(0),fSmodPar0(0),fSmodPar1(0),fSmodPar2(0),
    fInnerEdge(0),fCalFrame(0)
{
  //default ctor
  for(Int_t i = 0; i < 5 ; i++) fParEMOD[i]=0.0;
}

//______________________________________________________________________
AliEMCALv0::AliEMCALv0(const char *name, const char *title, 
                       const Bool_t checkGeoAndRun)
  : AliEMCAL(name,title,checkGeoAndRun),
    fShishKebabModules(),fEnvelop1(0),fIdRotm(0),fIdTmedArr(0),
    fSampleWidth(0),fSmodPar0(0),fSmodPar1(0),fSmodPar2(0),
    fInnerEdge(0),fCalFrame(0)
{
  // ctor : title is used to identify the layout
  // Apr 25, 2006
  // Nov 22, 2006 - case of 1X1  
  
  for(Int_t i = 0; i < 5 ; i++) fParEMOD[i]=0.0;
  
  TString ntmp(GetTitle());
  ntmp.ToUpper();

  fGeometry = GetGeometry() ; 
    
  TString gn(fGeometry->GetName()); 
  gn.ToUpper();
  
  fShishKebabModules = fGeometry->GetShishKebabTrd1Modules(); 

  fSampleWidth = Double_t(fGeometry->GetECPbRadThick() + fGeometry->GetECScintThick());
  
  if(gn.Contains("V1")) fSampleWidth += 2.*fGeometry->GetTrd1BondPaperThick();
  
  AliDebug(2,Form("fGeometry %p : TVirtualMC::GetMC() %p : fSampleWidth %5.4f\n", 
	 fGeometry, TVirtualMC::GetMC(), fSampleWidth));
  
  // Set geometry name again, in case it was changed during the initialization of the geometry.
  SetTitle(fGeometry->GetEMCGeometry()->GetName());
}

//______________________________________________________________________
void AliEMCALv0::CreateGeometry()
{
  // Create the EMCAL geometry for Geant
  // Geometry of a tower

  AliEMCALGeometry * geom = GetGeometry() ; 
  TString gn(geom->GetName());
  gn.ToUpper(); 
  
  if(!(geom->IsInitialized())){
    Error("CreateGeometry","EMCAL Geometry class has not been set up.");
  } // end if

  // Get pointer to the array containing media indices
  fIdTmedArr = fIdtmed->GetArray() - 1599 ;
  
  fIdRotm = 1;
  //  TVirtualMC::GetMC()->Matrix(nmat, theta1, phi1, theta2, phi2, theta3, phi3) - see AliModule
  AliMatrix(fIdRotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ; 
  
  // Create the EMCAL Mother Volume (a polygone) within which to place the Detector and named XEN1 
  
  Float_t envelopA[10];
  if(gn.Contains("WSUC") ) { // TRD1 for WSUC facility
    // Nov 25,2010
    envelopA[0] = 30.;
    envelopA[1] = 30;
    envelopA[2] = 20;
    TVirtualMC::GetMC()->Gsvolu("XEN1", "BOX", fIdTmedArr[kIdSC], envelopA, 3) ;
    fEnvelop1.Set(3);
    for(Int_t i=0; i<3; i++) fEnvelop1[i] = envelopA[i]; // 23-may-05  
    // Position the EMCAL Mother Volume (XEN1) in WSUC.
    // Look to AliEMCALWsuCosmicRaySetUp.  
    TVirtualMC::GetMC()->Gspos("XEN1", 1, "WSUC", 0.0, 0.0, + 265., fIdRotm, "ONLY") ;
  } else { 
    envelopA[0] = geom->GetArm1PhiMin();                         // minimum phi angle
    envelopA[1] = geom->GetArm1PhiMax() - geom->GetArm1PhiMin(); // angular range in phi
    envelopA[2] = envelopA[1]/geom->GetEMCGeometry()->GetPhiSuperModule();	 // Section of that
    envelopA[3] = 2;                                             // 2: z coordinates
    envelopA[4] = -geom->GetEnvelop(2)/2.;                       // zmin - includes padding
    envelopA[5] = geom->GetEnvelop(0) ;                          // rmin at z1 - includes padding
    envelopA[6] = geom->GetEnvelop(1) ;                          // rmax at z1 - includes padding
    envelopA[7] = geom->GetEnvelop(2)/2.;                        // zmax includes padding
    
    envelopA[8] = envelopA[5] ;                                  // radii are the same.
    envelopA[9] = envelopA[6] ;                                  // radii are the same.

    TVirtualMC::GetMC()->Gsvolu("XEN1", "PGON", fIdTmedArr[kIdAIR], envelopA, 10) ;   // Polygone filled with air 
    fEnvelop1.Set(10, envelopA);
    if (gDebug==2) {
      printf("CreateGeometry: XEN1 = %f, %f\n", envelopA[5], envelopA[6]); 
      printf("CreateGeometry: XU0 = %f, %f\n", envelopA[5], envelopA[6]); 
    }
    // Position the EMCAL Mother Volume (XEN1) in Alice (ALIC)  
    TVirtualMC::GetMC()->Gspos(geom->GetNameOfEMCALEnvelope(), 1, "ALIC", 0.0, 0.0, 0.0, fIdRotm, "ONLY") ;
  }

  // COMPACT, TRD1
  AliDebug(2,Form("Shish-Kebab geometry : %s", GetTitle())); 
  CreateShishKebabGeometry();

  if(gn.Contains("WSUC")==0) { // Nov 24,2010 for TB
  //Space Frame
    AliDebug(2,"Creating EMCAL Space Frame");
    fCalFrame = new AliEMCALSpaceFrame();
    fCalFrame->CreateGeometry();
  }
  
  // Set the sampling fraction used at creation hit level
  // Previously called in AliEMCALEMCGeometry::Init(), put it here for proper initialization by Geant3/4
  geom->GetEMCGeometry()->DefineSamplingFraction(TVirtualMC::GetMC()->GetName(),TVirtualMC::GetMC()->GetTitle());
}

//______________________________________________________________________
void AliEMCALv0::Init(void)
{
    // Just prints an information message
  AliEMCAL::Init();
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
    printf("%s",message.Data() ) ; 
  }
}

// 24-aug-04 by PAI
//______________________________________________________________________
void AliEMCALv0::CreateShishKebabGeometry()
{  
  // Oct 26,2010
  // TRD1
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper(); 
  Double_t trd1Angle = g->GetTrd1Angle()*TMath::DegToRad(), tanTrd1 = TMath::Tan(trd1Angle/2.);
  // see AliModule::fFIdTmedArr
  //  fIdTmedArr = fIdtmed->GetArray() - 1599 ; // see AliEMCAL::::CreateMaterials()
  //  Int_t kIdAIR=1599, kIdPB = 1600, kIdSC = 1601, kIdSTEEL = 1603;
  //  idAL = 1602;
  Double_t par[10], xpos=0., ypos=0., zpos=0.;

  CreateSmod(g->GetNameOfEMCALEnvelope());

  Int_t * SMTypeList = g->GetEMCSystem();
  Int_t tmpType = -1;
  for(Int_t i = 0 ; i < g->GetNumberOfSuperModules(); i++) {
    if( SMTypeList[i] == tmpType) continue;
    else  tmpType = SMTypeList[i];

    if( tmpType == AliEMCALGeometry::kEMCAL_Standard ) CreateEmod("SMOD","EMOD");   // 18-may-05 
    else if( tmpType == AliEMCALGeometry::kEMCAL_Half     ) CreateEmod("SM10","EMOD");   // Nov 1,2006 1/2 SM
    else if( tmpType == AliEMCALGeometry::kEMCAL_3rd      ) CreateEmod("SM3rd","EMOD");  // Feb 1,2012 1/3 SM
    else if( tmpType == AliEMCALGeometry::kDCAL_Standard  ) CreateEmod("DCSM","EMOD");   // Mar 13, 2012, 6 or 10 DCSM
    else if( tmpType == AliEMCALGeometry::kDCAL_Ext       ) CreateEmod("DCEXT","EMOD");  // Mar 13, 2012, DCAL extension SM
    else AliError("Unkown SM Type!!");
  }

  // Sensitive SC  (2x2 tiles)
  Double_t parSCM0[5]={0,0,0,0}, *dummy = 0, parTRAP[11];

  if(!gn.Contains("V1")) {
    Double_t wallThickness = g->GetPhiModuleSize()/g->GetNPHIdiv() -  g->GetPhiTileSize();
    for(Int_t i=0; i<3; i++) parSCM0[i] = fParEMOD[i] - wallThickness;
    parSCM0[3] = fParEMOD[3];
    TVirtualMC::GetMC()->Gsvolu("SCM0", "TRD1", fIdTmedArr[kIdAIR], parSCM0, 4);
    TVirtualMC::GetMC()->Gspos("SCM0", 1, "EMOD", 0., 0., 0., 0, "ONLY") ;
  } else {
    Double_t wTh = g->GetLateralSteelStrip();
    parSCM0[0] = fParEMOD[0] - wTh + tanTrd1*g->GetTrd1AlFrontThick();
    parSCM0[1] = fParEMOD[1] - wTh;
    parSCM0[2] = fParEMOD[2] - wTh;
    parSCM0[3] = fParEMOD[3] - g->GetTrd1AlFrontThick()/2.;
    TVirtualMC::GetMC()->Gsvolu("SCM0", "TRD1", fIdTmedArr[kIdAIR], parSCM0, 4);
    Double_t zshift = g->GetTrd1AlFrontThick()/2.;
    TVirtualMC::GetMC()->Gspos("SCM0", 1, "EMOD", 0., 0., zshift, 0, "ONLY");
    // 
    CreateAlFrontPlate("EMOD","ALFP");
  }

  if(g->GetNPHIdiv()==2 && g->GetNETAdiv()==2) {
    // Division to tile size - 1-oct-04
    AliDebug(2,Form(" Divide SCM0 on y-axis %i\n", g->GetNETAdiv()));
    TVirtualMC::GetMC()->Gsdvn("SCMY","SCM0", g->GetNETAdiv(), 2); // y-axis
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
    AliDebug(2,Form(" ** TRAP ** \n"));
    for(Int_t i=0; i<11; i++) AliDebug(3, Form(" par[%2.2i] %9.4f\n", i, parTRAP[i]));
    
    TVirtualMC::GetMC()->Gsvolu("SCMX", "TRAP", fIdTmedArr[kIdSC], parTRAP, 11);
    xpos = +(parSCM0[1]+parSCM0[0])/4.;
    TVirtualMC::GetMC()->Gspos("SCMX", 1, "SCMY", xpos, 0.0, 0.0, 0, "ONLY") ;

    // Using rotation because SCMX should be the same due to Pb tiles
    xpos = -xpos; 
    AliMatrix(fIdRotm, 90.0,180., 90.0, 270.0, 0.0,0.0) ;
    TVirtualMC::GetMC()->Gspos("SCMX", 2, "SCMY", xpos, 0.0, 0.0, fIdRotm, "ONLY");
    // put LED to the SCM0 
    AliEMCALShishKebabTrd1Module *mod = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(0);
    Double_t tanBetta = mod->GetTanBetta();

    Int_t nr=0; 
    ypos = 0.0; 
    Double_t xCenterSCMX =  (parTRAP[4] +  parTRAP[8])/2.;
    if(!gn.Contains("V1")) {
      par[1] = parSCM0[2]/2;            // y 
      par[2] = g->GetECPbRadThick()/2.; // z
      TVirtualMC::GetMC()->Gsvolu("PBTI", "BOX", fIdTmedArr[kIdPB], dummy, 0);
      zpos = -fSampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;
      AliDebug(2,Form(" Pb tiles \n"));
      for(Int_t iz=0; iz<g->GetNECLayers(); iz++){
        par[0] = (parSCM0[0] + tanBetta*fSampleWidth*iz)/2.;
        xpos   = par[0] - xCenterSCMX;
        TVirtualMC::GetMC()->Gsposp("PBTI", ++nr, "SCMX", xpos, ypos, zpos, 0, "ONLY", par, 3) ;
        AliDebug(3,Form(" %i xpos %f zpos %f par[0] %f \n", iz+1, xpos, zpos, par[0]));
        zpos += fSampleWidth;
      } 
      AliDebug(2,Form(" Number of Pb tiles in SCMX %i \n", nr));
    } else {
      // Oct 26, 2010
      // First sheet of paper
      par[1] = parSCM0[2]/2.;                 // y 
      par[2] = g->GetTrd1BondPaperThick()/2.; // z
      par[0] = parSCM0[0]/2.;                 // x 
      TVirtualMC::GetMC()->Gsvolu("PAP1", "BOX", fIdTmedArr[kIdPAPER], par, 3);
      xpos = par[0] - xCenterSCMX;
      zpos = -parSCM0[3] + g->GetTrd1BondPaperThick()/2.;
      TVirtualMC::GetMC()->Gspos("PAP1", 1, "SCMX", xpos, ypos, zpos, 0, "ONLY");
      for(Int_t iz=0; iz<g->GetNECLayers()-1; iz++){
        nr = iz + 1;
        Double_t dz = g->GetECScintThick() + g->GetTrd1BondPaperThick() + fSampleWidth*iz;
        // PB + 2 paper sheets
        par[2] = g->GetECPbRadThick()/2. + g->GetTrd1BondPaperThick(); // z
        par[0] = (parSCM0[0] + tanBetta*dz)/2.;
        TString pa(Form("PA%2.2i",nr));
        TVirtualMC::GetMC()->Gsvolu(pa.Data(), "BOX", fIdTmedArr[kIdPAPER], par, 3);
        xpos   = par[0] - xCenterSCMX;
        zpos   = -parSCM0[3] + dz + par[2];
        TVirtualMC::GetMC()->Gspos(pa.Data(), 1, "SCMX", xpos, ypos, zpos, 0, "ONLY") ;
        // Pb
        TString pb(Form("PB%2.2i",nr));
        par[2] = g->GetECPbRadThick()/2.; // z
        TVirtualMC::GetMC()->Gsvolu(pb.Data(), "BOX", fIdTmedArr[kIdPB], par, 3);
        TVirtualMC::GetMC()->Gspos(pb.Data(), 1, pa.Data(), 0.0, 0.0, 0.0, 0, "ONLY") ;
      } 
    }
    
  } else if(g->GetNPHIdiv()==3 && g->GetNETAdiv()==3) {
    printf(" before AliEMCALv0::Trd1Tower3X3() : parSCM0");
    for(Int_t i=0; i<4; i++) printf(" %7.4f ", parSCM0[i]);
    printf("\n"); 
    Trd1Tower3X3(parSCM0);
  } else if(g->GetNPHIdiv()==1 && g->GetNETAdiv()==1) {
    // no division in SCM0
    Trd1Tower1X1(parSCM0);
  }
  
}

//______________________________________________________________________
void AliEMCALv0::CreateSmod(const char* mother)
{ 
  // 18-may-05; mother="XEN1";
  // child="SMOD" from first to 10th, "SM10" (11th and 12th)
  // "DCSM" from 13th to 18/22th (TRD1 case), "DCEXT"(18th and 19th)  adapted for DCAL, Oct-23-2012
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper();

  Double_t par[3], xpos=0., ypos=0., zpos=0., rpos=0., dphi=0., phi=0.0, phiRad=0.;
  Double_t parC[3] = {0};
  TString smName;
  Int_t tmpType = -1;

  //  ===== define Super Module from air - 14x30 module ==== ;
  AliDebug(2,Form("\n ## Super Module | fSampleWidth %5.3f ## %s \n", fSampleWidth, gn.Data()));
  par[0] = g->GetShellThickness()/2.;
  par[1] = g->GetPhiModuleSize()*g->GetNPhi()/2.; 
  par[2] = g->GetEtaModuleSize()*g->GetNEta()/2.;
  fIdRotm=0;
  Int_t nSMod = g->GetNumberOfSuperModules(); 
  Int_t nphism = nSMod/2; // 20-may-05
  if(nphism > 0) {
    dphi = g->GetEMCGeometry()->GetPhiSuperModule();
    rpos = (g->GetEnvelop(0) + g->GetEnvelop(1))/2.;
    AliDebug(2,Form(" rpos %8.2f : dphi %6.1f degree \n", rpos, dphi));
  }

  if(gn.Contains("WSUC")) {
    Int_t nr=0;
    par[0] = g->GetPhiModuleSize()*g->GetNPhi()/2.; 
    par[1] = g->GetShellThickness()/2.;
    par[2] = g->GetEtaModuleSize()*g->GetNZ()/2. + 5; 

    TVirtualMC::GetMC()->Gsvolu("SMOD", "BOX", fIdTmedArr[kIdAIR], par, 3);

    AliDebug(2,Form("SMOD in WSUC : tmed %i | dx %7.2f dy %7.2f dz %7.2f (SMOD, BOX)\n", 
		    fIdTmedArr[kIdAIR], par[0],par[1],par[2]));
    fSmodPar0 = par[0]; 
    fSmodPar1 = par[1];
    fSmodPar2 = par[2];
    nphism   =  g->GetNumberOfSuperModules();
    for(Int_t i=0; i<nphism; i++) {
      xpos = ypos = zpos = 0.0;
      fIdRotm = 0;
      TVirtualMC::GetMC()->Gspos("SMOD", 1, mother, xpos, ypos, zpos, fIdRotm, "ONLY") ;
      printf(" fIdRotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f \n", 
      fIdRotm, phi, phiRad, xpos, ypos, zpos);
      nr++;
    }
  } else {// ALICE
    AliDebug(2,Form(" par[0] %7.2f (old) \n",  par[0]));
    for(Int_t i=0; i<3; i++) par[i] = g->GetSuperModulesPar(i);
    fSmodPar0 = par[0]; 
    fSmodPar2 = par[2];
    Int_t SMOrder = -1;
    tmpType = -1;
    for (Int_t smodnum = 0; smodnum < nSMod; ++smodnum) {
      for(Int_t i=0; i<3; i++) parC[i] = par[i];
      if(g->GetSMType(smodnum) == tmpType) {
        SMOrder++;
      } else {
        tmpType = g->GetSMType(smodnum);
        SMOrder = 1;
      }

      phiRad = g->GetPhiCenterOfSMSec(smodnum); // NEED  phi= 90, 110, 130, 150, 170, 190(not center)... 
      phi    = phiRad *180./TMath::Pi();
      Double_t phiy = 90. + phi;
      Double_t phiz = 0.;

      xpos = rpos * TMath::Cos(phiRad);
      ypos = rpos * TMath::Sin(phiRad);
      zpos = fSmodPar2; // 21-sep-04
      if(        tmpType == AliEMCALGeometry::kEMCAL_Standard ) {
        smName="SMOD";
      } else if( tmpType == AliEMCALGeometry::kEMCAL_Half ) {
        smName="SM10";
        parC[1] /= 2.;
        xpos += (par[1]/2. * TMath::Sin(phiRad));
        ypos -= (par[1]/2. * TMath::Cos(phiRad));
      } else if( tmpType == AliEMCALGeometry::kEMCAL_3rd ) {
        smName="SM3rd";
        parC[1] /= 3.;
        xpos += (2.*par[1]/3. * TMath::Sin(phiRad));
        ypos -= (2.*par[1]/3. * TMath::Cos(phiRad));
      } else if( tmpType == AliEMCALGeometry::kDCAL_Standard ) {
        smName="DCSM";
        parC[2] *= 2./3.;
        zpos = fSmodPar2 + g->GetDCALInnerEdge()/2.; // 21-sep-04
      } else if( tmpType == AliEMCALGeometry::kDCAL_Ext ) {
        smName="DCEXT";
        parC[1] /= 3.;
        xpos += (2.*par[1]/3. * TMath::Sin(phiRad));
        ypos -= (2.*par[1]/3. * TMath::Cos(phiRad));
      } else AliError("Unkown SM Type!!");

      if(SMOrder == 1) {//first time, create the SM
        TVirtualMC::GetMC()->Gsvolu(smName.Data(), "BOX", fIdTmedArr[kIdAIR], parC, 3);
        AliDebug(2,Form(" Super module with name \"%s\" was created in \"box\" with: par[0] = %f, par[1] = %f, par[2] = %f\n", smName.Data(), parC[0], parC[1], parC[2]));
      }

      if( smodnum%2 == 1) {
        phiy += 180.;
        if(phiy>=360.) phiy -= 360.;
        phiz = 180.;
        zpos *= -1.;
      }
      AliMatrix(fIdRotm, 90.0, phi, 90.0, phiy, phiz, 0.0);
      TVirtualMC::GetMC()->Gspos(smName.Data(), SMOrder, mother, xpos, ypos, zpos, fIdRotm, "ONLY") ;
      AliDebug(3, Form(" %s : %2i, fIdRotm %3i phi %6.1f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f : i %i \n",            
                       smName.Data(), SMOrder, fIdRotm, phi, phiRad, xpos, ypos, zpos, smodnum));
    }
  }
  AliDebug(2,Form(" Number of Super Modules %i \n", nSMod));

  // Steel plate
  if(g->GetSteelFrontThickness() > 0.0) { // 28-mar-05
    par[0] = g->GetSteelFrontThickness()/2.;
    TVirtualMC::GetMC()->Gsvolu("STPL", "BOX", fIdTmedArr[kIdSTEEL], par, 3);
    printf("tmed %i | dx %7.2f dy %7.2f dz %7.2f (STPL) \n", fIdTmedArr[kIdSTEEL], par[0],par[1],par[2]);
    xpos = -(g->GetShellThickness() - g->GetSteelFrontThickness())/2.;
    TVirtualMC::GetMC()->Gspos("STPL", 1, "SMOD", xpos, 0.0, 0.0, 0, "ONLY") ;
  }
}

//______________________________________________________________________
void AliEMCALv0::CreateEmod(const char* mother, const char* child)
{ 
  // 17-may-05; mother="SMOD"; child="EMOD"
  // Oct 26,2010 
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper(); 
  // Module definition
  Double_t xpos=0., ypos=0., zpos=0.;
  //Double_t trd1Angle = g->GetTrd1Angle()*TMath::DegToRad();tanTrd1 = TMath::Tan(trd1Angle/2.);

  if(strcmp(mother,"SMOD")==0) {
    fParEMOD[0] = g->GetEtaModuleSize()/2.;   // dx1
    fParEMOD[1] = g->Get2Trd1Dx2()/2.;        // dx2
    fParEMOD[2] = g->GetPhiModuleSize()/2.;;  // dy
    fParEMOD[3] = g->GetLongModuleSize()/2.;  // dz
    TVirtualMC::GetMC()->Gsvolu(child, "TRD1", fIdTmedArr[kIdSTEEL], fParEMOD, 4);
  }
  
  Int_t nr=0;
  fIdRotm=0;
  // X->Z(0, 0); Y->Y(90, 90); Z->X(90, 0)
  AliEMCALShishKebabTrd1Module *mod=0; // current module

  for(Int_t iz=0; iz<g->GetNZ(); iz++) {
    Double_t  angle=90., phiOK=0;
    mod = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(iz);
    angle = mod->GetThetaInDegree();
    if(!gn.Contains("WSUC")) { // ALICE 
      AliMatrix(fIdRotm, 90.-angle,180., 90.0,90.0, angle, 0.);
      phiOK = mod->GetCenterOfModule().Phi()*180./TMath::Pi(); 
      AliDebug(4,Form(" %2i | angle | %6.3f - %6.3f = %6.3f(eta %5.3f)\n",
                      iz+1, angle, phiOK, angle-phiOK, mod->GetEtaOfCenterOfModule()));
      xpos = mod->GetPosXfromR() + g->GetSteelFrontThickness() - fSmodPar0;
      zpos = mod->GetPosZ() - fSmodPar2;
      
      Int_t iyMax = g->GetNPhi();
      if(strcmp(mother,"SM10") == 0 ) {
         iyMax /= 2;
      } else if(strcmp(mother,"SM3rd") == 0 ) {
         iyMax /= 3;
      } else if(strcmp(mother,"DCEXT") == 0 ) {
         iyMax /= 3;
      } else if(strcmp(mother,"DCSM") == 0 ) {
        if(iz < 8 ) continue;//!!!DCSM from 8th to 23th
        zpos = mod->GetPosZ() - fSmodPar2 - g->GetDCALInnerEdge()/2.;
      } else if(strcmp(mother,"SMOD") != 0 ) 
        AliError("Unknown super module Type!!");
      for(Int_t iy=0; iy<iyMax; iy++) { // flat in phi
        ypos = g->GetPhiModuleSize()*(2*iy+1 - iyMax)/2.;
        TVirtualMC::GetMC()->Gspos(child, ++nr, mother, xpos, ypos, zpos, fIdRotm, "ONLY") ;
        //
        //printf(" %2i xpos %7.2f ypos %7.2f zpos %7.2f fIdRotm %i\n", nr, xpos, ypos, zpos, fIdRotm);
        AliDebug(3,Form("%3.3i(%2.2i,%2.2i) ", nr,iy+1,iz+1));
      }
      //PH          printf("\n");
    } else { //WSUC
      if(iz == 0) AliMatrix(fIdRotm, 0.,0., 90.,0., 90.,90.); // (x')z; y'(x); z'(y)
      else        AliMatrix(fIdRotm, 90-angle,270., 90.0,0.0, angle,90.);
      phiOK = mod->GetCenterOfModule().Phi()*180./TMath::Pi(); 
      AliDebug(4,Form(" %2i | angle -phiOK | %6.3f - %6.3f = %6.3f(eta %5.3f)\n",
                      iz+1, angle, phiOK, angle-phiOK, mod->GetEtaOfCenterOfModule()));
      zpos = mod->GetPosZ()      - fSmodPar2;
      ypos = mod->GetPosXfromR() - fSmodPar1;
      //printf(" zpos %7.2f ypos %7.2f fIdRotm %i\n xpos ", zpos, xpos, fIdRotm);
      for(Int_t ix=0; ix<g->GetNPhi(); ix++) 
      { // flat in phi
        xpos = g->GetPhiModuleSize()*(2*ix+1 - g->GetNPhi())/2.;
        TVirtualMC::GetMC()->Gspos(child, ++nr, mother, xpos, ypos, zpos, fIdRotm, "ONLY") ;
        //printf(" %7.2f ", xpos);
      }
      //printf("\n");
    }
  }
  AliDebug(2,Form(" Number of modules in Super Module(%s) %i \n", mother, nr));
}

void AliEMCALv0::CreateAlFrontPlate(const char* mother, const char* child)
{
  // Oct 26,2010 : Al front plate : ALFP
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()); gn.ToUpper(); 
  Double_t trd1Angle = g->GetTrd1Angle()*TMath::DegToRad(), tanTrd1 = TMath::Tan(trd1Angle/2.);
  Double_t parALFP[5], zposALFP=0.;

  parALFP[0] = g->GetEtaModuleSize()/2. - g->GetLateralSteelStrip();    // dx1
  parALFP[1] = parALFP[0]               + tanTrd1*g->GetTrd1AlFrontThick(); // dx2
  parALFP[2] = g->GetPhiModuleSize()/2. - g->GetLateralSteelStrip();    // dy
  parALFP[3] = g->GetTrd1AlFrontThick()/2.; // dz
  TVirtualMC::GetMC()->Gsvolu(child, "TRD1", fIdTmedArr[kIdAL], parALFP, 4);
  zposALFP   = -fParEMOD[3] + g->GetTrd1AlFrontThick()/2.;
  TVirtualMC::GetMC()->Gspos (child, 1, mother, 0.0, 0.0, zposALFP, 0, "ONLY");
}

//______________________________________________________________________
void AliEMCALv0::Trd1Tower3X3(const Double_t *parSCM0)
{
  // Started Dec 8,2004 by PAI
  // Fixed Nov 13,2006
  printf(" AliEMCALv0::Trd1Tower3X3() : parSCM0");
  for(Int_t i=0; i<4; i++) printf(" %7.4f ", parSCM0[i]);
  printf("\n"); 
  // Nov 10, 2006 - different name of SCMX
  Double_t parTRAP[11], *dummy=0;
  AliEMCALGeometry * g = GetGeometry(); 
  TString gn(g->GetName()), scmx; 
  gn.ToUpper(); 
 // Division to tile size 
  AliDebug(2,Form("Trd1Tower3X3() : Divide SCM0 on y-axis %i", g->GetNETAdiv()));
  TVirtualMC::GetMC()->Gsdvn("SCMY","SCM0", g->GetNETAdiv(), 2); // y-axis
  Double_t dx1=parSCM0[0], dx2=parSCM0[1], dy=parSCM0[2], dz=parSCM0[3];
  Double_t ndiv=3., xpos=0.0;
  // should be defined once
  TVirtualMC::GetMC()->Gsvolu("PBTI", "BOX", fIdTmedArr[kIdPB], dummy, 0);
  for(Int_t ix=1; ix<=3; ix++) { // 3X3
    scmx = "SCX"; // Nov 10,2006 
    // ix=1
    parTRAP[0] = dz;
    Double_t xCentBot = 2.*dx1/3.;
    Double_t xCentTop = 2.*(dx2/4. + dx1/12.);
    parTRAP[1] = TMath::ATan2((xCentTop-xCentBot),2.*dz)*TMath::RadToDeg(); // theta
    parTRAP[2] = 0.;           // phi
    // bottom
    parTRAP[3] = dy/ndiv;      // H1
    parTRAP[4] = dx1/ndiv;     // BL1
    parTRAP[5] = parTRAP[4];   // TL1
    parTRAP[6] = 0.0;          // ALP1
    // top
    parTRAP[7] = dy/ndiv;      // H2
    parTRAP[8] = dx2/2 - dx1/6.;// BL2
    parTRAP[9] = parTRAP[8];   // TL2
    parTRAP[10]= 0.0;          // ALP2
    xpos = (xCentBot+xCentTop)/2.;

    if      (ix==3) {
      parTRAP[1] = -parTRAP[1];
      xpos = -xpos;
    } else if(ix==2) { // central part is box but we treat as trapesoid due to numbering
      parTRAP[1] = 0.;
      parTRAP[8] = dx1/ndiv;     // BL2
      parTRAP[9] = parTRAP[8];   // TL2
      xpos = 0.0;
    }
    AliDebug(2,Form(" ** TRAP ** xpos %9.3f\n", xpos));
    for(Int_t i=0; i<11; i++) AliDebug(2,Form(" par[%2.2i] %9.4f\n", i, parTRAP[i]));

    scmx += ix;
    TVirtualMC::GetMC()->Gsvolu(scmx.Data(), "TRAP", fIdTmedArr[kIdSC], parTRAP, 11);
    TVirtualMC::GetMC()->Gspos(scmx.Data(), 1, "SCMY", xpos, 0.0, 0.0, 0, "ONLY") ;

    PbInTrap(parTRAP, scmx);
  }
  AliDebug(2,"Trd1Tower3X3 - Ver. 1.0 : was tested.");
}

// 8-dec-04 by PAI
//______________________________________________________________________
void AliEMCALv0::PbInTrap(const Double_t parTRAP[11], TString n)
{
 // see for example CreateShishKebabGeometry(); just for case TRD1
  static Int_t nr=0;
  AliDebug(2,Form(" Pb tiles : nrstart %i\n", nr));
  AliEMCALGeometry * g = GetGeometry(); 

  Double_t par[3];
  Double_t xpos = 0.0, ypos = 0.0;
  Double_t zpos = -fSampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;

  Double_t coef = (parTRAP[8] -  parTRAP[4]) / (2.*parTRAP[0]);
  Double_t xCenterSCMX =  (parTRAP[4] +  parTRAP[8])/2.; // ??
  //  Double_t tan = TMath::Tan(parTRAP[1]*TMath::DegToRad());

  par[1] = parTRAP[3];              // y 
  par[2] = g->GetECPbRadThick()/2.; // z
  for(Int_t iz=0; iz<g->GetNECLayers(); iz++){
    par[0] = parTRAP[4] + coef*fSampleWidth*iz;
    xpos   = par[0] - xCenterSCMX;
    if(parTRAP[1] < 0.) xpos = -xpos;
    TVirtualMC::GetMC()->Gsposp("PBTI", ++nr, n.Data(), xpos, ypos, zpos, 0, "ONLY", par, 3) ;
    AliDebug(2,Form(" %i xpos %9.3f zpos %9.3f par[0] %9.3f |", iz+1, xpos, zpos, par[0]));
    zpos += fSampleWidth;
    if(iz%2>0) printf("\n");
  } 
  AliDebug(2,Form(" Number of Pb tiles in SCMX %i coef %9.7f \n", nr, coef));
  AliDebug(2,Form(" par[1] %9.3f  par[2] %9.3f ypos %9.3f \n", par[1], par[2], ypos)); 
  AliDebug(2,Form(" PbInTrap Ver. 1.0 : was tested."));
}

//______________________________________________________________________
void AliEMCALv0::Trd1Tower1X1(Double_t *parSCM0)
{
  // Started Nov 22,2006 by PAI
  AliDebug(1," AliEMCALv0::Trd1Tower1X1() : parSCM0");
  for(Int_t i=0; i<4; i++) printf(" %7.4f ", parSCM0[i]);
  printf("\n"); 

  // No division - keeping the same volume logic 
  // and as consequence the same abs is scheme
  AliDebug(2,"Trd1Tower1X1() : Create SCMX(SCMY) as SCM0");

  TVirtualMC::GetMC()->Gsvolu("SCMY", "TRD1", fIdTmedArr[kIdAIR], parSCM0, 4);
  TVirtualMC::GetMC()->Gspos("SCMY", 1, "SCM0", 0.0, 0.0, 0.0, 0, "ONLY");
  TVirtualMC::GetMC()->Gsvolu("SCMX", "TRD1", fIdTmedArr[kIdSC], parSCM0, 4);
  TVirtualMC::GetMC()->Gspos("SCMX", 1, "SCMY", 0.0, 0.0, 0.0, 0, "ONLY");

  // should be defined once
  Double_t *dummy=0;
  TVirtualMC::GetMC()->Gsvolu("PBTI", "BOX", fIdTmedArr[kIdPB], dummy, 0);

  PbInTrd1(parSCM0, "SCMX");

  AliDebug(1,"Trd1Tower1X1() : Ver. 0.1 : was tested.");
}

//______________________________________________________________________
void AliEMCALv0::PbInTrd1(const Double_t *parTrd1, TString n)
{
 // see PbInTrap(const Double_t parTrd1[11], TString n)
  static Int_t nr=0, ndeb=2;
  AliDebug(ndeb,Form(" Pb tiles : nrstart %i\n", nr));
  AliEMCALGeometry * g = GetGeometry(); 

  Double_t par[3];
  Double_t xpos = 0.0, ypos = 0.0;
  Double_t zpos = -fSampleWidth*g->GetNECLayers()/2. + g->GetECPbRadThick()/2.;
  Double_t coef = (parTrd1[1] -  parTrd1[0]) / (2.*parTrd1[3]);

  par[1] = parTrd1[2];              // y 
  par[2] = g->GetECPbRadThick()/2.; // z

  for(Int_t iz=0; iz<g->GetNECLayers(); iz++){
    par[0] = parTrd1[0] + coef*fSampleWidth*iz;
    TVirtualMC::GetMC()->Gsposp("PBTI", ++nr, n.Data(), xpos, ypos, zpos, 0, "ONLY", par, 3) ;
    AliDebug(2,Form(" %i xpos %9.3f zpos %9.3f par[0] %9.3f |", iz+1, xpos, zpos, par[0]));
    zpos += fSampleWidth;
    if(iz%2>0) printf("\n");
  } 
  AliDebug(ndeb,Form(" Number of Pb tiles in SCMX %i coef %9.7f ", nr, coef));
  AliDebug(ndeb,Form(" PbInTrd1 Ver. 0.1 : was tested."));
}

//______________________________________________________________________
AliEMCALShishKebabTrd1Module* AliEMCALv0::GetShishKebabModule(Int_t neta)
{ 
  // 28-oct-05
  AliEMCALShishKebabTrd1Module* trd1=0;
  if(fShishKebabModules && neta>=0 && neta<fShishKebabModules->GetSize()) {
    trd1 = (AliEMCALShishKebabTrd1Module*)fShishKebabModules->At(neta);
  }
  return trd1;
}

//_____________________________________________________________________________
void AliEMCALv0::AddAlignableVolumes() const
{
  //Add volumes which are alignable (?)
  TString ntmp(GetTitle()); // name of EMCAL geometry

  if(ntmp.Contains("WSUC")) {
    AddAlignableVolumesInWSUC();  // WSUC case
  } else {
    AddAlignableVolumesInALICE(); // ALICE case
  }
}

//______________________________________________________________________
void AliEMCALv0::AddAlignableVolumesInALICE() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be synchronized with
  // eventual changes in the geometry.
  //

  Float_t pars[] = {GetGeometry()->GetSuperModulesPar(0),GetGeometry()->GetSuperModulesPar(1),GetGeometry()->GetSuperModulesPar(2)};
  Double_t rpos = (GetGeometry()->GetEnvelop(0) + GetGeometry()->GetEnvelop(1))/2.;
  Double_t phi, phiRad, xpos, ypos, zpos;

  AliGeomManager::ELayerID idEMCAL = AliGeomManager::kEMCAL;
  Int_t modUID, modnum = 0;
  TString volpath, symname;

  AliEMCALGeometry * geom = GetGeometry();
  Int_t nSMod = geom->GetNumberOfSuperModules(); 
  TString SMPathName;
  TString SMName;
  Int_t tmpType = -1;
  Int_t SMOrder = 0;

  for (Int_t smodnum = 0; smodnum < nSMod; ++smodnum) {
    modUID = AliGeomManager::LayerToVolUID(idEMCAL,modnum++);
    if(geom->GetSMType(smodnum) == AliEMCALGeometry::kEMCAL_Standard )      { SMPathName = "SMOD";  SMName = "FullSupermodule";}
    else if(geom->GetSMType(smodnum) == AliEMCALGeometry::kEMCAL_Half )     { SMPathName = "SM10";  SMName = "HalfSupermodule";}
    else if(geom->GetSMType(smodnum) == AliEMCALGeometry::kEMCAL_3rd )      { SMPathName = "SM3rd"; SMName = "OneThrdSupermodule";}
    else if( geom->GetSMType(smodnum) == AliEMCALGeometry::kDCAL_Standard ) { SMPathName = "DCSM";  SMName = "DCALSupermodule";}
    else if( geom->GetSMType(smodnum) == AliEMCALGeometry::kDCAL_Ext )      { SMPathName = "DCEXT"; SMName = "DCALExtensionSM";}
    else AliError("Unkown SM Type!!");

    if(geom->GetSMType(smodnum) == tmpType) {
      SMOrder++;
    } else {
      tmpType = geom->GetSMType(smodnum);
      SMOrder = 1;
    }

    volpath.Form("ALIC_1/XEN1_1/%s_%d",SMPathName.Data(), SMOrder);
    symname.Form("EMCAL/%s%d",SMName.Data(), SMOrder);

    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data(),modUID))
      AliFatal(Form("AliEMCALv0::Unable to set alignable entry!!\nName: %s\t Path: %s\t ModuleID: %d\n",symname.Data(),volpath.Data(), modUID));

    // Creates the Tracking to Local transformation matrix for EMCAL
    // modules                         
    TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(modUID) ;

    phiRad = GetGeometry()->GetPhiCenterOfSM(smodnum);  //comes in radians, not degrees
    phi = phiRad*180./TMath::Pi();             //need degrees for rot. matrix
    xpos = rpos * TMath::Cos(phiRad);
    ypos = rpos * TMath::Sin(phiRad);
    zpos = pars[2];
    if( geom->GetSMType(smodnum) == AliEMCALGeometry::kEMCAL_Half ) {
      xpos += (pars[1]/2. * TMath::Sin(phiRad));  //  half SM!
      ypos -= (pars[1]/2. * TMath::Cos(phiRad));
    } else if ( geom->GetSMType(smodnum) == AliEMCALGeometry::kEMCAL_3rd || geom->GetSMType(smodnum) == AliEMCALGeometry::kDCAL_Ext ) {
      xpos += (pars[1]/3. * TMath::Sin(phiRad));  // one_third SM !
      ypos -= (pars[1]/3. * TMath::Cos(phiRad));
    } else if( geom->GetSMType(smodnum) == AliEMCALGeometry::kDCAL_Standard ) {
      zpos = pars[2]*2./3. + GetGeometry()->GetDCALInnerEdge()/2.;
    }
    
    AliDebug(3, Form("  fIdRotm %3i phi %6.13f(%5.3f) xpos %7.2f ypos %7.2f zpos %7.2f : smodnum %i \n", 
		       fIdRotm, phi, phiRad, xpos, ypos, zpos, smodnum));
   
    TGeoHMatrix *matTtoL;
    TGeoHMatrix *globMatrix = alignableEntry->GetGlobalOrig();

    if(smodnum%2 == 0) {
      // pozitive z 
      TGeoTranslation geoTran0(xpos, ypos, zpos); 
      TGeoRotation geoRot0("geoRot0", 90.0, phi, 90.0, 90.0+phi, 0.0, 0.0);
      TGeoCombiTrans mat0(geoTran0, geoRot0);
      matTtoL = new TGeoHMatrix(mat0);

      matTtoL->MultiplyLeft(&(globMatrix->Inverse()));
      alignableEntry->SetMatrix(matTtoL);

    } else {
      // negative z                                                                                 
      Double_t phiy = 90. + phi + 180.;
      if(phiy>=360.) phiy -= 360.;
      TGeoTranslation geoTran1(xpos,ypos,-zpos);
      TGeoRotation geoRot1("geoRot1", 90.0, phi, 90.0, phiy, 180.0, 0.0);
      TGeoCombiTrans mat1(geoTran1, geoRot1);
      matTtoL = new TGeoHMatrix(mat1);

      matTtoL->MultiplyLeft(&(globMatrix->Inverse()));
      alignableEntry->SetMatrix(matTtoL);

    }

  }

}

//______________________________________________________________________
void AliEMCALv0::AddAlignableVolumesInWSUC() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be synchronized with
  // eventual changes in the geometry.
  //

  TString vpstr1 = "WSUC_1/XEN1_1/SMOD_";
  TString snstr1 = "EMCAL/CosmicTestSupermodule";
  TString volpath, symname;

  // #SM is just one 
  for (Int_t smodnum=0; smodnum < 1; smodnum++) {
    symname = snstr1;
    symname += (smodnum+1);
    volpath = vpstr1;
    volpath += (smodnum+1);
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data()))
      AliFatal("Unable to set alignable entry!!");
  }
}
