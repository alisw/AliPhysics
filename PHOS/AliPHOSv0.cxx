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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.82  2006/09/27 19:55:57  kharlov
 * Alignment object with symbolic volume names are introduced
 *
 * Revision 1.81  2006/03/04 20:25:56  kharlov
 * Set geom parameters from CDB
 *
 * Revision 1.80  2005/06/17 07:39:07  hristov
 * Removing GetDebug and SetDebug from AliRun and AliModule. Using AliLog for the messages
 *
 * Revision 1.79  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// An object of this class does not produce hits nor digits
// It is the one to use if you do not want to produce outputs in TREEH or TREED
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---

#include <TBRIK.h>
#include <TFolder.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TTRD1.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>

// --- Standard library ---

#include <string.h>
#include <stdlib.h>

// --- AliRoot header files ---

#include "AliConst.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSLoader.h"
#include "AliPHOSv0.h"
#include "AliRun.h"
#include "AliLog.h"

ClassImp(AliPHOSv0)

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(const char *name, const char *title):
  AliPHOS(name,title)
{
  // ctor : title is used to identify the layout
  GetGeometry() ; 
}

//____________________________________________________________________________
void AliPHOSv0::BuildGeometry()
{
  // Build the PHOS geometry for the ROOT display
  //BEGIN_HTML
  /*
    <H2>
     PHOS in ALICE displayed by root
    </H2>
    <UL>
    <LI> All Views
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="All Views" SRC="../images/AliPHOSv0AllViews.gif"> 
    </CENTER></P></LI>
    <LI> Front View
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="Front View" SRC="../images/AliPHOSv0FrontView.gif"> 
    </CENTER></P></LI>
     <LI> 3D View 1
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="3D View 1" SRC="../images/AliPHOSv03DView1.gif"> 
    </CENTER></P></LI>
    <LI> 3D View 2
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="3D View 2" SRC="../images/AliPHOSv03DView2.gif"> 
    </CENTER></P></LI>
    </UL>
  */
  //END_HTML  
  
  this->BuildGeometryforEMC() ; 
  this->BuildGeometryforCPV() ;
  
}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforEMC(void)
{
  // Build the PHOS-EMC geometry for the ROOT display
  
  const Int_t kColorPHOS = kRed ;
  const Int_t kColorXTAL = kBlue ;
  
  Double_t const kRADDEG = 180.0 / TMath::Pi() ;
  
  AliPHOSGeometry * geom = GetGeometry() ; 
  AliPHOSEMCAGeometry * emcg = geom->GetEMCAGeometry() ;
  Float_t * boxparams = emcg->GetEMCParams() ;

  new TTRD1("OuterBox", "PHOS box", "void",boxparams[0],boxparams[1],boxparams[2], boxparams[3] );
  
  
  // Crystals Box
  
  Float_t * cribox = emcg->GetInnerThermoHalfSize() ;  
  new TBRIK( "CrystalsBox", "PHOS crystals box", "void", cribox[0], cribox[2], cribox[1] ) ;
  
  // position PHOS into ALICE
  
  Float_t r = geom->GetIPtoOuterCoverDistance() + boxparams[3] ;
  Int_t number = 988 ; 
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
  
  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  new TRotMatrix("cribox", "cribox", 90, 0, 90, 90, 0, 0);  

  for( Int_t i = 1; i <= geom->GetNModules(); i++ ) { 

    Float_t angle = geom->GetPHOSAngle(i) ;
    sprintf(rotname, "%s%d", "rot", number++) ;
    new TRotMatrix(rotname, rotname, 90, angle, 0,  0,  90,  270 + angle);

    top->cd();
    sprintf(nodename,"%s%d", "Module", i) ;    
    Float_t x =  r * TMath::Sin( angle / kRADDEG ) ;
    Float_t y = -r * TMath::Cos( angle / kRADDEG ) ;
    TNode * outerboxnode = new TNode(nodename, nodename, "OuterBox", x, y, 0, rotname ) ;
    outerboxnode->SetLineColor(kColorPHOS) ;
    fNodes->Add(outerboxnode) ;
    outerboxnode->cd() ; 

    Float_t z = -boxparams[3] - geom->GetIPtoOuterCoverDistance() + 
                 cribox[1] +  geom->GetIPtoCrystalSurface() ;
    TNode * crystalsboxnode = new TNode(nodename, nodename, "CrystalsBox", 0, 0, z) ;    
    crystalsboxnode->SetLineColor(kColorXTAL) ; 
    fNodes->Add(crystalsboxnode) ; 
  }

  delete[] rotname ;  
  delete[] nodename ;
}


//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforCPV(void)
{
  //  Build the PHOS-CPV geometry for the ROOT display
  //  Author: Yuri Kharlov 11 September 2000
  //
  //BEGIN_HTML
  /*
    <H2>
    CPV displayed by root
    </H2>
    <table width=700>

    <tr>
         <td>CPV perspective view</td>
         <td>CPV front view      </td>
    </tr>

    <tr>
         <td> <img height=300 width=290 src="../images/CPVRootPersp.gif"> </td>
         <td> <img height=300 width=290 src="../images/CPVRootFront.gif"> </td>
    </tr>

    </table>

  */
  //END_HTML  

  const Double_t kRADDEG         = 180.0 / TMath::Pi() ;
  const Int_t    kColorCPV       = kGreen ;
  const Int_t    kColorFrame     = kYellow ;
  const Int_t    kColorGassiplex = kRed;
  const Int_t    kColorPCB       = kCyan;

  AliPHOSGeometry * geom = GetGeometry() ; 

  // Box for a full PHOS module

  new TBRIK ("CPVBox", "CPV box", "void",                   geom->GetCPVBoxSize(0)/2,
                                                            geom->GetCPVBoxSize(1)/2,
	                                                    geom->GetCPVBoxSize(2)/2 );
  new TBRIK ("CPVFrameLR", "CPV frame Left-Right", "void",  geom->GetCPVFrameSize(0)/2,
                                                            geom->GetCPVFrameSize(1)/2,
	                                                    geom->GetCPVBoxSize(2)/2 );
  new TBRIK ("CPVFrameUD", "CPV frame Up-Down",    "void",  geom->GetCPVBoxSize(0)/2 - geom->GetCPVFrameSize(0),
                                                            geom->GetCPVFrameSize(1)/2,
	                                                    geom->GetCPVFrameSize(2)/2);
  new TBRIK ("CPVPCB",    "CPV PCB",               "void",  geom->GetCPVActiveSize(0)/2,
                                                            geom->GetCPVTextoliteThickness()/2,
	                                                    geom->GetCPVActiveSize(1)/2);
  new TBRIK ("CPVGassiplex", "CPV Gassiplex PCB",  "void",  geom->GetGassiplexChipSize(0)/2,
                                                            geom->GetGassiplexChipSize(1)/2,
	                                                    geom->GetGassiplexChipSize(2)/2);

  // position CPV into ALICE

  char * nodename = new char[25] ;
  char * rotname  = new char[25] ;
  
  Float_t r = geom->GetIPtoCPVDistance() + geom->GetCPVBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;

  Int_t lastModule = 0 ;
  lastModule = geom->GetNModules();
  
  for( Int_t i = 1; i <= lastModule; i++ ) { // the number of PHOS modules
    
    // One CPV module
    
    Float_t angle = geom->GetPHOSAngle(i) ;
    sprintf(rotname, "%s%d", "rotg", number+i) ;
    new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
    top->cd();
    sprintf(nodename, "%s%d", "CPVModule", i) ;    
    Float_t x =  r * TMath::Sin( angle / kRADDEG ) ;
    Float_t y = -r * TMath::Cos( angle / kRADDEG ) ;
    Float_t z;
    TNode * cpvBoxNode = new TNode(nodename , nodename ,"CPVBox", x, y, 0, rotname ) ;
    cpvBoxNode->SetLineColor(kColorCPV) ;
    fNodes->Add(cpvBoxNode) ;
    cpvBoxNode->cd() ;

    // inside each CPV box:

    // Frame around CPV
    Int_t j;
    for (j=0; j<=1; j++) {
      sprintf(nodename, "CPVModule%d Frame%d", i, j+1) ;
      x = TMath::Sign(1,2*j-1) * (geom->GetCPVBoxSize(0) - geom->GetCPVFrameSize(0)) / 2;
      TNode * cpvFrameNode = new TNode(nodename , nodename ,"CPVFrameLR", x, 0, 0) ;
      cpvFrameNode->SetLineColor(kColorFrame) ;
      fNodes->Add(cpvFrameNode) ;

      sprintf(nodename, "CPVModule%d Frame%d", i, j+3) ;
      z = TMath::Sign(1,2*j-1) * (geom->GetCPVBoxSize(2) - geom->GetCPVFrameSize(2)) / 2;
      cpvFrameNode = new TNode(nodename , nodename ,"CPVFrameUD", 0, 0, z) ;
      cpvFrameNode->SetLineColor(kColorFrame) ;
      fNodes->Add(cpvFrameNode) ;
    }

    // 4 printed circuit boards
    for (j=0; j<4; j++) {
      sprintf(nodename, "CPVModule%d PCB%d", i, j+1) ;
      y = geom->GetCPVFrameSize(1) / 2 - geom->GetFTPosition(j) + geom->GetCPVTextoliteThickness()/2;
      TNode * cpvPCBNode = new TNode(nodename , nodename ,"CPVPCB", 0, y, 0) ;
      cpvPCBNode->SetLineColor(kColorPCB) ;
      fNodes->Add(cpvPCBNode) ;
    }

    // Gassiplex chips
    Float_t xStep = geom->GetCPVActiveSize(0) / (geom->GetNumberOfCPVChipsPhi() + 1);
    Float_t zStep = geom->GetCPVActiveSize(1) / (geom->GetNumberOfCPVChipsZ()   + 1);
    y = geom->GetCPVFrameSize(1)/2           - geom->GetFTPosition(0) +
        geom->GetCPVTextoliteThickness() / 2 + geom->GetGassiplexChipSize(1) / 2 + 0.1;
    for (Int_t ix=0; ix<geom->GetNumberOfCPVChipsPhi(); ix++) {
      x = xStep * (ix+1) - geom->GetCPVActiveSize(0)/2;
      for (Int_t iz=0; iz<geom->GetNumberOfCPVChipsZ(); iz++) {
	z = zStep * (iz+1) - geom->GetCPVActiveSize(1)/2;
	sprintf(nodename, "CPVModule%d Chip(%dx%d)", i, ix+1,iz+1) ;
	TNode * cpvGassiplexNode = new TNode(nodename , nodename ,"CPVGassiplex", x, y, z) ;
	cpvGassiplexNode->SetLineColor(kColorGassiplex) ;
	fNodes->Add(cpvGassiplexNode) ;
      }
    }

  } // PHOS modules
 
  delete[] rotname ;  
  delete[] nodename ; 
}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometry()
{
  // Create the PHOS geometry for Geant

  AliPHOSv0 *phostmp = dynamic_cast<AliPHOSv0*>(gAlice->GetModule("PHOS")) ;

  if ( phostmp == NULL ) {
    
    fprintf(stderr, "PHOS detector not found!\n") ;
    return;
    
  }

  AliPHOSGeometry * geom = GetGeometry() ; 

  // Get pointer to the array containing media indeces
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  // Create a PHOS module.
  
  gMC->Gsvolu("PHOS", "TRD1", idtmed[798], geom->GetPHOSParams(), 4) ;        
  
  this->CreateGeometryforEMC() ; 

  this->CreateGeometryforCPV() ;
  
  this->CreateGeometryforSupport() ; 
  
  // --- Position  PHOS mdules in ALICE setup ---
  
  Int_t idrotm[99] ;
  Int_t iXYZ,iAngle;
  for (Int_t iModule = 0; iModule < geom->GetNModules(); iModule++ ) {
    
    Float_t angle[3][2];
    for (iXYZ=0; iXYZ<3; iXYZ++)
      for (iAngle=0; iAngle<2; iAngle++)
	angle[iXYZ][iAngle] = geom->GetModuleAngle(iModule,iXYZ, iAngle);
    AliMatrix(idrotm[iModule],
	      angle[0][0],angle[0][1],
	      angle[1][0],angle[1][1],
	      angle[2][0],angle[2][1]) ;
    
    Float_t pos[3];
    for (iXYZ=0; iXYZ<3; iXYZ++)
      pos[iXYZ] = geom->GetModuleCenter(iModule,iXYZ);
    gMC->Gspos("PHOS", iModule+1, "ALIC", pos[0], pos[1], pos[2],
	       idrotm[iModule], "ONLY") ;
  }

}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforEMC()
{
  // Create the PHOS-EMC geometry for GEANT
  // Author: Dmitri Peressounko August 2001
  // The used coordinate system: 
  //   1. in Module: X along longer side, Y out of beam, Z along shorter side (along beam)
  //   2. In Strip the same: X along longer side, Y out of beam, Z along shorter side (along beam)


    //BEGIN_HTML
  /*
    <H2>
    Geant3 geometry tree of PHOS-EMC in ALICE
    </H2>
    <P><CENTER>
    <IMG Align=BOTTOM ALT="EMC geant tree" SRC="../images/EMCinAlice.gif"> 
    </CENTER><P>
  */
  //END_HTML  
  
  // Get pointer to the array containing media indexes
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  AliPHOSGeometry * geom = GetGeometry() ; 
  AliPHOSEMCAGeometry * emcg = geom->GetEMCAGeometry() ;

  // ======= Define the strip ===============

  gMC->Gsvolu("PSTR", "BOX ", idtmed[716], emcg->GetStripHalfSize(), 3) ;  //Made of stell
   
      // --- define air volume (cell of the honeycomb)
      gMC->Gsvolu("PCEL", "BOX ", idtmed[798], emcg->GetAirCellHalfSize(), 3);

      // --- define wrapped crystal and put it into AirCell

      gMC->Gsvolu("PWRA", "BOX ", idtmed[702], emcg->GetWrappedHalfSize(), 3);
      Float_t * pin = emcg->GetAPDHalfSize() ; 
      Float_t * preamp = emcg->GetPreampHalfSize() ;
      Float_t y = (emcg->GetAirGapLed()-2*pin[1]-2*preamp[1])/2;
      gMC->Gspos("PWRA", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY") ;
    
      // --- Define crystall and put it into wrapped crystall ---
      gMC->Gsvolu("PXTL", "BOX ", idtmed[699], emcg->GetCrystalHalfSize(), 3) ;
      gMC->Gspos("PXTL", 1, "PWRA", 0.0, 0.0, 0.0, 0, "ONLY") ;
      
      // --- define APD/PIN preamp and put it into AirCell
 
      gMC->Gsvolu("PPIN", "BOX ", idtmed[705], emcg->GetAPDHalfSize(), 3) ;
      Float_t * crystal = emcg->GetCrystalHalfSize() ;
      y = crystal[1] + emcg->GetAirGapLed() /2 - preamp[1]; 
      gMC->Gspos("PPIN", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY") ;

      gMC->Gsvolu("PREA", "BOX ", idtmed[711], emcg->GetPreampHalfSize(), 3) ;   // Here I assumed preamp
                                                                                 // as a printed Circuit
      y = crystal[1] + emcg->GetAirGapLed() /2 + pin[1]  ;                  // May it should be changed
      gMC->Gspos("PREA", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY") ;                    // to ceramics?
   

      // --- Fill strip with wrapped cristalls in Air Cells

      Float_t* splate = emcg->GetSupportPlateHalfSize();  
      y = -splate[1] ;
      Float_t* acel = emcg->GetAirCellHalfSize() ;
      Int_t icel ;
      for(icel = 1; icel <= emcg->GetNCellsInStrip(); icel++){
	Float_t x = (2*icel - 1 - emcg->GetNCellsInStrip())* acel[0] ;
	gMC->Gspos("PCEL", icel, "PSTR", x, y, 0.0, 0, "ONLY") ;
      }

      // --- define the support plate, hole in it and position it in strip ----
      gMC->Gsvolu("PSUP", "BOX ", idtmed[701], emcg->GetSupportPlateHalfSize(), 3) ;

      gMC->Gsvolu("PSHO", "BOX ", idtmed[798], emcg->GetSupportPlateInHalfSize(), 3) ;
      Float_t z = emcg->GetSupportPlateThickness()/2 ;
      gMC->Gspos("PSHO", 1, "PSUP", 0.0, 0.0, z, 0, "ONLY") ;

      y = acel[1] ;
      gMC->Gspos("PSUP", 1, "PSTR", 0.0, y, 0.0, 0, "ONLY") ;


    // ========== Fill module with strips and put them into inner thermoinsulation=============
      gMC->Gsvolu("PTII", "BOX ", idtmed[706], emcg->GetInnerThermoHalfSize(), 3) ;     

      Float_t * inthermo = emcg->GetInnerThermoHalfSize() ;
      Float_t * strip = emcg->GetStripHalfSize() ;
      y = inthermo[1] - strip[1] ;
      Int_t irow;
      Int_t nr = 1 ;
      Int_t icol ;

      for(irow = 0; irow < emcg->GetNStripX(); irow ++){
	Float_t x = (2*irow + 1 - emcg->GetNStripX())* strip[0] ;
	for(icol = 0; icol < emcg->GetNStripZ(); icol ++){
	  z = (2*icol + 1 - emcg->GetNStripZ()) * strip[2] ;
	  gMC->Gspos("PSTR", nr, "PTII", x, y, z, 0, "ONLY") ;
	  nr++ ;
	}
      }
	  

   // ------- define the air gap between thermoinsulation and cooler
      gMC->Gsvolu("PAGA", "BOX ", idtmed[798], emcg->GetAirGapHalfSize(), 3) ;   
      Float_t * agap = emcg->GetAirGapHalfSize() ;
      y = agap[1] - inthermo[1]  ;
      
      gMC->Gspos("PTII", 1, "PAGA", 0.0, y, 0.0, 0, "ONLY") ;



   // ------- define the Al passive cooler 
      gMC->Gsvolu("PCOR", "BOX ", idtmed[701], emcg->GetCoolerHalfSize(), 3) ;   
      Float_t * cooler = emcg->GetCoolerHalfSize() ;
      y = cooler[1] - agap[1]  ;
      
      gMC->Gspos("PAGA", 1, "PCOR", 0.0, y, 0.0, 0, "ONLY") ;

   // ------- define the outer thermoinsulating cover
      gMC->Gsvolu("PTIO", "TRD1", idtmed[706], emcg->GetOuterThermoParams(), 4) ;        
      Float_t * outparams = emcg->GetOuterThermoParams() ; 

      Int_t idrotm[99] ;
      AliMatrix(idrotm[1], 90.0, 0.0, 0.0, 0.0, 90.0, 270.0) ;
      // Frame in outer thermoinsulation and so on: z out of beam, y along beam, x across beam
 
      z = outparams[3] - cooler[1] ;
      gMC->Gspos("PCOR", 1, "PTIO", 0., 0.0, z, idrotm[1], "ONLY") ;
       
  // -------- Define the outer Aluminium cover -----
      gMC->Gsvolu("PCOL", "TRD1", idtmed[701], emcg->GetAlCoverParams(), 4) ;        
      Float_t * covparams = emcg->GetAlCoverParams() ; 
      z = covparams[3] - outparams[3] ;
      gMC->Gspos("PTIO", 1, "PCOL", 0., 0.0, z, 0, "ONLY") ;

 // --------- Define front fiberglass cover -----------
      gMC->Gsvolu("PFGC", "BOX ", idtmed[717], emcg->GetFiberGlassHalfSize(), 3) ;  
      z = - outparams[3] ;
      gMC->Gspos("PFGC", 1, "PCOL", 0., 0.0, z, 0, "ONLY") ;

 //=============This is all with cold section==============


      //------ Warm Section --------------
      gMC->Gsvolu("PWAR", "BOX ", idtmed[701], emcg->GetWarmAlCoverHalfSize(), 3) ; 
      Float_t * warmcov = emcg->GetWarmAlCoverHalfSize() ;

      // --- Define the outer thermoinsulation ---
      gMC->Gsvolu("PWTI", "BOX ", idtmed[706], emcg->GetWarmThermoHalfSize(), 3) ; 
      Float_t * warmthermo = emcg->GetWarmThermoHalfSize() ;
      z = -warmcov[2] + warmthermo[2] ;

      gMC->Gspos("PWTI", 1, "PWAR", 0., 0.0, z, 0, "ONLY") ;     

      // --- Define cables area and put in it T-supports ---- 
      gMC->Gsvolu("PCA1", "BOX ", idtmed[718], emcg->GetTCables1HalfSize(), 3) ; 
      Float_t * cbox = emcg->GetTCables1HalfSize() ;

      gMC->Gsvolu("PBE1", "BOX ", idtmed[701], emcg->GetTSupport1HalfSize(), 3) ;
      Float_t * beams = emcg->GetTSupport1HalfSize() ;
      Int_t isup ;
      for(isup = 0; isup < emcg->GetNTSuppots(); isup++){
	Float_t x = -cbox[0] + beams[0] + (2*beams[0]+emcg->GetTSupportDist())*isup ;
	gMC->Gspos("PBE1", isup, "PCA1", x, 0.0, 0.0, 0, "ONLY") ;
      }

      z = -warmthermo[2] + cbox[2] ;
      gMC->Gspos("PCA1", 1, "PWTI", 0.0, 0.0, z, 0, "ONLY") ;     

      gMC->Gsvolu("PCA2", "BOX ", idtmed[718], emcg->GetTCables2HalfSize(), 3) ; 
      Float_t * cbox2 = emcg->GetTCables2HalfSize() ;

      gMC->Gsvolu("PBE2", "BOX ", idtmed[701], emcg->GetTSupport2HalfSize(), 3) ;
      for(isup = 0; isup < emcg->GetNTSuppots(); isup++){
	Float_t x = -cbox[0] + beams[0] + (2*beams[0]+emcg->GetTSupportDist())*isup ;
	gMC->Gspos("PBE2", isup, "PCA2", x, 0.0, 0.0, 0, "ONLY") ;
      }

      z = -warmthermo[2] + 2*cbox[2] + cbox2[2];
      gMC->Gspos("PCA2", 1, "PWTI", 0.0, 0.0, z, 0, "ONLY") ;     


  // --- Define frame ---
      gMC->Gsvolu("PFRX", "BOX ", idtmed[716], emcg->GetFrameXHalfSize(), 3) ; 
      Float_t * posit = emcg->GetFrameXPosition() ;
      gMC->Gspos("PFRX", 1, "PWTI", posit[0],  posit[1], posit[2], 0, "ONLY") ;
      gMC->Gspos("PFRX", 2, "PWTI", posit[0], -posit[1], posit[2], 0, "ONLY") ;

      gMC->Gsvolu("PFRZ", "BOX ", idtmed[716], emcg->GetFrameZHalfSize(), 3) ; 
      posit = emcg->GetFrameZPosition() ;
      gMC->Gspos("PFRZ", 1, "PWTI", posit[0], posit[1],  posit[2], 0, "ONLY") ;
      gMC->Gspos("PFRZ", 2, "PWTI", -posit[0], posit[1], posit[2], 0, "ONLY") ;

 // --- Define Fiber Glass support ---
      gMC->Gsvolu("PFG1", "BOX ", idtmed[717], emcg->GetFGupXHalfSize(), 3) ; 
      posit = emcg->GetFGupXPosition() ;
      gMC->Gspos("PFG1", 1, "PWTI", posit[0],  posit[1], posit[2], 0, "ONLY") ;
      gMC->Gspos("PFG1", 2, "PWTI", posit[0], -posit[1], posit[2], 0, "ONLY") ;

      gMC->Gsvolu("PFG2", "BOX ", idtmed[717], emcg->GetFGupZHalfSize(), 3) ; 
      posit = emcg->GetFGupZPosition() ;
      gMC->Gspos("PFG2", 1, "PWTI",  posit[0], posit[1], posit[2], 0, "ONLY") ;
      gMC->Gspos("PFG2", 2, "PWTI", -posit[0], posit[1], posit[2], 0, "ONLY") ;

      gMC->Gsvolu("PFG3", "BOX ", idtmed[717], emcg->GetFGlowXHalfSize(), 3) ; 
      posit = emcg->GetFGlowXPosition() ;
      gMC->Gspos("PFG3", 1, "PWTI", posit[0],  posit[1], posit[2], 0, "ONLY") ;
      gMC->Gspos("PFG3", 2, "PWTI", posit[0], -posit[1], posit[2], 0, "ONLY") ;

      gMC->Gsvolu("PFG4", "BOX ", idtmed[717], emcg->GetFGlowZHalfSize(), 3) ; 
      posit = emcg->GetFGlowZPosition() ;
      gMC->Gspos("PFG4", 1, "PWTI",  posit[0], posit[1], posit[2], 0, "ONLY") ;
      gMC->Gspos("PFG4", 2, "PWTI", -posit[0], posit[1], posit[2], 0, "ONLY") ;

      // --- Define Air Gap for FEE electronics ----- 

      gMC->Gsvolu("PAFE", "BOX ", idtmed[798], emcg->GetFEEAirHalfSize(), 3) ; 
      posit = emcg->GetFEEAirPosition() ;
      gMC->Gspos("PAFE", 1, "PWTI",  posit[0], posit[1], posit[2], 0, "ONLY") ;

      // Define the EMC module volume and combine Cool and Warm sections

      gMC->Gsvolu("PEMC", "TRD1", idtmed[798], emcg->GetEMCParams(), 4) ;        

      z =  - warmcov[2] ;
      gMC->Gspos("PCOL", 1, "PEMC",  0., 0., z, 0, "ONLY") ;
      z = covparams[3] ;
      gMC->Gspos("PWAR", 1, "PEMC",  0., 0., z, 0, "ONLY") ;


      // Put created EMC geometry into PHOS volume
      
      z = geom->GetCPVBoxSize(1) / 2. ;
      gMC->Gspos("PEMC", 1, "PHOS", 0., 0., z, 0, "ONLY") ; 
            
}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforCPV()
{
  // Create the PHOS-CPV geometry for GEANT
  // Author: Yuri Kharlov 11 September 2000
  //BEGIN_HTML
  /*
    <H2>
    Geant3 geometry of PHOS-CPV in ALICE
    </H2>
    <table width=700>

    <tr>
         <td>CPV perspective view</td>
         <td>CPV front view      </td>
    </tr>

    <tr>
         <td> <img height=300 width=290 src="../images/CPVallPersp.gif"> </td>
         <td> <img height=300 width=290 src="../images/CPVallFront.gif"> </td>
    </tr>

    <tr>
         <td>One CPV module, perspective view                            </td>
         <td>One CPV module, front view (extended in vertical direction) </td>
    </tr>

    <tr>
         <td><img height=300 width=290 src="../images/CPVmodulePers.gif"></td>
         <td><img height=300 width=290 src="../images/CPVmoduleSide.gif"></td>
    </tr>

    </table>

    <H2>
    Geant3 geometry tree of PHOS-CPV in ALICE
    </H2>
    <center>
    <img height=300 width=290 src="../images/CPVtree.gif">
    </center>
  */
  //END_HTML  

  Float_t par[3], x,y,z;

  // Get pointer to the array containing media indexes
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  AliPHOSGeometry * geom = GetGeometry() ; 

  // The box containing all CPV for one PHOS module filled with air 
  par[0] = geom->GetCPVBoxSize(0) / 2.0 ;  
  par[1] = geom->GetCPVBoxSize(1) / 2.0 ; 
  par[2] = geom->GetCPVBoxSize(2) / 2.0 ;
  gMC->Gsvolu("PCPV", "BOX ", idtmed[798], par, 3) ;

  Float_t * emcParams = geom->GetEMCAGeometry()->GetEMCParams() ;
  z = - emcParams[3] ;
  Int_t rotm ;
  AliMatrix(rotm, 90.,0., 0., 0., 90., 90.) ;

  gMC->Gspos("PCPV", 1, "PHOS", 0.0, 0.0, z, rotm, "ONLY") ; 
  
  // Gassiplex board
  
  par[0] = geom->GetGassiplexChipSize(0)/2.;
  par[1] = geom->GetGassiplexChipSize(1)/2.;
  par[2] = geom->GetGassiplexChipSize(2)/2.;
  gMC->Gsvolu("PCPC","BOX ",idtmed[707],par,3);
  
  // Cu+Ni foil covers Gassiplex board

  par[1] = geom->GetCPVCuNiFoilThickness()/2;
  gMC->Gsvolu("PCPD","BOX ",idtmed[710],par,3);
  y      = -(geom->GetGassiplexChipSize(1)/2 - par[1]);
  gMC->Gspos("PCPD",1,"PCPC",0,y,0,0,"ONLY");

  // Position of the chip inside CPV

  Float_t xStep = geom->GetCPVActiveSize(0) / (geom->GetNumberOfCPVChipsPhi() + 1);
  Float_t zStep = geom->GetCPVActiveSize(1) / (geom->GetNumberOfCPVChipsZ()   + 1);
  Int_t   copy  = 0;
  y = geom->GetCPVFrameSize(1)/2           - geom->GetFTPosition(0) +
    geom->GetCPVTextoliteThickness() / 2 + geom->GetGassiplexChipSize(1) / 2 + 0.1;
  for (Int_t ix=0; ix<geom->GetNumberOfCPVChipsPhi(); ix++) {
    x = xStep * (ix+1) - geom->GetCPVActiveSize(0)/2;
    for (Int_t iz=0; iz<geom->GetNumberOfCPVChipsZ(); iz++) {
      copy++;
      z = zStep * (iz+1) - geom->GetCPVActiveSize(1)/2;
      gMC->Gspos("PCPC",copy,"PCPV",x,y,z,0,"ONLY");
    }
  }

  // Foiled textolite (1 mm of textolite + 50 mkm of Cu + 6 mkm of Ni)
  
  par[0] = geom->GetCPVActiveSize(0)        / 2;
  par[1] = geom->GetCPVTextoliteThickness() / 2;
  par[2] = geom->GetCPVActiveSize(1)        / 2;
  gMC->Gsvolu("PCPF","BOX ",idtmed[707],par,3);

  // Argon gas volume

  par[1] = (geom->GetFTPosition(2) - geom->GetFTPosition(1) - geom->GetCPVTextoliteThickness()) / 2;
  gMC->Gsvolu("PCPG","BOX ",idtmed[715],par,3);

  for (Int_t i=0; i<4; i++) {
    y = geom->GetCPVFrameSize(1) / 2 - geom->GetFTPosition(i) + geom->GetCPVTextoliteThickness()/2;
    gMC->Gspos("PCPF",i+1,"PCPV",0,y,0,0,"ONLY");
    if(i==1){
      y-= (geom->GetFTPosition(2) - geom->GetFTPosition(1)) / 2;
      gMC->Gspos("PCPG",1,"PCPV ",0,y,0,0,"ONLY");
    }
  }

  // Dummy sensitive plane in the middle of argone gas volume

  par[1]=0.001;
  gMC->Gsvolu("PCPQ","BOX ",idtmed[715],par,3);
  gMC->Gspos ("PCPQ",1,"PCPG",0,0,0,0,"ONLY");

  // Cu+Ni foil covers textolite

  par[1] = geom->GetCPVCuNiFoilThickness() / 2;
  gMC->Gsvolu("PCP1","BOX ",idtmed[710],par,3);
  y = geom->GetCPVTextoliteThickness()/2 - par[1];
  gMC->Gspos ("PCP1",1,"PCPF",0,y,0,0,"ONLY");

  // Aluminum frame around CPV

  par[0] = geom->GetCPVFrameSize(0)/2;
  par[1] = geom->GetCPVFrameSize(1)/2;
  par[2] = geom->GetCPVBoxSize(2)  /2;
  gMC->Gsvolu("PCF1","BOX ",idtmed[701],par,3);

  par[0] = geom->GetCPVBoxSize(0)/2 - geom->GetCPVFrameSize(0);
  par[1] = geom->GetCPVFrameSize(1)/2;
  par[2] = geom->GetCPVFrameSize(2)/2;
  gMC->Gsvolu("PCF2","BOX ",idtmed[701],par,3);

  for (Int_t j=0; j<=1; j++) {
    x = TMath::Sign(1,2*j-1) * (geom->GetCPVBoxSize(0) - geom->GetCPVFrameSize(0)) / 2;
    gMC->Gspos("PCF1",j+1,"PCPV", x,0,0,0,"ONLY");
    z = TMath::Sign(1,2*j-1) * (geom->GetCPVBoxSize(2) - geom->GetCPVFrameSize(2)) / 2;
    gMC->Gspos("PCF2",j+1,"PCPV",0, 0,z,0,"ONLY");
  }

}


//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforSupport()
{
  // Create the PHOS' support geometry for GEANT
    //BEGIN_HTML
  /*
    <H2>
    Geant3 geometry of the PHOS's support
    </H2>
    <P><CENTER>
    <IMG Align=BOTTOM ALT="EMC geant tree" SRC="../images/PHOS_support.gif"> 
    </CENTER><P>
  */
  //END_HTML  
  
  Float_t par[5], x0,y0,z0 ; 
  Int_t   i,j,copy;

  // Get pointer to the array containing media indexes
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  AliPHOSGeometry * geom = GetGeometry() ; 

  // --- Dummy box containing two rails on which PHOS support moves
  // --- Put these rails to the bottom of the L3 magnet

  par[0] =  geom->GetRailRoadSize(0) / 2.0 ;
  par[1] =  geom->GetRailRoadSize(1) / 2.0 ;
  par[2] =  geom->GetRailRoadSize(2) / 2.0 ;
  gMC->Gsvolu("PRRD", "BOX ", idtmed[798], par, 3) ;

  y0     = -(geom->GetRailsDistanceFromIP() - geom->GetRailRoadSize(1) / 2.0) ;
  gMC->Gspos("PRRD", 1, "ALIC", 0.0, y0, 0.0, 0, "ONLY") ; 

  // --- Dummy box containing one rail

  par[0] =  geom->GetRailOuterSize(0) / 2.0 ;
  par[1] =  geom->GetRailOuterSize(1) / 2.0 ;
  par[2] =  geom->GetRailOuterSize(2) / 2.0 ;
  gMC->Gsvolu("PRAI", "BOX ", idtmed[798], par, 3) ;

  for (i=0; i<2; i++) {
    x0     = (2*i-1) * geom->GetDistanceBetwRails()  / 2.0 ;
    gMC->Gspos("PRAI", i, "PRRD", x0, 0.0, 0.0, 0, "ONLY") ; 
  }

  // --- Upper and bottom steel parts of the rail

  par[0] =  geom->GetRailPart1(0) / 2.0 ;
  par[1] =  geom->GetRailPart1(1) / 2.0 ;
  par[2] =  geom->GetRailPart1(2) / 2.0 ;
  gMC->Gsvolu("PRP1", "BOX ", idtmed[716], par, 3) ;

  y0     = - (geom->GetRailOuterSize(1) - geom->GetRailPart1(1))  / 2.0 ;
  gMC->Gspos("PRP1", 1, "PRAI", 0.0, y0, 0.0, 0, "ONLY") ;
  y0     =   (geom->GetRailOuterSize(1) - geom->GetRailPart1(1))  / 2.0 - geom->GetRailPart3(1);
  gMC->Gspos("PRP1", 2, "PRAI", 0.0, y0, 0.0, 0, "ONLY") ;

  // --- The middle vertical steel parts of the rail

  par[0] =  geom->GetRailPart2(0) / 2.0 ;
  par[1] =  geom->GetRailPart2(1) / 2.0 ;
  par[2] =  geom->GetRailPart2(2) / 2.0 ;
  gMC->Gsvolu("PRP2", "BOX ", idtmed[716], par, 3) ;

  y0     =   - geom->GetRailPart3(1) / 2.0 ;
  gMC->Gspos("PRP2", 1, "PRAI", 0.0, y0, 0.0, 0, "ONLY") ; 

  // --- The most upper steel parts of the rail

  par[0] =  geom->GetRailPart3(0) / 2.0 ;
  par[1] =  geom->GetRailPart3(1) / 2.0 ;
  par[2] =  geom->GetRailPart3(2) / 2.0 ;
  gMC->Gsvolu("PRP3", "BOX ", idtmed[716], par, 3) ;

  y0     =   (geom->GetRailOuterSize(1) - geom->GetRailPart3(1))  / 2.0 ;
  gMC->Gspos("PRP3", 1, "PRAI", 0.0, y0, 0.0, 0, "ONLY") ; 

  // --- The wall of the cradle
  // --- The wall is empty: steel thin walls and air inside

  par[1] =  TMath::Sqrt(TMath::Power((geom->GetIPtoCPVDistance() + geom->GetOuterBoxSize(3)),2) +
			TMath::Power((geom->GetOuterBoxSize(1)/2),2))+10. ;
  par[0] =  par[1] - geom->GetCradleWall(1) ;
  par[2] =  geom->GetCradleWall(2) / 2.0 ;
  par[3] =  geom->GetCradleWall(3) ;
  par[4] =  geom->GetCradleWall(4) ;
  gMC->Gsvolu("PCRA", "TUBS", idtmed[716], par, 5) ;

  par[0] +=  geom->GetCradleWallThickness() ;
  par[1] -=  geom->GetCradleWallThickness() ;
  par[2] -=  geom->GetCradleWallThickness() ;
  gMC->Gsvolu("PCRE", "TUBS", idtmed[798], par, 5) ;
  gMC->Gspos ("PCRE", 1, "PCRA", 0.0, 0.0, 0.0, 0, "ONLY") ; 

  for (i=0; i<2; i++) {
    z0 = (2*i-1) * (geom->GetOuterBoxSize(2) + geom->GetCradleWall(2) )/ 2.0  ;
        gMC->Gspos("PCRA", i, "ALIC", 0.0, 0.0, z0, 0, "ONLY") ; 
  }

  // --- The "wheels" of the cradle
  
  par[0] = geom->GetCradleWheel(0) / 2;
  par[1] = geom->GetCradleWheel(1) / 2;
  par[2] = geom->GetCradleWheel(2) / 2;
  gMC->Gsvolu("PWHE", "BOX ", idtmed[716], par, 3) ;

  y0 = -(geom->GetRailsDistanceFromIP() - geom->GetRailRoadSize(1) -
	 geom->GetCradleWheel(1)/2) ;
  for (i=0; i<2; i++) {
    z0 = (2*i-1) * ((geom->GetOuterBoxSize(2) + geom->GetCradleWheel(2))/ 2.0 +
                    geom->GetCradleWall(2));
    for (j=0; j<2; j++) {
      copy = 2*i + j;
      x0 = (2*j-1) * geom->GetDistanceBetwRails()  / 2.0 ;
      gMC->Gspos("PWHE", copy, "ALIC", x0, y0, z0, 0, "ONLY") ; 
    }
  }

}

//_____________________________________________________________________________
void AliPHOSv0::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry
  // Alignable volumes are:
  // 1) PHOS modules as a whole
  // 2) Cradle
  // 3) Cradle wheels
  // 4) Strip units (group of 2x8 crystals)

  TString volpath, symname;

  // Alignable modules
  // Volume path /ALIC_1/PHOS_<i> => symbolic name /PHOS/Module<i>, <i>=1,2,3,4,5

  TString physModulePath="/ALIC_1/PHOS_";
  TString symbModuleName="PHOS/Module";
  Int_t nModules = GetGeometry()->GetNModules();

  for(Int_t iModule=1; iModule<=nModules; iModule++){
    volpath = physModulePath;
    volpath += iModule;
    symname = symbModuleName;
    symname += iModule;
    gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data());
  }

  // Alignable cradle walls
  // Volume path /ALIC_1/PCRA_<i> => symbolic name /PHOS/Cradle<i>, <i>=0,1

  TString physCradlePath="/ALIC_1/PCRA_";
  TString symbCradleName="PHOS/Cradle";
  Int_t nCradles = 2;

  for(Int_t iCradle=0; iCradle<nCradles; iCradle++){
    volpath = physCradlePath;
    volpath += iCradle;
    symname = symbCradleName;
    symname += iCradle;
    gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data());
  }

  // Alignable wheels
  // Volume path /ALIC_1/PWHE_<i> => symbolic name /PHOS/Wheel<i>, i=0,1,2,3

  TString physWheelPath="/ALIC_1/PWHE_";
  TString symbWheelName="PHOS/Wheel";
  Int_t nWheels = 4;

  for(Int_t iWheel=0; iWheel<nWheels; iWheel++){
    volpath = physWheelPath;
    volpath += iWheel;
    symname = symbWheelName;
    symname += iWheel;
    gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data());
  }

  // Alignable strip units are not implemented yet (27.09.2006)

}

//____________________________________________________________________________
Float_t AliPHOSv0::ZMin(void) const
{
  // Overall dimension of the PHOS (min)

  AliPHOSGeometry * geom = GetGeometry() ; 

  return -geom->GetOuterBoxSize(2)/2.;
}

//____________________________________________________________________________
Float_t AliPHOSv0::ZMax(void) const
{
  // Overall dimension of the PHOS (max)

  AliPHOSGeometry * geom = GetGeometry() ; 

  return  geom->GetOuterBoxSize(2)/2.;
}

//____________________________________________________________________________
void AliPHOSv0::Init(void)
{
  // Just prints an information message
  
  Int_t i;

  if(AliLog::GetGlobalDebugLevel()>0) {
    TString st ; 
    for(i=0;i<35;i++) 
      st += "*";
    Info("Init", "%s", st.Data()) ;  
    // Here the PHOS initialisation code (if any!)
    
    AliPHOSGeometry * geom = GetGeometry() ; 

    if (geom!=0)  
      Info("Init", "AliPHOS%s: PHOS geometry intialized for %s", Version().Data(), geom->GetName()) ;
    else
      Info("Init", "AliPHOS%s: PHOS geometry initialization failed !", Version().Data()) ;       

    Info("Init", "%s", st.Data()) ;  
  }
}
