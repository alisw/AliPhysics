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
// Implementation version v0 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
// Layout EMC + CPV  has name IHEP
// An object of this class does not produce hits nor digits
// It is the one to use if you do not want to produce outputs in TREEH or TREED
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"
#include "TRandom.h"
#include "TGeometry.h"
#include "TFolder.h"
#include "TROOT.h"
#include "TTree.h"


// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---

#include "AliPHOSv0.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMC.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h"

ClassImp(AliPHOSv0)

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(const char *name, const char *title):
  AliPHOS(name,title)
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)
  //        IHEP = 5 modules (EMC + CPV)
  //        MIXT = 4 modules (EMC + CPV) and 1 module (EMC + PPSD)
 
  // create the geometry parameters object  
  // and post it to a folder (Post retrieves the correct geometry)
  AliPHOSGetter::GetInstance(gDirectory->GetName(), 0)->PostGeometry() ; 

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
  
  AliPHOSGeometry * geom = GetGeometry() ; 

  this->BuildGeometryforPHOS() ; 
  if      (strcmp(geom->GetName(),"GPS2") == 0)
    this->BuildGeometryforPPSD() ;
  else if (strcmp(geom->GetName(),"IHEP") == 0)
    this->BuildGeometryforCPV() ;
  else if (strcmp(geom->GetName(),"MIXT") == 0) {
    this->BuildGeometryforPPSD() ;
    this->BuildGeometryforCPV() ;
  }
  else
    cout << "AliPHOSv0::BuildGeometry : no charged particle identification system installed: "
	 << "Geometry name = " << geom->GetName() << endl; 

}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforPHOS(void)
{
 // Build the PHOS-EMC geometry for the ROOT display

  const Int_t kColorPHOS = kRed ;
  const Int_t kColorXTAL = kBlue ;

  Double_t const kRADDEG = 180.0 / kPI ;

  AliPHOSGeometry * geom = GetGeometry() ; 

  new TBRIK( "OuterBox", "PHOS box", "void", geom->GetOuterBoxSize(0)/2, 
                                             geom->GetOuterBoxSize(1)/2, 
                                             geom->GetOuterBoxSize(2)/2 );

  // Textolit Wall box, position inside PHOS 
  
  new TBRIK( "TextolitBox", "PHOS Textolit box ", "void", geom->GetTextolitBoxSize(0)/2, 
                                                          geom->GetTextolitBoxSize(1)/2, 
                                                          geom->GetTextolitBoxSize(2)/2);

  // Polystyrene Foam Plate

  new TBRIK( "UpperFoamPlate", "PHOS Upper foam plate", "void", geom->GetTextolitBoxSize(0)/2, 
                                                                geom->GetSecondUpperPlateThickness()/2, 
                                                                geom->GetTextolitBoxSize(2)/2 ) ; 

  // Air Filled Box
 
  new TBRIK( "AirFilledBox", "PHOS air filled box", "void", geom->GetAirFilledBoxSize(0)/2, 
                                                            geom->GetAirFilledBoxSize(1)/2, 
                                                            geom->GetAirFilledBoxSize(2)/2 );

  // Crystals Box

  Float_t xtlX = geom->GetCrystalSize(0) ; 
  Float_t xtlY = geom->GetCrystalSize(1) ; 
  Float_t xtlZ = geom->GetCrystalSize(2) ; 

  Float_t xl =  geom->GetNPhi() * ( xtlX + 2 * geom->GetGapBetweenCrystals() ) / 2.0 + geom->GetModuleBoxThickness() ;
  Float_t yl =  ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() + geom->GetCrystalHolderThickness() ) / 2.0 
             + geom->GetModuleBoxThickness() / 2.0 ;
  Float_t zl =  geom->GetNZ() * ( xtlZ + 2 * geom->GetGapBetweenCrystals() ) / 2.0 +  geom->GetModuleBoxThickness() ;
  
  new TBRIK( "CrystalsBox", "PHOS crystals box", "void", xl, yl, zl ) ;

// position PHOS into ALICE

  Float_t r = geom->GetIPtoOuterCoverDistance() + geom->GetOuterBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  Float_t pphi =  TMath::ATan( geom->GetOuterBoxSize(0)  / ( 2.0 * geom->GetIPtoOuterCoverDistance() ) ) ;
  pphi *= kRADDEG ;
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
 
  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  for( Int_t i = 1; i <= geom->GetNModules(); i++ ) { 
   Float_t angle = pphi * 2 * ( i - geom->GetNModules() / 2.0 - 0.5 ) ;
   sprintf(rotname, "%s%d", "rot", number++) ;
   new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
   top->cd();
   sprintf(nodename,"%s%d", "Module", i) ;    
   Float_t x =  r * TMath::Sin( angle / kRADDEG ) ;
   Float_t y = -r * TMath::Cos( angle / kRADDEG ) ;
   TNode * outerboxnode = new TNode(nodename, nodename, "OuterBox", x, y, 0, rotname ) ;
   outerboxnode->SetLineColor(kColorPHOS) ;
   fNodes->Add(outerboxnode) ;
   outerboxnode->cd() ; 
   // now inside the outer box the textolit box
   y = ( geom->GetOuterBoxThickness(1) -  geom->GetUpperPlateThickness() ) / 2.  ;
   sprintf(nodename,"%s%d", "TexBox", i) ;  
   TNode * textolitboxnode = new TNode(nodename, nodename, "TextolitBox", 0, y, 0) ; 
   textolitboxnode->SetLineColor(kColorPHOS) ;
   fNodes->Add(textolitboxnode) ;
   // upper foam plate inside outre box
   outerboxnode->cd() ; 
   sprintf(nodename, "%s%d", "UFPlate", i) ;
   y =  ( geom->GetTextolitBoxSize(1) - geom->GetSecondUpperPlateThickness() ) / 2.0 ;
   TNode * upperfoamplatenode = new TNode(nodename, nodename, "UpperFoamPlate", 0, y, 0) ; 
   upperfoamplatenode->SetLineColor(kColorPHOS) ;
   fNodes->Add(upperfoamplatenode) ;  
   // air filled box inside textolit box (not drawn)
   textolitboxnode->cd();
   y = ( geom->GetTextolitBoxSize(1) - geom->GetAirFilledBoxSize(1) ) / 2.0 -  geom->GetSecondUpperPlateThickness() ;
   sprintf(nodename, "%s%d", "AFBox", i) ;
   TNode * airfilledboxnode = new TNode(nodename, nodename, "AirFilledBox", 0, y, 0) ; 
   fNodes->Add(airfilledboxnode) ; 
   // crystals box inside air filled box
   airfilledboxnode->cd() ; 
   y = geom->GetAirFilledBoxSize(1) / 2.0 - yl 
       - ( geom->GetIPtoCrystalSurface() - geom->GetIPtoOuterCoverDistance() - geom->GetModuleBoxThickness() 
       -  geom->GetUpperPlateThickness() -  geom->GetSecondUpperPlateThickness() ) ; 
   sprintf(nodename, "%s%d", "XTBox", i) ; 
   TNode * crystalsboxnode = new TNode(nodename, nodename, "CrystalsBox", 0, y, 0) ;    
   crystalsboxnode->SetLineColor(kColorXTAL) ; 
   fNodes->Add(crystalsboxnode) ; 
  }

  delete[] rotname ;  
  delete[] nodename ;
}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforPPSD(void)
{
 //  Build the PHOS-PPSD geometry for the ROOT display
 //BEGIN_HTML
  /*
    <H2>
     PPSD displayed by root
    </H2>
    <UL>
    <LI> Zoom on PPSD: Front View
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="PPSD Front View" SRC="../images/AliPHOSv0PPSDFrontView.gif"> 
    </CENTER></P></LI>
    <LI> Zoom on PPSD: Perspective View
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="PPSD Prespective View" SRC="../images/AliPHOSv0PPSDPerspectiveView.gif"> 
    </CENTER></P></LI>
    </UL>
  */
  //END_HTML  
  Double_t const kRADDEG = 180.0 / kPI ;

  const Int_t kColorPHOS = kRed ;
  const Int_t kColorPPSD = kGreen ;
  const Int_t kColorGas  = kBlue ;  
  const Int_t kColorAir  = kYellow ; 

  AliPHOSGeometry * geom = GetGeometry() ; 

  // Box for a full PHOS module

  new TBRIK( "PPSDBox", "PPSD box", "void",  geom->GetCPVBoxSize(0)/2, 
                                             geom->GetCPVBoxSize(1)/2, 
	                                     geom->GetCPVBoxSize(2)/2 );

  // Box containing one micromegas module 

  new TBRIK( "PPSDModule", "PPSD module", "void",  geom->GetPPSDModuleSize(0)/2, 
                                                   geom->GetPPSDModuleSize(1)/2, 
	                                           geom->GetPPSDModuleSize(2)/2 );
 // top lid

  new TBRIK ( "TopLid", "Micromegas top lid", "void",  geom->GetPPSDModuleSize(0)/2,
                                                       geom->GetLidThickness()/2,
                                                       geom->GetPPSDModuleSize(2)/2 ) ; 
 // composite panel (top and bottom)

  new TBRIK ( "TopPanel", "Composite top panel", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
                                                            geom->GetCompositeThickness()/2,
                                                          ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ;  
  
  new TBRIK ( "BottomPanel", "Composite bottom panel", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
                                                                  geom->GetCompositeThickness()/2,
                                                                ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ; 
 // gas gap (conversion and avalanche)

  new TBRIK ( "GasGap", "gas gap", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
	                                    ( geom->GetConversionGap() +  geom->GetAvalancheGap() )/2,
                                            ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ; 

 // anode and cathode 

  new TBRIK ( "Anode", "Anode", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
                                           geom->GetAnodeThickness()/2,
                                         ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ; 

  new TBRIK ( "Cathode", "Cathode", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
                                               geom->GetCathodeThickness()/2,
                                             ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ; 
 // PC  

  new TBRIK ( "PCBoard", "Printed Circuit", "void",  ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() )/2,
                                                       geom->GetPCThickness()/2,
                                                     ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() )/2 ) ; 
 // Gap between Lead and top micromegas

  new TBRIK ( "LeadToM", "Air Gap top", "void", geom->GetCPVBoxSize(0)/2,
                                                geom->GetMicro1ToLeadGap()/2,
                                                geom->GetCPVBoxSize(2)/2  ) ;  
 
// Gap between Lead and bottom micromegas

  new TBRIK ( "MToLead", "Air Gap bottom", "void", geom->GetCPVBoxSize(0)/2,
                                                   geom->GetLeadToMicro2Gap()/2,
                                                   geom->GetCPVBoxSize(2)/2  ) ; 
 // Lead converter
   
  new TBRIK ( "Lead", "Lead converter", "void", geom->GetCPVBoxSize(0)/2,
                                                geom->GetLeadConverterThickness()/2,
                                                geom->GetCPVBoxSize(2)/2  ) ; 

     // position PPSD into ALICE

  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  Float_t r = geom->GetIPtoTopLidDistance() + geom->GetCPVBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
 
  Int_t firstModule = 0 ; 
  if      (strcmp(geom->GetName(),"GPS2") == 0) 
    firstModule = 1;
  else if (strcmp(geom->GetName(),"MIXT") == 0) 
    firstModule = geom->GetNModules() - geom->GetNPPSDModules() + 1;
  
  for( Int_t i = firstModule; i <= geom->GetNModules(); i++ ) { // the number of PHOS modules
    Float_t angle = geom->GetPHOSAngle(i) ;
    sprintf(rotname, "%s%d", "rotg", number+i) ;
    new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
    top->cd();
    sprintf(nodename, "%s%d", "Moduleg", i) ;    
    Float_t x =  r * TMath::Sin( angle / kRADDEG ) ;
    Float_t y = -r * TMath::Cos( angle / kRADDEG ) ;
    TNode * ppsdboxnode = new TNode(nodename , nodename ,"PPSDBox", x, y, 0, rotname ) ;
    ppsdboxnode->SetLineColor(kColorPPSD) ;
    fNodes->Add(ppsdboxnode) ;
    ppsdboxnode->cd() ;
    // inside the PPSD box: 
    //   1.   fNumberOfModulesPhi x fNumberOfModulesZ top micromegas
    x = ( geom->GetCPVBoxSize(0) - geom->GetPPSDModuleSize(0) ) / 2. ;  
    {
      for ( Int_t iphi = 1; iphi <= geom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
	Float_t z = ( geom->GetCPVBoxSize(2) - geom->GetPPSDModuleSize(2) ) / 2. ;
	TNode * micro1node ; 
	for ( Int_t iz = 1; iz <= geom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	  y = ( geom->GetCPVBoxSize(1) - geom->GetMicromegas1Thickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Mic1", i, iphi, iz) ;
	  micro1node  = new TNode(nodename, nodename, "PPSDModule", x, y, z) ;
	  micro1node->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(micro1node) ; 
	  // inside top micromegas
	  micro1node->cd() ; 
	  //      a. top lid
	  y = ( geom->GetMicromegas1Thickness() - geom->GetLidThickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Lid", i, iphi, iz) ;
	  TNode * toplidnode = new TNode(nodename, nodename, "TopLid", 0, y, 0) ;
	  toplidnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(toplidnode) ; 
	  //      b. composite panel
	  y = y - geom->GetLidThickness() / 2. - geom->GetCompositeThickness() / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompU", i, iphi, iz) ;
	  TNode * compupnode = new TNode(nodename, nodename, "TopPanel", 0, y, 0) ;
	  compupnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compupnode) ; 
	  //      c. anode
	  y = y - geom->GetCompositeThickness() / 2. - geom->GetAnodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Ano", i, iphi, iz) ;
	  TNode * anodenode = new TNode(nodename, nodename, "Anode", 0, y, 0) ;
	  anodenode->SetLineColor(kColorPHOS) ;  
	  fNodes->Add(anodenode) ; 
	  //      d.  gas 
	  y = y - geom->GetAnodeThickness() / 2. - ( geom->GetConversionGap() +  geom->GetAvalancheGap() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "GGap", i, iphi, iz) ;
	  TNode * ggapnode = new TNode(nodename, nodename, "GasGap", 0, y, 0) ;
	  ggapnode->SetLineColor(kColorGas) ;  
	  fNodes->Add(ggapnode) ;          
	  //      f. cathode
	  y = y - ( geom->GetConversionGap() +  geom->GetAvalancheGap() ) / 2. - geom->GetCathodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Cathode", i, iphi, iz) ;
	  TNode * cathodenode = new TNode(nodename, nodename, "Cathode", 0, y, 0) ;
	  cathodenode->SetLineColor(kColorPHOS) ;  
	  fNodes->Add(cathodenode) ;        
	  //      g. printed circuit
	  y = y - geom->GetCathodeThickness() / 2. - geom->GetPCThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "PC", i, iphi, iz) ;
	  TNode * pcnode = new TNode(nodename, nodename, "PCBoard", 0, y, 0) ;
	  pcnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(pcnode) ;        
	  //      h. composite panel
	  y = y - geom->GetPCThickness() / 2. - geom->GetCompositeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompDown", i, iphi, iz) ;
	  TNode * compdownnode = new TNode(nodename, nodename, "BottomPanel", 0, y, 0) ;
	  compdownnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compdownnode) ;   
	  z = z - geom->GetPPSDModuleSize(2) ;
	  ppsdboxnode->cd() ;
	} // end of Z module loop     
	x = x -  geom->GetPPSDModuleSize(0) ; 
	ppsdboxnode->cd() ;
      } // end of phi module loop
    }
    //   2. air gap      
    ppsdboxnode->cd() ;
    y = ( geom->GetCPVBoxSize(1) - 2 * geom->GetMicromegas1Thickness() - geom->GetMicro1ToLeadGap() ) / 2. ; 
    sprintf(nodename, "%s%d", "GapUp", i) ;
    TNode * gapupnode = new TNode(nodename, nodename, "LeadToM", 0, y, 0) ;
    gapupnode->SetLineColor(kColorAir) ;  
    fNodes->Add(gapupnode) ;        
    //   3. lead converter
    y = y - geom->GetMicro1ToLeadGap() / 2. - geom->GetLeadConverterThickness() / 2. ; 
    sprintf(nodename, "%s%d", "LeadC", i) ;
    TNode * leadcnode = new TNode(nodename, nodename, "Lead", 0, y, 0) ;
    leadcnode->SetLineColor(kColorPPSD) ;  
    fNodes->Add(leadcnode) ;        
    //   4. air gap
    y = y - geom->GetLeadConverterThickness() / 2. - geom->GetLeadToMicro2Gap()  / 2. ; 
    sprintf(nodename, "%s%d", "GapDown", i) ;
    TNode * gapdownnode = new TNode(nodename, nodename, "MToLead", 0, y, 0) ;
    gapdownnode->SetLineColor(kColorAir) ;  
    fNodes->Add(gapdownnode) ;        
    //    5.  fNumberOfModulesPhi x fNumberOfModulesZ bottom micromegas
    x = ( geom->GetCPVBoxSize(0) - geom->GetPPSDModuleSize(0) ) / 2. - geom->GetPhiDisplacement() ;  
    {
      for ( Int_t iphi = 1; iphi <= geom->GetNumberOfModulesPhi(); iphi++ ) { 
	Float_t z = ( geom->GetCPVBoxSize(2) - geom->GetPPSDModuleSize(2) ) / 2.  - geom->GetZDisplacement() ;;
	TNode * micro2node ; 
	for ( Int_t iz = 1; iz <= geom->GetNumberOfModulesZ(); iz++ ) { 
	  y = - ( geom->GetCPVBoxSize(1) - geom->GetMicromegas2Thickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Mic2", i, iphi, iz) ;
	  micro2node  = new TNode(nodename, nodename, "PPSDModule", x, y, z) ;
	  micro2node->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(micro2node) ; 
	  // inside bottom micromegas
	  micro2node->cd() ; 
	  //      a. top lid

	  y = ( geom->GetMicromegas2Thickness() - geom->GetLidThickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Lidb", i, iphi, iz) ;

	  TNode * toplidbnode = new TNode(nodename, nodename, "TopLid", 0, y, 0) ;
	  toplidbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(toplidbnode) ; 
	  //      b. composite panel

	  y = y - geom->GetLidThickness() / 2. - geom->GetCompositeThickness() / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompUb", i, iphi, iz) ;

	  TNode * compupbnode = new TNode(nodename, nodename, "TopPanel", 0, y, 0) ;
	  compupbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compupbnode) ; 
	  //      c. anode

	  y = y - geom->GetCompositeThickness() / 2. - geom->GetAnodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Anob", i, iphi, iz) ;

	  TNode * anodebnode = new TNode(nodename, nodename, "Anode", 0, y, 0) ;
	  anodebnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(anodebnode) ; 
	  //      d. conversion gas

	  y = y - geom->GetAnodeThickness() / 2. - ( geom->GetConversionGap() +  geom->GetAvalancheGap() )  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "GGapb", i, iphi, iz) ;

	  TNode * ggapbnode = new TNode(nodename, nodename, "GasGap", 0, y, 0) ;
	  ggapbnode->SetLineColor(kColorGas) ;  
	  fNodes->Add(ggapbnode) ;           
	  //      f. cathode

	  y = y - ( geom->GetConversionGap() + geom->GetAvalancheGap() ) / 2. - geom->GetCathodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Cathodeb", i, iphi, iz) ;

	  TNode * cathodebnode = new TNode(nodename, nodename, "Cathode", 0, y, 0) ;
	  cathodebnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(cathodebnode) ;        
	  //      g. printed circuit
	  y = y - geom->GetCathodeThickness() / 2. - geom->GetPCThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "PCb", i, iphi, iz) ;
	  TNode * pcbnode = new TNode(nodename, nodename, "PCBoard", 0, y, 0) ;
	  pcbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(pcbnode) ;        
	  //      h. composite pane
	  y = y - geom->GetPCThickness() / 2. - geom->GetCompositeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompDownb", i, iphi, iz) ;
	  TNode * compdownbnode = new TNode(nodename, nodename, "BottomPanel", 0, y, 0) ;
	  compdownbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compdownbnode) ;        
       	  z = z - geom->GetPPSDModuleSize(2) ;
	  ppsdboxnode->cd() ;
	} // end of Z module loop     
	x = x -  geom->GetPPSDModuleSize(0) ; 
	ppsdboxnode->cd() ;
      } // end of phi module loop
    }
  } // PHOS modules
 
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

  const Double_t kRADDEG         = 180.0 / kPI ;
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
  if      (strcmp(geom->GetName(),"IHEP") == 0) 
    lastModule = geom->GetNModules();
  else if (strcmp(geom->GetName(),"MIXT") == 0) 
    lastModule = geom->GetNModules() - geom->GetNPPSDModules();
  
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

  AliPHOSv0 *phostmp = (AliPHOSv0*)gAlice->GetModule("PHOS") ;

  if ( phostmp == NULL ) {
    
    fprintf(stderr, "PHOS detector not found!\n") ;
    return;
    
  }

  AliPHOSGeometry * geom = GetGeometry() ; 

  // Get pointer to the array containing media indeces
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  // Create a box a PHOS module.
  // In case of MIXT geometry 2 different boxes are needed

  Float_t bigbox[3] ; 
  bigbox[0] =   geom->GetOuterBoxSize(0) / 2.0 ;
  bigbox[1] = ( geom->GetOuterBoxSize(1) + geom->GetCPVBoxSize(1) ) / 2.0 ;
  bigbox[2] =   geom->GetOuterBoxSize(2) / 2.0 ;
  
    gMC->Gsvolu("PHOS", "BOX ", idtmed[798], bigbox, 3) ;

  if ( strcmp( geom->GetName(),"MIXT") == 0 && geom->GetNPPSDModules() > 0) 
    gMC->Gsvolu("PHO1", "BOX ", idtmed[798], bigbox, 3) ;
  
    this->CreateGeometryforPHOS() ; 
  if      ( strcmp( geom->GetName(), "GPS2") == 0  ) 
    this->CreateGeometryforPPSD() ;
  else if ( strcmp( geom->GetName(), "IHEP") == 0  ) 
    this->CreateGeometryforCPV() ;
  else if ( strcmp( geom->GetName(), "MIXT") == 0  ) {
    this->CreateGeometryforPPSD() ;
    this->CreateGeometryforCPV() ;
  }
  else
    cout << "AliPHOSv0::CreateGeometry : no charged particle identification system installed" << endl; 

  this->CreateGeometryforSupport() ; 
  
  // --- Position  PHOS mdules in ALICE setup ---
  
  Int_t idrotm[99] ;
  Double_t const kRADDEG = 180.0 / kPI ;
  
  Int_t lastModule;
  if (strcmp(geom->GetName(),"MIXT") == 0) 
    lastModule = geom->GetNModules() - geom->GetNPPSDModules();
  else
    lastModule = geom->GetNModules();

  Int_t i;
  for( i = 1; i <= lastModule ; i++ ) {
    
    Float_t angle = geom->GetPHOSAngle(i) ;
    AliMatrix(idrotm[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0) ;
 
    Float_t r = geom->GetIPtoOuterCoverDistance() + ( geom->GetOuterBoxSize(1) + geom->GetCPVBoxSize(1) ) / 2.0 ;

    Float_t xP1 =  r * TMath::Sin( angle / kRADDEG ) ;
    Float_t yP1 = -r * TMath::Cos( angle / kRADDEG ) ;

    gMC->Gspos("PHOS", i, "ALIC", xP1, yP1, 0.0, idrotm[i-1], "ONLY") ;
 
  } // for GetNModules

  for( i = lastModule+1; i <= geom->GetNModules(); i++ ) {
    
    Float_t angle = geom->GetPHOSAngle(i) ;
    AliMatrix(idrotm[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0) ;
 
    Float_t r = geom->GetIPtoOuterCoverDistance() + ( geom->GetOuterBoxSize(1) + geom->GetCPVBoxSize(1) ) / 2.0 ;

    Float_t xP1 =  r * TMath::Sin( angle / kRADDEG ) ;
    Float_t yP1 = -r * TMath::Cos( angle / kRADDEG ) ;

    gMC->Gspos("PHO1", i-lastModule, "ALIC", xP1, yP1, 0.0, idrotm[i-1], "ONLY") ;
 
  } // for GetNModules

}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforPHOS()
{
  // Create the PHOS-EMC geometry for GEANT
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

  // ---
  // --- Define PHOS box volume, fPUFPill with thermo insulating foam ---
  // --- Foam Thermo Insulating outer cover dimensions ---
  // --- Put it in bigbox = PHOS

  Float_t dphos[3] ; 
  dphos[0] =  geom->GetOuterBoxSize(0) / 2.0 ;
  dphos[1] =  geom->GetOuterBoxSize(1) / 2.0 ;
  dphos[2] =  geom->GetOuterBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PEMC", "BOX ", idtmed[706], dphos, 3) ;

  Float_t yO =  - geom->GetCPVBoxSize(1)  / 2.0 ;

    gMC->Gspos("PEMC", 1, "PHOS", 0.0, yO, 0.0, 0, "ONLY") ; 
  if ( strcmp( geom->GetName(),"MIXT") == 0 && geom->GetNPPSDModules() > 0) 
    gMC->Gspos("PEMC", 1, "PHO1", 0.0, yO, 0.0, 0, "ONLY") ; 

  // ---
  // --- Define Textolit Wall box, position inside PEMC ---
  // --- Textolit Wall box dimentions ---
 
 
  Float_t dptxw[3];
  dptxw[0] = geom->GetTextolitBoxSize(0) / 2.0 ;
  dptxw[1] = geom->GetTextolitBoxSize(1) / 2.0 ;
  dptxw[2] = geom->GetTextolitBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXW", "BOX ", idtmed[707], dptxw, 3);

  yO =   (  geom->GetOuterBoxThickness(1) -   geom->GetUpperPlateThickness() ) / 2.  ;
   
  gMC->Gspos("PTXW", 1, "PEMC", 0.0, yO, 0.0, 0, "ONLY") ;

  // --- 
  // --- Define Upper Polystyrene Foam Plate, place inside PTXW ---
  // --- immediately below Foam Thermo Insulation Upper plate ---

  // --- Upper Polystyrene Foam plate thickness ---
 
  Float_t  dpufp[3] ;
  dpufp[0] = geom->GetTextolitBoxSize(0) / 2.0 ; 
  dpufp[1] = geom->GetSecondUpperPlateThickness() / 2. ;
  dpufp[2] = geom->GetTextolitBoxSize(2) /2.0 ; 

  gMC->Gsvolu("PUFP", "BOX ", idtmed[703], dpufp, 3) ;
  
  yO = ( geom->GetTextolitBoxSize(1) -  geom->GetSecondUpperPlateThickness() ) / 2.0 ;
  
  gMC->Gspos("PUFP", 1, "PTXW", 0.0, yO, 0.0, 0, "ONLY") ;
  
  // ---
  // --- Define air-filled box, place inside PTXW ---
  // --- Inner AIR volume dimensions ---
 

  Float_t  dpair[3] ;
  dpair[0] = geom->GetAirFilledBoxSize(0) / 2.0 ;
  dpair[1] = geom->GetAirFilledBoxSize(1) / 2.0 ;
  dpair[2] = geom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PAIR", "BOX ", idtmed[798], dpair, 3) ;
  
  yO = ( geom->GetTextolitBoxSize(1) -  geom->GetAirFilledBoxSize(1) ) / 2.0 -   geom->GetSecondUpperPlateThickness() ;
  
  gMC->Gspos("PAIR", 1, "PTXW", 0.0, yO, 0.0, 0, "ONLY") ;

// --- Dimensions of PbWO4 crystal ---

  Float_t xtlX =  geom->GetCrystalSize(0) ; 
  Float_t xtlY =  geom->GetCrystalSize(1) ; 
  Float_t xtlZ =  geom->GetCrystalSize(2) ; 

  Float_t dptcb[3] ;  
  dptcb[0] =  geom->GetNPhi() * ( xtlX + 2 *  geom->GetGapBetweenCrystals() ) / 2.0 + geom->GetModuleBoxThickness() ;
  dptcb[1] = ( xtlY +  geom->GetCrystalSupportHeight() +  geom->GetCrystalWrapThickness() + geom->GetCrystalHolderThickness() ) / 2.0 
             + geom->GetModuleBoxThickness() / 2.0 ;
  dptcb[2] = geom->GetNZ() * ( xtlZ + 2 * geom->GetGapBetweenCrystals() ) / 2.0 +  geom->GetModuleBoxThickness() ;
  
  gMC->Gsvolu("PTCB", "BOX ", idtmed[706], dptcb, 3) ;

  yO =  geom->GetAirFilledBoxSize(1) / 2.0 - dptcb[1] 
       - ( geom->GetIPtoCrystalSurface() - geom->GetIPtoOuterCoverDistance() - geom->GetModuleBoxThickness() 
       -  geom->GetUpperPlateThickness() -  geom->GetSecondUpperPlateThickness() ) ;
  
  gMC->Gspos("PTCB", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Crystal BLock filled with air, position it inside PTCB ---
  Float_t dpcbl[3] ; 
  
  dpcbl[0] = geom->GetNPhi() * ( xtlX + 2 * geom->GetGapBetweenCrystals() ) / 2.0 ;
  dpcbl[1] = ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() + geom->GetCrystalHolderThickness() ) / 2.0 ;
  dpcbl[2] = geom->GetNZ() * ( xtlZ + 2 * geom->GetGapBetweenCrystals() ) / 2.0 ;
  
  gMC->Gsvolu("PCBL", "BOX ", idtmed[798], dpcbl, 3) ;
  
  // --- Divide PCBL in X (phi) and Z directions --
  gMC->Gsdvn("PROW", "PCBL", Int_t (geom->GetNPhi()), 1) ;
  gMC->Gsdvn("PCEL", "PROW", Int_t (geom->GetNZ()), 3) ;

  yO = -geom->GetModuleBoxThickness() / 2.0 ;
  
  gMC->Gspos("PCBL", 1, "PTCB", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define STeel (actually, it's titanium) Cover volume, place inside PCEL
  Float_t  dpstc[3] ; 
  
  dpstc[0] = ( xtlX + 2 * geom->GetCrystalWrapThickness() ) / 2.0 ;
  dpstc[1] = ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() + geom->GetCrystalHolderThickness() ) / 2.0 ;
  dpstc[2] = ( xtlZ + 2 * geom->GetCrystalWrapThickness()  + 2 *  geom->GetCrystalHolderThickness() ) / 2.0 ;
  
  gMC->Gsvolu("PSTC", "BOX ", idtmed[704], dpstc, 3) ;

  gMC->Gspos("PSTC", 1, "PCEL", 0.0, 0.0, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Tyvek volume, place inside PSTC ---
  Float_t  dppap[3] ;

  dppap[0] = xtlX / 2.0 + geom->GetCrystalWrapThickness() ;
  dppap[1] = ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() ) / 2.0 ;
  dppap[2] = xtlZ / 2.0 + geom->GetCrystalWrapThickness() ;
  
  gMC->Gsvolu("PPAP", "BOX ", idtmed[702], dppap, 3) ;
  
  yO = ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() ) / 2.0 
              - ( xtlY +  geom->GetCrystalSupportHeight() +  geom->GetCrystalWrapThickness() + geom->GetCrystalHolderThickness() ) / 2.0 ;
   
  gMC->Gspos("PPAP", 1, "PSTC", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PbWO4 crystal volume, place inside PPAP ---
  Float_t  dpxtl[3] ; 

  dpxtl[0] = xtlX / 2.0 ;
  dpxtl[1] = xtlY / 2.0 ;
  dpxtl[2] = xtlZ / 2.0 ;
  
  gMC->Gsvolu("PXTL", "BOX ", idtmed[699], dpxtl, 3) ;

  yO = ( xtlY + geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() ) / 2.0 - xtlY / 2.0 - geom->GetCrystalWrapThickness() ;
  
  gMC->Gspos("PXTL", 1, "PPAP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define crystal support volume, place inside PPAP ---
  Float_t dpsup[3] ; 

  dpsup[0] = xtlX / 2.0 + geom->GetCrystalWrapThickness()  ;
  dpsup[1] = geom->GetCrystalSupportHeight() / 2.0 ;
  dpsup[2] = xtlZ / 2.0 +  geom->GetCrystalWrapThickness() ;

  gMC->Gsvolu("PSUP", "BOX ", idtmed[798], dpsup, 3) ;

  yO =  geom->GetCrystalSupportHeight() / 2.0 - ( xtlY +  geom->GetCrystalSupportHeight() + geom->GetCrystalWrapThickness() ) / 2.0 ;

  gMC->Gspos("PSUP", 1, "PPAP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PIN-diode volume and position it inside crystal support ---
  // --- right behind PbWO4 crystal

  // --- PIN-diode dimensions ---

 
  Float_t dppin[3] ;
  dppin[0] = geom->GetPinDiodeSize(0) / 2.0 ;
  dppin[1] = geom->GetPinDiodeSize(1) / 2.0 ;
  dppin[2] = geom->GetPinDiodeSize(2) / 2.0 ;
 
  gMC->Gsvolu("PPIN", "BOX ", idtmed[705], dppin, 3) ;
 
  yO = geom->GetCrystalSupportHeight() / 2.0 - geom->GetPinDiodeSize(1) / 2.0 ;
 
  gMC->Gspos("PPIN", 1, "PSUP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Upper Cooling Panel, place it on top of PTCB ---
  Float_t dpucp[3] ;
 // --- Upper Cooling Plate thickness ---
 
  dpucp[0] = dptcb[0] ;
  dpucp[1] = geom->GetUpperCoolingPlateThickness() ;
  dpucp[2] = dptcb[2] ;
  
  gMC->Gsvolu("PUCP", "BOX ", idtmed[701], dpucp,3) ;
  
  yO = geom->GetAirFilledBoxSize(1) / 2. 
    -( geom->GetIPtoCrystalSurface()  - geom->GetIPtoOuterCoverDistance()    - geom->GetModuleBoxThickness()
      -geom->GetUpperPlateThickness() - geom->GetSecondUpperPlateThickness() - geom->GetUpperCoolingPlateThickness() ) ; 
  
  gMC->Gspos("PUCP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Al Support Plate, position it inside PAIR ---
  // --- right beneath PTCB ---
 // --- Al Support Plate thickness ---
 
  Float_t dpasp[3] ;
  dpasp[0] =  geom->GetAirFilledBoxSize(0) / 2.0 ;
  dpasp[1] = geom->GetSupportPlateThickness() / 2.0 ;
  dpasp[2] =  geom->GetAirFilledBoxSize(2) / 2.0 ;
  
  gMC->Gsvolu("PASP", "BOX ", idtmed[701], dpasp, 3) ;
  
  yO = (  geom->GetAirFilledBoxSize(1) - geom->GetSupportPlateThickness() ) / 2. 
       -  ( geom->GetIPtoCrystalSurface() - geom->GetIPtoOuterCoverDistance()
           - geom->GetUpperPlateThickness() - geom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 ) ;
  
  gMC->Gspos("PASP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Thermo Insulating Plate, position it inside PAIR ---
  // --- right beneath PASP ---
  // --- Lower Thermo Insulating Plate thickness ---
  
  Float_t dptip[3] ;
  dptip[0] = geom->GetAirFilledBoxSize(0) / 2.0 ;
  dptip[1] = geom->GetLowerThermoPlateThickness() / 2.0 ;
  dptip[2] = geom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTIP", "BOX ", idtmed[706], dptip, 3) ;

  yO =  ( geom->GetAirFilledBoxSize(1) - geom->GetLowerThermoPlateThickness() ) / 2. 
       -  ( geom->GetIPtoCrystalSurface() - geom->GetIPtoOuterCoverDistance() - geom->GetUpperPlateThickness() 
            - geom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 + geom->GetSupportPlateThickness() ) ;

  gMC->Gspos("PTIP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Textolit Plate, position it inside PAIR ---
  // --- right beneath PTIP ---
  // --- Lower Textolit Plate thickness ---
 
  Float_t dptxp[3] ;
  dptxp[0] = geom->GetAirFilledBoxSize(0) / 2.0 ;
  dptxp[1] = geom->GetLowerTextolitPlateThickness() / 2.0 ;
  dptxp[2] = geom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXP", "BOX ", idtmed[707], dptxp, 3) ;

  yO =  ( geom->GetAirFilledBoxSize(1) - geom->GetLowerTextolitPlateThickness() ) / 2. 
       -  ( geom->GetIPtoCrystalSurface() - geom->GetIPtoOuterCoverDistance() - geom->GetUpperPlateThickness() 
            - geom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 + geom->GetSupportPlateThickness() 
            +  geom->GetLowerThermoPlateThickness() ) ;

  gMC->Gspos("PTXP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforPPSD()
{
  // Create the PHOS-PPSD geometry for GEANT
  //BEGIN_HTML
  /*
    <H2>
    Geant3 geometry tree of PHOS-PPSD in ALICE
    </H2>
    <P><CENTER>
    <IMG Align=BOTTOM ALT="PPSD geant tree" SRC="../images/PPSDinAlice.gif"> 
    </CENTER><P>
  */
  //END_HTML  

  // Get pointer to the array containing media indexes
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  AliPHOSGeometry * geom = GetGeometry() ; 

  // The box containing all ppsd's for one PHOS module filled with air 
  Float_t ppsd[3] ; 
  ppsd[0] = geom->GetCPVBoxSize(0) / 2.0 ;  
  ppsd[1] = geom->GetCPVBoxSize(1) / 2.0 ; 
  ppsd[2] = geom->GetCPVBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PPSD", "BOX ", idtmed[798], ppsd, 3) ;

  Float_t yO =  geom->GetOuterBoxSize(1) / 2.0 ;

  if ( strcmp( geom->GetName(),"MIXT") == 0 && geom->GetNPPSDModules() > 0) 
    gMC->Gspos("PPSD", 1, "PHO1", 0.0, yO, 0.0, 0, "ONLY") ; 
  else
    gMC->Gspos("PPSD", 1, "PHOS", 0.0, yO, 0.0, 0, "ONLY") ; 

  // Now we build a micromegas module
  // The box containing the whole module filled with epoxy (FR4)

  Float_t mppsd[3] ;  
  mppsd[0] = geom->GetPPSDModuleSize(0) / 2.0 ;  
  mppsd[1] = geom->GetPPSDModuleSize(1) / 2.0 ;  
  mppsd[2] = geom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("PMPP", "BOX ", idtmed[708], mppsd, 3) ;  
 
  // Inside mppsd :
  // 1. The Top Lid made of epoxy (FR4) 

  Float_t tlppsd[3] ; 
  tlppsd[0] = geom->GetPPSDModuleSize(0) / 2.0 ; 
  tlppsd[1] = geom->GetLidThickness() / 2.0 ;
  tlppsd[2] = geom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("PTLP", "BOX ", idtmed[708], tlppsd, 3) ; 

  Float_t  y0 = ( geom->GetMicromegas1Thickness() - geom->GetLidThickness() ) / 2. ; 

  gMC->Gspos("PTLP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 
 
  // 2. the upper panel made of composite material

  Float_t upppsd[3] ; 
  upppsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2.0 ;
  upppsd[1] = geom->GetCompositeThickness() / 2.0 ;
  upppsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0 ;
 
  gMC->Gsvolu("PUPP", "BOX ", idtmed[709], upppsd, 3) ; 
  
  y0 = y0 - geom->GetLidThickness() / 2. - geom->GetCompositeThickness() / 2. ; 

  gMC->Gspos("PUPP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 3. the anode made of Copper
  
  Float_t anppsd[3] ; 
  anppsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2.0 ; 
  anppsd[1] = geom->GetAnodeThickness() / 2.0 ; 
  anppsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0  ; 

  gMC->Gsvolu("PANP", "BOX ", idtmed[710], anppsd, 3) ; 
  
  y0 = y0 - geom->GetCompositeThickness() / 2. - geom->GetAnodeThickness()  / 2. ; 
  
  gMC->Gspos("PANP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 4. the conversion gap + avalanche gap filled with gas

  Float_t ggppsd[3] ; 
  ggppsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2.0 ;
  ggppsd[1] = ( geom->GetConversionGap() +  geom->GetAvalancheGap() ) / 2.0 ; 
  ggppsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("PGGP", "BOX ", idtmed[715], ggppsd, 3) ; 
  
  // --- Divide GGPP in X (phi) and Z directions --
  gMC->Gsdvn("PPRO", "PGGP", geom->GetNumberOfPadsPhi(), 1) ;
  gMC->Gsdvn("PPCE", "PPRO", geom->GetNumberOfPadsZ() ,  3) ;

  y0 = y0 - geom->GetAnodeThickness() / 2.  - ( geom->GetConversionGap() +  geom->GetAvalancheGap() ) / 2. ; 

  gMC->Gspos("PGGP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 


  // 6. the cathode made of Copper

  Float_t cappsd[3] ;
  cappsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2.0 ;
  cappsd[1] = geom->GetCathodeThickness() / 2.0 ; 
  cappsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0  ;

  gMC->Gsvolu("PCAP", "BOX ", idtmed[710], cappsd, 3) ; 

  y0 = y0 - ( geom->GetConversionGap() +  geom->GetAvalancheGap() ) / 2. - geom->GetCathodeThickness()  / 2. ; 

  gMC->Gspos("PCAP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 7. the printed circuit made of G10       

  Float_t pcppsd[3] ; 
  pcppsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2,.0 ; 
  pcppsd[1] = geom->GetPCThickness() / 2.0 ; 
  pcppsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("PCPS", "BOX ", idtmed[711], cappsd, 3) ; 

  y0 = y0 - geom->GetCathodeThickness() / 2. - geom->GetPCThickness()  / 2. ; 

  gMC->Gspos("PCPS", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 8. the lower panel made of composite material
						    
  Float_t lpppsd[3] ; 
  lpppsd[0] = ( geom->GetPPSDModuleSize(0) - geom->GetMicromegasWallThickness() ) / 2.0 ; 
  lpppsd[1] = geom->GetCompositeThickness() / 2.0 ; 
  lpppsd[2] = ( geom->GetPPSDModuleSize(2) - geom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("PLPP", "BOX ", idtmed[709], lpppsd, 3) ; 
 
  y0 = y0 - geom->GetPCThickness() / 2. - geom->GetCompositeThickness()  / 2. ; 

  gMC->Gspos("PLPP", 1, "PMPP", 0.0, y0, 0.0, 0, "ONLY") ; 

  // Position the  fNumberOfModulesPhi x fNumberOfModulesZ modules (mppsd) inside PPSD to cover a PHOS module
  // the top and bottom one's (which are assumed identical) :

   Float_t yt = ( geom->GetCPVBoxSize(1) - geom->GetMicromegas1Thickness() ) / 2. ; 
   Float_t yb = - ( geom->GetCPVBoxSize(1) - geom->GetMicromegas2Thickness() ) / 2. ; 

   Int_t copyNumbertop = 0 ; 
   Int_t copyNumberbot = geom->GetNumberOfModulesPhi() *  geom->GetNumberOfModulesZ() ; 

   Float_t x  = ( geom->GetCPVBoxSize(0) - geom->GetPPSDModuleSize(0) ) / 2. ;  

   for ( Int_t iphi = 1; iphi <= geom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
      Float_t z = ( geom->GetCPVBoxSize(2) - geom->GetPPSDModuleSize(2) ) / 2. ;

      for ( Int_t iz = 1; iz <= geom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	gMC->Gspos("PMPP", ++copyNumbertop, "PPSD", x, yt, z, 0, "ONLY") ;
	gMC->Gspos("PMPP", ++copyNumberbot, "PPSD", x, yb, z, 0, "ONLY") ; 
	z = z - geom->GetPPSDModuleSize(2) ;
      } // end of Z module loop   
      x = x -  geom->GetPPSDModuleSize(0) ; 
    } // end of phi module loop

   // The Lead converter between two air gaps
   // 1. Upper air gap

   Float_t uappsd[3] ;
   uappsd[0] = geom->GetCPVBoxSize(0) / 2.0 ;
   uappsd[1] = geom->GetMicro1ToLeadGap() / 2.0 ; 
   uappsd[2] = geom->GetCPVBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PUAPPS", "BOX ", idtmed[798], uappsd, 3) ; 

  y0 = ( geom->GetCPVBoxSize(1) - 2 * geom->GetMicromegas1Thickness() - geom->GetMicro1ToLeadGap() ) / 2. ; 

  gMC->Gspos("PUAPPS", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 

   // 2. Lead converter
 
  Float_t lcppsd[3] ; 
  lcppsd[0] = geom->GetCPVBoxSize(0) / 2.0 ;
  lcppsd[1] = geom->GetLeadConverterThickness() / 2.0 ; 
  lcppsd[2] = geom->GetCPVBoxSize(2) / 2.0 ;
 
  gMC->Gsvolu("PLCPPS", "BOX ", idtmed[712], lcppsd, 3) ; 
  
  y0 = y0 - geom->GetMicro1ToLeadGap() / 2. - geom->GetLeadConverterThickness() / 2. ; 

  gMC->Gspos("PLCPPS", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 3. Lower air gap

  Float_t lappsd[3] ; 
  lappsd[0] = geom->GetCPVBoxSize(0) / 2.0 ; 
  lappsd[1] = geom->GetLeadToMicro2Gap() / 2.0 ; 
  lappsd[2] = geom->GetCPVBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PLAPPS", "BOX ", idtmed[798], lappsd, 3) ; 
    
  y0 = y0 - geom->GetLeadConverterThickness() / 2. - geom->GetLeadToMicro2Gap()  / 2. ; 
  
  gMC->Gspos("PLAPPS", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 
   
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
  
  y = geom->GetOuterBoxSize(1) / 2.0 ;
  gMC->Gspos("PCPV", 1, "PHOS", 0.0, y, 0.0, 0, "ONLY") ; 
  
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

  par[1] =  TMath::Sqrt(
			TMath::Power((geom->GetIPtoOuterCoverDistance() + geom->GetOuterBoxSize(1)),2) +
			TMath::Power((geom->GetOuterBoxSize(0)/2),2)) + 10.;
  par[0] =  par[1] - geom->GetCradleWall(1) ;
  par[2] =  geom->GetCradleWall(2) / 2.0 ;
  par[3] =  geom->GetCradleWall(3) ;
  par[4] =  geom->GetCradleWall(4) ;
  gMC->Gsvolu("PCRA", "TUBS", idtmed[716], par, 5) ;

  par[0] -=  geom->GetCradleWallThickness() ;
  par[1] -=  geom->GetCradleWallThickness() ;
  par[2] -=  geom->GetCradleWallThickness() ;
  gMC->Gsvolu("PCRE", "TUBS", idtmed[798], par, 5) ;
  gMC->Gspos ("PCRE", 1, "PCRA", 0.0, 0.0, 0.0, 0, "ONLY") ; 

  for (i=0; i<2; i++) {
    z0 = (2*i-1) * (geom->GetOuterBoxSize(2) + geom->GetCradleWall(2)) / 2.0 ;
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
    z0 = (2*i-1) * ((geom->GetOuterBoxSize(2) + geom->GetCradleWheel(2)) / 2.0 +
                    geom->GetCradleWall(2));
    for (j=0; j<2; j++) {
      copy = 2*i + j;
      x0 = (2*j-1) * geom->GetDistanceBetwRails()  / 2.0 ;
      gMC->Gspos("PWHE", copy, "ALIC", x0, y0, z0, 0, "ONLY") ; 
    }
  }

}

//____________________________________________________________________________
Float_t AliPHOSv0::ZMin(void) const
{
  // Overall dimension of the PHOS (min)
  // Take it twice more than the PHOS module size

  AliPHOSGeometry * geom = GetGeometry() ; 

  return -geom->GetOuterBoxSize(2);
}

//____________________________________________________________________________
Float_t AliPHOSv0::ZMax(void) const
{
  // Overall dimension of the PHOS (max)
  // Take it twice more than the PHOS module size

  AliPHOSGeometry * geom = GetGeometry() ; 

  return  geom->GetOuterBoxSize(2);
}

//____________________________________________________________________________
void AliPHOSv0::Init(void)
{
  // Just prints an information message
  
  Int_t i;

  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" PHOS_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    
    
    // Here the PHOS initialisation code (if any!)

    AliPHOSGeometry * geom = GetGeometry() ; 

    if (geom!=0)  
      cout << "AliPHOS" << Version() << " : PHOS geometry intialized for " << geom->GetName() << endl ;
    else
      cout << "AliPHOS" << Version() << " : PHOS geometry initialization failed !" << endl ;   
    
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }  
}


