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
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"
#include "TRandom.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>

// --- AliRoot header files ---

#include "AliPHOSv0.h"
#include "AliPHOSHit.h"
#include "AliPHOSDigit.h"
#include "AliPHOSReconstructioner.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliPHOSv0)

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0()
{
  // ctor
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ; 
}

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(const char *name, const char *title):
  AliPHOS(name,title)
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)   
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fPinElectronicNoise = 0.010 ;
  fDigitThreshold      = 1. ;   // 1 GeV 

  // We do not want to save in TreeH the raw hits
  // fHits   = new TClonesArray("AliPHOSHit",100) ;
  // gAlice->AddHitList(fHits) ; 

  // But save the cumulated hits instead (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
  // and the TreeD at the end of the event (branch is set in FinishEvent() ).
  
  fTmpHits= new TClonesArray("AliPHOSHit",100) ;

  fNTmpHits = fNhits = 0 ;

  fDigits = new TClonesArray("AliPHOSDigit",100) ;


  fIshunt     =  1 ; // All hits are associated with primary particles
 
  // gets an instance of the geometry parameters class  
   
  fGeom =  AliPHOSGeometry::GetInstance(title, "") ; 

  if (fGeom->IsInitialized() ) 
    cout << "AliPHOSv0 : PHOS geometry intialized for " << fGeom->GetName() << endl ;
  else
   cout << "AliPHOSv0 : PHOS geometry initialization failed !" << endl ;   
}

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title):
  AliPHOS(name,title)
{
  // ctor : title is used to identify the layout
  //        GPS2 = 5 modules (EMC + PPSD)   
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fPinElectronicNoise = 0.010 ;

  // We do not want to save in TreeH the raw hits
  //fHits   = new TClonesArray("AliPHOSHit",100) ;

  fDigits = new TClonesArray("AliPHOSDigit",100) ;
  fTmpHits= new TClonesArray("AliPHOSHit",100) ;

  fNTmpHits = fNhits = 0 ;

  fIshunt     =  1 ; // All hits are associated with primary particles
 
  // gets an instance of the geometry parameters class  
  fGeom =  AliPHOSGeometry::GetInstance(title, "") ; 

  if (fGeom->IsInitialized() ) 
    cout << "AliPHOSv0 : PHOS geometry intialized for " << fGeom->GetName() << endl ;
  else
   cout << "AliPHOSv0 : PHOS geometry initialization failed !" << endl ;   

  // Defining the PHOS Reconstructioner
 
 fReconstructioner = Reconstructioner ;
}

//____________________________________________________________________________
AliPHOSv0::~AliPHOSv0()
{
  // dtor

  fTmpHits->Delete() ; 
  delete fTmpHits ;
  fTmpHits = 0 ; 

  fEmcRecPoints->Delete() ; 
  delete fEmcRecPoints ; 
  fEmcRecPoints = 0 ; 

  fPpsdRecPoints->Delete() ;
  delete fPpsdRecPoints ;
  fPpsdRecPoints = 0 ; 

  fTrackSegments->Delete() ; 
  delete fTrackSegments ;
  fTrackSegments = 0 ; 
}

//____________________________________________________________________________
void AliPHOSv0::AddHit(Int_t primary, Int_t Id, Float_t * hits)
{
  // Add a hit to the hit list.
  // A PHOS hit is the sum of all hits in a single crystal
  //   or in a single PPSD gas cell

  Int_t hitCounter ;
  TClonesArray &ltmphits = *fTmpHits ;
  AliPHOSHit *newHit ;
  AliPHOSHit *curHit ;
  //  AliPHOSHit *curHit2 ;
  Bool_t deja = kFALSE ;

  // In any case, fills the fTmpHit TClonesArray (with "accumulated hits")

  newHit = new AliPHOSHit(primary, Id, hits) ;

  // We do not want to save in TreeH the raw hits 
  //  TClonesArray &lhits = *fHits;

  for ( hitCounter = 0 ; hitCounter < fNTmpHits && !deja ; hitCounter++ ) {
    curHit = (AliPHOSHit*) ltmphits[hitCounter] ;
  if( *curHit == *newHit ) {
    *curHit = *curHit + *newHit ;
    deja = kTRUE ;
    }
  }
         
  if ( !deja ) {
    new(ltmphits[fNTmpHits]) AliPHOSHit(*newHit) ;
    fNTmpHits++ ;
  }

  // We do not want to save in TreeH the raw hits 
  //   new(lhits[fNhits]) AliPHOSHit(*newHit) ;    
  //   fNhits++ ;

  // Please note that the fTmpHits array must survive up to the
  // end of the events, so it does not appear e.g. in ResetHits() (
  // which is called at the end of each primary).  

  delete newHit;

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

  this->BuildGeometryforPHOS() ; 
  if ( ( strcmp(fGeom->GetName(), "GPS2" ) == 0 ) ) 
    this->BuildGeometryforPPSD() ;
  else
    cout << "AliPHOSv0::BuildGeometry : no charged particle identification system installed" << endl; 

}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforPHOS(void)
{
 // Build the PHOS-EMC geometry for the ROOT display

  const Int_t kColorPHOS = kRed ;
  const Int_t kColorXTAL = kBlue ;

  Double_t const kRADDEG = 180.0 / kPI ;
 
  new TBRIK( "OuterBox", "PHOS box", "void", fGeom->GetOuterBoxSize(0)/2, 
                                             fGeom->GetOuterBoxSize(1)/2, 
                                             fGeom->GetOuterBoxSize(2)/2 );

  // Textolit Wall box, position inside PHOS 
  
  new TBRIK( "TextolitBox", "PHOS Textolit box ", "void", fGeom->GetTextolitBoxSize(0)/2, 
                                                          fGeom->GetTextolitBoxSize(1)/2, 
                                                          fGeom->GetTextolitBoxSize(2)/2);

  // Polystyrene Foam Plate

  new TBRIK( "UpperFoamPlate", "PHOS Upper foam plate", "void", fGeom->GetTextolitBoxSize(0)/2, 
                                                                fGeom->GetSecondUpperPlateThickness()/2, 
                                                                fGeom->GetTextolitBoxSize(2)/2 ) ; 

  // Air Filled Box
 
  new TBRIK( "AirFilledBox", "PHOS air filled box", "void", fGeom->GetAirFilledBoxSize(0)/2, 
                                                            fGeom->GetAirFilledBoxSize(1)/2, 
                                                            fGeom->GetAirFilledBoxSize(2)/2 );

  // Crystals Box

  Float_t xtlX = fGeom->GetCrystalSize(0) ; 
  Float_t xtlY = fGeom->GetCrystalSize(1) ; 
  Float_t xtlZ = fGeom->GetCrystalSize(2) ; 

  Float_t xl =  fGeom->GetNPhi() * ( xtlX + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 + fGeom->GetModuleBoxThickness() ;
  Float_t yl =  ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 
             + fGeom->GetModuleBoxThickness() / 2.0 ;
  Float_t zl =  fGeom->GetNZ() * ( xtlZ + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 +  fGeom->GetModuleBoxThickness() ;
  
  new TBRIK( "CrystalsBox", "PHOS crystals box", "void", xl, yl, zl ) ;

// position PHOS into ALICE

  Float_t r = fGeom->GetIPtoOuterCoverDistance() + fGeom->GetOuterBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  Float_t pphi =  TMath::ATan( fGeom->GetOuterBoxSize(0)  / ( 2.0 * fGeom->GetIPtoOuterCoverDistance() ) ) ;
  pphi *= kRADDEG ;
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
 
  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) { 
   Float_t angle = pphi * 2 * ( i - fGeom->GetNModules() / 2.0 - 0.5 ) ;
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
   y = ( fGeom->GetOuterBoxThickness(1) -  fGeom->GetUpperPlateThickness() ) / 2.  ;
   sprintf(nodename,"%s%d", "TexBox", i) ;  
   TNode * textolitboxnode = new TNode(nodename, nodename, "TextolitBox", 0, y, 0) ; 
   textolitboxnode->SetLineColor(kColorPHOS) ;
   fNodes->Add(textolitboxnode) ;
   // upper foam plate inside outre box
   outerboxnode->cd() ; 
   sprintf(nodename, "%s%d", "UFPlate", i) ;
   y =  ( fGeom->GetTextolitBoxSize(1) - fGeom->GetSecondUpperPlateThickness() ) / 2.0 ;
   TNode * upperfoamplatenode = new TNode(nodename, nodename, "UpperFoamPlate", 0, y, 0) ; 
   upperfoamplatenode->SetLineColor(kColorPHOS) ;
   fNodes->Add(upperfoamplatenode) ;  
   // air filled box inside textolit box (not drawn)
   textolitboxnode->cd();
   y = ( fGeom->GetTextolitBoxSize(1) - fGeom->GetAirFilledBoxSize(1) ) / 2.0 -  fGeom->GetSecondUpperPlateThickness() ;
   sprintf(nodename, "%s%d", "AFBox", i) ;
   TNode * airfilledboxnode = new TNode(nodename, nodename, "AirFilledBox", 0, y, 0) ; 
   fNodes->Add(airfilledboxnode) ; 
   // crystals box inside air filled box
   airfilledboxnode->cd() ; 
   y = fGeom->GetAirFilledBoxSize(1) / 2.0 - yl 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness() 
       -  fGeom->GetUpperPlateThickness() -  fGeom->GetSecondUpperPlateThickness() ) ; 
   sprintf(nodename, "%s%d", "XTBox", i) ; 
   TNode * crystalsboxnode = new TNode(nodename, nodename, "CrystalsBox", 0, y, 0) ;    
   crystalsboxnode->SetLineColor(kColorXTAL) ; 
   fNodes->Add(crystalsboxnode) ; 
  }
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

  // Box for a full PHOS module

  new TBRIK( "PPSDBox", "PPSD box", "void",  fGeom->GetPPSDBoxSize(0)/2, 
                                             fGeom->GetPPSDBoxSize(1)/2, 
	                                     fGeom->GetPPSDBoxSize(2)/2 );

  // Box containing one micromegas module 

  new TBRIK( "PPSDModule", "PPSD module", "void",  fGeom->GetPPSDModuleSize(0)/2, 
                                                   fGeom->GetPPSDModuleSize(1)/2, 
	                                           fGeom->GetPPSDModuleSize(2)/2 );
 // top lid

  new TBRIK ( "TopLid", "Micromegas top lid", "void",  fGeom->GetPPSDModuleSize(0)/2,
                                                       fGeom->GetLidThickness()/2,
                                                       fGeom->GetPPSDModuleSize(2)/2 ) ; 
 // composite panel (top and bottom)

  new TBRIK ( "TopPanel", "Composite top panel", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
                                                            fGeom->GetCompositeThickness()/2,
                                                          ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ;  
  
  new TBRIK ( "BottomPanel", "Composite bottom panel", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
                                                                  fGeom->GetCompositeThickness()/2,
                                                                ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ; 
 // gas gap (conversion and avalanche)

  new TBRIK ( "GasGap", "gas gap", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
	                                    ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() )/2,
                                            ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ; 

 // anode and cathode 

  new TBRIK ( "Anode", "Anode", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
                                           fGeom->GetAnodeThickness()/2,
                                         ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ; 

  new TBRIK ( "Cathode", "Cathode", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
                                               fGeom->GetCathodeThickness()/2,
                                             ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ; 
 // PC  

  new TBRIK ( "PCBoard", "Printed Circuit", "void",  ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() )/2,
                                                       fGeom->GetPCThickness()/2,
                                                     ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() )/2 ) ; 
 // Gap between Lead and top micromegas

  new TBRIK ( "LeadToM", "Air Gap top", "void", fGeom->GetPPSDBoxSize(0)/2,
                                                fGeom->GetMicro1ToLeadGap()/2,
                                                fGeom->GetPPSDBoxSize(2)/2  ) ;  
 
// Gap between Lead and bottom micromegas

  new TBRIK ( "MToLead", "Air Gap bottom", "void", fGeom->GetPPSDBoxSize(0)/2,
                                                   fGeom->GetLeadToMicro2Gap()/2,
                                                   fGeom->GetPPSDBoxSize(2)/2  ) ; 
 // Lead converter
   
  new TBRIK ( "Lead", "Lead converter", "void", fGeom->GetPPSDBoxSize(0)/2,
                                                fGeom->GetLeadConverterThickness()/2,
                                                fGeom->GetPPSDBoxSize(2)/2  ) ; 

     // position PPSD into ALICE

  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  Float_t r = fGeom->GetIPtoTopLidDistance() + fGeom->GetPPSDBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
 
  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) { // the number of PHOS modules
    Float_t angle = fGeom->GetPHOSAngle(i) ;
    sprintf(rotname, "%s%d", "rotg", number++) ;
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
    x = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. ;  
    {
      for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
	Float_t z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2. ;
	TNode * micro1node ; 
	for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	  y = ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas1Thickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Mic1", i, iphi, iz) ;
	  micro1node  = new TNode(nodename, nodename, "PPSDModule", x, y, z) ;
	  micro1node->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(micro1node) ; 
	  // inside top micromegas
	  micro1node->cd() ; 
	  //      a. top lid
	  y = ( fGeom->GetMicromegas1Thickness() - fGeom->GetLidThickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Lid", i, iphi, iz) ;
	  TNode * toplidnode = new TNode(nodename, nodename, "TopLid", 0, y, 0) ;
	  toplidnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(toplidnode) ; 
	  //      b. composite panel
	  y = y - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompU", i, iphi, iz) ;
	  TNode * compupnode = new TNode(nodename, nodename, "TopPanel", 0, y, 0) ;
	  compupnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compupnode) ; 
	  //      c. anode
	  y = y - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Ano", i, iphi, iz) ;
	  TNode * anodenode = new TNode(nodename, nodename, "Anode", 0, y, 0) ;
	  anodenode->SetLineColor(kColorPHOS) ;  
	  fNodes->Add(anodenode) ; 
	  //      d.  gas 
	  y = y - fGeom->GetAnodeThickness() / 2. - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "GGap", i, iphi, iz) ;
	  TNode * ggapnode = new TNode(nodename, nodename, "GasGap", 0, y, 0) ;
	  ggapnode->SetLineColor(kColorGas) ;  
	  fNodes->Add(ggapnode) ;          
	  //      f. cathode
	  y = y - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Cathode", i, iphi, iz) ;
	  TNode * cathodenode = new TNode(nodename, nodename, "Cathode", 0, y, 0) ;
	  cathodenode->SetLineColor(kColorPHOS) ;  
	  fNodes->Add(cathodenode) ;        
	  //      g. printed circuit
	  y = y - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "PC", i, iphi, iz) ;
	  TNode * pcnode = new TNode(nodename, nodename, "PCBoard", 0, y, 0) ;
	  pcnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(pcnode) ;        
	  //      h. composite panel
	  y = y - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "CompDown", i, iphi, iz) ;
	  TNode * compdownnode = new TNode(nodename, nodename, "BottomPanel", 0, y, 0) ;
	  compdownnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compdownnode) ;   
	  z = z - fGeom->GetPPSDModuleSize(2) ;
	  ppsdboxnode->cd() ;
	} // end of Z module loop     
	x = x -  fGeom->GetPPSDModuleSize(0) ; 
	ppsdboxnode->cd() ;
      } // end of phi module loop
    }
    //   2. air gap      
    ppsdboxnode->cd() ;
    y = ( fGeom->GetPPSDBoxSize(1) - 2 * fGeom->GetMicromegas1Thickness() - fGeom->GetMicro1ToLeadGap() ) / 2. ; 
    sprintf(nodename, "%s%d", "GapUp", i) ;
    TNode * gapupnode = new TNode(nodename, nodename, "LeadToM", 0, y, 0) ;
    gapupnode->SetLineColor(kColorAir) ;  
    fNodes->Add(gapupnode) ;        
    //   3. lead converter
    y = y - fGeom->GetMicro1ToLeadGap() / 2. - fGeom->GetLeadConverterThickness() / 2. ; 
    sprintf(nodename, "%s%d", "LeadC", i) ;
    TNode * leadcnode = new TNode(nodename, nodename, "Lead", 0, y, 0) ;
    leadcnode->SetLineColor(kColorPPSD) ;  
    fNodes->Add(leadcnode) ;        
    //   4. air gap
    y = y - fGeom->GetLeadConverterThickness() / 2. - fGeom->GetLeadToMicro2Gap()  / 2. ; 
    sprintf(nodename, "%s%d", "GapDown", i) ;
    TNode * gapdownnode = new TNode(nodename, nodename, "MToLead", 0, y, 0) ;
    gapdownnode->SetLineColor(kColorAir) ;  
    fNodes->Add(gapdownnode) ;        
    //    5.  fNumberOfModulesPhi x fNumberOfModulesZ bottom micromegas
    x = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. - fGeom->GetPhiDisplacement() ;  
    {
      for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { 
	Float_t z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2.  - fGeom->GetZDisplacement() ;;
	TNode * micro2node ; 
	for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { 
	  y = - ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas2Thickness() ) / 2. ; 
	  sprintf(nodename, "%s%d%d%d", "Mic2", i, iphi, iz) ;
	  micro2node  = new TNode(nodename, nodename, "PPSDModule", x, y, z) ;
	  micro2node->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(micro2node) ; 
	  // inside bottom micromegas
	  micro2node->cd() ; 
	  //      a. top lid
	  y = ( fGeom->GetMicromegas2Thickness() - fGeom->GetLidThickness() ) / 2. ; 
	  sprintf(nodename, "%s%d", "Lidb", i) ;
	  TNode * toplidbnode = new TNode(nodename, nodename, "TopLid", 0, y, 0) ;
	  toplidbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(toplidbnode) ; 
	  //      b. composite panel
	  y = y - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 
	  sprintf(nodename, "%s%d", "CompUb", i) ;
	  TNode * compupbnode = new TNode(nodename, nodename, "TopPanel", 0, y, 0) ;
	  compupbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compupbnode) ; 
	  //      c. anode
	  y = y - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "Anob", i) ;
	  TNode * anodebnode = new TNode(nodename, nodename, "Anode", 0, y, 0) ;
	  anodebnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(anodebnode) ; 
	  //      d. conversion gas
	  y = y - fGeom->GetAnodeThickness() / 2. - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() )  / 2. ; 
	  sprintf(nodename, "%s%d", "GGapb", i) ;
	  TNode * ggapbnode = new TNode(nodename, nodename, "GasGap", 0, y, 0) ;
	  ggapbnode->SetLineColor(kColorGas) ;  
	  fNodes->Add(ggapbnode) ;           
	  //      f. cathode
	  y = y - ( fGeom->GetConversionGap() + fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "Cathodeb", i) ;
	  TNode * cathodebnode = new TNode(nodename, nodename, "Cathode", 0, y, 0) ;
	  cathodebnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(cathodebnode) ;        
	  //      g. printed circuit
	  y = y - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "PCb", i) ;
	  TNode * pcbnode = new TNode(nodename, nodename, "PCBoard", 0, y, 0) ;
	  pcbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(pcbnode) ;        
	  //      h. composite pane
	  y = y - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "CompDownb", i) ;
	  TNode * compdownbnode = new TNode(nodename, nodename, "BottomPanel", 0, y, 0) ;
	  compdownbnode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(compdownbnode) ;        
       	  z = z - fGeom->GetPPSDModuleSize(2) ;
	  ppsdboxnode->cd() ;
	} // end of Z module loop     
	x = x -  fGeom->GetPPSDModuleSize(0) ; 
	ppsdboxnode->cd() ;
      } // end of phi module loop
    }
  } // PHOS modules
 
  delete rotname ;  
  delete nodename ; 

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
  // Get pointer to the array containing media indeces
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;

  Float_t bigbox[3] ; 
  bigbox[0] =   fGeom->GetOuterBoxSize(0) / 2.0 ;
  bigbox[1] = ( fGeom->GetOuterBoxSize(1) + fGeom->GetPPSDBoxSize(1) ) / 2.0 ;
  bigbox[2] =   fGeom->GetOuterBoxSize(2) / 2.0 ;
  
  gMC->Gsvolu("PHOS", "BOX ", idtmed[798], bigbox, 3) ;
  
  this->CreateGeometryforPHOS() ; 
  if ( strcmp( fGeom->GetName(), "GPS2") == 0  ) 
    this->CreateGeometryforPPSD() ;
  else
    cout << "AliPHOSv0::CreateGeometry : no charged particle identification system installed" << endl; 
  
  // --- Position  PHOS mdules in ALICE setup ---
  
  Int_t idrotm[99] ;
  Double_t const kRADDEG = 180.0 / kPI ;
  
  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) {
    
    Float_t angle = fGeom->GetPHOSAngle(i) ;
    AliMatrix(idrotm[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0) ;
 
    Float_t r = fGeom->GetIPtoOuterCoverDistance() + ( fGeom->GetOuterBoxSize(1) + fGeom->GetPPSDBoxSize(1) ) / 2.0 ;

    Float_t xP1 = r * TMath::Sin( angle / kRADDEG ) ;
    Float_t yP1 = -r * TMath::Cos( angle / kRADDEG ) ;

    gMC->Gspos("PHOS", i, "ALIC", xP1, yP1, 0.0, idrotm[i-1], "ONLY") ;
 
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

  // ---
  // --- Define PHOS box volume, fPUFPill with thermo insulating foam ---
  // --- Foam Thermo Insulating outer cover dimensions ---
  // --- Put it in bigbox = PHOS

  Float_t dphos[3] ; 
  dphos[0] =  fGeom->GetOuterBoxSize(0) / 2.0 ;
  dphos[1] =  fGeom->GetOuterBoxSize(1) / 2.0 ;
  dphos[2] =  fGeom->GetOuterBoxSize(2) / 2.0 ;

  gMC->Gsvolu("EMCA", "BOX ", idtmed[706], dphos, 3) ;

  Float_t yO =  - fGeom->GetPPSDBoxSize(1)  / 2.0 ;

  gMC->Gspos("EMCA", 1, "PHOS", 0.0, yO, 0.0, 0, "ONLY") ; 

  // ---
  // --- Define Textolit Wall box, position inside EMCA ---
  // --- Textolit Wall box dimentions ---
 
 
  Float_t dptxw[3];
  dptxw[0] = fGeom->GetTextolitBoxSize(0) / 2.0 ;
  dptxw[1] = fGeom->GetTextolitBoxSize(1) / 2.0 ;
  dptxw[2] = fGeom->GetTextolitBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXW", "BOX ", idtmed[707], dptxw, 3);

  yO =   (  fGeom->GetOuterBoxThickness(1) -   fGeom->GetUpperPlateThickness() ) / 2.  ;
   
  gMC->Gspos("PTXW", 1, "EMCA", 0.0, yO, 0.0, 0, "ONLY") ;

  // --- 
  // --- Define Upper Polystyrene Foam Plate, place inside PTXW ---
  // --- immediately below Foam Thermo Insulation Upper plate ---

  // --- Upper Polystyrene Foam plate thickness ---
 
  Float_t  dpufp[3] ;
  dpufp[0] = fGeom->GetTextolitBoxSize(0) / 2.0 ; 
  dpufp[1] = fGeom->GetSecondUpperPlateThickness() / 2. ;
  dpufp[2] = fGeom->GetTextolitBoxSize(2) /2.0 ; 

  gMC->Gsvolu("PUFP", "BOX ", idtmed[703], dpufp, 3) ;
  
  yO = ( fGeom->GetTextolitBoxSize(1) -  fGeom->GetSecondUpperPlateThickness() ) / 2.0 ;
  
  gMC->Gspos("PUFP", 1, "PTXW", 0.0, yO, 0.0, 0, "ONLY") ;
  
  // ---
  // --- Define air-filled box, place inside PTXW ---
  // --- Inner AIR volume dimensions ---
 

  Float_t  dpair[3] ;
  dpair[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  dpair[1] = fGeom->GetAirFilledBoxSize(1) / 2.0 ;
  dpair[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PAIR", "BOX ", idtmed[798], dpair, 3) ;
  
  yO = ( fGeom->GetTextolitBoxSize(1) -  fGeom->GetAirFilledBoxSize(1) ) / 2.0 -   fGeom->GetSecondUpperPlateThickness() ;
  
  gMC->Gspos("PAIR", 1, "PTXW", 0.0, yO, 0.0, 0, "ONLY") ;

// --- Dimensions of PbWO4 crystal ---

  Float_t xtlX =  fGeom->GetCrystalSize(0) ; 
  Float_t xtlY =  fGeom->GetCrystalSize(1) ; 
  Float_t xtlZ =  fGeom->GetCrystalSize(2) ; 

  Float_t dptcb[3] ;  
  dptcb[0] =  fGeom->GetNPhi() * ( xtlX + 2 *  fGeom->GetGapBetweenCrystals() ) / 2.0 + fGeom->GetModuleBoxThickness() ;
  dptcb[1] = ( xtlY +  fGeom->GetCrystalSupportHeight() +  fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 
             + fGeom->GetModuleBoxThickness() / 2.0 ;
  dptcb[2] = fGeom->GetNZ() * ( xtlZ + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 +  fGeom->GetModuleBoxThickness() ;
  
  gMC->Gsvolu("PTCB", "BOX ", idtmed[706], dptcb, 3) ;

  yO =  fGeom->GetAirFilledBoxSize(1) / 2.0 - dptcb[1] 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness() 
       -  fGeom->GetUpperPlateThickness() -  fGeom->GetSecondUpperPlateThickness() ) ;
  
  gMC->Gspos("PTCB", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Crystal BLock filled with air, position it inside PTCB ---
  Float_t dpcbl[3] ; 
  
  dpcbl[0] = fGeom->GetNPhi() * ( xtlX + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 ;
  dpcbl[1] = ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  dpcbl[2] = fGeom->GetNZ() * ( xtlZ + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 ;
  
  gMC->Gsvolu("PCBL", "BOX ", idtmed[798], dpcbl, 3) ;
  
  // --- Divide PCBL in X (phi) and Z directions --
  gMC->Gsdvn("PROW", "PCBL", Int_t (fGeom->GetNPhi()), 1) ;
  gMC->Gsdvn("PCEL", "PROW", Int_t (fGeom->GetNZ()), 3) ;

  yO = -fGeom->GetModuleBoxThickness() / 2.0 ;
  
  gMC->Gspos("PCBL", 1, "PTCB", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define STeel (actually, it's titanium) Cover volume, place inside PCEL
  Float_t  dpstc[3] ; 
  
  dpstc[0] = ( xtlX + 2 * fGeom->GetCrystalWrapThickness() ) / 2.0 ;
  dpstc[1] = ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  dpstc[2] = ( xtlZ + 2 * fGeom->GetCrystalWrapThickness()  + 2 *  fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  
  gMC->Gsvolu("PSTC", "BOX ", idtmed[704], dpstc, 3) ;

  gMC->Gspos("PSTC", 1, "PCEL", 0.0, 0.0, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Tyvek volume, place inside PSTC ---
  Float_t  dppap[3] ;

  dppap[0] = xtlX / 2.0 + fGeom->GetCrystalWrapThickness() ;
  dppap[1] = ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 ;
  dppap[2] = xtlZ / 2.0 + fGeom->GetCrystalWrapThickness() ;
  
  gMC->Gsvolu("PPAP", "BOX ", idtmed[702], dppap, 3) ;
  
  yO = ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 
              - ( xtlY +  fGeom->GetCrystalSupportHeight() +  fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
   
  gMC->Gspos("PPAP", 1, "PSTC", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PbWO4 crystal volume, place inside PPAP ---
  Float_t  dpxtl[3] ; 

  dpxtl[0] = xtlX / 2.0 ;
  dpxtl[1] = xtlY / 2.0 ;
  dpxtl[2] = xtlZ / 2.0 ;
  
  gMC->Gsvolu("PXTL", "BOX ", idtmed[699], dpxtl, 3) ;

  yO = ( xtlY + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 - xtlY / 2.0 - fGeom->GetCrystalWrapThickness() ;
  
  gMC->Gspos("PXTL", 1, "PPAP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define crystal support volume, place inside PPAP ---
  Float_t dpsup[3] ; 

  dpsup[0] = xtlX / 2.0 + fGeom->GetCrystalWrapThickness()  ;
  dpsup[1] = fGeom->GetCrystalSupportHeight() / 2.0 ;
  dpsup[2] = xtlZ / 2.0 +  fGeom->GetCrystalWrapThickness() ;

  gMC->Gsvolu("PSUP", "BOX ", idtmed[798], dpsup, 3) ;

  yO =  fGeom->GetCrystalSupportHeight() / 2.0 - ( xtlY +  fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 ;

  gMC->Gspos("PSUP", 1, "PPAP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PIN-diode volume and position it inside crystal support ---
  // --- right behind PbWO4 crystal

  // --- PIN-diode dimensions ---

 
  Float_t dppin[3] ;
  dppin[0] = fGeom->GetPinDiodeSize(0) / 2.0 ;
  dppin[1] = fGeom->GetPinDiodeSize(1) / 2.0 ;
  dppin[2] = fGeom->GetPinDiodeSize(2) / 2.0 ;
 
  gMC->Gsvolu("PPIN", "BOX ", idtmed[705], dppin, 3) ;
 
  yO = fGeom->GetCrystalSupportHeight() / 2.0 - fGeom->GetPinDiodeSize(1) / 2.0 ;
 
  gMC->Gspos("PPIN", 1, "PSUP", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Upper Cooling Panel, place it on top of PTCB ---
  Float_t dpucp[3] ;
 // --- Upper Cooling Plate thickness ---
 
  dpucp[0] = dptcb[0] ;
  dpucp[1] = fGeom->GetUpperCoolingPlateThickness() ;
  dpucp[2] = dptcb[2] ;
  
  gMC->Gsvolu("PUCP", "BOX ", idtmed[701], dpucp,3) ;
  
  yO = (  fGeom->GetAirFilledBoxSize(1) -  fGeom->GetUpperCoolingPlateThickness() ) / 2. 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness()
           - fGeom->GetUpperPlateThickness() - fGeom->GetSecondUpperPlateThickness() - fGeom->GetUpperCoolingPlateThickness() ) ; 
  
  gMC->Gspos("PUCP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Al Support Plate, position it inside PAIR ---
  // --- right beneath PTCB ---
 // --- Al Support Plate thickness ---
 
  Float_t dpasp[3] ;
  dpasp[0] =  fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  dpasp[1] = fGeom->GetSupportPlateThickness() / 2.0 ;
  dpasp[2] =  fGeom->GetAirFilledBoxSize(2) / 2.0 ;
  
  gMC->Gsvolu("PASP", "BOX ", idtmed[701], dpasp, 3) ;
  
  yO = (  fGeom->GetAirFilledBoxSize(1) - fGeom->GetSupportPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance()
           - fGeom->GetUpperPlateThickness() - fGeom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 ) ;
  
  gMC->Gspos("PASP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Thermo Insulating Plate, position it inside PAIR ---
  // --- right beneath PASP ---
  // --- Lower Thermo Insulating Plate thickness ---
  
  Float_t dptip[3] ;
  dptip[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  dptip[1] = fGeom->GetLowerThermoPlateThickness() / 2.0 ;
  dptip[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTIP", "BOX ", idtmed[706], dptip, 3) ;

  yO =  ( fGeom->GetAirFilledBoxSize(1) - fGeom->GetLowerThermoPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetUpperPlateThickness() 
            - fGeom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 + fGeom->GetSupportPlateThickness() ) ;

  gMC->Gspos("PTIP", 1, "PAIR", 0.0, yO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Textolit Plate, position it inside PAIR ---
  // --- right beneath PTIP ---
  // --- Lower Textolit Plate thickness ---
 
  Float_t dptxp[3] ;
  dptxp[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  dptxp[1] = fGeom->GetLowerTextolitPlateThickness() / 2.0 ;
  dptxp[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXP", "BOX ", idtmed[707], dptxp, 3) ;

  yO =  ( fGeom->GetAirFilledBoxSize(1) - fGeom->GetLowerTextolitPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetUpperPlateThickness() 
            - fGeom->GetSecondUpperPlateThickness() + dpcbl[1] * 2 + fGeom->GetSupportPlateThickness() 
            +  fGeom->GetLowerThermoPlateThickness() ) ;

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
  
  // The box containing all ppsd's for one PHOS module filled with air 
  Float_t ppsd[3] ; 
  ppsd[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;  
  ppsd[1] = fGeom->GetPPSDBoxSize(1) / 2.0 ; 
  ppsd[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PPSD", "BOX ", idtmed[798], ppsd, 3) ;

  Float_t yO =  fGeom->GetOuterBoxSize(1) / 2.0 ;

  gMC->Gspos("PPSD", 1, "PHOS", 0.0, yO, 0.0, 0, "ONLY") ; 

  // Now we build a micromegas module
  // The box containing the whole module filled with epoxy (FR4)

  Float_t mppsd[3] ;  
  mppsd[0] = fGeom->GetPPSDModuleSize(0) / 2.0 ;  
  mppsd[1] = fGeom->GetPPSDModuleSize(1) / 2.0 ;  
  mppsd[2] = fGeom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("MPPS", "BOX ", idtmed[708], mppsd, 3) ;  
 
  // Inside mppsd :
  // 1. The Top Lid made of epoxy (FR4) 

  Float_t tlppsd[3] ; 
  tlppsd[0] = fGeom->GetPPSDModuleSize(0) / 2.0 ; 
  tlppsd[1] = fGeom->GetLidThickness() / 2.0 ;
  tlppsd[2] = fGeom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("TLPS", "BOX ", idtmed[708], tlppsd, 3) ; 

  Float_t  y0 = ( fGeom->GetMicromegas1Thickness() - fGeom->GetLidThickness() ) / 2. ; 

  gMC->Gspos("TLPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 
 
  // 2. the upper panel made of composite material

  Float_t upppsd[3] ; 
  upppsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  upppsd[1] = fGeom->GetCompositeThickness() / 2.0 ;
  upppsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
 
  gMC->Gsvolu("UPPS", "BOX ", idtmed[709], upppsd, 3) ; 
  
  y0 = y0 - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 

  gMC->Gspos("UPPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 3. the anode made of Copper
  
  Float_t anppsd[3] ; 
  anppsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ; 
  anppsd[1] = fGeom->GetAnodeThickness() / 2.0 ; 
  anppsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0  ; 

  gMC->Gsvolu("ANPS", "BOX ", idtmed[710], anppsd, 3) ; 
  
  y0 = y0 - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
  
  gMC->Gspos("ANPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 4. the conversion gap + avalanche gap filled with gas

  Float_t ggppsd[3] ; 
  ggppsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  ggppsd[1] = ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2.0 ; 
  ggppsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("GGPS", "BOX ", idtmed[715], ggppsd, 3) ; 
  
  // --- Divide GGPP in X (phi) and Z directions --
  gMC->Gsdvn("GROW", "GGPS", fGeom->GetNumberOfPadsPhi(), 1) ;
  gMC->Gsdvn("GCEL", "GROW", fGeom->GetNumberOfPadsZ() , 3) ;

  y0 = y0 - fGeom->GetAnodeThickness() / 2.  - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. ; 

  gMC->Gspos("GGPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 


  // 6. the cathode made of Copper

  Float_t cappsd[3] ;
  cappsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  cappsd[1] = fGeom->GetCathodeThickness() / 2.0 ; 
  cappsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0  ;

  gMC->Gsvolu("CAPS", "BOX ", idtmed[710], cappsd, 3) ; 

  y0 = y0 - ( fGeom->GetAvalancheGap() +  fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 

  gMC->Gspos("CAPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 7. the printed circuit made of G10       

  Float_t pcppsd[3] ; 
  pcppsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2,.0 ; 
  pcppsd[1] = fGeom->GetPCThickness() / 2.0 ; 
  pcppsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("PCPS", "BOX ", idtmed[711], cappsd, 3) ; 

  y0 = y0 - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 

  gMC->Gspos("PCPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 8. the lower panel made of composite material
						    
  Float_t lpppsd[3] ; 
  lpppsd[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ; 
  lpppsd[1] = fGeom->GetCompositeThickness() / 2.0 ; 
  lpppsd[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("LPPS", "BOX ", idtmed[709], lpppsd, 3) ; 
 
  y0 = y0 - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 

  gMC->Gspos("LPPS", 1, "MPPS", 0.0, y0, 0.0, 0, "ONLY") ; 

  // Position the  fNumberOfModulesPhi x fNumberOfModulesZ modules (mppsd) inside PPSD to cover a PHOS module
  // the top and bottom one's (which are assumed identical) :

   Float_t yt = ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas1Thickness() ) / 2. ; 
   Float_t yb = - ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas2Thickness() ) / 2. ; 

   Int_t copyNumbertop = 0 ; 
   Int_t copyNumberbot = fGeom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() ; 

   Float_t x  = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. ;  

   for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
      Float_t z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2. ;

      for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	gMC->Gspos("MPPS", ++copyNumbertop, "PPSD", x, yt, z, 0, "ONLY") ;
	gMC->Gspos("MPPS", ++copyNumberbot, "PPSD", x, yb, z, 0, "ONLY") ; 
	z = z - fGeom->GetPPSDModuleSize(2) ;
      } // end of Z module loop   
      x = x -  fGeom->GetPPSDModuleSize(0) ; 
    } // end of phi module loop

   // The Lead converter between two air gaps
   // 1. Upper air gap

   Float_t uappsd[3] ;
   uappsd[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;
   uappsd[1] = fGeom->GetMicro1ToLeadGap() / 2.0 ; 
   uappsd[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("UAPPSD", "BOX ", idtmed[798], uappsd, 3) ; 

  y0 = ( fGeom->GetPPSDBoxSize(1) - 2 * fGeom->GetMicromegas1Thickness() - fGeom->GetMicro1ToLeadGap() ) / 2. ; 

  gMC->Gspos("UAPPSD", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 

   // 2. Lead converter
 
  Float_t lcppsd[3] ; 
  lcppsd[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;
  lcppsd[1] = fGeom->GetLeadConverterThickness() / 2.0 ; 
  lcppsd[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;
 
  gMC->Gsvolu("LCPPSD", "BOX ", idtmed[712], lcppsd, 3) ; 
  
  y0 = y0 - fGeom->GetMicro1ToLeadGap() / 2. - fGeom->GetLeadConverterThickness() / 2. ; 

  gMC->Gspos("LCPPSD", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 

  // 3. Lower air gap

  Float_t lappsd[3] ; 
  lappsd[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ; 
  lappsd[1] = fGeom->GetLeadToMicro2Gap() / 2.0 ; 
  lappsd[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("LAPPSD", "BOX ", idtmed[798], lappsd, 3) ; 
    
  y0 = y0 - fGeom->GetLeadConverterThickness() / 2. - fGeom->GetLeadToMicro2Gap()  / 2. ; 
  
  gMC->Gspos("LAPPSD", 1, "PPSD", 0.0, y0, 0.0, 0, "ONLY") ; 
   
}

//___________________________________________________________________________
Int_t AliPHOSv0::Digitize(Float_t Energy)
{
  // Applies the energy calibration
  
  Float_t fB = 100000000. ;
  Float_t fA = 0. ;
  Int_t chan = Int_t(fA + Energy*fB ) ;
  return chan ;
}

//___________________________________________________________________________
void AliPHOSv0::FinishEvent()
{
  // Makes the digits from the sum of summed hit in a single crystal or PPSD gas cell
  // Adds to the energy the electronic noise
  // Keeps digits with energy above fDigitThreshold

  // Save the cumulated hits instead of raw hits (need to create the branch myself)
  // It is put in the Digit Tree because the TreeH is filled after each primary
  // and the TreeD at the end of the event.
  if ( fTmpHits && gAlice->TreeD() ) {
    char branchname[10] ;
    sprintf(branchname, "%sCH", GetName()) ;
    gAlice->TreeD()->Branch(branchname, &fTmpHits, fBufferSize) ; 
  } else 
    cout << "AliPHOSv0::AliPHOSv0: Failed to create branch PHOSCH in TreeD " << endl ;  
  
  Int_t i ;
  Int_t relid[4];
  Int_t j ; 
  TClonesArray &lDigits = *fDigits ;
  AliPHOSHit  * hit ;
  AliPHOSDigit * newdigit ;
  AliPHOSDigit * curdigit ;
  Bool_t deja = kFALSE ; 
  
  for ( i = 0 ; i < fNTmpHits ; i++ ) {
    hit = (AliPHOSHit*)fTmpHits->At(i) ;
    newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
    deja =kFALSE ;
    for ( j = 0 ; j < fNdigits ;  j++) { 
      curdigit = (AliPHOSDigit*) lDigits[j] ;
      if ( *curdigit == *newdigit) {
	*curdigit = *curdigit + *newdigit ; 
	deja = kTRUE ; 
      }
    }
    if ( !deja ) {
      new(lDigits[fNdigits]) AliPHOSDigit(* newdigit) ;
      fNdigits++ ;  
    }
 
    delete newdigit ;    
  } 
  
  // Noise induced by the PIN diode of the PbWO crystals

  Float_t energyandnoise ;
  for ( i = 0 ; i < fNdigits ; i++ ) {
    newdigit =  (AliPHOSDigit * ) fDigits->At(i) ;
    fGeom->AbsToRelNumbering(newdigit->GetId(), relid) ;

    if (relid[1]==0){   // Digits belong to EMC (PbW0_4 crystals)
      energyandnoise = newdigit->GetAmp() + Digitize(gRandom->Gaus(0., fPinElectronicNoise)) ;

      if (energyandnoise < 0 ) 
	energyandnoise = 0 ;

      if ( newdigit->GetAmp() < fDigitThreshold ) // if threshold not surpassed, remove digit from list
	fDigits->RemoveAt(i) ; 
    }
  }
  
  fDigits->Compress() ;  

  fNdigits =  fDigits->GetEntries() ; 
  for (i = 0 ; i < fNdigits ; i++) { 
    newdigit = (AliPHOSDigit *) fDigits->At(i) ; 
    newdigit->SetIndexInList(i) ; 
  }

}

//____________________________________________________________________________
void AliPHOSv0::Init(void)
{
  // Just prints an information message
  
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

//___________________________________________________________________________
void AliPHOSv0::MakeBranch(Option_t* opt)
{  
  // Create new branche in the current Root Tree in the digit Tree

  AliDetector::MakeBranch(opt) ;
  
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  char *cdD = strstr(opt,"D");
  
  if (fDigits && gAlice->TreeD() && cdD) {
    gAlice->TreeD()->Branch(branchname, &fDigits, fBufferSize);
  }
}

//____________________________________________________________________________
RecPointsList * AliPHOSv0::PpsdRecPoints(Int_t evt) 
{
  // returns the pointer to the PPSD RecPoints list
  // if the list is empty, get it from TreeR on the disk file

  RecPointsList * rv = 0 ; 

  if ( fPpsdRecPoints ) 
    rv = fPpsdRecPoints ; 

  else {
    fPpsdRecPoints = new TClonesArray("AliPHOSPpsdRecPoint", 100) ; 
    gAlice->GetEvent(evt) ; 
    TTree * fReconstruct = gAlice->TreeR() ; 
    fReconstruct->SetBranchAddress( "PHOSPpsdRP", &fPpsdRecPoints) ;
    fReconstruct->GetEvent(0) ;
    rv =  fPpsdRecPoints ;
  }
  
  fPpsdRecPoints->Expand( fPpsdRecPoints->GetEntries() ) ; 
    
  return rv ; 
  
}

//_____________________________________________________________________________
void AliPHOSv0::Reconstruction(AliPHOSReconstructioner * Reconstructioner)
{ 
  // 1. Reinitializes the existing RecPoint, TrackSegment, and RecParticles Lists and 
  // 2. Creates TreeR wit a branch for each list
  // 3. Steers the reconstruction processes
  // 4. Saves the 3 lists in TreeR
  // 5. Write the Tree to File
  
  fReconstructioner = Reconstructioner ;
  
  char branchname[10] ;

  // 1.

  gAlice->MakeTree("R") ; 
  Int_t splitlevel = 0 ; 
  
  if (fEmcRecPoints) { 
    fEmcRecPoints->Delete() ; 
    delete fEmcRecPoints ;
    fEmcRecPoints = 0 ; 
  }

  //  fEmcRecPoints= new RecPointsList("AliPHOSEmcRecPoint", 100) ; if TClonesArray
  fEmcRecPoints= new RecPointsList(100) ; 

  if ( fEmcRecPoints && gAlice->TreeR() ) {
    sprintf(branchname,"%sEmcRP",GetName()) ;
    
    // gAlice->TreeR()->Branch(branchname, &fEmcRecPoints, fBufferSize); if TClonesArray
    gAlice->TreeR()->Branch(branchname, "TObjArray", &fEmcRecPoints, fBufferSize, splitlevel) ; 
  }

  if (fPpsdRecPoints) { 
    fPpsdRecPoints->Delete() ; 
    delete fPpsdRecPoints ; 
    fPpsdRecPoints = 0 ; 
  }

  //  fPpsdRecPoints = new RecPointsList("AliPHOSPpsdRecPoint", 100) ; if TClonesArray
  fPpsdRecPoints = new RecPointsList(100) ;

  if ( fPpsdRecPoints && gAlice->TreeR() ) {
    sprintf(branchname,"%sPpsdRP",GetName()) ;
     
     // gAlice->TreeR()->Branch(branchname, &fPpsdRecPoints, fBufferSize); if TClonesArray
    gAlice->TreeR()->Branch(branchname, "TObjArray", &fPpsdRecPoints, fBufferSize, splitlevel) ;
  }

  if (fTrackSegments) { 
   fTrackSegments->Delete() ; 
    delete fTrackSegments ; 
    fTrackSegments = 0 ; 
  }

  fTrackSegments = new TrackSegmentsList("AliPHOSTrackSegment", 100) ;
  if ( fTrackSegments && gAlice->TreeR() ) { 
    sprintf(branchname,"%sTS",GetName()) ;
    gAlice->TreeR()->Branch(branchname, &fTrackSegments, fBufferSize) ;
  }

  if (fRecParticles) {  
    fRecParticles->Delete() ; 
    delete fRecParticles ; 
    fRecParticles = 0 ; 
  }
  fRecParticles = new RecParticlesList("AliPHOSRecParticle", 100) ;
  if ( fRecParticles && gAlice->TreeR() ) { 
     sprintf(branchname,"%sRP",GetName()) ;
     gAlice->TreeR()->Branch(branchname, &fRecParticles, fBufferSize) ;
  }
  
  // 3.

  fReconstructioner->Make(fDigits, fEmcRecPoints, fPpsdRecPoints, fTrackSegments, fRecParticles);

  // 4. Expand or Shrink the arrays to the proper size
  
  Int_t size ;
  
  size = fEmcRecPoints->GetEntries() ;
  fEmcRecPoints->Expand(size) ;
 
  size = fPpsdRecPoints->GetEntries() ;
  fPpsdRecPoints->Expand(size) ;

  size = fTrackSegments->GetEntries() ;
  fTrackSegments->Expand(size) ;

  size = fRecParticles->GetEntries() ;
  fRecParticles->Expand(size) ;

  gAlice->TreeR()->Fill() ;
 
  // 5.

  gAlice->TreeR()->Write() ;
   
}

//____________________________________________________________________________
void AliPHOSv0::ResetDigits()
{
  // May sound strange, but cumulative hits are store in digits Tree

  if(  fTmpHits ) {
    fTmpHits->Delete();
    fNTmpHits = 0 ;
  }
}
  
//____________________________________________________________________________
void AliPHOSv0::StepManager(void)
{
  // Accumulates hits as long as the track stays in a single crystal or PPSD gas Cell

  Int_t          relid[4] ;      // (box, layer, row, column) indices
  Float_t        xyze[4] ;       // position wrt MRS and energy deposited
  TLorentzVector pos ;
  Int_t copy ;

  Int_t primary =  gAlice->GetPrimary( gAlice->CurrentTrack() ); 
  TString name = fGeom->GetName() ; 
  if ( name == "GPS2" ) { // the CPV is a PPSD
    if( gMC->CurrentVolID(copy) == gMC->VolId("GCEL") ) // We are inside a gas cell 
    {
      gMC->TrackPosition(pos) ;
      xyze[0] = pos[0] ;
      xyze[1] = pos[1] ;
      xyze[2] = pos[2] ;
      xyze[3] = gMC->Edep() ; 

      if ( xyze[3] != 0 ) { // there is deposited energy 
       	gMC->CurrentVolOffID(5, relid[0]) ;  // get the PHOS Module number
       	gMC->CurrentVolOffID(3, relid[1]) ;  // get the Micromegas Module number 
      // 1-> Geom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() upper                         
      //  >  fGeom->GetNumberOfModulesPhi()  *  fGeom->GetNumberOfModulesZ() lower
       	gMC->CurrentVolOffID(1, relid[2]) ;  // get the row number of the cell
        gMC->CurrentVolID(relid[3]) ;        // get the column number 

	// get the absolute Id number

	Int_t absid ; 
       	fGeom->RelToAbsNumbering(relid, absid) ; 

	// add current hit to the hit list      
	AddHit(primary, absid, xyze);

      } // there is deposited energy 
     } // We are inside the gas of the CPV  
   } // GPS2 configuration
  
   if(gMC->CurrentVolID(copy) == gMC->VolId("PXTL") )  //  We are inside a PBWO crystal
     {
       gMC->TrackPosition(pos) ;
       xyze[0] = pos[0] ;
       xyze[1] = pos[1] ;
       xyze[2] = pos[2] ;
       xyze[3] = gMC->Edep() ;

       if ( xyze[3] != 0 ) {
          gMC->CurrentVolOffID(10, relid[0]) ; // get the PHOS module number ;
          relid[1] = 0   ;                    // means PBW04
          gMC->CurrentVolOffID(4, relid[2]) ; // get the row number inside the module
          gMC->CurrentVolOffID(3, relid[3]) ; // get the cell number inside the module

      // get the absolute Id number

          Int_t absid ; 
          fGeom->RelToAbsNumbering(relid, absid) ; 
 
      // add current hit to the hit list

          AddHit(primary, absid, xyze);
    
       } // there is deposited energy
    } // we are inside a PHOS Xtal
}

