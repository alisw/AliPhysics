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

//_________________________________________________________________________
// Manager class for PHOS version SUBATECH
//*-- Author : Y. Schutz SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TBRIK.h"
#include "TNode.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>
#include <assert.h>

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
  fNTmpHits = 0 ; 
  fTmpHits  = 0 ; 
}

//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(const char *name, const char *title):
  AliPHOS(name,title)
{
  
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fHits   = new TClonesArray("AliPHOSHit",100) ;
  fDigits = new TClonesArray("AliPHOSDigit",100) ;
  fTmpHits= new TClonesArray("AliPHOSHit",100) ;

  assert ( fHits != 0 ) ;
  assert ( fDigits != 0 ) ;
  assert ( fTmpHits != 0 ) ;

  fNTmpHits = fNhits = 0 ;

  fIshunt     =  1 ; // All hits are associated with primary particles
 
  // gets an instance of the geometry parameters class  
  fGeom =  AliPHOSGeometry::GetInstance(title, "") ; 

  if (fGeom->IsInitialized() ) 
    cout << "AliPHOSv0 : PHOS geometry intialized for " << fGeom->GetName() << endl ;
  else
   cout << "AliPHOSv0 : PHOS geometry initialization failed !" << endl ;   
}
//____________________________________________________________________________
AliPHOSv0::AliPHOSv0(AliPHOSReconstructioner&  Reconstructioner, const char *name, const char *title):
  AliPHOS(name,title)
{
  
  // We use 2 arrays of hits :
  //
  //   - fHits (the "normal" one), which retains the hits associated with
  //     the current primary particle being tracked
  //     (this array is reset after each primary has been tracked).
  //
  //   - fTmpHits, which retains all the hits of the current event. It 
  //     is used for the digitization part.

  fHits   = new TClonesArray("AliPHOSHit",100) ;
  fDigits = new TClonesArray("AliPHOSDigit",100) ;
  fTmpHits= new TClonesArray("AliPHOSHit",100) ;

  assert ( fHits != 0 ) ;
  assert ( fDigits != 0 ) ;
  assert ( fTmpHits != 0 ) ;

  fNTmpHits = fNhits = 0 ;

  fIshunt     =  1 ; // All hits are associated with primary particles
 
  // gets an instance of the geometry parameters class  
  fGeom =  AliPHOSGeometry::GetInstance(title, "") ; 

  if (fGeom->IsInitialized() ) 
    cout << "AliPHOSv0 : PHOS geometry intialized for " << fGeom->GetName() << endl ;
  else
   cout << "AliPHOSv0 : PHOS geometry initialization failed !" << endl ;   

  // Defining the PHOS Reconstructioner
 
 fReconstructioner = &Reconstructioner;
}

//____________________________________________________________________________
AliPHOSv0::~AliPHOSv0()
{
  delete fHits ;
  delete fTmpHits ;
  delete fDigits ;
}

//____________________________________________________________________________
void AliPHOSv0::AddHit(Int_t track, Int_t Id, Float_t * hits)
{
  Int_t hitCounter ;
  TClonesArray &ltmphits = *fTmpHits;
  AliPHOSHit *newHit ;
  AliPHOSHit *curHit;
  bool already = false ;

  // In any case, fills the fTmpHit TClonesArray (with "accumulated hits")

  newHit = new AliPHOSHit(fIshunt, track, Id, hits) ;

  for ( hitCounter = 0 ; hitCounter < fNTmpHits && !already ; hitCounter++ ) {
    curHit = (AliPHOSHit*) ltmphits[hitCounter] ;
    if( *curHit == *newHit ) {
      *curHit = *curHit + *newHit ;
      already = true ;
    }
  }
       
  if ( !already ) {
    new(ltmphits[fNTmpHits]) AliPHOSHit(*newHit) ;
    fNTmpHits++ ;
  }

  // Please note that the fTmpHits array must survive up to the
  // end of the events, so it does not appear e.g. in ResetHits() (
  // which is called at the end of each primary).  

  //  if (IsTreeSelected('H')) {
    // And, if we really want raw hits tree, have the fHits array filled also
  //    TClonesArray &lhits = *fHits;
  //    new(lhits[fNhits]) AliPHOSHit(*newHit) ;
  //    fNhits++ ;
  //  }

   delete newHit;

}


//____________________________________________________________________________
void AliPHOSv0::BuildGeometry()
{

  this->BuildGeometryforPHOS() ; 
  if ( ( strcmp(fGeom->GetName(), "GPS2" ) == 0 ) ) 
    this->BuildGeometryforPPSD() ;
  else
    cout << "AliPHOSv0::BuildGeometry : no charged particle identification system installed" << endl; 

}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforPHOS(void)
{
 // Build the PHOS geometry for the ROOT display

  const Int_t kColorPHOS = kRed ;
  const Int_t kColorXTAL = kBlue ;

  Double_t const RADDEG = 180.0 / kPI ;
 
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

  Float_t XTL_X = fGeom->GetCrystalSize(0) ; 
  Float_t XTL_Y = fGeom->GetCrystalSize(1) ; 
  Float_t XTL_Z = fGeom->GetCrystalSize(2) ; 

  Float_t XL =  fGeom->GetNPhi() * ( XTL_X + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 + fGeom->GetModuleBoxThickness() ;
  Float_t YL =  ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 
             + fGeom->GetModuleBoxThickness() / 2.0 ;
  Float_t ZL =  fGeom->GetNZ() * ( XTL_Z + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 +  fGeom->GetModuleBoxThickness() ;
  
  new TBRIK( "CrystalsBox", "PHOS crystals box", "void", XL, YL, ZL ) ;

// position PHOS into ALICE

  Float_t R = fGeom->GetIPtoOuterCoverDistance() + fGeom->GetOuterBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  Float_t pphi =  TMath::ATan( fGeom->GetOuterBoxSize(0)  / ( 2.0 * fGeom->GetIPtoOuterCoverDistance() ) ) ;
  pphi *= RADDEG ;
  TNode * Top = gAlice->GetGeometry()->GetNode("alice") ;
 
  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) { 
   Float_t angle = pphi * 2 * ( i - fGeom->GetNModules() / 2.0 - 0.5 ) ;
   sprintf(rotname, "%s%d", "rot", number++) ;
   new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
   Top->cd();
   sprintf(nodename,"%s%d", "Module", i) ;    
   Float_t X =  R * TMath::Sin( angle / RADDEG ) ;
   Float_t Y = -R * TMath::Cos( angle / RADDEG ) ;
   TNode * OuterBoxNode = new TNode(nodename, nodename, "OuterBox", X, Y, 0, rotname ) ;
   OuterBoxNode->SetLineColor(kColorPHOS) ;
   fNodes->Add(OuterBoxNode) ;
   OuterBoxNode->cd() ; 
   // now inside the outer box the textolit box
   Y = ( fGeom->GetOuterBoxThickness(1) -  fGeom->GetUpperPlateThickness() ) / 2.  ;
   sprintf(nodename,"%s%d", "TexBox", i) ;  
   TNode * TextolitBoxNode = new TNode(nodename, nodename, "TextolitBox", 0, Y, 0) ; 
   TextolitBoxNode->SetLineColor(kColorPHOS) ;
   fNodes->Add(TextolitBoxNode) ;
   // upper foam plate inside outre box
   OuterBoxNode->cd() ; 
   sprintf(nodename, "%s%d", "UFPlate", i) ;
   Y =  ( fGeom->GetTextolitBoxSize(1) - fGeom->GetSecondUpperPlateThickness() ) / 2.0 ;
   TNode * UpperFoamPlateNode = new TNode(nodename, nodename, "UpperFoamPlate", 0, Y, 0) ; 
   UpperFoamPlateNode->SetLineColor(kColorPHOS) ;
   fNodes->Add(UpperFoamPlateNode) ;  
   // air filled box inside textolit box (not drawn)
   TextolitBoxNode->cd();
   Y = ( fGeom->GetTextolitBoxSize(1) - fGeom->GetAirFilledBoxSize(1) ) / 2.0 -  fGeom->GetSecondUpperPlateThickness() ;
   sprintf(nodename, "%s%d", "AFBox", i) ;
   TNode * AirFilledBoxNode = new TNode(nodename, nodename, "AirFilledBox", 0, Y, 0) ; 
   fNodes->Add(AirFilledBoxNode) ; 
   // crystals box inside air filled box
   AirFilledBoxNode->cd() ; 
   Y = fGeom->GetAirFilledBoxSize(1) / 2.0 - YL 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness() 
       -  fGeom->GetUpperPlateThickness() -  fGeom->GetSecondUpperPlateThickness() ) ; 
   sprintf(nodename, "%s%d", "XTBox", i) ; 
   TNode * CrystalsBoxNode = new TNode(nodename, nodename, "CrystalsBox", 0, Y, 0) ;    
   CrystalsBoxNode->SetLineColor(kColorXTAL) ; 
   fNodes->Add(CrystalsBoxNode) ; 
  }
}

//____________________________________________________________________________
void AliPHOSv0:: BuildGeometryforPPSD(void)
{
 //  Build the PPSD geometry for the ROOT display

  Double_t const RADDEG = 180.0 / kPI ;

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

  Float_t R = fGeom->GetIPtoTopLidDistance() + fGeom->GetPPSDBoxSize(1) / 2.0 ;
  Int_t number = 988 ; 
  TNode * Top = gAlice->GetGeometry()->GetNode("alice") ;
 
  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) { // the number of PHOS modules
    Float_t angle = fGeom->GetPHOSAngle(i) ;
    sprintf(rotname, "%s%d", "rotg", number++) ;
    new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
    Top->cd();
    sprintf(nodename, "%s%d", "Moduleg", i) ;    
    Float_t X =  R * TMath::Sin( angle / RADDEG ) ;
    Float_t Y = -R * TMath::Cos( angle / RADDEG ) ;
    TNode * PPSDBoxNode = new TNode(nodename , nodename ,"PPSDBox", X, Y, 0, rotname ) ;
    PPSDBoxNode->SetLineColor(kColorPPSD) ;
    fNodes->Add(PPSDBoxNode) ;
    PPSDBoxNode->cd() ;
    // inside the PPSD box: 
    //   1.   fNumberOfModulesPhi x fNumberOfModulesZ top micromegas
    X = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. ;  
    for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
      Float_t Z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2. ;
      TNode * Micro1Node ; 
      for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	Y = ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas1Thickness() ) / 2. ; 
	sprintf(nodename, "%s%d%d%d", "Mic1", i, iphi, iz) ;
	Micro1Node  = new TNode(nodename, nodename, "PPSDModule", X, Y, Z) ;
	Micro1Node->SetLineColor(kColorPPSD) ;  
	fNodes->Add(Micro1Node) ; 
	// inside top micromegas
	Micro1Node->cd() ; 
	//      a. top lid
	Y = ( fGeom->GetMicromegas1Thickness() - fGeom->GetLidThickness() ) / 2. ; 
	sprintf(nodename, "%s%d%d%d", "Lid", i, iphi, iz) ;
	TNode * TopLidNode = new TNode(nodename, nodename, "TopLid", 0, Y, 0) ;
	TopLidNode->SetLineColor(kColorPPSD) ;  
	fNodes->Add(TopLidNode) ; 
	//      b. composite panel
	Y = Y - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 
	sprintf(nodename, "%s%d%d%d", "CompU", i, iphi, iz) ;
	TNode * CompUpNode = new TNode(nodename, nodename, "TopPanel", 0, Y, 0) ;
	CompUpNode->SetLineColor(kColorPPSD) ;  
	fNodes->Add(CompUpNode) ; 
	//      c. anode
	Y = Y - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
	sprintf(nodename, "%s%d%d%d", "Ano", i, iphi, iz) ;
	TNode * AnodeNode = new TNode(nodename, nodename, "Anode", 0, Y, 0) ;
	AnodeNode->SetLineColor(kColorPHOS) ;  
	fNodes->Add(AnodeNode) ; 
	//      d.  gas 
	Y = Y - fGeom->GetAnodeThickness() / 2. - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. ; 
	sprintf(nodename, "%s%d%d%d", "GGap", i, iphi, iz) ;
	TNode * GGapNode = new TNode(nodename, nodename, "GasGap", 0, Y, 0) ;
	GGapNode->SetLineColor(kColorGas) ;  
	fNodes->Add(GGapNode) ;          
	  //      f. cathode
	Y = Y - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 
	sprintf(nodename, "%s%d%d%d", "Cathode", i, iphi, iz) ;
	TNode * CathodeNode = new TNode(nodename, nodename, "Cathode", 0, Y, 0) ;
	CathodeNode->SetLineColor(kColorPHOS) ;  
	fNodes->Add(CathodeNode) ;        
	//      g. printed circuit
	Y = Y - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 
	sprintf(nodename, "%s%d%d%d", "PC", i, iphi, iz) ;
	TNode * PCNode = new TNode(nodename, nodename, "PCBoard", 0, Y, 0) ;
	PCNode->SetLineColor(kColorPPSD) ;  
	fNodes->Add(PCNode) ;        
	//      h. composite panel
	Y = Y - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 
	sprintf(nodename, "%s%d%d%d", "CompDown", i, iphi, iz) ;
	TNode * CompDownNode = new TNode(nodename, nodename, "BottomPanel", 0, Y, 0) ;
	CompDownNode->SetLineColor(kColorPPSD) ;  
	fNodes->Add(CompDownNode) ;   
	Z = Z - fGeom->GetPPSDModuleSize(2) ;
	PPSDBoxNode->cd() ;
      } // end of Z module loop     
      X = X -  fGeom->GetPPSDModuleSize(0) ; 
      PPSDBoxNode->cd() ;
    } // end of phi module loop
    //   2. air gap      
    PPSDBoxNode->cd() ;
    Y = ( fGeom->GetPPSDBoxSize(1) - 2 * fGeom->GetMicromegas1Thickness() - fGeom->GetMicro1ToLeadGap() ) / 2. ; 
    sprintf(nodename, "%s%d", "GapUp", i) ;
    TNode * GapUpNode = new TNode(nodename, nodename, "LeadToM", 0, Y, 0) ;
    GapUpNode->SetLineColor(kColorAir) ;  
    fNodes->Add(GapUpNode) ;        
    //   3. lead converter
    Y = Y - fGeom->GetMicro1ToLeadGap() / 2. - fGeom->GetLeadConverterThickness() / 2. ; 
    sprintf(nodename, "%s%d", "LeadC", i) ;
    TNode * LeadCNode = new TNode(nodename, nodename, "Lead", 0, Y, 0) ;
    LeadCNode->SetLineColor(kColorPPSD) ;  
    fNodes->Add(LeadCNode) ;        
    //   4. air gap
    Y = Y - fGeom->GetLeadConverterThickness() / 2. - fGeom->GetLeadToMicro2Gap()  / 2. ; 
    sprintf(nodename, "%s%d", "GapDown", i) ;
    TNode * GapDownNode = new TNode(nodename, nodename, "MToLead", 0, Y, 0) ;
    GapDownNode->SetLineColor(kColorAir) ;  
    fNodes->Add(GapDownNode) ;        
    //    5.  fNumberOfModulesPhi x fNumberOfModulesZ bottom micromegas
    X = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. - fGeom->GetPhiDisplacement() ;  
    for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { 
      Float_t Z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2.  - fGeom->GetZDisplacement() ;;
      TNode * Micro2Node ; 
      for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { 
	Y = - ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas2Thickness() ) / 2. ; 
	sprintf(nodename, "%s%d%d%d", "Mic2", i, iphi, iz) ;
	Micro2Node  = new TNode(nodename, nodename, "PPSDModule", X, Y, Z) ;
	Micro2Node->SetLineColor(kColorPPSD) ;  
	fNodes->Add(Micro2Node) ; 
	// inside bottom micromegas
	Micro2Node->cd() ; 
	  //      a. top lid
	  Y = ( fGeom->GetMicromegas2Thickness() - fGeom->GetLidThickness() ) / 2. ; 
	  sprintf(nodename, "%s%d", "Lidb", i) ;
	  TNode * TopLidbNode = new TNode(nodename, nodename, "TopLid", 0, Y, 0) ;
	  TopLidbNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(TopLidbNode) ; 
	  //      b. composite panel
	  Y = Y - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 
	  sprintf(nodename, "%s%d", "CompUb", i) ;
	  TNode * CompUpbNode = new TNode(nodename, nodename, "TopPanel", 0, Y, 0) ;
	  CompUpbNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(CompUpbNode) ; 
	  //      c. anode
	  Y = Y - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "Anob", i) ;
	  TNode * AnodebNode = new TNode(nodename, nodename, "Anode", 0, Y, 0) ;
	  AnodebNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(AnodebNode) ; 
	  //      d. conversion gas
	  Y = Y - fGeom->GetAnodeThickness() / 2. - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() )  / 2. ; 
	  sprintf(nodename, "%s%d", "GGapb", i) ;
	  TNode * GGapbNode = new TNode(nodename, nodename, "GasGap", 0, Y, 0) ;
	  GGapbNode->SetLineColor(kColorGas) ;  
	  fNodes->Add(GGapbNode) ;           
	  //      f. cathode
	  Y = Y - ( fGeom->GetConversionGap() + fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "Cathodeb", i) ;
	  TNode * CathodebNode = new TNode(nodename, nodename, "Cathode", 0, Y, 0) ;
	  CathodebNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(CathodebNode) ;        
	  //      g. printed circuit
	  Y = Y - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "PCb", i) ;
	  TNode * PCbNode = new TNode(nodename, nodename, "PCBoard", 0, Y, 0) ;
	  PCbNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(PCbNode) ;        
	  //      h. composite pane
	  Y = Y - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 
	  sprintf(nodename, "%s%d", "CompDownb", i) ;
	  TNode * CompDownbNode = new TNode(nodename, nodename, "BottomPanel", 0, Y, 0) ;
	  CompDownbNode->SetLineColor(kColorPPSD) ;  
	  fNodes->Add(CompDownbNode) ;        
       	  Z = Z - fGeom->GetPPSDModuleSize(2) ;
	  PPSDBoxNode->cd() ;
	} // end of Z module loop     
	X = X -  fGeom->GetPPSDModuleSize(0) ; 
	PPSDBoxNode->cd() ;
       } // end of phi module loop
     } // PHOS modules
 delete rotname ; 
 delete nodename ; 
}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometry()
{

  AliPHOSv0 *PHOS_tmp = (AliPHOSv0*)gAlice->GetModule("PHOS") ;

  if ( PHOS_tmp == NULL ) {
    
    fprintf(stderr, "PHOS detector not found!\n") ;
    return;
    
  }

  // Get pointer to the array containing media indeces
  Int_t *IDTMED = fIdtmed->GetArray() - 699 ;

  Float_t BigBox[3] ; 
  BigBox[0] =   fGeom->GetOuterBoxSize(0) / 2.0 ;
  BigBox[1] = ( fGeom->GetOuterBoxSize(1) + fGeom->GetPPSDBoxSize(1) ) / 2.0 ;
  BigBox[2] =   fGeom->GetOuterBoxSize(2) / 2.0 ;
  
  gMC->Gsvolu("PHOS", "BOX ", IDTMED[798], BigBox, 3) ;
  
  this->CreateGeometryforPHOS() ; 
  if ( strcmp( fGeom->GetName(), "GPS2") == 0  ) 
    this->CreateGeometryforPPSD() ;
  else
    cout << "AliPHOSv0::CreateGeometry : no charged particle identification system installed" << endl; 
  
  // --- Position  PHOS mdules in ALICE setup ---
  
  Int_t IDROTM[99] ;
  Double_t const RADDEG = 180.0 / kPI ;
  
  for( Int_t i = 1; i <= fGeom->GetNModules(); i++ ) {
    
    Float_t angle = fGeom->GetPHOSAngle(i) ;
    AliMatrix(IDROTM[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0) ;
 
    Float_t R = fGeom->GetIPtoOuterCoverDistance() + ( fGeom->GetOuterBoxSize(1) + fGeom->GetPPSDBoxSize(1) ) / 2.0 ;

    Float_t XP1 = R * TMath::Sin( angle / RADDEG ) ;
    Float_t YP1 = -R * TMath::Cos( angle / RADDEG ) ;

    gMC->Gspos("PHOS", i, "ALIC", XP1, YP1, 0.0, IDROTM[i-1], "ONLY") ;
 
  } // for GetNModules

}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforPHOS()
{
  // Get pointer to the array containing media indeces
  Int_t *IDTMED = fIdtmed->GetArray() - 699 ;

  // ---
  // --- Define PHOS box volume, fPUFPill with thermo insulating foam ---
  // --- Foam Thermo Insulating outer cover dimensions ---
  // --- Put it in BigBox = PHOS

  Float_t DPHOS[3] ; 
  DPHOS[0] =  fGeom->GetOuterBoxSize(0) / 2.0 ;
  DPHOS[1] =  fGeom->GetOuterBoxSize(1) / 2.0 ;
  DPHOS[2] =  fGeom->GetOuterBoxSize(2) / 2.0 ;

  gMC->Gsvolu("EMCA", "BOX ", IDTMED[706], DPHOS, 3) ;

  Float_t YO =  - fGeom->GetPPSDBoxSize(1)  / 2.0 ;

  gMC->Gspos("EMCA", 1, "PHOS", 0.0, YO, 0.0, 0, "ONLY") ; 

  // ---
  // --- Define Textolit Wall box, position inside EMCA ---
  // --- Textolit Wall box dimentions ---
 
 
  Float_t DPTXW[3];
  DPTXW[0] = fGeom->GetTextolitBoxSize(0) / 2.0 ;
  DPTXW[1] = fGeom->GetTextolitBoxSize(1) / 2.0 ;
  DPTXW[2] = fGeom->GetTextolitBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXW", "BOX ", IDTMED[707], DPTXW, 3);

  YO =   (  fGeom->GetOuterBoxThickness(1) -   fGeom->GetUpperPlateThickness() ) / 2.  ;
   
  gMC->Gspos("PTXW", 1, "EMCA", 0.0, YO, 0.0, 0, "ONLY") ;

  // --- 
  // --- Define Upper Polystyrene Foam Plate, place inside PTXW ---
  // --- immediately below Foam Thermo Insulation Upper plate ---

  // --- Upper Polystyrene Foam plate thickness ---
 
  Float_t  DPUFP[3] ;
  DPUFP[0] = fGeom->GetTextolitBoxSize(0) / 2.0 ; 
  DPUFP[1] = fGeom->GetSecondUpperPlateThickness() / 2. ;
  DPUFP[2] = fGeom->GetTextolitBoxSize(2) /2.0 ; 

  gMC->Gsvolu("PUFP", "BOX ", IDTMED[703], DPUFP, 3) ;
  
  YO = ( fGeom->GetTextolitBoxSize(1) -  fGeom->GetSecondUpperPlateThickness() ) / 2.0 ;
  
  gMC->Gspos("PUFP", 1, "PTXW", 0.0, YO, 0.0, 0, "ONLY") ;
  
  // ---
  // --- Define air-filled box, place inside PTXW ---
  // --- Inner AIR volume dimensions ---
 

  Float_t  DPAIR[3] ;
  DPAIR[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  DPAIR[1] = fGeom->GetAirFilledBoxSize(1) / 2.0 ;
  DPAIR[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PAIR", "BOX ", IDTMED[798], DPAIR, 3) ;
  
  YO = ( fGeom->GetTextolitBoxSize(1) -  fGeom->GetAirFilledBoxSize(1) ) / 2.0 -   fGeom->GetSecondUpperPlateThickness() ;
  
  gMC->Gspos("PAIR", 1, "PTXW", 0.0, YO, 0.0, 0, "ONLY") ;

// --- Dimensions of PbWO4 crystal ---

  Float_t XTL_X =  fGeom->GetCrystalSize(0) ; 
  Float_t XTL_Y =  fGeom->GetCrystalSize(1) ; 
  Float_t XTL_Z =  fGeom->GetCrystalSize(2) ; 

  Float_t DPTCB[3] ;  
  DPTCB[0] =  fGeom->GetNPhi() * ( XTL_X + 2 *  fGeom->GetGapBetweenCrystals() ) / 2.0 + fGeom->GetModuleBoxThickness() ;
  DPTCB[1] = ( XTL_Y +  fGeom->GetCrystalSupportHeight() +  fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 
             + fGeom->GetModuleBoxThickness() / 2.0 ;
  DPTCB[2] = fGeom->GetNZ() * ( XTL_Z + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 +  fGeom->GetModuleBoxThickness() ;
  
  gMC->Gsvolu("PTCB", "BOX ", IDTMED[706], DPTCB, 3) ;

  YO =  fGeom->GetAirFilledBoxSize(1) / 2.0 - DPTCB[1] 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness() 
       -  fGeom->GetUpperPlateThickness() -  fGeom->GetSecondUpperPlateThickness() ) ;
  
  gMC->Gspos("PTCB", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Crystal BLock filled with air, position it inside PTCB ---
  Float_t DPCBL[3] ; 
  
  DPCBL[0] = fGeom->GetNPhi() * ( XTL_X + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 ;
  DPCBL[1] = ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  DPCBL[2] = fGeom->GetNZ() * ( XTL_Z + 2 * fGeom->GetGapBetweenCrystals() ) / 2.0 ;
  
  gMC->Gsvolu("PCBL", "BOX ", IDTMED[798], DPCBL, 3) ;
  
  // --- Divide PCBL in X (phi) and Z directions --
  gMC->Gsdvn("PROW", "PCBL", Int_t (fGeom->GetNPhi()), 1) ;
  gMC->Gsdvn("PCEL", "PROW", Int_t (fGeom->GetNZ()), 3) ;

  YO = -fGeom->GetModuleBoxThickness() / 2.0 ;
  
  gMC->Gspos("PCBL", 1, "PTCB", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define STeel (actually, it's titanium) Cover volume, place inside PCEL
  Float_t  DPSTC[3] ; 
  
  DPSTC[0] = ( XTL_X + 2 * fGeom->GetCrystalWrapThickness() ) / 2.0 ;
  DPSTC[1] = ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  DPSTC[2] = ( XTL_Z + 2 * fGeom->GetCrystalWrapThickness()  + 2 *  fGeom->GetCrystalHolderThickness() ) / 2.0 ;
  
  gMC->Gsvolu("PSTC", "BOX ", IDTMED[704], DPSTC, 3) ;

  gMC->Gspos("PSTC", 1, "PCEL", 0.0, 0.0, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Tyvek volume, place inside PSTC ---
  Float_t  DPPAP[3] ;

  DPPAP[0] = XTL_X / 2.0 + fGeom->GetCrystalWrapThickness() ;
  DPPAP[1] = ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 ;
  DPPAP[2] = XTL_Z / 2.0 + fGeom->GetCrystalWrapThickness() ;
  
  gMC->Gsvolu("PPAP", "BOX ", IDTMED[702], DPPAP, 3) ;
  
  YO = ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 
              - ( XTL_Y +  fGeom->GetCrystalSupportHeight() +  fGeom->GetCrystalWrapThickness() + fGeom->GetCrystalHolderThickness() ) / 2.0 ;
   
  gMC->Gspos("PPAP", 1, "PSTC", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PbWO4 crystal volume, place inside PPAP ---
  Float_t  DPXTL[3] ; 

  DPXTL[0] = XTL_X / 2.0 ;
  DPXTL[1] = XTL_Y / 2.0 ;
  DPXTL[2] = XTL_Z / 2.0 ;
  
  gMC->Gsvolu("PXTL", "BOX ", IDTMED[699], DPXTL, 3) ;

  YO = ( XTL_Y + fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 - XTL_Y / 2.0 - fGeom->GetCrystalWrapThickness() ;
  
  gMC->Gspos("PXTL", 1, "PPAP", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define crystal support volume, place inside PPAP ---
  Float_t DPSUP[3] ; 

  DPSUP[0] = XTL_X / 2.0 + fGeom->GetCrystalWrapThickness()  ;
  DPSUP[1] = fGeom->GetCrystalSupportHeight() / 2.0 ;
  DPSUP[2] = XTL_Z / 2.0 +  fGeom->GetCrystalWrapThickness() ;

  gMC->Gsvolu("PSUP", "BOX ", IDTMED[798], DPSUP, 3) ;

  YO =  fGeom->GetCrystalSupportHeight() / 2.0 - ( XTL_Y +  fGeom->GetCrystalSupportHeight() + fGeom->GetCrystalWrapThickness() ) / 2.0 ;

  gMC->Gspos("PSUP", 1, "PPAP", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define PIN-diode volume and position it inside crystal support ---
  // --- right behind PbWO4 crystal

  // --- PIN-diode dimensions ---

 
  Float_t DPPIN[3] ;
  DPPIN[0] = fGeom->GetPinDiodeSize(0) / 2.0 ;
  DPPIN[1] = fGeom->GetPinDiodeSize(1) / 2.0 ;
  DPPIN[2] = fGeom->GetPinDiodeSize(2) / 2.0 ;
 
  gMC->Gsvolu("PPIN", "BOX ", IDTMED[705], DPPIN, 3) ;
 
  YO = fGeom->GetCrystalSupportHeight() / 2.0 - fGeom->GetPinDiodeSize(1) / 2.0 ;
 
  gMC->Gspos("PPIN", 1, "PSUP", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Upper Cooling Panel, place it on top of PTCB ---
  Float_t DPUCP[3] ;
 // --- Upper Cooling Plate thickness ---
 
  DPUCP[0] = DPTCB[0] ;
  DPUCP[1] = fGeom->GetUpperCoolingPlateThickness() ;
  DPUCP[2] = DPTCB[2] ;
  
  gMC->Gsvolu("PUCP", "BOX ", IDTMED[701], DPUCP,3) ;
  
  YO = (  fGeom->GetAirFilledBoxSize(1) -  fGeom->GetUpperCoolingPlateThickness() ) / 2. 
       - ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetModuleBoxThickness()
           - fGeom->GetUpperPlateThickness() - fGeom->GetSecondUpperPlateThickness() - fGeom->GetUpperCoolingPlateThickness() ) ; 
  
  gMC->Gspos("PUCP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Al Support Plate, position it inside PAIR ---
  // --- right beneath PTCB ---
 // --- Al Support Plate thickness ---
 
  Float_t DPASP[3] ;
  DPASP[0] =  fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  DPASP[1] = fGeom->GetSupportPlateThickness() / 2.0 ;
  DPASP[2] =  fGeom->GetAirFilledBoxSize(2) / 2.0 ;
  
  gMC->Gsvolu("PASP", "BOX ", IDTMED[701], DPASP, 3) ;
  
  YO = (  fGeom->GetAirFilledBoxSize(1) - fGeom->GetSupportPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance()
           - fGeom->GetUpperPlateThickness() - fGeom->GetSecondUpperPlateThickness() + DPCBL[1] * 2 ) ;
  
  gMC->Gspos("PASP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Thermo Insulating Plate, position it inside PAIR ---
  // --- right beneath PASP ---
  // --- Lower Thermo Insulating Plate thickness ---
  
  Float_t DPTIP[3] ;
  DPTIP[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  DPTIP[1] = fGeom->GetLowerThermoPlateThickness() / 2.0 ;
  DPTIP[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTIP", "BOX ", IDTMED[706], DPTIP, 3) ;

  YO =  ( fGeom->GetAirFilledBoxSize(1) - fGeom->GetLowerThermoPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetUpperPlateThickness() 
            - fGeom->GetSecondUpperPlateThickness() + DPCBL[1] * 2 + fGeom->GetSupportPlateThickness() ) ;

  gMC->Gspos("PTIP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY") ;

  // ---
  // --- Define Textolit Plate, position it inside PAIR ---
  // --- right beneath PTIP ---
  // --- Lower Textolit Plate thickness ---
 
  Float_t DPTXP[3] ;
  DPTXP[0] = fGeom->GetAirFilledBoxSize(0) / 2.0 ;
  DPTXP[1] = fGeom->GetLowerTextolitPlateThickness() / 2.0 ;
  DPTXP[2] = fGeom->GetAirFilledBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PTXP", "BOX ", IDTMED[707], DPTXP, 3) ;

  YO =  ( fGeom->GetAirFilledBoxSize(1) - fGeom->GetLowerTextolitPlateThickness() ) / 2. 
       -  ( fGeom->GetIPtoCrystalSurface() - fGeom->GetIPtoOuterCoverDistance() - fGeom->GetUpperPlateThickness() 
            - fGeom->GetSecondUpperPlateThickness() + DPCBL[1] * 2 + fGeom->GetSupportPlateThickness() 
            +  fGeom->GetLowerThermoPlateThickness() ) ;

  gMC->Gspos("PTXP", 1, "PAIR", 0.0, YO, 0.0, 0, "ONLY") ;

}

//____________________________________________________________________________
void AliPHOSv0::CreateGeometryforPPSD()
{
  // Get pointer to the array containing media indeces
  Int_t *IDTMED = fIdtmed->GetArray() - 699 ;
  
  // The box containing all PPSD's for one PHOS module filled with air 
  Float_t PPSD[3] ; 
  PPSD[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;  
  PPSD[1] = fGeom->GetPPSDBoxSize(1) / 2.0 ; 
  PPSD[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("PPSD", "BOX ", IDTMED[798], PPSD, 3) ;

  Float_t YO =  fGeom->GetOuterBoxSize(1) / 2.0 ;

  gMC->Gspos("PPSD", 1, "PHOS", 0.0, YO, 0.0, 0, "ONLY") ; 

  // Now we build a micromegas module
  // The box containing the whole module filled with epoxy (FR4)

  Float_t MPPSD[3] ;  
  MPPSD[0] = fGeom->GetPPSDModuleSize(0) / 2.0 ;  
  MPPSD[1] = fGeom->GetPPSDModuleSize(1) / 2.0 ;  
  MPPSD[2] = fGeom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("MPPS", "BOX ", IDTMED[708], MPPSD, 3) ;  
 
  // Inside MPPSD :
  // 1. The Top Lid made of epoxy (FR4) 

  Float_t TLPPSD[3] ; 
  TLPPSD[0] = fGeom->GetPPSDModuleSize(0) / 2.0 ; 
  TLPPSD[1] = fGeom->GetLidThickness() / 2.0 ;
  TLPPSD[2] = fGeom->GetPPSDModuleSize(2) / 2.0 ;

  gMC->Gsvolu("TLPS", "BOX ", IDTMED[708], TLPPSD, 3) ; 

  Float_t  Y0 = ( fGeom->GetMicromegas1Thickness() - fGeom->GetLidThickness() ) / 2. ; 

  gMC->Gspos("TLPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 
 
  // 2. the upper panel made of composite material

  Float_t UPPPSD[3] ; 
  UPPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  UPPPSD[1] = fGeom->GetCompositeThickness() / 2.0 ;
  UPPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
 
  gMC->Gsvolu("UPPS", "BOX ", IDTMED[709], UPPPSD, 3) ; 
  
  Y0 = Y0 - fGeom->GetLidThickness() / 2. - fGeom->GetCompositeThickness() / 2. ; 

  gMC->Gspos("UPPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // 3. the anode made of Copper
  
  Float_t ANPPSD[3] ; 
  ANPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ; 
  ANPPSD[1] = fGeom->GetAnodeThickness() / 2.0 ; 
  ANPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0  ; 

  gMC->Gsvolu("ANPS", "BOX ", IDTMED[710], ANPPSD, 3) ; 
  
  Y0 = Y0 - fGeom->GetCompositeThickness() / 2. - fGeom->GetAnodeThickness()  / 2. ; 
  
  gMC->Gspos("ANPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // 4. the conversion gap + avalanche gap filled with gas

  Float_t GGPPSD[3] ; 
  GGPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  GGPPSD[1] = ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2.0 ; 
  GGPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("GGPS", "BOX ", IDTMED[715], GGPPSD, 3) ; 
  
  // --- Divide GGPP in X (phi) and Z directions --
  gMC->Gsdvn("GROW", "GGPS", fGeom->GetNumberOfPadsPhi(), 1) ;
  gMC->Gsdvn("GCEL", "GROW", fGeom->GetNumberOfPadsZ() , 3) ;

  Y0 = Y0 - fGeom->GetAnodeThickness() / 2.  - ( fGeom->GetConversionGap() +  fGeom->GetAvalancheGap() ) / 2. ; 

  gMC->Gspos("GGPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 


  // 6. the cathode made of Copper

  Float_t CAPPSD[3] ;
  CAPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;
  CAPPSD[1] = fGeom->GetCathodeThickness() / 2.0 ; 
  CAPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0  ;

  gMC->Gsvolu("CAPS", "BOX ", IDTMED[710], CAPPSD, 3) ; 

  Y0 = Y0 - ( fGeom->GetAvalancheGap() +  fGeom->GetAvalancheGap() ) / 2. - fGeom->GetCathodeThickness()  / 2. ; 

  gMC->Gspos("CAPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // 7. the printed circuit made of G10       

  Float_t PCPPSD[3] ; 
  PCPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2,.0 ; 
  PCPPSD[1] = fGeom->GetPCThickness() / 2.0 ; 
  PCPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("PCPS", "BOX ", IDTMED[711], CAPPSD, 3) ; 

  Y0 = Y0 - fGeom->GetCathodeThickness() / 2. - fGeom->GetPCThickness()  / 2. ; 

  gMC->Gspos("PCPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // 8. the lower panel made of composite material
						    
  Float_t LPPPSD[3] ; 
  LPPPSD[0] = ( fGeom->GetPPSDModuleSize(0) - fGeom->GetMicromegasWallThickness() ) / 2.0 ; 
  LPPPSD[1] = fGeom->GetCompositeThickness() / 2.0 ; 
  LPPPSD[2] = ( fGeom->GetPPSDModuleSize(2) - fGeom->GetMicromegasWallThickness() ) / 2.0 ;

  gMC->Gsvolu("LPPS", "BOX ", IDTMED[709], LPPPSD, 3) ; 
 
  Y0 = Y0 - fGeom->GetPCThickness() / 2. - fGeom->GetCompositeThickness()  / 2. ; 

  gMC->Gspos("LPPS", 1, "MPPS", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // Position the  fNumberOfModulesPhi x fNumberOfModulesZ modules (MPPSD) inside PPSD to cover a PHOS module
  // the top and bottom one's (which are assumed identical) :

   Float_t Yt = ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas1Thickness() ) / 2. ; 
   Float_t Yb = - ( fGeom->GetPPSDBoxSize(1) - fGeom->GetMicromegas2Thickness() ) / 2. ; 

   Int_t CopyNumbertop = 0 ; 
   Int_t CopyNumberbot = fGeom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() ; 

   Float_t X  = ( fGeom->GetPPSDBoxSize(0) - fGeom->GetPPSDModuleSize(0) ) / 2. ;  

   for ( Int_t iphi = 1; iphi <= fGeom->GetNumberOfModulesPhi(); iphi++ ) { // the number of micromegas modules in phi per PHOS module
      Float_t Z = ( fGeom->GetPPSDBoxSize(2) - fGeom->GetPPSDModuleSize(2) ) / 2. ;

      for ( Int_t iz = 1; iz <= fGeom->GetNumberOfModulesZ(); iz++ ) { // the number of micromegas modules in z per PHOS module
	gMC->Gspos("MPPS", ++CopyNumbertop, "PPSD", X, Yt, Z, 0, "ONLY") ;
	gMC->Gspos("MPPS", ++CopyNumberbot, "PPSD", X, Yb, Z, 0, "ONLY") ; 
	Z = Z - fGeom->GetPPSDModuleSize(2) ;
      } // end of Z module loop   
      X = X -  fGeom->GetPPSDModuleSize(0) ; 
    } // end of phi module loop

   // The Lead converter between two air gaps
   // 1. Upper air gap

   Float_t UAPPSD[3] ;
   UAPPSD[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;
   UAPPSD[1] = fGeom->GetMicro1ToLeadGap() / 2.0 ; 
   UAPPSD[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("UAPPSD", "BOX ", IDTMED[798], UAPPSD, 3) ; 

  Y0 = ( fGeom->GetPPSDBoxSize(1) - 2 * fGeom->GetMicromegas1Thickness() - fGeom->GetMicro1ToLeadGap() ) / 2. ; 

  gMC->Gspos("UAPPSD", 1, "PPSD", 0.0, Y0, 0.0, 0, "ONLY") ; 

   // 2. Lead converter
 
  Float_t LCPPSD[3] ; 
  LCPPSD[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ;
  LCPPSD[1] = fGeom->GetLeadConverterThickness() / 2.0 ; 
  LCPPSD[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;
 
  gMC->Gsvolu("LCPPSD", "BOX ", IDTMED[712], LCPPSD, 3) ; 
  
  Y0 = Y0 - fGeom->GetMicro1ToLeadGap() / 2. - fGeom->GetLeadConverterThickness() / 2. ; 

  gMC->Gspos("LCPPSD", 1, "PPSD", 0.0, Y0, 0.0, 0, "ONLY") ; 

  // 3. Lower air gap

  Float_t LAPPSD[3] ; 
  LAPPSD[0] = fGeom->GetPPSDBoxSize(0) / 2.0 ; 
  LAPPSD[1] = fGeom->GetLeadToMicro2Gap() / 2.0 ; 
  LAPPSD[2] = fGeom->GetPPSDBoxSize(2) / 2.0 ;

  gMC->Gsvolu("LAPPSD", "BOX ", IDTMED[798], LAPPSD, 3) ; 
    
  Y0 = Y0 - fGeom->GetLeadConverterThickness() / 2. - fGeom->GetLeadToMicro2Gap()  / 2. ; 
  
  gMC->Gspos("LAPPSD", 1, "PPSD", 0.0, Y0, 0.0, 0, "ONLY") ; 
   
}

//___________________________________________________________________________
Int_t AliPHOSv0::Digitize(Float_t Energy){
  Float_t fB = 10000000. ;
  Float_t fA = 0. ;
  Int_t chan = Int_t(fA + Energy*fB ) ;
  return chan ;
}
//___________________________________________________________________________
void AliPHOSv0::FinishEvent()
{
  cout << "//_____________________________________________________" << endl ;
  cout << "<I> AliPHOSv0::FinishEvent() -- Starting digitalization" << endl ;
  Int_t i ;
  TClonesArray &lDigits = *fDigits ;
  AliPHOSHit  * Hit ;
  AliPHOSDigit * Digit ;

  for ( i = 0 ; i < fNTmpHits ; i++ ) {
    Hit = (AliPHOSHit*)fTmpHits->At(i) ;
    assert (Hit!=0) ;
    Digit = new AliPHOSDigit(Hit->GetId(),Digitize(Hit->GetEnergy())) ;
    new(lDigits[fNdigits]) AliPHOSDigit(* Digit) ;
    fNdigits++;  delete Digit ;    
  }

  // Reset the array of all the "accumulated hits" of this event.
  fNTmpHits = 0 ;
  fTmpHits->Delete();
}

//____________________________________________________________________________
void AliPHOSv0::Init(void)
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

//___________________________________________________________________________
void AliPHOSv0::MakeBranch(Option_t* opt)
{  
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  //
  AliDetector::MakeBranch(opt) ;
  
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  char *D = strstr(opt,"D");
  
  if (fDigits && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, fBufferSize);
    printf("* AliPHOS::MakeBranch * Making Branch %s for digits\n",branchname);
  }
}
//_____________________________________________________________________________

void AliPHOSv0::Reconstruction(AliPHOSReconstructioner& Reconstructioner)
{ 
  fReconstructioner = &Reconstructioner;
  if (fEmcClusters) 
    { fEmcClusters->Delete();}
  else
    { fEmcClusters= new TClonesArray("AliPHOSEmcRecPoint", 100); } ;

  if (fPpsdClusters)
    { fPpsdClusters->Delete(); }
  else
    { fPpsdClusters = new TClonesArray("AliPHOSPpsdRecPoint", 100) ;}

  if (fTrackSegments)
    { fTrackSegments->Delete(); }
  else
    { fTrackSegments = new TObjArray(100) ;}
  
  fReconstructioner->Make(fDigits, fEmcClusters, fPpsdClusters, fTrackSegments);
  
}

//____________________________________________________________________________
void AliPHOSv0::StepManager(void)
{
  Int_t          RelId[4] ;      // (box, layer, row, column) indices
  Float_t        xyze[4] ;       // position wrt MRS and energy deposited
  TLorentzVector pos ;
  Int_t copy;

  TString name = fGeom->GetName() ; 

  if ( name == "GPS2" ) { // the CPV is a PPSD
    if( gMC->CurrentVolID(copy) == gMC->VolId("GCEL") )
    //     if( strcmp ( gMC->CurrentVolName(), "GCEL" ) == 0 )  // We are inside a gas cell 
    {
      gMC->TrackPosition(pos) ;
      xyze[0] = pos[0] ;
      xyze[1] = pos[1] ;
      xyze[2] = pos[2] ;
      xyze[3] = gMC->Edep() ; 

      if ( xyze[3] != 0 ) { // there is deposited energy 
       	gMC->CurrentVolOffID(5, RelId[0]) ;  // get the PHOS Module number
       	gMC->CurrentVolOffID(3, RelId[1]) ;  // get the Micromegas Module number 
      // 1-> Geom->GetNumberOfModulesPhi() *  fGeom->GetNumberOfModulesZ() upper                         
      //  >  fGeom->GetNumberOfModulesPhi()  *  fGeom->GetNumberOfModulesZ() lower
       	gMC->CurrentVolOffID(1, RelId[2]) ;  // get the row number of the cell
        gMC->CurrentVolID(RelId[3]) ;        // get the column number 

	// get the absolute Id number

	Int_t AbsId ; 
       	fGeom->RelToAbsNumbering(RelId,AbsId) ; 

	// add current hit to the hit list      
	AddHit(gAlice->CurrentTrack(), AbsId, xyze);

      } // there is deposited energy 
     } // We are inside the gas of the CPV  
   } // GPS2 configuration
  
   if(gMC->CurrentVolID(copy) == gMC->VolId("PXTL") ) 
  //      if( strcmp ( gMC->CurrentVolName(), "PXTL" ) == 0 ) { //  We are inside a PWO crystal
     {
       gMC->TrackPosition(pos) ;
       xyze[0] = pos[0] ;
       xyze[1] = pos[1] ;
       xyze[2] = pos[2] ;
       xyze[3] = gMC->Edep() ;

       if ( xyze[3] != 0 ) {
          gMC->CurrentVolOffID(10, RelId[0]) ; // get the PHOS module number ;
          RelId[1] = 0   ;                    // means PW04
          gMC->CurrentVolOffID(4, RelId[2]) ; // get the row number inside the module
          gMC->CurrentVolOffID(3, RelId[3]) ; // get the cell number inside the module

      // get the absolute Id number

          Int_t AbsId ; 
          fGeom->RelToAbsNumbering(RelId,AbsId) ; 
 
      // add current hit to the hit list

          AddHit(gAlice->CurrentTrack(), AbsId, xyze);
    
       } // there is deposited energy
    } // we are inside a PHOS Xtal
}

