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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                          V-Zero   Detector                            //
//  This class contains the base procedures for the VZERO  detector      //
//  Geometry of November 2003 : V0R box is now 4.4 cm thick              //
//                              scintillators are 2 cm thick             //
//  All comments should be sent to Brigitte CHEYNIS :                    //
//                                 b.cheynis@ipnl.in2p3.fr               //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


// --- Standard libraries ---
#include <Riostream.h>

// --- ROOT libraries ---
#include <TNamed.h>
#include <TTree.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliVZERO.h"
#include "AliVZEROLoader.h"

ClassImp(AliVZERO)
 
//_____________________________________________________________________________
AliVZERO::AliVZERO(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for VZERO Detector
  //
  
//  fIshunt       =  1;  // All hits are associated with primary particles  
   
  fHits         =  new TClonesArray("AliVZEROhit", 400);
  fDigits       =  new TClonesArray("AliVZEROdigit",400); 
   
  gAlice->GetMCApp()->AddHitList(fHits);

  fThickness    =  4.4;   // total thickness of the V0R box in cm
  fThickness1   =  2.0;   // thickness of scintillating cells in cm
  
  fMaxStepQua   =  0.05; 
  fMaxStepAlu   =  0.01; 
  
  fMaxDestepQua =  -1.0;
  fMaxDestepAlu =  -1.0;
  
  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliVZERO::~AliVZERO()
{
  //
  // Default destructor for VZERO Detector
  //
  
    if (fHits) {
        fHits->Delete();
        delete fHits;
	fHits=0;
    }
    
    if (fDigits) {
        fDigits->Delete();
        delete fDigits;
        fDigits=0;
    }
}

//_____________________________________________________________________________
void AliVZERO::BuildGeometry()
{
  //
  // Builds simple ROOT TNode geometry for event display
  //
}
 
//_____________________________________________________________________________
void AliVZERO::CreateGeometry()
{
  //
  // Builds simple Geant3 geometry 
  //
}
//_____________________________________________________________________________
void AliVZERO::CreateMaterials()
{
  //
  // Creates materials used for Geant3 geometry 
  //
}

//_____________________________________________________________________________
Int_t AliVZERO::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
  //
  // Calculates the distance from the mouse to the VZERO on the screen
  // Dummy routine
  //
  
  return 9999;
}
 
//_____________________________________________________________________________
void AliVZERO::Init()
{
  //
  // Initialises the VZERO  class after it has been built
  //
}


//_____________________________________________________________________________
void AliVZERO::SetMaxStepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials
  //
     fMaxStepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxStepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum foils (not used in 
  // version v2)
  //
    fMaxStepAlu = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials (quartz)
  //
    fMaxDestepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum (not used in 
  // version v2)
  //
    fMaxDestepAlu = p1;
}

//_____________________________________________________________________________
AliLoader* AliVZERO::MakeLoader(const char* topfoldername)
{ 
  //
  // Builds VZEROgetter (AliLoader type)
  // if detector wants to use customized getter, it must overload this method
  //

  Info("MakeLoader","Creating AliVZEROLoader. Top folder is %s.",topfoldername);
  fLoader = new AliVZEROLoader(GetName(),topfoldername);
  return fLoader;
}

//_____________________________________________________________________________
void AliVZERO::SetTreeAddress()
{
  // 
  // Sets tree address for hits.
  //

  if (fLoader->TreeH() && (fHits == 0x0))
    fHits = new  TClonesArray("AliVZEROhit", 400);

  AliDetector::SetTreeAddress();
}


