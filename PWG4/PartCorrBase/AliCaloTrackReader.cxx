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
/* $Id:  $ */

//_________________________________________________________________________
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors (CTS).
// Not all MC particles/tracks/clusters are kept, some kinematical/fidutial restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader: Fills Kinematics data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TClonesArrays (PHOS, EMCAL, CTS) 
//                
//-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TMath.h>
#include <TLorentzVector.h>
#include <TString.h>

//---- ANALYSIS system ----
#include "AliCaloTrackReader.h"
#include "AliLog.h"
#include "AliStack.h"  
#include "AliHeader.h"  
#include "AliGenEventHeader.h"  

ClassImp(AliCaloTrackReader)
  
  
//____________________________________________________________________________
  AliCaloTrackReader::AliCaloTrackReader() : 
    TObject(), fEventNumber(-1), fDataType(0), fDebug(0), 
	fFidutialCut(0x0),
    fCTSPtMin(0), fEMCALPtMin(0),fPHOSPtMin(0),
    fAODCTS(0x0), fAODEMCAL(0x0), fAODPHOS(0x0),
    fEMCALCells(0x0), fPHOSCells(0x0),
    fESD(0x0), fAOD(0x0),fMC(0x0),
    fFillCTS(0),fFillEMCAL(0),fFillPHOS(0),
    fFillEMCALCells(0),fFillPHOSCells(0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliCaloTrackReader::AliCaloTrackReader(const AliCaloTrackReader & g) :   
  TObject(g), fEventNumber(g.fEventNumber), fDataType(g.fDataType), fDebug(g.fDebug),
  fFidutialCut(g.fFidutialCut),
  fCTSPtMin(g.fCTSPtMin), fEMCALPtMin(g.fEMCALPtMin),fPHOSPtMin(g.fPHOSPtMin), 
  fAODCTS(new TClonesArray(*g.fAODCTS)),  
  fAODEMCAL(new TClonesArray(*g.fAODEMCAL)),
  fAODPHOS(new TClonesArray(*g.fAODPHOS)),
  fEMCALCells(new TNamed(*g.fEMCALCells)),
  fPHOSCells(new TNamed(*g.fPHOSCells)),
  fESD(g.fESD), fAOD(g.fAOD), fMC(g.fMC),
  fFillCTS(g.fFillCTS),fFillEMCAL(g.fFillEMCAL),fFillPHOS(g.fFillPHOS),
  fFillEMCALCells(g.fFillEMCALCells),fFillPHOSCells(g.fFillPHOSCells)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliCaloTrackReader & AliCaloTrackReader::operator = (const AliCaloTrackReader & source)
{
  // assignment operator
  
  if(&source == this) return *this;
  
  fDataType    = source.fDataType ;
  fDebug       = source.fDebug ;
  fEventNumber = source.fEventNumber ;
  fFidutialCut = source.fFidutialCut;
  
  fCTSPtMin   = source.fCTSPtMin ;
  fEMCALPtMin = source.fEMCALPtMin ;
  fPHOSPtMin  = source.fPHOSPtMin ; 
  
  fAODCTS     = new TClonesArray(*source.fAODCTS) ;
  fAODEMCAL   = new TClonesArray(*source.fAODEMCAL) ;
  fAODPHOS    = new TClonesArray(*source.fAODPHOS) ;
  fEMCALCells = new TNamed(*source.fEMCALCells) ;
  fPHOSCells  = new TNamed(*source.fPHOSCells) ;

  fESD = source.fESD;
  fAOD = source.fAOD;
  fMC  = source.fMC;
  
  fFillCTS        = source.fFillCTS;
  fFillEMCAL      = source.fFillEMCAL;
  fFillPHOS       = source.fFillPHOS;
  fFillEMCALCells = source.fFillEMCALCells;
  fFillPHOSCells  = source.fFillPHOSCells;

  return *this;
  
}

//_________________________________
AliCaloTrackReader::~AliCaloTrackReader() {
  //Dtor
  
  if(fFidutialCut) delete fFidutialCut ;
  
  if(fAODCTS){
    fAODCTS->Clear() ; 
    delete fAODCTS ;
  }
  
  if(fAODEMCAL){
    fAODEMCAL->Clear() ; 
    delete fAODEMCAL ;
  }
  
  if(fAODPHOS){
    fAODPHOS->Clear() ; 
    delete fAODPHOS ;
  }
  
  if(fEMCALCells){
    fEMCALCells->Clear() ; 
    delete fEMCALCells ;
  }
  
  if(fPHOSCells){
    fPHOSCells->Clear() ; 
    delete fPHOSCells ;
  }

  if(fESD) delete fESD ;
  if(fAOD) delete fAOD ;
  if(fMC) delete fMC ;
}

//____________________________________________________________________________
AliStack* AliCaloTrackReader::GetStack() const {
  //Return pointer to stack
  if(fMC)
    return fMC->Stack();
  else{
    printf("stack is not available"); 
    return 0x0 ;
  }
}

//____________________________________________________________________________
AliHeader* AliCaloTrackReader::GetHeader() const {
  //Return pointer to header
  if(fMC)
    return fMC->Header();
  else{
    printf("Header is not available\n"); 
    return 0x0 ;
  }
}
//____________________________________________________________________________
AliGenEventHeader* AliCaloTrackReader::GetGenEventHeader() const {
  //Return pointer to Generated event header
  if(fMC)
    return fMC->GenEventHeader();
  else{
    printf("GenEventHeader is not available\n"); 
    return 0x0 ;
  }
}

//_______________________________________________________________
void AliCaloTrackReader::InitParameters()
{
 
  //Initialize the parameters of the analysis.
  fDataType = kESD ;

  fCTSPtMin   = 0.2 ;
  fEMCALPtMin   = 0.5 ;
  fPHOSPtMin   = 0.5 ;

  fFillEMCAL = kTRUE;
  fFillPHOS = kTRUE;
  fFillCTS = kTRUE;
  fFillEMCALCells = kFALSE;
  fFillPHOSCells = kFALSE;

  fFidutialCut = new AliFidutialCut();

}


//________________________________________________________________
void AliCaloTrackReader::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  printf("Data type      : %d\n", fDataType) ;
  printf("CTS Min pT     : %2.1f GeV/c\n", fCTSPtMin) ;
  printf("EMCAL Min pT   : %2.1f GeV/c\n", fEMCALPtMin) ;
  printf("PHOS Min pT    : %2.1f GeV/c\n", fPHOSPtMin) ;
  printf("Use CTS         =     %d\n", fFillCTS) ;
  printf("Use EMCAL       =     %d\n", fFillEMCAL) ;
  printf("Use PHOS        =     %d\n", fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n", fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n", fFillPHOSCells) ;
  printf("    \n") ;
} 

//___________________________________________________
void AliCaloTrackReader::FillInputEvent(Int_t iEntry) {
  //Fill the event counter and input lists that are needed, called by the analysis maker.
  
  fEventNumber = iEntry;
  if(fFillCTS) FillInputCTS();
  if(fFillEMCAL) FillInputEMCAL();
  if(fFillPHOS) FillInputPHOS();
  if(fFillEMCALCells) FillInputEMCALCells();
  if(fFillPHOSCells) FillInputPHOSCells();

}

//__________________________________________________
void AliCaloTrackReader::ResetLists() {
  //  Reset lists, called by the analysis maker 

  if(fAODCTS) fAODCTS -> Clear();
  if(fAODEMCAL) fAODEMCAL -> Clear();
  if(fAODPHOS) fAODPHOS -> Clear();
  if(fEMCALCells) fEMCALCells -> Clear();
  if(fPHOSCells) fPHOSCells -> Clear();

}
