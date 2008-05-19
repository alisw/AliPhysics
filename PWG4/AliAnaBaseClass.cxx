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
/* $Id: $ */

/* History of cvs commits:
 *
 * $Log$
 *
 */
//_________________________________________________________________________
// Base class for analysis algorithms
//-- Author: Gustavo Conesa (LNF-INFN) 
//_________________________________________________________________________
  

// --- ROOT system ---
#include <TClonesArray.h>
#include <Riostream.h>

//---- AliRoot system ----
#include "AliAODParticleCorrelation.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliAODTrack.h"
#include "AliAnaBaseClass.h"
#include "AliCaloTrackReader.h"
#include "AliFidutialCut.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliLog.h"
// #include "AliStack.h"
// #include "AliHeader.h"
// #include "AliGenEventHeader.h"

ClassImp(AliAnaBaseClass)
  
  
//_______________________________________________
  AliAnaBaseClass::AliAnaBaseClass() : 
    TObject(), fDataMC(0), fDebug(0), fCheckFidCut(0),
    fCheckCaloPID(0), fRecalculateCaloPID(0), fReader(0x0), 
    fAODBranch(0x0),  fAODCaloClusters(0x0), fAODCaloCells(0x0), 
    fCaloPID(0x0), fFidCut(0x0), fIC(0x0),fNMS(0x0) 
{
  //Default Ctor
  
  fReader = new AliCaloTrackReader();
  fCaloPID = new AliCaloPID();
  fFidCut = new AliFidutialCut();
  fIC = new AliIsolationCut();
  
  //Initialize parameters
  InitParameters();
}

//___________________________________________________________
AliAnaBaseClass::AliAnaBaseClass(const AliAnaBaseClass & abc) :   
  TObject(), fDataMC(abc.fDataMC), fDebug(abc.fDebug),
  fCheckFidCut(abc.fCheckFidCut),  fCheckCaloPID(abc. fCheckCaloPID),
  fRecalculateCaloPID(abc.fRecalculateCaloPID),fReader(abc.fReader),  
  fAODBranch(new TClonesArray(*abc.fAODBranch)),
  fAODCaloClusters(new TClonesArray(*abc.fAODCaloClusters)),
  fAODCaloCells(new AliAODCaloCells(*abc.fAODCaloCells)),
  fCaloPID(abc.fCaloPID), fFidCut(abc.fFidCut), fIC(abc.fIC),fNMS(abc.fNMS)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaBaseClass & AliAnaBaseClass::operator = (const AliAnaBaseClass & abc)
{
  // assignment operator
  
  if(this == &abc) return *this;
  ((TObject *)this)->operator=(abc);
  
  fDataMC = abc.fDataMC;
  fDebug = abc.fDebug ;
  fRecalculateCaloPID =  abc.fRecalculateCaloPID ;
  fCheckCaloPID = abc. fCheckCaloPID ;
  fCheckFidCut = abc.fCheckFidCut ; 

  fReader = abc.fReader ;
  
  fAODBranch = new TClonesArray(*abc.fAODBranch) ;
  fAODCaloClusters = new TClonesArray(*abc.fAODCaloClusters) ;
  fAODCaloCells = new AliAODCaloCells(*abc.fAODCaloCells) ;

  fCaloPID = abc.fCaloPID;  
  fFidCut = abc.fFidCut;
  fIC = abc.fIC;
  fNMS = abc.fNMS;

  return *this;
  
}

//____________________________________________________________________________
AliAnaBaseClass::~AliAnaBaseClass() 
{
  // Remove all pointers except analysis output pointers.
  
  if(fAODBranch){
    fAODBranch->Clear() ; 
    delete fAODBranch ;
  }
  
  if(fAODCaloClusters){
    fAODCaloClusters->Clear() ; 
    delete fAODCaloClusters ;
  }

  if(fAODCaloCells){
    fAODCaloCells->Clear() ; 
    delete fAODCaloCells ;
  }
  
  if(fReader) delete fReader ;
  if(fCaloPID) delete fCaloPID ;
  if(fFidCut) delete fFidCut ;
  if(fIC) delete fIC ;
  if(fNMS) delete fNMS ;

}

//____________________________________________________________________________
void AliAnaBaseClass::AddAODCaloCluster(AliAODCaloCluster calo) {
  //Put AOD calo cluster in the CaloClusters array

  Int_t i = fAODCaloClusters->GetEntries();
  new((*fAODCaloClusters)[i])  AliAODCaloCluster(calo);

}

//____________________________________________________________________________
void AliAnaBaseClass::AddAODParticleCorrelation(AliAODParticleCorrelation pc) {
  //Put AOD calo cluster in the AODParticleCorrelation array

  Int_t i = fAODBranch->GetEntries();
  new((*fAODBranch)[i])  AliAODParticleCorrelation(pc);

}

//___________________________________________________
void AliAnaBaseClass::ConnectAODCaloClusters() {
  //Recover the list of AODCaloClusters

  fAODCaloClusters = fReader->GetAOD()->GetCaloClusters();

}

//___________________________________________________
void AliAnaBaseClass::ConnectAODPHOSCells() {
  //Recover the list of PHOS AODCaloCells 

  fAODCaloCells = fReader->GetAOD()->GetPHOSCells();

}

//___________________________________________________
void AliAnaBaseClass::ConnectAODEMCALCells() {
  //Recover the list of EMCAL AODCaloCells 

  fAODCaloCells = fReader->GetAOD()->GetEMCALCells();

}

//__________________________________________________
TClonesArray *  AliAnaBaseClass::GetAODCTS() const {
  //Get list of tracks from reader

  return fReader->GetAODCTS(); 

}

//__________________________________________________
TClonesArray *  AliAnaBaseClass::GetAODPHOS() const {
  //Get list of PHOS calo clusters from reader

  return fReader->GetAODPHOS(); 

}


//__________________________________________________
TClonesArray *  AliAnaBaseClass::GetAODEMCAL() const {
  //Get list of emcal caloclusters from reader

  return fReader->GetAODEMCAL(); 

}

//__________________________________________________
TNamed *  AliAnaBaseClass::GetPHOSCells() const {
  //Get list of PHOS calo cells (ESD or AOD) from reader
  
  return fReader->GetPHOSCells(); 
  
}


//__________________________________________________
TNamed *  AliAnaBaseClass::GetEMCALCells() const {
  //Get list of emcal calo cells (ESD or AOD) from reader
  
  return fReader->GetEMCALCells(); 

}

//__________________________________________________
AliStack *  AliAnaBaseClass::GetMCStack() const {
  //Get stack pointer from reader

  return fReader->GetStack(); 

}
//__________________________________________________
AliHeader *  AliAnaBaseClass::GetMCHeader() const {
  //Get header pointer from reader

  return fReader->GetHeader(); 

}

//__________________________________________________
AliGenEventHeader *  AliAnaBaseClass::GetMCGenEventHeader() const {
  //Get GenEventHeader pointer from reader

  return fReader->GetGenEventHeader(); 

}


void AliAnaBaseClass::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  fDataMC = kFALSE;
  fDebug = -1;
  fCheckCaloPID = kTRUE ;
  fCheckFidCut = kFALSE ;
  fRecalculateCaloPID = kFALSE ;
  fMinPt = 2. ; //Min pt in particle analysis
  fMaxPt = 300. ; //Max pt in particle analysis

  fCaloPID = new AliCaloPID ;  
  fFidCut = new AliFidutialCut;
  fIC = new AliIsolationCut;
  fNMS = new AliNeutralMesonSelection;

}

