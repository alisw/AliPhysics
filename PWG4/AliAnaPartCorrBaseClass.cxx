/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *GetEntriesFast(
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

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
#include "AliAnaPartCorrBaseClass.h"
#include "AliCaloTrackReader.h"
#include "AliFidutialCut.h"
#include "AliIsolationCut.h"
#include "AliNeutralMesonSelection.h"
#include "AliLog.h"
// #include "AliStack.h"
// #include "AliHeader.h"
// #include "AliGenEventHeader.h"

ClassImp(AliAnaPartCorrBaseClass)
  
  
//_______________________________________________
  AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass() : 
    TObject(), fDataMC(0), fDebug(0), fCheckFidCut(0),
    fCheckCaloPID(0), fRecalculateCaloPID(0), fMinPt(0), fMaxPt(0),
    fReader(0x0), fAODBranch(0x0),  fAODCaloClusters(0x0), fAODCaloCells(0x0), 
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
AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & abc) :   
  TObject(), fDataMC(abc.fDataMC), fDebug(abc.fDebug),
  fCheckFidCut(abc.fCheckFidCut),  fCheckCaloPID(abc. fCheckCaloPID),
  fRecalculateCaloPID(abc.fRecalculateCaloPID),
  fMinPt(abc.fMinPt), fMaxPt(abc.fMaxPt), fReader(abc.fReader),  
  fAODBranch(new TClonesArray(*abc.fAODBranch)),
  fAODCaloClusters(new TClonesArray(*abc.fAODCaloClusters)),
  fAODCaloCells(new AliAODCaloCells(*abc.fAODCaloCells)),
  fCaloPID(abc.fCaloPID), fFidCut(abc.fFidCut), fIC(abc.fIC),fNMS(abc.fNMS)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaPartCorrBaseClass & AliAnaPartCorrBaseClass::operator = (const AliAnaPartCorrBaseClass & abc)
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

  fMinPt = abc.fMinPt;
  fMaxPt = abc.fMaxPt;

  fCaloPID = abc.fCaloPID;  
  fFidCut = abc.fFidCut;
  fIC = abc.fIC;
  fNMS = abc.fNMS;

  return *this;
  
}

//____________________________________________________________________________
AliAnaPartCorrBaseClass::~AliAnaPartCorrBaseClass() 
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
void AliAnaPartCorrBaseClass::AddAODCaloCluster(AliAODCaloCluster calo) {
  //Put AOD calo cluster in the CaloClusters array

  Int_t i = fAODCaloClusters->GetEntriesFast();
  new((*fAODCaloClusters)[i])  AliAODCaloCluster(calo);

}

//____________________________________________________________________________
void AliAnaPartCorrBaseClass::AddAODParticleCorrelation(AliAODParticleCorrelation pc) {
  //Put AOD calo cluster in the AODParticleCorrelation array

  Int_t i = fAODBranch->GetEntriesFast();
  new((*fAODBranch)[i])  AliAODParticleCorrelation(pc);

}

//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectAODCaloClusters() {
  //Recover the list of AODCaloClusters

  fAODCaloClusters = fReader->GetAOD()->GetCaloClusters();

}

//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectAODPHOSCells() {
  //Recover the list of PHOS AODCaloCells 

  fAODCaloCells = fReader->GetAOD()->GetPHOSCells();

}

//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectAODEMCALCells() {
  //Recover the list of EMCAL AODCaloCells 

  fAODCaloCells = fReader->GetAOD()->GetEMCALCells();

}

//__________________________________________________
TClonesArray *  AliAnaPartCorrBaseClass::GetAODCTS() const {
  //Get list of tracks from reader

  return fReader->GetAODCTS(); 

}

//__________________________________________________
TClonesArray *  AliAnaPartCorrBaseClass::GetAODPHOS() const {
  //Get list of PHOS calo clusters from reader

  return fReader->GetAODPHOS(); 

}


//__________________________________________________
TClonesArray *  AliAnaPartCorrBaseClass::GetAODEMCAL() const {
  //Get list of emcal caloclusters from reader

  return fReader->GetAODEMCAL(); 

}

//__________________________________________________
TNamed *  AliAnaPartCorrBaseClass::GetPHOSCells() const {
  //Get list of PHOS calo cells (ESD or AOD) from reader
  
  return fReader->GetPHOSCells(); 
  
}


//__________________________________________________
TNamed *  AliAnaPartCorrBaseClass::GetEMCALCells() const {
  //Get list of emcal calo cells (ESD or AOD) from reader
  
  return fReader->GetEMCALCells(); 

}

//__________________________________________________
AliStack *  AliAnaPartCorrBaseClass::GetMCStack() const {
  //Get stack pointer from reader

  return fReader->GetStack(); 

}
//__________________________________________________
AliHeader *  AliAnaPartCorrBaseClass::GetMCHeader() const {
  //Get header pointer from reader

  return fReader->GetHeader(); 

}

//__________________________________________________
AliGenEventHeader *  AliAnaPartCorrBaseClass::GetMCGenEventHeader() const {
  //Get GenEventHeader pointer from reader

  return fReader->GetGenEventHeader(); 

}


void AliAnaPartCorrBaseClass::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  fDataMC = kFALSE;
  fDebug = -1;
  fCheckCaloPID = kTRUE ;
  fCheckFidCut = kFALSE ;
  fRecalculateCaloPID = kFALSE ;
  fMinPt = 0.2 ; //Min pt in particle analysis
  fMaxPt = 300. ; //Max pt in particle analysis

  fCaloPID = new AliCaloPID ;  
  fFidCut = new AliFidutialCut;
  fIC = new AliIsolationCut;
  fNMS = new AliNeutralMesonSelection;

}

