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

//_________________________________________________________________________
// Base class for analysis algorithms
//-- Author: Gustavo Conesa (LNF-INFN) 
//_________________________________________________________________________
  

// --- ROOT system ---
#include <TClonesArray.h>
//#include <Riostream.h>

//---- AliRoot system ----
#include "AliAnaPartCorrBaseClass.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliFidutialCut.h"
#include "AliIsolationCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliAODCaloCells.h" 
#include "AliAODEvent.h"

ClassImp(AliAnaPartCorrBaseClass)
  
  
//_______________________________________________
  AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass() : 
    TObject(), fDataMC(0), fDebug(0), fCheckFidCut(0),
    fCheckCaloPID(0), fRecalculateCaloPID(0), fMinPt(0), fMaxPt(0),
    fReader(0x0), fInputAODBranch(0x0), fInputAODName(""),
    fOutputAODBranch(0x0), fNewAOD(kFALSE),
    fOutputAODName(""), fOutputAODClassName(""),
    fAODObjArrayName(""), fAddToHistogramsName(""),
    fAODCaloCells(0x0),//fAODCaloClusters(0x0),  
    fCaloPID(0x0), fFidCut(0x0), fIC(0x0),fMCUtils(0x0), fNMS(0x0),
    fHistoNPtBins(0),  fHistoPtMax(0.),  fHistoPtMin(0.),
    fHistoNPhiBins(0), fHistoPhiMax(0.), fHistoPhiMin(0.),
    fHistoNEtaBins(0), fHistoEtaMax(0.), fHistoEtaMin(0.)
{
  //Default Ctor
  
  fReader  = new AliCaloTrackReader();
  fCaloPID = new AliCaloPID();
  fFidCut  = new AliFidutialCut();
  fIC      = new AliIsolationCut();
  fMCUtils = new AliMCAnalysisUtils();
  
  //Initialize parameters
  InitParameters();
}

//___________________________________________________________
AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & abc) :   
  TObject(), fDataMC(abc.fDataMC), fDebug(abc.fDebug),
  fCheckFidCut(abc.fCheckFidCut),  fCheckCaloPID(abc. fCheckCaloPID),
  fRecalculateCaloPID(abc.fRecalculateCaloPID),
  fMinPt(abc.fMinPt), fMaxPt(abc.fMaxPt), fReader(abc.fReader),  
  fInputAODBranch(new TClonesArray(*abc.fInputAODBranch)), fInputAODName(abc.fInputAODName),
  fOutputAODBranch(new TClonesArray(*abc.fOutputAODBranch)),fNewAOD(abc.fNewAOD), 
  fOutputAODName(abc.fOutputAODName), fOutputAODClassName(abc.fOutputAODClassName),
  fAODObjArrayName(abc.fAODObjArrayName),
  fAddToHistogramsName(abc.fAddToHistogramsName),
  //fAODCaloClusters(new TClonesArray(*abc.fAODCaloClusters)),
  fAODCaloCells(new AliAODCaloCells(*abc.fAODCaloCells)),
  fCaloPID(abc.fCaloPID), fFidCut(abc.fFidCut), fIC(abc.fIC),fMCUtils(abc.fMCUtils), fNMS(abc.fNMS),
  fHistoNPtBins(abc.fHistoNPtBins),   fHistoPtMax(abc.fHistoPtMax),   fHistoPtMin(abc.fHistoPtMin),
  fHistoNPhiBins(abc.fHistoNPhiBins), fHistoPhiMax(abc.fHistoPhiMax), fHistoPhiMin(abc.fHistoPhiMin),
  fHistoNEtaBins(abc.fHistoNEtaBins), fHistoEtaMax(abc.fHistoEtaMax), fHistoEtaMin(abc.fHistoEtaMin)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaPartCorrBaseClass & AliAnaPartCorrBaseClass::operator = (const AliAnaPartCorrBaseClass & abc)
{
  // assignment operator
  
  if(this == &abc) return *this;
  ((TObject *)this)->operator=(abc);
  
  fDataMC             = abc.fDataMC;
  fDebug              = abc.fDebug ;
  fRecalculateCaloPID = abc.fRecalculateCaloPID ;
  fCheckCaloPID       = abc.fCheckCaloPID ;
  fCheckFidCut        = abc.fCheckFidCut ; 
	
  fReader             = abc.fReader ;
  //fAODCaloClusters   = new TClonesArray(*abc.fAODCaloClusters) ;
  fAODCaloCells      = new AliAODCaloCells(*abc.fAODCaloCells) ;
  
  fMinPt   = abc.fMinPt;
  fMaxPt   = abc.fMaxPt;
  fCaloPID = abc.fCaloPID;  
  fFidCut  = abc.fFidCut;
  fIC      = abc.fIC;
  fMCUtils = abc.fMCUtils;
  fNMS     = abc.fNMS;
  
  fInputAODBranch      = new TClonesArray(*abc.fInputAODBranch) ;
  fInputAODName        = abc.fInputAODName;
  fOutputAODBranch     = new TClonesArray(*abc.fOutputAODBranch) ;
  fNewAOD              = abc.fNewAOD ; 
  fOutputAODName       = abc.fOutputAODName; 
  fOutputAODClassName  = abc.fOutputAODClassName;
  fAddToHistogramsName = abc.fAddToHistogramsName;
  fAODObjArrayName     = abc.fAODObjArrayName;

  fHistoNPtBins  = abc.fHistoNPtBins;  fHistoPtMax  = abc.fHistoPtMax;  fHistoPtMin  = abc.fHistoPtMin;
  fHistoNPhiBins = abc.fHistoNPhiBins; fHistoPhiMax = abc.fHistoPhiMax; fHistoPhiMin = abc.fHistoPhiMin;
  fHistoNEtaBins = abc.fHistoNEtaBins; fHistoEtaMax = abc.fHistoEtaMax; fHistoEtaMin = abc.fHistoEtaMin;
  
  return *this;
  
}

//____________________________________________________________________________
AliAnaPartCorrBaseClass::~AliAnaPartCorrBaseClass() 
{
  // Remove all pointers except analysis output pointers.
  	
  if(fOutputAODBranch){
    fOutputAODBranch->Clear() ; 
    delete fOutputAODBranch ;
  }
  
  if(fInputAODBranch){
    fInputAODBranch->Clear() ; 
    delete fInputAODBranch ;
  }
  
//	if(fAODCaloClusters){
//		fAODCaloClusters->Clear() ; 
//		delete fAODCaloClusters ;
//	}
	
  if(fAODCaloCells){
    fAODCaloCells->Clear() ; 
    delete fAODCaloCells ;
  }
  
  if(fReader)  delete fReader ;
  if(fCaloPID) delete fCaloPID ;
  if(fFidCut)  delete fFidCut ;
  if(fIC)      delete fIC ;
  if(fMCUtils) delete fMCUtils ;
  if(fNMS)     delete fNMS ;
  
}

////____________________________________________________________________________
//void AliAnaPartCorrBaseClass::AddAODCaloCluster(AliAODCaloCluster calo) {
//  //Put AOD calo cluster in the CaloClusters array
//
//  Int_t i = fAODCaloClusters->GetEntriesFast();
//  new((*fAODCaloClusters)[i])  AliAODCaloCluster(calo);
//
//}


//____________________________________________________________________________
void AliAnaPartCorrBaseClass::AddAODParticle(AliAODPWG4Particle pc) {
  //Put AOD calo cluster in the AODParticleCorrelation array
  
  if(fOutputAODBranch){
    
    Int_t i = fOutputAODBranch->GetEntriesFast();
    //new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
    if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4Particle")==0)
      new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
    else   if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4ParticleCorrelation")==0)
      new((*fOutputAODBranch)[i])  AliAODPWG4ParticleCorrelation(pc);
    else {
      printf("AliAnaPartCorrBaseClass::AddAODParticle() - Cannot add an object of type < %s >, to the AOD TClonesArray \n",  
	     fOutputAODBranch->GetClass()->GetName());
      abort();    
    }
  }
  else {
    printf(" AliAnaPartCorrBaseClass::AddAODParticle() - No AOD branch available!!!\n");
    abort();
  }

}	


//___________________________________________________
//void AliAnaPartCorrBaseClass::ConnectAODCaloClusters() {
//  //Recover the list of AODCaloClusters
//
//  fAODCaloClusters = fReader->GetOutputEvent()->GetCaloClusters();
//
//}
//
//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectAODPHOSCells() {
  //Recover the list of PHOS AODCaloCells 

  fAODCaloCells = fReader->GetOutputEvent()->GetPHOSCells();

}

//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectAODEMCALCells() {
  //Recover the list of EMCAL AODCaloCells 

  fAODCaloCells = fReader->GetOutputEvent()->GetEMCALCells();

}

//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectInputOutputAODBranches() {
  //Recover ouput and input AOD pointers for each event in the maker

  fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);		

}

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetAODCTS() const {
  //Get list of referenced tracks from reader

  return fReader->GetAODCTS(); 

}

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetAODPHOS() const {
  //Get list of PHOS reference caloclusters from reader

  return fReader->GetAODPHOS(); 

}

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetAODEMCAL() const {
  //Get list of emcal referenced caloclusters from reader

  return fReader->GetAODEMCAL(); 

}

//__________________________________________________
TClonesArray *  AliAnaPartCorrBaseClass::GetAODCaloClusters() const {
  //Get list of all caloclusters in AOD output file 

  return fReader->GetOutputEvent()->GetCaloClusters(); 

}

//__________________________________________________
TClonesArray *  AliAnaPartCorrBaseClass::GetAODTracks() const {
  //Get list of all tracks in AOD output file 

  return fReader->GetOutputEvent()->GetTracks(); 

}

//__________________________________________________
TString  AliAnaPartCorrBaseClass::GetBaseParametersList()  {
  //Put data member values in string to keep in output container

  TString parList ; //this will be list of parameters used for this analysis.
  char onePar[255] ;
  sprintf(onePar,"--- AliAnaPartCorrBaseClass ---\n") ;
  parList+=onePar ;	
  sprintf(onePar,"Minimal P_t: %2.2f ; Max\n", fMinPt) ;
  parList+=onePar ;
  sprintf(onePar,"Minimal P_t: %2.2f ; Max\n", fMaxPt) ;
  parList+=onePar ;
  sprintf(onePar,"fDataMC =%d (Check MC information, on/off) \n",fDataMC) ;
  parList+=onePar ;
  sprintf(onePar,"fCheckFidCut=%d (Check Fidutial cut selection on/off) \n",fCheckFidCut) ;
  parList+=onePar ;
  sprintf(onePar,"fCheckCaloPID =%d (Use Bayesian PID in calorimetes, on/off) \n",fCheckCaloPID) ;
  parList+=onePar ;
  sprintf(onePar,"fRecalculateCaloPID  =%d (Calculate PID from shower/tof/tracking parameters, on/off) \n",fRecalculateCaloPID) ;
  parList+=onePar ;
  sprintf(onePar,"fInputAODName  =%s Input AOD name \n",fInputAODName.Data()) ;
  parList+=onePar ;	
  if(fNewAOD){
     sprintf(onePar,"fOutputAODName  =%s Output AOD name \n",fOutputAODName.Data()) ;
     parList+=onePar ;	
	 sprintf(onePar,"fOutputAODClassName  =%s Output AOD class name \n",fOutputAODClassName.Data()) ;
	 parList+=onePar ;	
  }
  sprintf(onePar,"fAODObjArrayName  =%s Reference arrays in AOD name \n",fAODObjArrayName.Data()) ;
  parList+=onePar ;	
  sprintf(onePar,"fAddToHistogramsName  =%s String added to beginning of histograms name \n",fAddToHistogramsName.Data()) ;
  parList+=onePar ;	
	
  return parList; 

}

//__________________________________________________
 TClonesArray * AliAnaPartCorrBaseClass::GetCreateOutputAODBranch() {
   //Create AOD branch filled in the analysis
   
   printf("Create AOD branch of %s objects and with name < %s >\n",
	  fOutputAODClassName.Data(),fOutputAODName.Data()) ;
   
   TClonesArray * aodBranch = new TClonesArray(fOutputAODClassName, 0);
   aodBranch->SetName(fOutputAODName);
   return aodBranch ;
   
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
Int_t AliAnaPartCorrBaseClass::GetEventNumber() const {
	//Get current event number
	
	return fReader->GetEventNumber() ; 
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

//__________________________________________________
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
  fNewAOD = kFALSE ;
  fOutputAODName = "PartCorr";
  fOutputAODClassName = "AliAODPWG4Particle";
  fInputAODName = "PartCorr";
  fAddToHistogramsName = "";
  fAODObjArrayName="Ref";
	  
  //Histogrammes settings
  fHistoNPtBins = 240 ;
  fHistoPtMax   = 120 ;
  fHistoPtMin   = 0.  ;

  fHistoNPhiBins = 120 ;
  fHistoPhiMax   = TMath::TwoPi();
  fHistoPhiMin   = 0.  ;

  fHistoNEtaBins = 100 ;
  fHistoEtaMax   =  1  ;
  fHistoEtaMin   = -1  ;

}

//__________________________________________________________________
void AliAnaPartCorrBaseClass::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
	
  printf("New AOD:            =     %d\n",fNewAOD);
  printf("Input AOD name:     =     %s\n",fInputAODName.Data());
  printf("Output AOD name:    =     %s\n",fOutputAODName.Data());
  printf("Output AOD Class name: =  %s\n",fOutputAODClassName.Data());
  printf("Min Photon pT       =     %2.2f\n",  fMinPt) ;
  printf("Max Photon pT       =     %3.2f\n",  fMaxPt) ;
  printf("Check PID           =     %d\n",     fCheckCaloPID) ;
  printf("Recalculate PID     =     %d\n",     fRecalculateCaloPID) ;
  printf("Check Fidutial cut  =     %d\n",     fCheckFidCut) ;
  printf("Check MC labels     =     %d\n",     fDataMC);
  printf("Debug Level         =     %d\n",     fDebug);
  printf("Histograms: %3.1f < pT < %3.1f,  Nbin = %d\n", fHistoPtMin,  fHistoPtMax,  fHistoNPtBins);
  printf("Histograms: %3.1f < phi < %3.1f, Nbin = %d\n", fHistoPhiMin, fHistoPhiMax, fHistoNPhiBins);
  printf("Histograms: %3.1f < eta < %3.1f, Nbin = %d\n", fHistoEtaMin, fHistoEtaMax, fHistoNEtaBins);
  printf("Name of reference array      : %s\n", fAODObjArrayName.Data());	
  printf("String added histograms name : %s\n",fAddToHistogramsName.Data());
	
  printf("    \n") ;
  
} 
