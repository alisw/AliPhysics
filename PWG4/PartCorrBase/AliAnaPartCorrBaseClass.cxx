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
#include "AliCalorimeterUtils.h"
#include "AliCaloPID.h"
#include "AliFiducialCut.h"
#include "AliIsolationCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliNeutralMesonSelection.h"
#include "AliVCaloCells.h" 
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliMixedEvent.h"
#include "AliAODPWG4Particle.h"

ClassImp(AliAnaPartCorrBaseClass)
  
  
//_______________________________________________
  AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass() : 
    TObject(), fDataMC(0), fDebug(0), fCheckFidCut(0),
    fCheckCaloPID(0), fRecalculateCaloPID(0), fMinPt(0), fMaxPt(0),
    fReader(0x0), fInputAODBranch(0x0), fInputAODName(""),
    fOutputAODBranch(0x0), fNewAOD(kFALSE),
    fOutputAODName(""), fOutputAODClassName(""),
    fAODObjArrayName(""), fAddToHistogramsName(""),
    fCaloPID(0x0), fFidCut(0x0), fIC(0x0),fMCUtils(0x0), fNMS(0x0),
    fCaloUtils(0x0),
    fHistoPtBins(0),   fHistoPtMax(0.),   fHistoPtMin(0.),
    fHistoPhiBins(0),  fHistoPhiMax(0.),  fHistoPhiMin(0.),
    fHistoEtaBins(0),  fHistoEtaMax(0.),  fHistoEtaMin(0.),
    fHistoMassBins(0), fHistoMassMax(0.), fHistoMassMin(0.),
    fHistoAsymBins(0), fHistoAsymMax(0.), fHistoAsymMin(0.), fMixedEvent(NULL), fNMixedEvent(1), fVertex(NULL)
{
  //Default Ctor
    
  //Initialize parameters
  InitParameters();
}
/*
//___________________________________________________________
AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass(const AliAnaPartCorrBaseClass & abc) :   
  TObject(), fDataMC(abc.fDataMC), fDebug(abc.fDebug),
  fCheckFidCut(abc.fCheckFidCut),  fCheckCaloPID(abc. fCheckCaloPID),
  fRecalculateCaloPID(abc.fRecalculateCaloPID),
  fMinPt(abc.fMinPt), fMaxPt(abc.fMaxPt), fReader(new AliCaloTrackReader(*abc.fReader)),  
  fInputAODBranch(new TClonesArray(*abc.fInputAODBranch)), fInputAODName(abc.fInputAODName),
  fOutputAODBranch(new TClonesArray(*abc.fOutputAODBranch)),fNewAOD(abc.fNewAOD), 
  fOutputAODName(abc.fOutputAODName), fOutputAODClassName(abc.fOutputAODClassName),
  fAODObjArrayName(abc.fAODObjArrayName),
  fAddToHistogramsName(abc.fAddToHistogramsName),
  fCaloPID(new AliCaloPID(*abc.fCaloPID)), fFidCut(new AliFiducialCut(*abc.fFidCut)), fIC(new AliIsolationCut(*abc.fIC)),
  fMCUtils(new AliMCAnalysisUtils(*abc.fMCUtils)), fNMS(new AliNeutralMesonSelection(*abc.fNMS)),
  fCaloUtils(new AliCalorimeterUtils(*abc.fCaloUtils)),
  fHistoPtBins(abc.fHistoPtBins),     fHistoPtMax(abc.fHistoPtMax),     fHistoPtMin(abc.fHistoPtMin),
  fHistoPhiBins(abc.fHistoPhiBins),   fHistoPhiMax(abc.fHistoPhiMax),   fHistoPhiMin(abc.fHistoPhiMin),
  fHistoEtaBins(abc.fHistoEtaBins),   fHistoEtaMax(abc.fHistoEtaMax),   fHistoEtaMin(abc.fHistoEtaMin),
  fHistoMassBins(abc.fHistoMassBins), fHistoMassMax(abc.fHistoMassMax), fHistoMassMin(abc.fHistoMassMin),
  fHistoAsymBins(abc.fHistoAsymBins), fHistoAsymMax(abc.fHistoAsymMax), fHistoAsymMin(abc.fHistoAsymMin)
{
  // cpy ctor

}
*/
  /*
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
	
  fMinPt   = abc.fMinPt;
  fMaxPt   = abc.fMaxPt;
	
  delete fCaloPID;   fCaloPID   = new AliCaloPID        (*abc.fCaloPID);
  delete fFidCut;    fFidCut    = new AliFiducialCut    (*abc.fFidCut);
  delete fMCUtils;   fMCUtils   = new AliMCAnalysisUtils(*abc.fMCUtils);
  delete fIC;        fIC        = new AliIsolationCut   (*abc.fIC);
  delete fNMS;       fNMS       = new AliNeutralMesonSelection(*abc.fNMS);
  delete fCaloUtils; fCaloUtils = new AliCalorimeterUtils(*abc.fCaloUtils);
  delete fReader;    fReader    = new AliCaloTrackReader(*abc.fReader) ;
	
  delete fInputAODBranch;  fInputAODBranch      = new TClonesArray(*abc.fInputAODBranch) ;
  fInputAODName        = abc.fInputAODName;
  delete fOutputAODBranch; fOutputAODBranch     = new TClonesArray(*abc.fOutputAODBranch) ;
  fNewAOD              = abc.fNewAOD ; 
  fOutputAODName       = abc.fOutputAODName; 
  fOutputAODClassName  = abc.fOutputAODClassName;
  fAddToHistogramsName = abc.fAddToHistogramsName;
  fAODObjArrayName     = abc.fAODObjArrayName;

  fHistoPtBins  = abc.fHistoPtBins;  fHistoPtMax  = abc.fHistoPtMax;  fHistoPtMin  = abc.fHistoPtMin;
  fHistoPhiBins = abc.fHistoPhiBins; fHistoPhiMax = abc.fHistoPhiMax; fHistoPhiMin = abc.fHistoPhiMin;
  fHistoEtaBins = abc.fHistoEtaBins; fHistoEtaMax = abc.fHistoEtaMax; fHistoEtaMin = abc.fHistoEtaMin;
  
  return *this;
  
}
  */
//____________________________________________________________________________
AliAnaPartCorrBaseClass::~AliAnaPartCorrBaseClass() 
{
  // Remove all pointers except analysis output pointers.
//printf("---Delete analysis %s\n", fAODObjArrayName.Data());
  // Not sure if it should be me who deletes the delta AOD branches.
//  if(fOutputAODBranch){
//    fOutputAODBranch->Clear() ; 
//    delete fOutputAODBranch ;
//  }
//  
//  if(fInputAODBranch){
//    fInputAODBranch->Clear() ; 
//    delete fInputAODBranch ;
//  }
  
  //if(fCaloUtils)    delete fCaloUtils ; //Already deleted in maker
  //if(fReader)       delete fReader ;    //Already deleted in maker
	
  if(fCaloPID)   delete fCaloPID ;
  if(fFidCut)    delete fFidCut ;
  if(fIC)        delete fIC ;
  if(fMCUtils)   delete fMCUtils ;
  if(fNMS)       delete fNMS ;
  
  for (Int_t i = 0; i < fNMixedEvent; i++) {
     delete [] fVertex[i] ; 
  }
  delete [] fVertex ;  
  
//   printf("--- analysis deleted \n");
}

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

//___________________________________________________________________________
TClonesArray * AliAnaPartCorrBaseClass::GetAODBranch(TString aodName) const {
	//Recover ouput and input AOD pointers for each event in the maker
	
	//Delta AODs
	if(fDebug > 3) printf("AliAnaPartCorrBaseClass::GetAODBranch() - Get Input Branch with name: <%s>; \n",aodName.Data());
	
  //Get the AOD handler, if output AOD is created use it, if not get the branches from the input which should be deltaAODs
  AliAODHandler* aodHandler = 0x0;
  Bool_t outAOD = kFALSE;
  if((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler()) outAOD = kTRUE;
  if(outAOD) aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler()); 
  else 	     aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  
  if(!GetReader()->WriteDeltaAODToFile())
  {
    return  (TClonesArray *) (fReader->GetAODBranchList())->FindObject(aodName);
  }
  else if (aodHandler->GetExtensions())
  { 
    AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject(GetReader()->GetDeltaAODFileName()); 
    if(ext){
      AliAODEvent *aodEvent = ext->GetAOD(); 
      TClonesArray * aodbranch =  (TClonesArray*) aodEvent->FindListObject(aodName);
      if(aodbranch) return aodbranch;
      else {
        if(outAOD) return  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(aodName);
        else       return  (TClonesArray *) fReader->GetInputEvent() ->FindListObject(aodName);
      }
    }
    else{//If no Delta AODs, kept in standard branch, to revise. 
      if(outAOD) return (TClonesArray *) fReader->GetOutputEvent()->FindListObject(aodName);
      else       return (TClonesArray *) fReader->GetInputEvent() ->FindListObject(aodName);
    }
  }
  else{ //If no Delta AODs, kept in standard branch, to revise. 
    if(outAOD) return (TClonesArray *) fReader->GetOutputEvent()->FindListObject(aodName);
    else       return  (TClonesArray *) fReader->GetInputEvent() ->FindListObject(aodName);
  }
  
}



//___________________________________________________
void AliAnaPartCorrBaseClass::ConnectInputOutputAODBranches() {
  //Recover ouput and input AOD pointers for each event in the maker
	
  //Delta AODs
  if(fDebug > 3) printf("AliAnaPartCorrBaseClass::ConnectInputOutputAODBranches() - Connect Input with name: <%s>; Connect output with name <%s>\n",fInputAODName.Data(),fOutputAODName.Data());
  
 	//Get the AOD handler, if output AOD is created use it, if not get the branches from the input which should be deltaAODs
	AliAODHandler* aodHandler = 0x0;
	Bool_t outAOD = kFALSE;
	if((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler()) outAOD = kTRUE;
	if(outAOD) aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler()); 
	else 	     aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
	
  if(!GetReader()->WriteDeltaAODToFile())
  {
    fOutputAODBranch =  (TClonesArray *) (fReader->GetAODBranchList())->FindObject(fOutputAODName);
    fInputAODBranch  =  (TClonesArray *) (fReader->GetAODBranchList())->FindObject(fInputAODName);	
  }
  else if (aodHandler->GetExtensions()) { 
    
	  AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject(GetReader()->GetDeltaAODFileName()); 
	  if(ext){
      AliAODEvent *aodEvent = ext->GetAOD(); 
      if(fNewAOD)fOutputAODBranch = (TClonesArray*) aodEvent->FindListObject(fOutputAODName);
      fInputAODBranch = (TClonesArray*) aodEvent->FindListObject(fInputAODName); 	  
      if(!fOutputAODBranch && fNewAOD) fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
      if(!fInputAODBranch)  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);
	  }
	  else{//If no Delta AODs, kept in standard branch, to revise. 
		  if(fNewAOD && fReader->GetOutputEvent()) {
			  fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
			  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);	
		  }
		  else {
			  fInputAODBranch  =  (TClonesArray *) fReader->GetInputEvent()->FindListObject(fInputAODName);	
			  if(!fInputAODBranch && fReader->GetOutputEvent() ) 
				  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);//Try the output event.
		  }
	  }
  }
  else{ //If no Delta AODs, kept in standard branch
	  if(fNewAOD && fReader->GetOutputEvent()) {
		  fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
		  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);	
	  }
	  else{ 
		  fInputAODBranch  =  (TClonesArray *) fReader->GetInputEvent()->FindListObject(fInputAODName);
		  if(!fInputAODBranch && fReader->GetOutputEvent())  
			  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);//Try the output event.
	  }
  }
	
  if(GetDebug() > 1){
    if(fNewAOD && !fOutputAODBranch) 
      printf(" AliAnaPartCorrBaseClass::ConnectInputOutputAODBranches() - Output Branch <%s>, not found!\n",fOutputAODName.Data());
    if(!fNewAOD && !fInputAODBranch) 
      printf(" AliAnaPartCorrBaseClass::ConnectInputOutputAODBranches() - Input Branch  <%s>, not found!\n",fInputAODName.Data());
  }
}

//__________________________________________________________________________
Bool_t AliAnaPartCorrBaseClass::IsTrackMatched(AliVCluster* cluster) const {
  //Check if there is any track attached to this cluster
  
  Int_t nMatches = cluster->GetNTracksMatched();
//  printf("N matches %d, first match %d\n",nMatches,cluster->GetTrackMatchedIndex());
//  if     (cluster->GetTrackMatched(0))        printf("\t matched track id %d\n",((AliVTrack*)cluster->GetTrackMatched(0))->GetID()) ;
//  else if(cluster->GetTrackMatchedIndex()>=0) printf("\t matched track id %d\n",((AliVTrack*) GetReader()->GetInputEvent()->GetTrack(cluster->GetTrackMatchedIndex()))->GetID()) ;

  if(fReader->GetDataType()==AliCaloTrackReader::kESD)
  {
    
    if (nMatches > 0) {
      if (nMatches == 1 ) {
        Int_t iESDtrack = cluster->GetTrackMatchedIndex();
        //printf("Track Matched index %d\n",iESDtrack);
        if(iESDtrack==-1) return kFALSE ;// Default value of array, there is no match
        else              return kTRUE;
      }//Just one, check
      else return kTRUE ;//More than one, there is a match.
    }// > 0
    else return kFALSE; //It does not happen, but in case
      
  }//ESDs
  else
  {
    //AODs
    if(nMatches > 0) return kTRUE; //There is at least one match.
    else             return kFALSE;
    
  }//AODs or MC (copy into AOD)
  
  return kFALSE;
  
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
  sprintf(onePar,"fCheckFidCut=%d (Check Fiducial cut selection on/off) \n",fCheckFidCut) ;
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
  fMinPt = 0.1  ; //Min pt in particle analysis
  fMaxPt = 300. ; //Max pt in particle analysis

  //fReader    = new AliCaloTrackReader(); //Initialized in maker
  //fCaloUtils = new AliCalorimeterUtils();//Initialized in maker
  	
  fNewAOD              = kFALSE ;
  fOutputAODName       = "PartCorr";
  fOutputAODClassName  = "AliAODPWG4Particle";
  fInputAODName        = "PartCorr";
  fAddToHistogramsName = "";
  fAODObjArrayName     = "Ref";
	  
  //Histogram settings
  fHistoPtBins    = 240 ;
  fHistoPtMax     = 120 ;
  fHistoPtMin     = 0.  ;

  fHistoPhiBins   = 120 ;
  fHistoPhiMax    = TMath::TwoPi();
  fHistoPhiMin    = 0.  ;

  fHistoEtaBins   = 100 ;
  fHistoEtaMax    =  1  ;
  fHistoEtaMin    = -1  ;

  fHistoMassBins  = 200;
  fHistoMassMax   = 1. ;
  fHistoMassMin   = 0. ;
	
  fHistoAsymBins  = 10 ;
  fHistoAsymMax   = 1. ;
  fHistoAsymMin   = 0. ;
  	
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
  printf("Check Fiducial cut  =     %d\n",     fCheckFidCut) ;
  printf("Check MC labels     =     %d\n",     fDataMC);
  printf("Debug Level         =     %d\n",     fDebug);
  printf("Histograms: %3.1f < pT  < %3.1f,  Nbin = %d\n", fHistoPtMin,  fHistoPtMax,  fHistoPtBins);
  printf("Histograms: %3.1f < phi < %3.1f, Nbin = %d\n", fHistoPhiMin, fHistoPhiMax, fHistoPhiBins);
  printf("Histograms: %3.1f < eta < %3.1f, Nbin = %d\n", fHistoEtaMin, fHistoEtaMax, fHistoEtaBins);
  printf("Histograms: %3.1f < mass < %3.1f, Nbin = %d\n", fHistoMassMin, fHistoMassMax, fHistoMassBins);
  printf("Histograms: %3.1f < asymmetry < %3.1f, Nbin = %d\n", fHistoAsymMin, fHistoAsymMax, fHistoAsymBins);
  printf("Name of reference array      : %s\n", fAODObjArrayName.Data());	
  printf("String added histograms name : %s\n",fAddToHistogramsName.Data());
	
  printf("    \n") ;
  
} 

//__________________________________________________________________
AliMixedEvent * AliAnaPartCorrBaseClass::GetMixedEvent() 
{ 
  //gets the mixed event objects and does some setting
  if (!fMixedEvent) {
    fMixedEvent = dynamic_cast<AliMixedEvent*>(GetReader()->GetInputEvent()) ; 
    if (fMixedEvent) {
      fNMixedEvent = fMixedEvent->GetNumberOfEvents() ; 
    }
    
    //Delete previous vertex
    if(fVertex){
      for (Int_t i = 0; i < fNMixedEvent; i++) {
	delete [] fVertex[i] ; 
      }
      delete [] fVertex ;  
    }
    
    fVertex = new Double_t*[fNMixedEvent] ; 
    for (Int_t i = 0; i < fNMixedEvent; i++) {
      fVertex[i] = new Double_t[3] ; 
      fVertex[i][0] = 0.0 ; 
      fVertex[i][1] = 0.0 ; 
      fVertex[i][2] = 0.0 ; 
    }          
  }  
  return fMixedEvent ; 
}



