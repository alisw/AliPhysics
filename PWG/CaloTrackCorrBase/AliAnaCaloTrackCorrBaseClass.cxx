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
// Base class for CaloTrackCorr analysis algorithms
//-- Author: Gustavo Conesa (LNF-INFN, LPSC-Grenoble) 
//
//
//_________________________________________________________________________


// --- ROOT system ---
#include <TClonesArray.h>
//#include <Riostream.h>

//---- AliRoot system ----
#include "AliAnaCaloTrackCorrBaseClass.h"
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
#include "AliAODPWG4Particle.h"

ClassImp(AliAnaCaloTrackCorrBaseClass)


//__________________________________________________________
AliAnaCaloTrackCorrBaseClass::AliAnaCaloTrackCorrBaseClass() : 
TObject(), 
fDataMC(0),                   fDebug(0),                   fCheckFidCut(0),
fCheckCaloPID(0),             fRecalculateCaloPID(0), 
fMinPt(0),                    fMaxPt(0),                   fPairTimeCut(200), 
fMultiBin(0),                 fNZvertBin(0),
fNrpBin(0),                   fNCentrBin(0),
fNmaxMixEv(0),                fMaxMulti(0),                fMinMulti(0),
fUseSelectEvent(kFALSE),      fMakePlots(kFALSE),
fInputAODBranch(0x0),         fInputAODName(""),
fOutputAODBranch(0x0),        fNewAOD(kFALSE),
fOutputAODName(""),           fOutputAODClassName(""),
fAODObjArrayName(""),         fAddToHistogramsName(""),
fCaloPID(0x0),                fCaloUtils(0x0),
fFidCut(0x0),                 fHisto(0x0),
fIC(0x0),                     fMCUtils(0x0),                
fNMS(0x0),                    fReader(0x0)
{
  //Default Ctor
  
  //Initialize parameters
  InitParameters();
}

//___________________________________________________________
AliAnaCaloTrackCorrBaseClass::~AliAnaCaloTrackCorrBaseClass() 
{
  // Remove all pointers except analysis output pointers.
  
  //delete fCaloUtils ; //Already deleted in maker
  //delete fReader ;    //Already deleted in maker
	
  delete fCaloPID ; 
  delete fFidCut  ;  
  delete fIC      ;      
  delete fMCUtils ; 
  delete fNMS     ;     
  delete fHisto   ;    
}

//______________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::AddAODParticle(AliAODPWG4Particle pc)
{
  //Put AOD calo cluster in the AODParticleCorrelation array
  
  if(fOutputAODBranch){
    
    Int_t i = fOutputAODBranch->GetEntriesFast();
    //new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
    if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4Particle")==0)
      new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
    else   if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4ParticleCorrelation")==0)
      new((*fOutputAODBranch)[i])  AliAODPWG4ParticleCorrelation(pc);
    else {
      printf("AliAnaCaloTrackCorrBaseClass::AddAODParticle() - Cannot add an object of type < %s >, to the AOD TClonesArray \n",  
             fOutputAODBranch->GetClass()->GetName());
      abort();    
    }
  }
  else {
    printf(" AliAnaCaloTrackCorrBaseClass::AddAODParticle() - No AOD branch available!!!\n");
    abort();
  }
  
}	

//_______________________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::CheckMixedEventVertex(const Int_t caloLabel, 
                                                          const Int_t trackLabel)
{
  // Check vertex in mixed events
  
  if (!GetMixedEvent()) return 1; // Not mixed event continue normal processing
  
  Int_t evt = -1;
  
  if     (caloLabel  >= 0 )
  {
    evt = GetMixedEvent()->EventIndexForCaloCluster(caloLabel) ;
  }
  else if(trackLabel >= 0 )
  {
    evt = GetMixedEvent()->EventIndex(trackLabel) ;
  }
  else  
    return 0; // go to next entry in the particle list
  
  if(evt == -1) 
    return 0 ; // to content coverity
  
  if (TMath::Abs(GetVertex(evt)[2]) > GetZvertexCut())  return -1; // Vertex out of range process next event
  
  return 1 ; // continue processing normally
  
}

//________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() 
{
  //Recover ouput and input AOD pointers for each event in the maker
	
  //Delta AODs
  if(fDebug > 3) printf("AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() - Connect Input with name: <%s>; Connect output with name <%s>\n",fInputAODName.Data(),fOutputAODName.Data());
  
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
      printf(" AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() - Output Branch <%s>, not found!\n",fOutputAODName.Data());
    if(!fNewAOD && !fInputAODBranch) 
      printf(" AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() - Input Branch  <%s>, not found!\n",fInputAODName.Data());
  }
}

//__________________________________________________________________________________________
AliVCluster * AliAnaCaloTrackCorrBaseClass::FindCluster(TObjArray* clusters, const Int_t id, 
                                                        Int_t & iclus, const Int_t first) 
{
  // Given the cluster ID stored in AliAODPWG4Particle, get the originator cluster and its index in the array
  
  if(!clusters) return 0x0;
  
  for(iclus = first; iclus < clusters->GetEntriesFast(); iclus++){
    AliVCluster *cluster= dynamic_cast<AliVCluster*> (clusters->At(iclus));
    if(cluster){
      if     (cluster->GetID()==id) {
        return cluster;
      }
    }      
  }// calorimeter clusters loop
  
  return 0x0;
  
}

//______________________________________________________________________________
TClonesArray * AliAnaCaloTrackCorrBaseClass::GetAODBranch(TString aodName) const 
{
	//Recover ouput and input AOD pointers for each event in the maker
	
	//Delta AODs
	if(fDebug > 3) printf("AliAnaCaloTrackCorrBaseClass::GetAODBranch() - Get Input Branch with name: <%s>; \n",aodName.Data());
	
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
    if(outAOD) return (TClonesArray *)  fReader->GetOutputEvent()->FindListObject(aodName);
    else       return  (TClonesArray *) fReader->GetInputEvent() ->FindListObject(aodName);
  }
  
}

//_____________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetCTSTracks() const 
{
  //Get list of referenced tracks from reader
  
  return fReader->GetCTSTracks(); 
  
}

//________________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetPHOSClusters() const 
{
  //Get list of PHOS reference caloclusters from reader
  
  return fReader->GetPHOSClusters(); 
  
}

//_________________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetEMCALClusters() const 
{
  //Get list of emcal referenced caloclusters from reader
  
  return fReader->GetEMCALClusters(); 
  
}

//______________________________________________________________________
TClonesArray *  AliAnaCaloTrackCorrBaseClass::GetAODCaloClusters() const 
{
  //Get list of all caloclusters in AOD output file 
  
  return fReader->GetOutputEvent()->GetCaloClusters(); 
  
}

//________________________________________________________________
TClonesArray *  AliAnaCaloTrackCorrBaseClass::GetAODTracks() const 
{
  //Get list of all tracks in AOD output file 
  
  return fReader->GetOutputEvent()->GetTracks(); 
  
}

//____________________________________________________________
TString  AliAnaCaloTrackCorrBaseClass::GetBaseParametersList()  
{
  //Put data member values in string to keep in output container
  
  TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliAnaCaloTrackCorrBaseClass ---\n") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Minimal P_t: %2.2f ; Max\n", fMinPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Minimal P_t: %2.2f ; Max\n", fMaxPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"|t_{1}-t_{2}| < %2.2f ; Max\n", fPairTimeCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fDataMC =%d (Check MC information, on/off) \n",fDataMC) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckFidCut=%d (Check Fiducial cut selection on/off) \n",fCheckFidCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckCaloPID =%d (Use Bayesian PID in calorimetes, on/off) \n",fCheckCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRecalculateCaloPID  =%d (Calculate PID from shower/tof/tracking parameters, on/off) \n",fRecalculateCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fInputAODName  =%s Input AOD name \n",fInputAODName.Data()) ;
  parList+=onePar ;	
  if(fNewAOD){
    snprintf(onePar,buffersize,"fOutputAODName  =%s Output AOD name \n",fOutputAODName.Data()) ;
    parList+=onePar ;	
    snprintf(onePar,buffersize,"fOutputAODClassName  =%s Output AOD class name \n",fOutputAODClassName.Data()) ;
    parList+=onePar ;	
  }
  snprintf(onePar,buffersize,"fAODObjArrayName  =%s Reference arrays in AOD name \n",fAODObjArrayName.Data()) ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fAddToHistogramsName  =%s String added to beginning of histograms name \n",fAddToHistogramsName.Data()) ;
  parList+=onePar ;	
	
  return parList; 
  
}

//_____________________________________________________________________
TClonesArray * AliAnaCaloTrackCorrBaseClass::GetCreateOutputAODBranch() 
{
  //Create AOD branch filled in the analysis
  
  printf("Create AOD branch of %s objects and with name < %s >\n",
         fOutputAODClassName.Data(),fOutputAODName.Data()) ;
  
  TClonesArray * aodBranch = new TClonesArray(fOutputAODClassName, 0);
  aodBranch->SetName(fOutputAODName);
  return aodBranch ;
  
}

//________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventNumber() const 
{
  //Get current event number
  
  return fReader->GetEventNumber() ; 
  
}


//__________________________________________________________
AliStack *  AliAnaCaloTrackCorrBaseClass::GetMCStack() const 
{
  //Get stack pointer from reader
  
  return fReader->GetStack(); 
  
}
//____________________________________________________________
AliHeader *  AliAnaCaloTrackCorrBaseClass::GetMCHeader() const
{
  //Get header pointer from reader
  
  return fReader->GetHeader(); 
  
}

//____________________________________________________________________________
AliGenEventHeader *  AliAnaCaloTrackCorrBaseClass::GetMCGenEventHeader() const
{
  //Get GenEventHeader pointer from reader
  
  return fReader->GetGenEventHeader(); 
  
}

//_________________________________________________
void AliAnaCaloTrackCorrBaseClass::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  fDataMC              = kFALSE;
  fDebug               = -1;
  fCheckCaloPID        = kTRUE ;
  fCheckFidCut         = kFALSE ;
  fRecalculateCaloPID  = kFALSE ;
  fMinPt               = 0.1  ; //Min pt in particle analysis
  fMaxPt               = 300. ; //Max pt in particle analysis
  fMultiBin            = 1;
  fNZvertBin           = 1;
  fNrpBin              = 1;
  fMaxMulti            = 1000;
  fMinMulti            = 0;
  fUseSelectEvent      = kFALSE ;
  
  //fReader    = new AliCaloTrackReader(); //Initialized in maker
  //fCaloUtils = new AliCalorimeterUtils();//Initialized in maker
  
  fNewAOD              = kFALSE ;
  fOutputAODName       = "CaloTrackCorr";
  fOutputAODClassName  = "AliAODPWG4Particle";
  fInputAODName        = "CaloTrackCorr";
  fAddToHistogramsName = "";
  fAODObjArrayName     = "Ref";
  
}

//__________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::Print(const Option_t * opt) const
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
  printf("Make plots?         =     %d \n",    fMakePlots); 	
  printf("Debug Level         =     %d\n",     fDebug);
  
  printf("Name of reference array      : %s\n", fAODObjArrayName.Data());	
  printf("String added histograms name : %s\n", fAddToHistogramsName.Data());
	
  printf("    \n") ;
  
} 



