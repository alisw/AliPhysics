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
// --Yaxian Mao: Add the possibality for event selection analysis based on vertex and multiplicity bins (10/10/2010)
  

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
#include "AliAODPWG4Particle.h"

ClassImp(AliAnaPartCorrBaseClass)
  
  
//_______________________________________________
  AliAnaPartCorrBaseClass::AliAnaPartCorrBaseClass() : 
    TObject(), 
    fDataMC(0),                   fDebug(0),                   fCheckFidCut(0),
    fCheckCaloPID(0),             fRecalculateCaloPID(0), 
    fMinPt(0),                    fMaxPt(0),                   fPairTimeCut(200), 
    fMultiBin(0),                 fNZvertBin(0),
    fNrpBin(0),                   fNCentrBin(0),
    fNmaxMixEv(0),                fMaxMulti(0),                fMinMulti(0),
    fUseSelectEvent(kFALSE),      fMakePlots(kFALSE),
    fReader(0x0),                 fInputAODBranch(0x0),        fInputAODName(""),
    fOutputAODBranch(0x0),        fNewAOD(kFALSE),
    fOutputAODName(""),           fOutputAODClassName(""),
    fAODObjArrayName(""),         fAddToHistogramsName(""),
    fCaloPID(0x0),                fFidCut(0x0),                fIC(0x0),
    fMCUtils(0x0),                fNMS(0x0),                   fCaloUtils(0x0),
    fHistoPtBins(0),              fHistoPtMax(0.),             fHistoPtMin(0.),
    fHistoPhiBins(0),             fHistoPhiMax(0.),            fHistoPhiMin(0.),
    fHistoEtaBins(0),             fHistoEtaMax(0.),            fHistoEtaMin(0.),
    fHistoMassBins(0),            fHistoMassMax(0.),           fHistoMassMin(0.),
    fHistoAsymBins(0),            fHistoAsymMax(0.),           fHistoAsymMin(0.),
    fHistoV0SBins(0),             fHistoV0SMax(0),             fHistoV0SMin(0),
    fHistoV0MBins(0),             fHistoV0MMax(0),             fHistoV0MMin(0),
    fHistoTrMBins(0),             fHistoTrMMax(0),             fHistoTrMMin(0),
    fHistoFinePtBins(1000),       fHistoFinePtMax(5.),         fHistoFinePtMin(0.),
    fHistoPOverEBins(100),        fHistoPOverEMax(100.),       fHistoPOverEMin(0.),
    fHistodEdxBins(100),          fHistodEdxMax(100.),         fHistodEdxMin(0.),
    fHistodRBins(100),            fHistodRMax(100.),           fHistodRMin(0.),
    fHistoTimeBins(100),          fHistoTimeMax(100.),         fHistoTimeMin(0.),
    fHistoNClusCellBins(100),     fHistoNClusCellMax(100),     fHistoNClusCellMin(0),
    fHistoNCellsBins(100),        fHistoNCellsMax(100),        fHistoNCellsMin(0),
    fHistoNClustersBins(100),     fHistoNClustersMax(100),     fHistoNClustersMin(0),
    fHistoRatioBins(100),         fHistoRatioMax(100.),        fHistoRatioMin(0.),
    fHistoVertexDistBins(100),    fHistoVertexDistMax(100.),   fHistoVertexDistMin(0.),
    fHistoRBins(100),             fHistoRMax(1000),            fHistoRMin(-1000),
    fHistoXBins(100),             fHistoXMax(1000),            fHistoXMin(-1000),
    fHistoYBins(100),             fHistoYMax(1000),            fHistoYMin(-1000),
    fHistoZBins(100),             fHistoZMax(1000),            fHistoZMin(-1000),
    fHistoSSBins(0),              fHistoSSMax(0),              fHistoSSMin(0),
    fHistoDiffTimeBins(0),        fHistoDiffTimeMax(0),        fHistoDiffTimeMin(0)
{
  //Default Ctor
    
  //Initialize parameters
  InitParameters();
}

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
    if(outAOD) return (TClonesArray *)  fReader->GetOutputEvent()->FindListObject(aodName);
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

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetCTSTracks() const {
  //Get list of referenced tracks from reader

  return fReader->GetCTSTracks(); 

}

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetPHOSClusters() const {
  //Get list of PHOS reference caloclusters from reader

  return fReader->GetPHOSClusters(); 

}

//__________________________________________________
TObjArray *  AliAnaPartCorrBaseClass::GetEMCALClusters() const {
  //Get list of emcal referenced caloclusters from reader

  return fReader->GetEMCALClusters(); 

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
  const Int_t buffersize = 255;
  char onePar[buffersize] ;
  snprintf(onePar,buffersize,"--- AliAnaPartCorrBaseClass ---\n") ;
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
  
  fHistoV0SBins   = 100 ;
  fHistoV0SMax    = 10000  ;
  fHistoV0SMin    = 0  ;
  
  fHistoV0MBins  = 100;
  fHistoV0MMax   = 10000 ;
  fHistoV0MMin   = 0 ;
	
  fHistoTrMBins  = 100 ;
  fHistoTrMMax   = 2000 ;
  fHistoTrMMin   = 0 ;
  
  fHistoPOverEBins       = 100 ;  fHistoPOverEMax       = 10.  ;  fHistoPOverEMin       = 0. ;
  fHistodEdxBins         = 200 ;  fHistodEdxMax         = 400. ;  fHistodEdxMin         = 0. ;  
  fHistodRBins           = 300 ;  fHistodRMax           = 3.15 ;  fHistodRMin           = 0. ;
  fHistoTimeBins         = 1000;  fHistoTimeMax         = 1.e3 ;  fHistoTimeMin         = 0. ;//ns
  fHistoNClusCellBins    = 200 ;  fHistoNClusCellMax    = 200  ;  fHistoNClusCellMin    = 0  ;
  fHistoNCellsBins       = 300 ;  fHistoNCellsMax       = 300  ;  fHistoNCellsMin       = 0  ;
  fHistoNClustersBins    = 50  ;  fHistoNClustersMax    = 50   ;  fHistoNClustersMin    = 0  ;
  fHistoRatioBins        = 200 ;  fHistoRatioMax        = 2    ;  fHistoRatioMin        = 0. ;
  fHistoVertexDistBins   = 100 ;  fHistoVertexDistMax   = 500. ;  fHistoVertexDistMin   = 0. ;
  fHistoRBins            = 100 ;  fHistoRMax            = 500  ;  fHistoRMin            = -500  ;//cm
  fHistoXBins            = 100 ;  fHistoXMax            = 500  ;  fHistoXMin            = -500  ;//cm
  fHistoYBins            = 100 ;  fHistoYMax            = 500  ;  fHistoYMin            = -500  ;//cm
  fHistoZBins            = 100 ;  fHistoZMax            = 600  ;  fHistoZMin            = -500  ;//cm
  fHistoSSBins           = 500 ;  fHistoSSMax           = 5    ;  fHistoSSMin           = 0  ;  
  fHistoDiffTimeBins     = 800 ;  fHistoDiffTimeMax     = 400  ;  fHistoDiffTimeMin     = -400; // ns

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
  printf("Make plots?         =     %d \n",    fMakePlots); 	
  printf("Debug Level         =     %d\n",     fDebug);
  printf("Histograms: %3.1f < pT  < %3.1f,  Nbin = %d\n", fHistoPtMin,  fHistoPtMax,  fHistoPtBins);
  printf("Histograms: %3.1f < phi < %3.1f, Nbin = %d\n", fHistoPhiMin, fHistoPhiMax, fHistoPhiBins);
  printf("Histograms: %3.1f < eta < %3.1f, Nbin = %d\n", fHistoEtaMin, fHistoEtaMax, fHistoEtaBins);
  printf("Histograms: %3.1f < mass < %3.1f, Nbin = %d\n", fHistoMassMin, fHistoMassMax, fHistoMassBins);
  printf("Histograms: %3.1f < asymmetry < %3.1f, Nbin = %d\n", fHistoAsymMin, fHistoAsymMax, fHistoAsymBins);
  printf("Histograms: %d < V0 Signal < %d, Nbin = %d\n", fHistoV0SMin, fHistoV0SMax, fHistoV0SBins);
  printf("Histograms: %d < V0 Mult < %d, Nbin = %d\n", fHistoV0MMin, fHistoV0MMax, fHistoV0MBins);
  printf("Histograms: %d < Track Mult < %d, Nbin = %d\n", fHistoTrMMin, fHistoTrMMax, fHistoTrMBins);
  printf("Histograms: %3.1f < p/E  < %3.1f, Nbin = %d\n", fHistoPOverEMin, fHistoPOverEMax, fHistoPOverEBins);
  printf("Histograms: %3.1f < dEdx < %3.1f, Nbin = %d\n", fHistodEdxMin,   fHistodEdxMax,   fHistodEdxBins);
  printf("Histograms: %3.1f < dR (track cluster)   < %3.1f, Nbin = %d\n", fHistodRMin,     fHistodRMax,     fHistodRBins);
  printf("Histograms: %3.1f < R=sqrt{x^2+y^2}    < %3.1f, Nbin = %d\n", fHistoRMin,      fHistoRMax,      fHistoRBins);
  printf("Histograms: %3.1f < X    < %3.1f, Nbin = %d\n", fHistoXMin,      fHistoXMax,      fHistoXBins);
  printf("Histograms: %3.1f < Y    < %3.1f, Nbin = %d\n", fHistoYMin,      fHistoYMax,      fHistoYBins);
  printf("Histograms: %3.1f < Z    < %3.1f, Nbin = %d\n", fHistoZMin,      fHistoZMax,      fHistoZBins);
  printf("Histograms: %g < Time < %g, Nbin = %d\n"      , fHistoTimeMin,   fHistoTimeMax,   fHistoTimeBins);
  printf("Histograms: %d < N cells per cluster    < %d, Nbin = %d\n", fHistoNClusCellMin,   fHistoNClusCellMax,   fHistoNClusCellBins);
  printf("Histograms: %d < N cells   < %d, Nbin = %d\n"             , fHistoNCellsMin,      fHistoNCellsMax,      fHistoNCellsBins);
  printf("Histograms: %d < N clusters   < %d, Nbin = %d\n"          , fHistoNClustersMin,   fHistoNClustersMax,   fHistoNClustersBins);
  printf("Histograms: %3.1f < Ratio< %3.1f, Nbin = %d\n", fHistoRatioMin,  fHistoRatioMax,  fHistoRatioBins);
  printf("Histograms: %3.1f < Vertex Distance < %3.1f, Nbin = %d\n", fHistoVertexDistMin, fHistoVertexDistMax, fHistoVertexDistBins);
    
  printf("Name of reference array      : %s\n", fAODObjArrayName.Data());	
  printf("String added histograms name : %s\n",fAddToHistogramsName.Data());
	
  printf("    \n") ;
  
} 



