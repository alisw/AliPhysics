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
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODPWG4Particle.h"

/// \cond CLASSIMP
ClassImp(AliAnaCaloTrackCorrBaseClass) ;
/// \endcond

//__________________________________________________________
/// Default constructor.
/// Initialize parameters.
//__________________________________________________________
AliAnaCaloTrackCorrBaseClass::AliAnaCaloTrackCorrBaseClass() : 
TObject(), 
fNModules(20),                fNRCU(2),        
fFirstModule(0),              fLastModule(19),
fNMaxCols(48),                fNMaxRows(24),  
fNMaxColsFull(48),            fNMaxRowsFull(24),  
fNMaxRowsFullMin(0),          fNMaxRowsFullMax(24),  
fDataMC(0),                   fDebug(0),
fCalorimeter(-1),             fCalorimeterString(""),
fCheckFidCut(0),              fCheckRealCaloAcc(0),
fCheckCaloPID(0),             fRecalculateCaloPID(0), 
fMinPt(0),                    fMaxPt(0),
fPairTimeCut(200),            fTRDSMCovered(-1),
fNZvertBin(0),                fNrpBin(0),
fNCentrBin(0),                fNmaxMixEv(0),
fDoOwnMix(0),                 fUseTrackMultBins(0),
fFillPileUpHistograms(0),     fFillHighMultHistograms(0),
fMakePlots(kFALSE),
fInputAODBranch(0x0),         fInputAODName(""),
fOutputAODBranch(0x0),        fNewAOD(kFALSE),
fOutputAODName(""),           fOutputAODClassName(""),
fAODObjArrayName(""),         fAddToHistogramsName(""),
fCaloPID(0x0),                fCaloUtils(0x0),
fFidCut(0x0),                 fHisto(0x0),
fIC(0x0),                     fMCUtils(0x0),                
fNMS(0x0),                    fReader(0x0),
fStudyClusterOverlapsPerGenerator(0),
fNCocktailGenNames(0)
{
  InitParameters();
}

//___________________________________________________________
/// Destructor.
//___________________________________________________________
AliAnaCaloTrackCorrBaseClass::~AliAnaCaloTrackCorrBaseClass() 
{  
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
/// Put cluster/track or created particle object
/// in the AODParticleCorrelation array.
//______________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::AddAODParticle(AliAODPWG4Particle pc)
{  
  if(!fOutputAODBranch)
  {
    AliFatal("No AOD branch available!!!\n");
    return; // coverity
  }
  
  Int_t i = fOutputAODBranch->GetEntriesFast();
  //new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
  if     (strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4Particle")==0)
  {
    new((*fOutputAODBranch)[i])  AliAODPWG4Particle(pc);
  }
  else if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliAODPWG4ParticleCorrelation")==0)
  {
    new((*fOutputAODBranch)[i])  AliAODPWG4ParticleCorrelation(pc);
  }
  else
  {
    AliFatal(Form("Cannot add an object of type < %s >, to the AOD TClonesArray \n", fOutputAODBranch->GetClass()->GetName()));
  }
}

//__________________________________________________________________________________________
/// Check vertex in mixed events. Input:
/// \param caloLabel: Index of cluster.
/// \param trackLabel: Index of track.
/// \return bool with quality of vertex, kFALSE not ok.
//__________________________________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::CheckMixedEventVertex(Int_t caloLabel, Int_t trackLabel)
{  
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
/// Recover ouput and input AOD pointers for each event in the maker.
/// put them in the corresponding pointer.
//________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() 
{	
  AliDebug(3,Form("AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() - Connect Input with name: <%s>; Connect output with name <%s>\n",fInputAODName.Data(),fOutputAODName.Data()));
  
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
  else if (aodHandler->GetExtensions())
  {
    AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject(GetReader()->GetDeltaAODFileName()); 
    if(ext)
    {
      AliAODEvent *aodEvent = ext->GetAOD(); 
      if(fNewAOD)fOutputAODBranch = (TClonesArray*) aodEvent->FindListObject(fOutputAODName);
      fInputAODBranch = (TClonesArray*) aodEvent->FindListObject(fInputAODName); 	  
      if(!fOutputAODBranch && fNewAOD) fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
      if(!fInputAODBranch)  fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);
    }
    else
    { // If no Delta AODs, kept in standard branch, to revise.
      if(fNewAOD && fReader->GetOutputEvent())
      {
        fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
        fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);	
      }
      else
      {
        fInputAODBranch  =  (TClonesArray *) fReader->GetInputEvent()->FindListObject(fInputAODName);	
        if(!fInputAODBranch && fReader->GetOutputEvent() ) 
          fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);//Try the output event.
      }
    }
  }
  else
  { // If no Delta AODs, kept in standard branch
    if(fNewAOD && fReader->GetOutputEvent())
    {
      fOutputAODBranch =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fOutputAODName);
      fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);	
    }
    else
    {
      fInputAODBranch  =  (TClonesArray *) fReader->GetInputEvent()->FindListObject(fInputAODName);
      if(!fInputAODBranch && fReader->GetOutputEvent())  
        fInputAODBranch  =  (TClonesArray *) fReader->GetOutputEvent()->FindListObject(fInputAODName);//Try the output event.
    }
  }
  
//  if(GetDebug() > 1)
//  {
//    if(fNewAOD && !fOutputAODBranch) 
//      AliInfo(Form("Output Branch <%s>, not found!\n",fOutputAODName.Data()));
//    if(!fNewAOD && !fInputAODBranch) 
//      AliInfo(Form("Input Branch  <%s>, not found!\n",fInputAODName.Data()));
//  }
}

//_____________________________________________________________________________________
/// Given the cluster ID stored in AliAODPWG4Particle, get the originator cluster 
/// and its index in the array. Input parameters:
/// \param clusters: Full TObjarray of clusters.
/// \param clId: Integer with the searched cluster ID.
/// \param iclus: Integer with index of cluster in cluster array.
/// \param first: Integer with the first index to consider in the array of clusters.
/// \return pointer to cluster.
//_____________________________________________________________________________________
AliVCluster * AliAnaCaloTrackCorrBaseClass::FindCluster(TObjArray* clusters, Int_t clId,
                                                        Int_t & iclus, Int_t first)
{  
  if ( !clusters ) return 0x0;
  
  for(iclus = first; iclus < clusters->GetEntriesFast(); iclus++)
  {
    AliVCluster *cluster= dynamic_cast<AliVCluster*> (clusters->At(iclus));
    if ( cluster )
    {
      if ( cluster->GetID() == clId )
      {
        return cluster;
      }
    }      
  }// calorimeter clusters loop
  
  return 0x0;
}

//______________________________________________________________________________________
/// Recover ouput and input AOD pointers for each event in AliCaloTrackMaker.
//______________________________________________________________________________________
TClonesArray * AliAnaCaloTrackCorrBaseClass::GetAODBranch(const TString & aodName) const
{
	AliDebug(3,Form("AliAnaCaloTrackCorrBaseClass::GetAODBranch() - Get Input Branch with name: <%s>; \n",aodName.Data()));
	
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
/// \return list of filtered Tracks from AliCaloTrackReader.
//_____________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetCTSTracks() const 
{  
  return fReader->GetCTSTracks(); 
}

//________________________________________________________________
/// \return list of PHOS filtered caloClusters from AliCaloTrackReader.
//________________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetPHOSClusters() const 
{  
  return fReader->GetPHOSClusters(); 
}

//_________________________________________________________________
/// \return list of EMCAL filtered caloClusters from AliCaloTrackReader.
//_________________________________________________________________
TObjArray *  AliAnaCaloTrackCorrBaseClass::GetEMCALClusters() const 
{  
  return fReader->GetEMCALClusters(); 
}

//______________________________________________________________________
/// \return list of all caloClusters in AOD output file.
//______________________________________________________________________
TClonesArray *  AliAnaCaloTrackCorrBaseClass::GetAODCaloClusters() const 
{  
  return fReader->GetOutputEvent()->GetCaloClusters(); 
}

//________________________________________________________________
/// \return list of all tracks in AOD output file. 
//________________________________________________________________
TClonesArray *  AliAnaCaloTrackCorrBaseClass::GetAODTracks() const 
{  
  return fReader->GetOutputEvent()->GetTracks(); 
}

//____________________________________________________________
/// Put data member values in string to keep in output container.
/// \return String with list of parameters.
//____________________________________________________________
TString  AliAnaCaloTrackCorrBaseClass::GetBaseParametersList()  
{
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
  snprintf(onePar,buffersize,"fCheckRealCaloAcc=%d (Check Real Calo Acceptance on/off) \n",fCheckRealCaloAcc) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckCaloPID =%d (Use Bayesian PID in calorimetes, on/off) \n",fCheckCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRecalculateCaloPID  =%d (Calculate PID from shower/tof/tracking parameters, on/off) \n",fRecalculateCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fInputAODName  =%s Input AOD name \n",fInputAODName.Data()) ;
  parList+=onePar ;	
  if(fNewAOD)
  {
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
/// Check the content of the cluster other than the generator that 
/// deposited more energy. Assign a tag.
///
/// \param cluster: pointer to AliVCluster
/// \param mctag: MC tag
/// \param genName: name of generator main contributor to the cluster
/// \param index: label
/// \param genNameBkg: name of generators overlapped to the cluster
/// \param indexBkg: bkg label
/// \return tag : 0-no other generator, 1 hijing, 2 other generator, 3 both hijing and other
//_____________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetCocktailGeneratorBackgroundTag(AliVCluster * cluster, Int_t mctag,
                                                                      TString & genName   , Int_t & index,
                                                                      TString & genNameBkg, Int_t & indexBkg)
{
  if(cluster->GetNLabels() == 0 || cluster->GetLabel() < 0 ) return -1;

  //(GetReader()->GetMC())->
  index = GetReader()->GetCocktailGeneratorAndIndex(cluster->GetLabel(), genName);
  //TString genNameOrg;
  //(GetReader()->GetMC())->GetCocktailGenerator(cluster->GetLabel(), genNameOrg);
 
  //printf("Generator?: %s, index %d\n",genName.Data(), index);
  
  Bool_t overlapGener       = kFALSE;
  Bool_t overlapGenerHIJING = kFALSE;
  Bool_t overlapGenerOther  = kFALSE;
  
  genNameBkg = "";
  indexBkg = -1;
  
  const UInt_t nlabels = cluster->GetNLabels();
  //Int_t noverlapsGen = 0;
  //TString genName2Prev = genName;
  for(UInt_t ilabel = 1; ilabel < nlabels; ilabel++)
  {
    Int_t label2 = cluster->GetLabels()[ilabel];
    TString genName2;
    //(GetReader()->GetMC())->GetCocktailGenerator(label2,genName2);
    Int_t index2 = GetReader()->GetCocktailGeneratorAndIndex(label2, genName2);

    //if(genName2 != genName2Prev) noverlapsGen++;
    
    if(genName2 != genName) 
    {
      //genName2Prev = genName2;
      
      if(!genNameBkg.Contains(genName2) && indexBkg != index2)
      {
        genNameBkg = Form("%s_%s",genNameBkg.Data(), genName2.Data());
        indexBkg = index2;
      }
      
      overlapGener = kTRUE; 
      
      if( genName2.Contains("ijing") && !genName.Contains("ijing")) 
        overlapGenerHIJING = kTRUE;
      
      if(!genName2.Contains("ijing"))
        overlapGenerOther  = kTRUE;
    }
  }

  Int_t genBkgTag = -1;
  
  if      ( !overlapGener )                            genBkgTag = 0; // Pure
  else if (  overlapGenerHIJING && !overlapGenerOther) genBkgTag = 1; // Gen+Hij
  else if ( !overlapGenerHIJING &&  overlapGenerOther) genBkgTag = 2; // GenX+GenY
  else if (  overlapGenerHIJING &&  overlapGenerOther) genBkgTag = 3; // GenX+GenY+Hij
  else                                                 genBkgTag = 4; 
  
  // check overlap with same generator, but not hijing
  Int_t overpdg[nlabels];
  Int_t overlab[nlabels];
  Int_t noverlaps = GetMCAnalysisUtils()->GetNOverlaps(cluster->GetLabels(), nlabels,mctag,-1,GetMC(),overpdg,overlab);
  Bool_t sameGenOverlap   = kFALSE;
  Bool_t sameGenOverlapHI = kFALSE;
  for(Int_t iover = 0; iover < noverlaps; iover++)
  {
    TString genName2;
    //(GetReader()->GetMC())->GetCocktailGenerator(overlab[iover],genName2);
    Int_t index2 = GetReader()->GetCocktailGeneratorAndIndex(overlab[iover], genName2);

    if ( genName2==genName  && index==index2)
    {
      if ( !genName.Contains("ijing") ) sameGenOverlap   = kTRUE;
      else                              sameGenOverlapHI = kTRUE;  
    }
  }
  
  //printf("bkg tag %d, noverlaps %d; same gen overlap %d\n",genBkgTag,noverlaps,sameGenOverlap);
  if(sameGenOverlap)
  {
    if(genBkgTag == 0) genBkgTag = 2; // GenX+GenX
    if(genBkgTag == 1) genBkgTag = 3; // GenX+GenX+Hij
  }

  // Logic a bit different for hijing main particles
  if(genName.Contains("ijing"))
  {
    
    if(sameGenOverlapHI)
    {
      if(!overlapGener) genBkgTag = 1; // Hij+Hij
      else              genBkgTag = 3; // Hij+Gen+Hij   
    }
    else
    {
      if(!overlapGener) genBkgTag = 0; // Pure
      else              genBkgTag = 2; // Hij+Gen  
    }
  }
  
  return genBkgTag;
}

//_____________________________________________________________________
/// Create AOD branch filled in the analysis.
//_____________________________________________________________________
TClonesArray * AliAnaCaloTrackCorrBaseClass::GetCreateOutputAODBranch() 
{  
  AliInfo(Form("Create AOD branch of %s objects and with name < %s >\n",
          fOutputAODClassName.Data(),fOutputAODName.Data())) ;
  
  TClonesArray * aodBranch = new TClonesArray(fOutputAODClassName, 0);
  
  aodBranch->SetName(fOutputAODName);
  
  return aodBranch ;
}

//________________________________________________________
/// \return Current event number.
//________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventNumber() const 
{  
  return fReader->GetEventNumber() ; 
}

//__________________________________________________________
/// \return  AliMCEvent pointer from AliCaloTrackReader.
//__________________________________________________________
AliMCEvent *  AliAnaCaloTrackCorrBaseClass::GetMC() const 
{  
  return fReader->GetMC(); 
}

//____________________________________________________________
/// \return Header pointer from AliCaloTrackReader.
//____________________________________________________________
AliHeader *  AliAnaCaloTrackCorrBaseClass::GetMCHeader() const
{  
  return fReader->GetHeader(); 
}

//____________________________________________________________________________
/// \return AliGenEventHeader pointer from AliCaloTrackReader.
//____________________________________________________________________________
AliGenEventHeader *  AliAnaCaloTrackCorrBaseClass::GetMCGenEventHeader() const
{  
  return fReader->GetGenEventHeader(); 
}

//_________________________________________________________________
/// Tag the event within a given track multiplicity bin.
/// Bin range defined in multiplicity limit array fTrackMultBins.
/// Multiplicity got from AliCaloTrackReader.
/// \return event track multiplicity bin.
//_________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetTrackMultiplicityBin() const
{  
  //curCentrBin = (GetTrackMultiplicity()-1)/5;
  //if(curCentrBin > GetNCentrBin()-1) curCentrBin=GetNCentrBin()-1;
  Int_t trackMult = GetReader()->GetTrackMultiplicity();
  
  for(Int_t ibin = 0; ibin < GetNTrackMultBin()-1; ibin++)
  {
    if(trackMult >= fTrackMultBins[ibin] && trackMult < fTrackMultBins[ibin+1]) return ibin;
  }
  
  AliWarning(Form("Bin not found for track multiplicity %d",trackMult));
  
  return -1;
}

//________________________________________________________________
/// Define the centrality bin for mixing.
/// In pp collisions analysis hardcoded track multiplicities.
/// \return centrality bin
//________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventCentralityBin() const
{
  Int_t curCentrBin = 0;
  
  if(fUseTrackMultBins) // pp collisions
  {
    return GetTrackMultiplicityBin();
  }
  else // Set centrality based on centrality task, PbPb collisions
  {
    Float_t minCent = GetReader()->GetCentralityBin(0);
    Float_t maxCent = GetReader()->GetCentralityBin(1);
    
    if((minCent< 0 && maxCent< 0) || minCent>=maxCent)
    {
      curCentrBin = GetEventCentrality() * GetNCentrBin() / GetReader()->GetCentralityOpt();
      if(curCentrBin==GetNCentrBin())
      {
        curCentrBin = GetNCentrBin()-1;
        AliDebug(1,Form("Centrality = %d, put it in last bin \n",GetEventCentrality()));
      }
    }
    else
    {
      curCentrBin = (Int_t)((GetEventCentrality()-minCent) * GetNCentrBin() / (maxCent-minCent));
      if(curCentrBin==GetNCentrBin()) curCentrBin = GetNCentrBin()-1;
    }
    
    AliDebug(1,Form("Current CentrBin %d, centrality %d, n bins %d, max bin from centrality %d",
                    curCentrBin, GetEventCentrality(), GetNCentrBin(), GetReader()->GetCentralityOpt()));
  }
  
  return curCentrBin;
}

//_______________________________________________________
/// \return Reaction plane bin.
//_______________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventRPBin() const
{
  Int_t curRPBin  = 0 ;
  
  if(GetNRPBin() > 1 && GetEventPlane())
  {
    Float_t epAngle = GetEventPlaneAngle();//->GetEventplane(GetEventPlaneMethod(),fReader->GetInputEvent());
    
    if(epAngle < 0 || epAngle >TMath::Pi())
    { 
      AliWarning(Form("Wrong event plane angle : %f \n",epAngle));
      return -1;
    }
    
    curRPBin = TMath::Nint(epAngle*(GetNRPBin()-1)/TMath::Pi());
    if(curRPBin >= GetNRPBin()) printf("RP Bin %d out of range %d",curRPBin,GetNRPBin());
    
    AliDebug(1,Form("Current RP bin %d, bin float %f, angle %f, n bins %d",
                    curRPBin,epAngle*(GetNRPBin()-1)/TMath::Pi(),epAngle,GetNRPBin()));
  }  
  
  return curRPBin ;
}

//_______________________________________________________
/// \return Vz bin, divide vertex in GetNZvertBin() bins, 
/// depending on the vertex cut
//_______________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventVzBin() const
{
  Double_t v[3] = {0,0,0}; //vertex 
  GetReader()->GetVertex(v);
  
  Int_t curZvertBin = (Int_t)(0.5*GetNZvertBin()*(v[2]+GetZvertexCut())/GetZvertexCut());
  
  AliDebug(1,Form("AliAnaCaloTrackCorrBaseClass::GetEventVzBin() - %d, vz %2.2f, n bins %d",
                  curZvertBin, v[2], GetNZvertBin()));
  
  return curZvertBin;
}

//________________________________________________________________________________________
/// \return  Event mixing bin, combination of vz, centrality and reaction plane bins.
//________________________________________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventMixBin(Int_t iCen, Int_t iVz, Int_t iRP) const
{  
  if(iCen<0 || iVz < 0 || iRP < 0)
    return -1;
  else
    return iCen*GetNZvertBin()*GetNRPBin()+iVz*GetNRPBin()+iRP;
}

//________________________________________________________
/// \return  Event mixing bin, combination of vz, centrality and reaction plane bins.
//________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventMixBin() const
{  
  //Get vertex z bin
  Int_t iVz =  GetEventVzBin();
  
  // centrality (PbPb) or tracks multiplicity (pp) bin
  Int_t iCen = GetEventCentralityBin();
  
  // reaction plane bin (PbPb)
  Int_t iRP = GetEventRPBin();  
  
  Int_t eventBin = GetEventMixBin(iCen, iVz, iRP);
  
  AliDebug(1,Form("Bins : cent %d, vz %d, RP %d, event %d/%d",
                  iCen,iVz, iRP, eventBin, GetNZvertBin()*GetNRPBin()*GetNCentrBin()));
  
  return eventBin;
}

//____________________________________________
/// Init once the debugging level, if requested.
/// Activate debug level in analysis.
//____________________________________________
void AliAnaCaloTrackCorrBaseClass::InitDebug()
{
  if( fDebug >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(this->ClassName(),fDebug);
  
  if( GetMCAnalysisUtils()->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetMCAnalysisUtils()->ClassName(),GetMCAnalysisUtils()->GetDebug());
  
  if( GetIsolationCut()->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetIsolationCut()->ClassName(),GetIsolationCut()->GetDebug());

  if( GetNeutralMesonSelection()->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetNeutralMesonSelection()->ClassName(),GetNeutralMesonSelection()->GetDebug());
    
  //printf("Debug levels: Ana %d, MC %d, Iso %d\n",fDebug,GetMCAnalysisUtils()->GetDebug(),GetIsolationCut()->GetDebug());
}

//_________________________________________________
/// Initialize the parameters of the analysis.
//_________________________________________________
void AliAnaCaloTrackCorrBaseClass::InitParameters()
{ 
  fDataMC              = kFALSE;
  fDebug               = 0;
  fCheckCaloPID        = kTRUE ;
  fCheckFidCut         = kFALSE ;
  fCheckRealCaloAcc    = kFALSE ;
  fRecalculateCaloPID  = kFALSE ;
  fMinPt               = 0.2  ; //Min pt in particle analysis
  fMaxPt               = 300. ; //Max pt in particle analysis
  fNZvertBin           = 1;
  fNrpBin              = 1;
  
  fCalorimeterString   = "EMCAL";
  fCalorimeter         = kEMCAL ;
  
  fTrackMultBins[0] =  0;  fTrackMultBins[1] =  5;  fTrackMultBins[2] = 10;
  fTrackMultBins[3] = 15;  fTrackMultBins[4] = 20;  fTrackMultBins[5] = 30;
  fTrackMultBins[6] = 40;  fTrackMultBins[7] = 55;  fTrackMultBins[8] = 70;
  for(Int_t ibin=9; ibin < 20; ibin++) fTrackMultBins[ibin] = 10000;
  
  //fReader    = new AliCaloTrackReader(); //Initialized in maker
  //fCaloUtils = new AliCalorimeterUtils();//Initialized in maker
  
  fNewAOD              = kFALSE ;
  fOutputAODName       = "CaloTrackCorr";
  fOutputAODClassName  = "AliAODPWG4Particle";
  fInputAODName        = "CaloTrackCorr";
  fAddToHistogramsName = "";
  fAODObjArrayName     = "Ref";
  
  fNCocktailGenNames = 7;
  // Order matters, here cocktail of MC LHC14a1a
  fCocktailGenNames[0] = ""; // First must be always empty
  fCocktailGenNames[1] = "pi0EMC";
  fCocktailGenNames[2] = "pi0";
  fCocktailGenNames[3] = "etaEMC";
  fCocktailGenNames[4] = "eta";
  fCocktailGenNames[5] = "hijing";
  fCocktailGenNames[6] = "other";  
  
  for(Int_t igen = 7; igen < 10; igen++)
    fCocktailGenNames[igen] = "";
  
  for(Int_t igen = 0; igen < 10; igen++)
    fCocktailGenIndeces[igen] = -1;
}

//_________________________________________________
/// Initialize the parameters related to the calorimeters
/// number of modules and modules to analyze.
/// Mainly used to set histogram ranges. Add the line at 
/// the Init() or GetCreateOutputObjects() methods of derived classes
//_________________________________________________
void AliAnaCaloTrackCorrBaseClass::InitCaloParameters()
{
  fNModules = GetCaloUtils()->GetNumberOfSuperModulesUsed();
  if(GetCalorimeter()==kPHOS && fNModules > 4) fNModules = 4;
  
  fFirstModule = 0; 
  fLastModule  = fNModules-1;
  
  // Set First/Last SM depending on CaloUtils or fiducial cut settings
   
  if ( IsFiducialCutOn() )
  {
    //printf("Get SM range from FiducialCut\n");

    if(GetCalorimeter() != kPHOS)
    {
      Int_t nSections = GetFiducialCut()->GetEMCALFidCutMaxPhiArray()->GetSize();
      if( nSections == 1 )
      {
        Float_t minPhi = GetFiducialCut()->GetEMCALFidCutMinPhiArray()->At(0);
        Float_t maxPhi = GetFiducialCut()->GetEMCALFidCutMaxPhiArray()->At(0);
        //printf("sections %d, min %f, max %f\n",nSections,minPhi,maxPhi);
        
        if     ( minPhi > 70  && maxPhi < 190) // EMCal
        {
          fFirstModule = 0;
          fLastModule  = 11;
        }
        else if( minPhi > 250 && maxPhi < 330) // DCal
        {
          fFirstModule = 12;
          fLastModule  = 19;
        }
      }
    }
  }

  // Overwrite what used in FidCut, if more strict on CaloUtils
  // Needed for special case in QA analysis train
  if ( GetCaloUtils()->GetFirstSuperModuleUsed() >= 0 )
  {
    if(fFirstModule < GetCaloUtils()->GetFirstSuperModuleUsed() || 
       fLastModule  > GetCaloUtils()->GetLastSuperModuleUsed())
    {
      //printf("Get SM range from CaloUtils\n");
      
      fFirstModule = GetCaloUtils()->GetFirstSuperModuleUsed();
      fLastModule  = GetCaloUtils()->GetLastSuperModuleUsed();
    }
  }
  
  // EMCAL
  fNMaxCols = 48;
  fNMaxRows = 24;
  fNRCU     = 2 ;
  // PHOS
  if(GetCalorimeter()==kPHOS)
  {
    fNMaxCols = 56;
    fNMaxRows = 64;
    fNRCU     = 4 ;
  }
  
  fNMaxColsFull = fNMaxCols;
  fNMaxRowsFull = fNMaxRows;
  
  fNMaxRowsFullMax = fNMaxRowsFull;
  fNMaxRowsFullMin = 0;
  
  if(GetCalorimeter()==kEMCAL)
  {
    fNMaxColsFull=2*fNMaxCols;
    
    fNMaxRowsFull=Int_t(fNModules/2)*fNMaxRows;
    if(fNMaxRowsFull > 208) 
      fNMaxRowsFull = 208; //  8*24+2/3.*24, reduce since 1/3 SM should not be counted full.
    
    fNMaxRowsFullMin = 0;
    fNMaxRowsFullMax = fNMaxRowsFull; 

    if(fLastModule < 12)
    {
      fNMaxRowsFullMax = Int_t((fLastModule-fFirstModule+1)/2)*fNMaxRows;
      if(fNMaxRowsFullMax > 127) fNMaxRowsFullMax = 127; // 24*5+8
      fNMaxRowsFullMin = 0; 
    }
    else if (fFirstModule > 11)
    {
      fNMaxRowsFullMax = fNMaxRowsFull; 
      fNMaxRowsFullMin = Int_t(fFirstModule/2)*fNMaxRows-Int_t(2./3.*fNMaxRows); // remove 2/3*24
    }    
  }
  else
  {
    fNMaxRowsFull=fNModules*fNMaxRows;
  }
  
//  printf("%s: N SM %d, first SM %d, last SM %d, SM col-row (%d,%d), Full detector col-row (%d, %d), partial calo row min-max(%d,%d) \n",
//                  GetName(),fNModules,fFirstModule,fLastModule, fNMaxCols,fNMaxRows, 
//                  fNMaxColsFull,fNMaxRowsFull, fNMaxRowsFullMin,fNMaxRowsFullMax);

  AliDebug(1,Form("N SM %d, first SM %d, last SM %d, SM col-row (%d,%d), Full detector col-row (%d, %d), partial calo row min-max(%d,%d)",
                  fNModules,fFirstModule,fLastModule, fNMaxCols,fNMaxRows, 
                  fNMaxColsFull,fNMaxRowsFull, fNMaxRowsFullMin,fNMaxRowsFullMax));

  
}

//__________________________________________________________________
/// Print some relevant parameters set for the analysis.
//__________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::Print(const Option_t * opt) const
{
  if(! opt)
    return;
	
  printf("New AOD:            =     %d\n",      fNewAOD);
  printf("Input AOD name:     =     %s\n",      fInputAODName.Data());
  printf("Output AOD name:    =     %s\n",      fOutputAODName.Data());
  printf("Output AOD Class name: =  %s\n",      fOutputAODClassName.Data());
  printf("Name of reference array      : %s\n", fAODObjArrayName.Data());
  printf("String added histograms name : %s\n", fAddToHistogramsName.Data());

  printf("Min Photon pT       =     %2.2f\n", fMinPt) ;
  printf("Max Photon pT       =     %3.2f\n", fMaxPt) ;
  printf("Check PID           =     %d\n",    fCheckCaloPID) ;
  printf("Recalculate PID     =     %d\n",    fRecalculateCaloPID) ;
  printf("Check Fiducial cut  =     %d\n",    fCheckFidCut) ;
  printf("Check Real Calo Acc =     %d\n",    fCheckRealCaloAcc) ;
  printf("Check MC labels     =     %d\n",    fDataMC);
  printf("Make plots?         =     %d\n",    fMakePlots);
  printf("Debug Level         =     %d\n",    fDebug);
  
  printf("    \n") ;
} 

//_______________________________________________________________
/// Set the calorimeter for the analysis. A string.
//_______________________________________________________________
void AliAnaCaloTrackCorrBaseClass::SetCalorimeter(TString & calo)
{
  fCalorimeterString = calo;
  
  if     (calo=="EMCAL") fCalorimeter = kEMCAL;
  else if(calo=="PHOS" ) fCalorimeter = kPHOS;
  else if(calo=="CTS")   fCalorimeter = kCTS;
  else if(calo=="DCAL")  fCalorimeter = kDCAL;
  else if(calo.Contains("DCAL") && calo.Contains("PHOS")) fCalorimeter = kDCALPHOS;
  else AliFatal(Form("Detector < %s > not known!", calo.Data()));

}

//_______________________________________________________________
/// Set the calorimeter for the analysis. An integer.
//_______________________________________________________________
void AliAnaCaloTrackCorrBaseClass::SetCalorimeter(Int_t calo)
{  
  fCalorimeter = calo;
  
  if     (calo==kEMCAL)    fCalorimeterString = "EMCAL";
  else if(calo==kPHOS )    fCalorimeterString = "PHOS";
  else if(calo==kCTS)      fCalorimeterString = "CTS";
  else if(calo==kDCAL)     fCalorimeterString = "DCAL";
  else if(calo==kDCALPHOS) fCalorimeterString = "DCAL_PHOS";
  else AliFatal(Form("Detector < %d > not known!", calo));
}


