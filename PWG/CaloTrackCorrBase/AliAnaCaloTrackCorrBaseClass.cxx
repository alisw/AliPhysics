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
#include <TCustomBinning.h>

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
#include "AliCaloTrackParticle.h"

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
fNSectors(10),       
fFirstSector(0),              fLastSector(9),
fNMaxCols(48),                fNMaxRows(24),  
fNMaxColsFull(48),            fNMaxRowsFull(24),  
fNMaxRowsFullMin(0),          fNMaxRowsFullMax(24),  
fTotalUsedSM(20),             fHistoSMArr(0),
fHistoNColumns(50),           fHistoColumnArr(0),
fHistoColumnMin(-1.5),        fHistoColumnMax(48.5),
fHistoNRows(402),             fHistoRowArr(0),
fHistoRowMin(-1.5),           fHistoRowMax(120.5),
fHistoPtBinNonConstantInArray(0),
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
fFillGenPartHisto(1),         fMakePlots(kFALSE),
fFillEmbedHistograms(0),      fSelectEmbededSignal(0),

fInputAODBranch(0x0),         fInputAODName(""),
fOutputAODBranch(0x0),        fNewAOD(kFALSE),
fOutputAODName(""),           fOutputAODClassName(""),
fAODObjArrayName(""),         fAddToHistogramsName(""),
fCaloPID(0x0),                fCaloUtils(0x0),
fFidCut(0x0),                 fHisto(0x0),
fIC(0x0),                                    
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
  delete fNMS     ;     
  delete fHisto   ;    
}

//______________________________________________________________________
/// Put cluster/track or created particle object
/// in the AODParticleCorrelation array.
//______________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::AddAODParticle(AliCaloTrackParticle pc)
{  
  if(!fOutputAODBranch)
  {
    AliFatal("No AOD branch available!!!");
    return; // coverity
  }
  
  Int_t i = fOutputAODBranch->GetEntriesFast();
  //new((*fOutputAODBranch)[i])  AliCaloTrackParticle(pc);
  if     (strcmp(fOutputAODBranch->GetClass()->GetName(),"AliCaloTrackParticle")==0)
  {
    new((*fOutputAODBranch)[i])  AliCaloTrackParticle(pc);
  }
  else if(strcmp(fOutputAODBranch->GetClass()->GetName(),"AliCaloTrackParticleCorrelation")==0)
  {
    new((*fOutputAODBranch)[i])  AliCaloTrackParticleCorrelation(pc);
  }
  else
  {
    AliFatal(Form("Cannot add an object of type < %s >, to the AOD TClonesArray", fOutputAODBranch->GetClass()->GetName()));
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
  AliDebug(3,Form("AliAnaCaloTrackCorrBaseClass::ConnectInputOutputAODBranches() - Connect Input with name: <%s>; Connect output with name <%s>",fInputAODName.Data(),fOutputAODName.Data()));
  
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
//      AliInfo(Form("Output Branch <%s>, not found!",fOutputAODName.Data()));
//    if(!fNewAOD && !fInputAODBranch) 
//      AliInfo(Form("Input Branch  <%s>, not found!",fInputAODName.Data()));
//  }
}

//_____________________________________________________________________________________
/// Given the cluster ID stored in AliCaloTrackParticle, get the originator cluster 
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
TClonesArray * AliAnaCaloTrackCorrBaseClass::GetAODBranch(TString & aodName) const
{
	AliDebug(3,Form("AliAnaCaloTrackCorrBaseClass::GetAODBranch() - Get Input Branch with name: <%s>;",aodName.Data()));
	
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
  snprintf(onePar,buffersize,"--- AliAnaCaloTrackCorrBaseClass ---") ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"Minimal P_t: %2.2f ;", fMinPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"Minimal P_t: %2.2f ;", fMaxPt) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"|t_{1}-t_{2}| < %2.2f;", fPairTimeCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fDataMC =%d (Check MC information, on/off);",fDataMC) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckFidCut=%d (Check Fiducial cut selection on/off);",fCheckFidCut) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckRealCaloAcc=%d (Check Real Calo Acceptance on/off);",fCheckRealCaloAcc) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fCheckCaloPID =%d (Use Bayesian PID in calorimetes, on/off);",fCheckCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fRecalculateCaloPID  =%d (Calculate PID from shower/tof/tracking parameters, on/off);",fRecalculateCaloPID) ;
  parList+=onePar ;
  snprintf(onePar,buffersize,"fInputAODName  =%s Input AOD name;",fInputAODName.Data()) ;
  parList+=onePar ;	
  if(fNewAOD)
  {
    snprintf(onePar,buffersize,"fOutputAODName  =%s Output AOD name",fOutputAODName.Data()) ;
    parList+=onePar ;	
    snprintf(onePar,buffersize,"fOutputAODClassName  =%s Output AOD class name;",fOutputAODClassName.Data()) ;
    parList+=onePar ;	
  }
  snprintf(onePar,buffersize,"fAODObjArrayName  =%s Reference arrays in AOD name;",fAODObjArrayName.Data()) ;
  parList+=onePar ;	
  snprintf(onePar,buffersize,"fAddToHistogramsName  =%s String added to beginning of histograms name;",fAddToHistogramsName.Data()) ;
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
  AliInfo(Form("Create AOD branch of %s objects and with name < %s >",
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
Int_t AliAnaCaloTrackCorrBaseClass::GetEventCentralityBin() 
{
  Int_t curCentrBin = 0;
  
  if ( fNCentrBin <= 1)
  {
    return 0;
  }
  else if ( fUseTrackMultBins ) // pp collisions
  {
    return GetTrackMultiplicityBin();
  }
  else // Set centrality based on centrality task, PbPb collisions
  {
    TArrayD cenArr = GetHistogramRanges()->GetHistoCentralityArr();

    if ( cenArr.GetSize() == 0 )
       SetEventCentralityBins();
    
    curCentrBin = TMath::BinarySearch(cenArr.GetSize(),cenArr.GetArray(), (Double_t) GetEventCentrality());
    
    AliDebug(1,Form("Current CentrBin %d, centrality %d, n bins %d, max bin from centrality %d",
                    curCentrBin, GetEventCentrality(), GetNCentrBin(), GetReader()->GetCentralityOpt()));
  }
  
  return curCentrBin;
}

//________________________________________________________________
/// Initialize centrality bins used for histogramming
//________________________________________________________________
void AliAnaCaloTrackCorrBaseClass::SetEventCentralityBins()
{
  if ( GetHistogramRanges()->GetHistoCentralityArr().GetSize() != 0 ) return;
  
  Float_t max  = GetHistogramRanges()->GetHistoCentralityMax();
  Float_t min  = GetHistogramRanges()->GetHistoCentralityMin();
  Int_t   nbin = GetHistogramRanges()->GetHistoCentralityBins();
  
  Float_t minCent = GetReader()->GetCentralityBin(0);
  Float_t maxCent = GetReader()->GetCentralityBin(1);
  
  if ( minCent >=0 && TMath::Floor(minCent) != TMath::Floor(min) ) min = minCent;
  if ( maxCent >=0 && TMath::Floor(maxCent) != TMath::Floor(max) ) max = maxCent;
  
  TCustomBinning cenBinning;
  cenBinning.SetMinimum(min);
  
  Float_t binWidth = 1;
  if      ( fNCentrBin > 0 ) // coarser binning
    binWidth = (max-min) / fNCentrBin ; 
  else if ( nbin       > 0 ) 
    binWidth = (max-min) / nbin;
  
  cenBinning.AddStep(max, binWidth); 
  
  TArrayD cenBinsArray;
  cenBinning.CreateBinEdges(cenBinsArray);
  GetHistogramRanges()->SetHistoCentralityArr(cenBinsArray);
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
      AliWarning(Form("Wrong event plane angle : %f",epAngle));
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
Int_t AliAnaCaloTrackCorrBaseClass::GetEventMixBin(Int_t iCen, Int_t iVz, Int_t iRP) 
{  
  if(iCen<0 || iVz < 0 || iRP < 0)
    return -1;
  else
    return iCen*GetNZvertBin()*GetNRPBin()+iVz*GetNRPBin()+iRP;
}

//________________________________________________________
/// \return  Event mixing bin, combination of vz, centrality and reaction plane bins.
//________________________________________________________
Int_t AliAnaCaloTrackCorrBaseClass::GetEventMixBin() 
{  
  //Get vertex z bin
  Int_t iVz =  GetEventVzBin();
  
  // centrality (PbPb) or tracks multiplicity (pp) bin
  Int_t iCen = GetEventCentralityBin();
  
  // reaction plane bin (PbPb)
  Int_t iRP = GetEventRPBin();  
  
  Int_t eventBin = GetEventMixBin(iCen, iVz, iRP);
  
  AliDebug(1,Form("Bins : cent %d (<%d), vz %d (<%d), RP %d (<%d), event %d/%d",
                  iCen, fNCentrBin,
                  iVz, GetNZvertBin(),
                  iRP, GetNRPBin(),
                  eventBin, GetNZvertBin()*GetNRPBin()*GetNCentrBin()));
  
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
  
  if( GetIsolationCut()->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetIsolationCut()->ClassName(),GetIsolationCut()->GetDebug());

  if( GetNeutralMesonSelection()->GetDebug() >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(GetNeutralMesonSelection()->ClassName(),GetNeutralMesonSelection()->GetDebug());
    
  //printf("Debug levels: Ana %d, Neutral Sel %d, Iso %d\n",fDebug,GetNeutralMesonSelection()->GetDebug(),GetIsolationCut()->GetDebug());
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
  
  fHistoPtBinNonConstantInArray = kTRUE;
  
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
  fOutputAODClassName  = "AliCaloTrackParticle";
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
  
  fFirstSector = 0;
  fLastSector  = fNModules/2 - 1;
  
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
          fFirstModule =  0;
          fLastModule  = 11;
          fFirstSector =  0;
          fLastSector  =  5;
        }
        else if( minPhi > 250 && maxPhi < 330) // DCal
        {
          fFirstModule = 12;
          fLastModule  = 19;
          fFirstSector =  6;
          fLastSector  =  9;
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
      
      fFirstSector =  fFirstModule / 2;
      fLastSector  =  fLastModule  / 2;
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
  
  // For histograms:
  fTotalUsedSM    = fLastModule-fFirstModule+1;
  
  TCustomBinning smBinning;
  smBinning.SetMinimum(fFirstModule-0.5);
  smBinning.AddStep(fLastModule+0.5, 1); 
  smBinning.CreateBinEdges(fHistoSMArr);
  
  fHistoNColumns  = fNMaxColsFull+2;
  fHistoColumnMin = -1.5;
  fHistoColumnMax = fNMaxColsFull+0.5;
  
  fHistoNRows     = fNMaxRowsFullMax-fNMaxRowsFullMin+2;
  fHistoRowMin    = fNMaxRowsFullMin-1.5;
  fHistoRowMax    = fNMaxRowsFullMax+0.5;
  
  // Cell column-row histograms, see base class for data members setting
  TCustomBinning rowBinning;
  rowBinning.SetMinimum(fHistoRowMin-1.5);
  rowBinning.AddStep(fHistoRowMax+0.5,1);   
  rowBinning.CreateBinEdges(fHistoRowArr);
  //
  TCustomBinning colBinning;
  colBinning.SetMinimum(fHistoColumnMin-1.5);
  colBinning.AddStep(fHistoColumnMax+0.5,1);   
  colBinning.CreateBinEdges(fHistoColumnArr);  
  
//  printf("%s: N SM %d, first SM %d, last SM %d, SM col-row (%d,%d), Full detector col-row (%d, %d), partial calo row min-max(%d,%d) \n",
//                  GetName(),fNModules,fFirstModule,fLastModule, fNMaxCols,fNMaxRows, 
//                  fNMaxColsFull,fNMaxRowsFull, fNMaxRowsFullMin,fNMaxRowsFullMax);

  AliDebug(1,Form("N SM %d, first SM %d, last SM %d, N sector %d, first sector %d, last sector %d,"
                  " SM col-row (%d,%d), Full detector col-row (%d, %d), partial calo row min-max(%d,%d)",
                  fNModules,fFirstModule,fLastModule, 
                  fNSectors,fFirstSector,fLastSector, 
                  fNMaxCols,fNMaxRows, 
                  fNMaxColsFull,fNMaxRowsFull, fNMaxRowsFullMin,fNMaxRowsFullMax));
}

///
/// In case of using the array instead of constant bins and limits and array not passed 
/// set here the array ranges, depending on analysis settings for some.
/// Call it in the GetCreateOutputObjects() for each analysis class.
/// Mainly used on TH3 histograms
///
void AliAnaCaloTrackCorrBaseClass::InitHistoRangeArrays()
{
  // Momentum and energy
  if ( GetHistogramRanges()->GetHistoPtArr().GetSize() == 0 )
  {
    Float_t ptmaxhi = GetHistogramRanges()->GetHistoPtMax();
    Float_t ptminhi = GetHistogramRanges()->GetHistoPtMin();
    
    TCustomBinning ptBinning;
    
    if ( fHistoPtBinNonConstantInArray )
    {
      Float_t min = TMath::Floor(GetMinPt()); // if min 0, 0.5 or 0.7 then start at 0, if  5 or 5.5 then 5
      if ( ptminhi > min ) min = ptminhi;
      ptBinning.SetMinimum(min); 
      
      Float_t max = GetMaxPt();
      if ( ptmaxhi < GetMaxPt() ) max = ptmaxhi;
      
      if ( min <   5 && max >   5 ) ptBinning.AddStep(5  , 0.5);
      if ( min <  12 && max >  12 ) ptBinning.AddStep(12 , 1.0);
      if ( min <  20 && max >  20 ) ptBinning.AddStep(20 , 2.0); 
      if ( min <  50 && max >  50 ) ptBinning.AddStep(50 , 5.0); 
      if ( min < 100 && max > 100 ) ptBinning.AddStep(100, 10.); 
      if ( min < 200 && max > 200 ) ptBinning.AddStep(200, 20.); 
      
      if      ( max <=   5 ) ptBinning.AddStep(max, 0.5);
      if      ( max <=  12 ) ptBinning.AddStep(max, 1.0);
      else if ( max <=  20 ) ptBinning.AddStep(max, 2.0); 
      else if ( max <=  50 ) ptBinning.AddStep(max, 5.0); 
      else if ( max <= 100 ) ptBinning.AddStep(max, 10.); 
      else if ( max <= 200 ) ptBinning.AddStep(max, 20.); 
      else                   ptBinning.AddStep(max, 50.); 
    }
    else 
    {
      ptBinning.SetMinimum(ptminhi);
      ptBinning.AddStep(ptmaxhi, (ptmaxhi-ptminhi) / GetHistogramRanges()->GetHistoPtBins()); 
    }
    
    TArrayD ptBinsArray;
    ptBinning.CreateBinEdges(ptBinsArray);
    GetHistogramRanges()->SetHistoPtArr(ptBinsArray);
  }

  // Cell, in cone energy/pt histograms, finer binning at low energy (same as GetHistoFinePt)
  if ( GetHistogramRanges()->GetHistoCellEnArr().GetSize() == 0 )
  {
    Float_t eCellMaxhi = GetHistogramRanges()->GetHistoCellEnMax();
    Float_t eCellMinhi = GetHistogramRanges()->GetHistoCellEnMin();
    TCustomBinning eCellBinning;
    
    if ( eCellMinhi <  0 ) eCellBinning.SetMinimum(0.1);
    else                   eCellBinning.SetMinimum(eCellMinhi);
    
    if ( eCellMaxhi >   1 ) eCellBinning.AddStep(  1, 0.05);
    if ( eCellMaxhi >   3 ) eCellBinning.AddStep(  3, 0.10); 
    if ( eCellMaxhi >   5 ) eCellBinning.AddStep(  5, 0.50); 
    if ( eCellMaxhi >  10 ) eCellBinning.AddStep( 10, 1.00); 
    if ( eCellMaxhi >  50 ) eCellBinning.AddStep( 50, 5.00); 
    if ( eCellMaxhi > 100 ) eCellBinning.AddStep(100,10.00); 
    if ( eCellMaxhi > 200 ) eCellBinning.AddStep(200,20.00); 
    
    if      ( eCellMaxhi <=   1 ) eCellBinning.AddStep(eCellMaxhi, 0.05);
    else if ( eCellMaxhi <=   3 ) eCellBinning.AddStep(eCellMaxhi, 0.10); 
    else if ( eCellMaxhi <=   5 ) eCellBinning.AddStep(eCellMaxhi, 0.50); 
    else if ( eCellMaxhi <=  10 ) eCellBinning.AddStep(eCellMaxhi, 1.00); 
    else if ( eCellMaxhi <=  50 ) eCellBinning.AddStep(eCellMaxhi, 5.00); 
    else if ( eCellMaxhi <= 100 ) eCellBinning.AddStep(eCellMaxhi,10.00); 
    else if ( eCellMaxhi <= 200 ) eCellBinning.AddStep(eCellMaxhi,20.00); 
    else                          eCellBinning.AddStep(eCellMaxhi,50.00); 
    
    TArrayD eCellBinsArray;
    eCellBinning.CreateBinEdges(eCellBinsArray);
    GetHistogramRanges()->SetHistoCellEnArr(eCellBinsArray);
  }
  
  // Cluster pT/energy histograms, coarser binning
  if ( GetHistogramRanges()->GetHistoWidePtArr().GetSize() == 0 )
  {
    Float_t ptWideMaxhi = GetHistogramRanges()->GetHistoWidePtMax();
    Float_t ptWideMinhi = GetHistogramRanges()->GetHistoWidePtMin();
    
    TCustomBinning ptWideBinning;
    
    Float_t min = TMath::Floor(GetMinPt());// if min 0, 0.5 or 0.7 then start at 0, if  5 or 5.5 then 5
    if ( ptWideMinhi > min ) min = ptWideMinhi;
    ptWideBinning.SetMinimum(min); 
          
    Float_t max     = GetMaxPt();
    if ( ptWideMaxhi < GetMaxPt() ) max = ptWideMaxhi;

    if ( max >   5 ) ptWideBinning.AddStep(  5,  1.);
    if ( max >  25 ) ptWideBinning.AddStep( 25,  5.);
    if ( max > 100 ) ptWideBinning.AddStep(100, 25.); 
    if ( max > 200 ) ptWideBinning.AddStep(200, 50.); 
    
    if      ( max <=   5 ) ptWideBinning.AddStep(  5,  1.);
    else if ( max <=  25 ) ptWideBinning.AddStep( 25,  5.);
    else if ( max <= 100 ) ptWideBinning.AddStep(100, 25.); 
    else if ( max <= 200 ) ptWideBinning.AddStep(200, 50.); 
    else                   ptWideBinning.AddStep(max,100.); 
    
    TArrayD ptWideBinsArray;
    ptWideBinning.CreateBinEdges(ptWideBinsArray);
    GetHistogramRanges()->SetHistoWidePtArr(ptWideBinsArray);
  }
  
  // Cluster/Track pt in isolation cone or UE bands
  if ( GetHistogramRanges()->GetHistoPtInConeArr().GetSize() == 0 )
  {
    Float_t ptmax = GetHistogramRanges()->GetHistoPtInConeMax();
    Float_t ptmin = GetHistogramRanges()->GetHistoPtInConeMin();
    
    TCustomBinning ptCBinning;    
    ptCBinning.SetMinimum(ptmin);   
    if ( ptmax > 2  ) ptCBinning.AddStep(  2, 0.2);                            
    if ( ptmax > 15 ) ptCBinning.AddStep( 15, 0.5);                            
    if ( ptmax > 30 ) ptCBinning.AddStep( 30, 1.0); 
    if ( ptmax > 60 ) ptCBinning.AddStep( 60, 2.5); 
    if ( ptmax > 100) ptCBinning.AddStep(100, 5.0);   
    if ( ptmax > 200) ptCBinning.AddStep(200,10.0);  

    if      ( ptmax <= 2  ) ptCBinning.AddStep(ptmax, 0.2);                            
    else if ( ptmax <= 15 ) ptCBinning.AddStep(ptmax, 0.5);                            
    else if ( ptmax <= 30 ) ptCBinning.AddStep(ptmax, 1.0); 
    else if ( ptmax <= 60 ) ptCBinning.AddStep(ptmax, 2.5);  
    else if ( ptmax <= 100) ptCBinning.AddStep(ptmax, 5.0);
    else if ( ptmax <= 200) ptCBinning.AddStep(ptmax,10.0); 
    else                    ptCBinning.AddStep(ptmax,20.0); 
    
    TArrayD ptCBinsArray;  
    ptCBinning.CreateBinEdges(ptCBinsArray);
    GetHistogramRanges()->SetHistoPtInConeArr(ptCBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoShowerShapeArr().GetSize() == 0 )
  {
    Float_t max = GetHistogramRanges()->GetHistoShowerShapeMax();
    Float_t min = GetHistogramRanges()->GetHistoShowerShapeMin();

    TCustomBinning ssBinning;
    if (min <= 0 ) ssBinning.SetMinimum(-0.01);
    else           ssBinning.SetMinimum(min);
    
    if ( fHistoPtBinNonConstantInArray )
    {
      if ( min < 0.5 && max > 0.5 ) ssBinning.AddStep(0.50,0.01);  // 51
      if ( min < 1.0 && max > 1.0 ) ssBinning.AddStep(1.00,0.05);  // 10
      if ( min < 3.0 && max > 3.0 ) ssBinning.AddStep(3.00,0.10);  // 20
      if ( min < 5.0 && max > 5.0 ) ssBinning.AddStep(5.00,0.25);  // 20

      if      ( max <= 0.5 ) ssBinning.AddStep(max,0.01);  // 51
      else if ( max <= 1.0 ) ssBinning.AddStep(max,0.05);  // 10
      else if ( max <= 3.0 ) ssBinning.AddStep(max,0.10);  // 20
      else if ( max <= 5.0 ) ssBinning.AddStep(max,0.25);  // 20
      else                   ssBinning.AddStep(max,0.50);
    }
    else
    {
      ssBinning.AddStep(max, (max-min) / GetHistogramRanges()->GetHistoShowerShapeBins());
    }
    
    TArrayD ssBinsArray;
    ssBinning.CreateBinEdges(ssBinsArray);
    GetHistogramRanges()->SetHistoShowerShapeArr(ssBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoEtaArr().GetSize() == 0 )
  {
    TCustomBinning etaBinning;
    etaBinning.SetMinimum(GetHistogramRanges()->GetHistoEtaMin());
    Float_t binWidth = ( GetHistogramRanges()->GetHistoEtaMax() - GetHistogramRanges()->GetHistoEtaMin() ) / GetHistogramRanges()->GetHistoEtaBins();
    etaBinning.AddStep(GetHistogramRanges()->GetHistoEtaMax(), binWidth); 
    
    TArrayD etaBinsArray;
    etaBinning.CreateBinEdges(etaBinsArray);
    GetHistogramRanges()->SetHistoEtaArr(etaBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoPhiArr().GetSize() == 0 )
  {
    TCustomBinning phiBinning;
    phiBinning.SetMinimum(GetHistogramRanges()->GetHistoPhiMin());
    Float_t binWidth = ( GetHistogramRanges()->GetHistoPhiMax() - GetHistogramRanges()->GetHistoPhiMin() ) / GetHistogramRanges()->GetHistoPhiBins();
    phiBinning.AddStep(GetHistogramRanges()->GetHistoPhiMax(), binWidth); 
    
    TArrayD phiBinsArray;
    phiBinning.CreateBinEdges(phiBinsArray);
    GetHistogramRanges()->SetHistoPhiArr(phiBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoTrackResidualEtaArr().GetSize() == 0 )
  {
    TCustomBinning resBinning; // track matching residuals
    resBinning.SetMinimum(-0.100);
    resBinning.AddStep(-0.050,0.010);
    resBinning.AddStep(-0.025,0.005);
    resBinning.AddStep( 0.025,0.001);
    resBinning.AddStep( 0.050,0.005);
    resBinning.AddStep( 0.100,0.010);
    
    TArrayD resBinsArray;
    resBinning.CreateBinEdges(resBinsArray);
    GetHistogramRanges()->SetHistoTrackResidualEtaArr(resBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoTrackResidualPhiArr().GetSize() == 0 )
  {
    TCustomBinning resBinning; // track matching residuals
    resBinning.SetMinimum(-0.100);
    resBinning.AddStep(-0.050,0.010);                            
    resBinning.AddStep(-0.025,0.005);                            
    resBinning.AddStep( 0.025,0.001);                            
    resBinning.AddStep( 0.050,0.005);                            
    resBinning.AddStep( 0.100,0.010);                            
    
    TArrayD resBinsArray;
    resBinning.CreateBinEdges(resBinsArray);
    GetHistogramRanges()->SetHistoTrackResidualPhiArr(resBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoEOverPArr().GetSize() == 0 )
  {
    TCustomBinning eopBinning;
    
    eopBinning.SetMinimum(GetHistogramRanges()->GetHistoEOverPMin());
    Float_t binWidth = ( GetHistogramRanges()->GetHistoEOverPMax() - GetHistogramRanges()->GetHistoEOverPMin() ) / GetHistogramRanges()->GetHistoEOverPBins();
    eopBinning.AddStep(GetHistogramRanges()->GetHistoEOverPMax(), binWidth);  
        
    TArrayD eopBinsArray;
    eopBinning.CreateBinEdges(eopBinsArray);
    GetHistogramRanges()->SetHistoEOverPArr(eopBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoNClusterCellArr().GetSize() == 0 )
  {
    TCustomBinning nceBinning;
    Int_t min = GetHistogramRanges()->GetHistoNClusterCellMin();
    if ( min <= 0 ) min = 1;
    nceBinning.SetMinimum(min);
    
    Int_t max = GetHistogramRanges()->GetHistoNClusterCellMax();

    if ( max >  30 ) nceBinning.AddStep(30 , 1);
    if ( max >  50 ) nceBinning.AddStep(50 , 2);
    if ( max > 100 ) nceBinning.AddStep(100, 5); 
    if ( max > 200 ) nceBinning.AddStep(200, 10); 
    
    if      ( max <=  30 ) nceBinning.AddStep(max, 1);
    else if ( max <=  50 ) nceBinning.AddStep(max, 2);
    else if ( max <= 100 ) nceBinning.AddStep(max, 5); 
    else if ( max <= 200 ) nceBinning.AddStep(max, 10); 
    else                   nceBinning.AddStep(max, 20);
    
    TArrayD nceBinsArray;
    nceBinning.CreateBinEdges(nceBinsArray);
    GetHistogramRanges()->SetHistoNClusterCellArr(nceBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoEDiffArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoEDiffMin();
    Float_t max = GetHistogramRanges()->GetHistoEDiffMax();
    Float_t binWidth = (max - min) / GetHistogramRanges()->GetHistoEDiffBins();

    TCustomBinning diffBinning;
    diffBinning.SetMinimum(min);
    diffBinning.AddStep(max, binWidth);  
  
    TArrayD diffBinsArray;
    diffBinning.CreateBinEdges(diffBinsArray);
    GetHistogramRanges()->SetHistoEDiffArr(diffBinsArray);
  }
  
  // Ratio, max 1
  if ( GetHistogramRanges()->GetHistoRatio1Arr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoRatio1Min();
    //Float_t max = GetHistogramRanges()->GetHistoRatio1Max();
    //Float_t binWidth = (max - min) / GetHistogramRanges()->GetHistoRatio1Bins();
    TCustomBinning fr1Binning;
    fr1Binning.SetMinimum(min);
    fr1Binning.AddStep(1, 0.05); 
    
    TArrayD fr1BinsArray;
    fr1Binning.CreateBinEdges(fr1BinsArray);
    GetHistogramRanges()->SetHistoRatio1Arr(fr1BinsArray);
  }
  
  // Ratio
  if ( GetHistogramRanges()->GetHistoRatioArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoRatioMin();
    Float_t max = GetHistogramRanges()->GetHistoRatioMax();
    //Float_t binWidth = (max - min) / GetHistogramRanges()->GetHistoRatioBins();
    TCustomBinning frBinning;
    frBinning.SetMinimum(min);
    frBinning.AddStep(max, 0.05); 
    
    TArrayD frBinsArray;
    frBinning.CreateBinEdges(frBinsArray);
    GetHistogramRanges()->SetHistoRatioArr(frBinsArray);
  }
  
  // Number of local maxima
  if ( GetHistogramRanges()->GetHistoNLMArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoNLMMin();
    if ( min <=0 ) min = 1;
    Float_t max = GetHistogramRanges()->GetHistoNLMMax();
    TCustomBinning nlmBinning;
    nlmBinning.SetMinimum(min);
    nlmBinning.AddStep(max, 1); 
    
    TArrayD nlmBinsArray;
    nlmBinning.CreateBinEdges(nlmBinsArray);
    GetHistogramRanges()->SetHistoNLMArr(nlmBinsArray);
  }
  
  // Number of MC overlaps
  if ( GetHistogramRanges()->GetHistoNoverlapArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoNoverlapMin();
    Float_t max = GetHistogramRanges()->GetHistoNoverlapMax();
    TCustomBinning novBinning;
    novBinning.SetMinimum(min);
    novBinning.AddStep(max, 1); 
    
    TArrayD novBinsArray;
    novBinning.CreateBinEdges(novBinsArray);
    GetHistogramRanges()->SetHistoNoverlapArr(novBinsArray);
  }
  
  // Track multiplicity
  if ( GetHistogramRanges()->GetHistoTrackMultiplicityArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoTrackMultiplicityMin();
    Float_t max = GetHistogramRanges()->GetHistoTrackMultiplicityMax();
    TCustomBinning mulBinning;
    mulBinning.SetMinimum(min);
    
    if ( !fFillHighMultHistograms )
    {
      mulBinning.AddStep(max, 1); 
    }
    else
    {
      if ( max >   50 ) mulBinning.AddStep(  50,  1);
      if ( max >  100 ) mulBinning.AddStep( 100,  2);
      if ( max >  200 ) mulBinning.AddStep( 200,  5);
      if ( max >  400 ) mulBinning.AddStep( 400, 10);
      if ( max > 1000 ) mulBinning.AddStep(1000, 20);
      if ( max > 2000 ) mulBinning.AddStep(2000, 50);
      if ( max > 5000 ) mulBinning.AddStep(5000,100);
      
      if      ( max <=   50 ) mulBinning.AddStep(max,  1);
      else if ( max <=  100 ) mulBinning.AddStep(max,  2);
      else if ( max <=  200 ) mulBinning.AddStep(max,  5);
      else if ( max <=  400 ) mulBinning.AddStep(max, 10);
      else if ( max <= 1000 ) mulBinning.AddStep(max, 20);
      else if ( max <= 2000 ) mulBinning.AddStep(max, 50);
      else if ( max <= 5000 ) mulBinning.AddStep(max,100);
      else                    mulBinning.AddStep(max,200);
    }
    TArrayD mulBinsArray;
    mulBinning.CreateBinEdges(mulBinsArray);
    GetHistogramRanges()->SetHistoTrackMultiplicityArr(mulBinsArray);
  }
  
  // Cluster multiplicity
  if ( GetHistogramRanges()->GetHistoNClustersArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoNClustersMin();
    Float_t max = GetHistogramRanges()->GetHistoNClustersMax();
    TCustomBinning mulBinning;
    mulBinning.SetMinimum(min);
    mulBinning.AddStep(max, 1); 
    
    TArrayD mulBinsArray;
    mulBinning.CreateBinEdges(mulBinsArray);
    GetHistogramRanges()->SetHistoNClustersArr(mulBinsArray);
  }
  
  // Isolation sum pT in cone
  if ( GetHistogramRanges()->GetHistoPtSumArr().GetSize() == 0 )
  {
    Float_t max = GetHistogramRanges()->GetHistoPtSumMax();
    TCustomBinning sumBinning;
    sumBinning.SetMinimum(0);
    
    if ( max >   2 ) sumBinning.AddStep(  2, 0.10); // 20
    if ( max >   4 ) sumBinning.AddStep(  4, 0.20); // 20
    if ( max >  10 ) sumBinning.AddStep( 10, 0.50); // 12
    if ( max >  25 ) sumBinning.AddStep( 25, 1.00); // 15
    if ( max >  50 ) sumBinning.AddStep( 50, 2.50); // 10
    if ( max > 100 ) sumBinning.AddStep(100, 5.00); // 10
    if ( max > 200 ) sumBinning.AddStep(200,10.00); // 10
 
    if      ( max <=   2 ) sumBinning.AddStep(max, 0.10); 
    else if ( max <=   4 ) sumBinning.AddStep(max, 0.20); 
    else if ( max <=  10 ) sumBinning.AddStep(max, 0.50); 
    else if ( max <=  25 ) sumBinning.AddStep(max, 1.00); 
    else if ( max <=  50 ) sumBinning.AddStep(max, 2.50); 
    else if ( max <= 100 ) sumBinning.AddStep(max, 5.00); 
    else if ( max <= 200 ) sumBinning.AddStep(max,10.00); 
    else                   sumBinning.AddStep(max,20.00); 
    
    TArrayD sumBinsArray;
    sumBinning.CreateBinEdges(sumBinsArray);
    GetHistogramRanges()->SetHistoPtSumArr(sumBinsArray);
  }
  
  // Isolation sum pT in cone after UE subtraction
  if ( GetHistogramRanges()->GetHistoPtSumSubArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoPtSumSubMin();
    Float_t max = GetHistogramRanges()->GetHistoPtSumSubMax();

    TCustomBinning sueBinning;
    sueBinning.SetMinimum(min);
    
    if (min <-200 && max >-200 ) sueBinning.AddStep(-200,20.); 
    if (min <-100 && max >-100 ) sueBinning.AddStep(-100,10.);
    if (min < -50 && max > -50 ) sueBinning.AddStep(-50,  5.); // 10
    if (min < -25 && max > -25 ) sueBinning.AddStep(-25, 2.5); // 10
    if (min < -10 && max > -10 ) sueBinning.AddStep(-10, 1.0); // 15
    if (min < -4  && max > -4  ) sueBinning.AddStep(-4 , 0.5); // 12
    if (min < -2  && max > -2  ) sueBinning.AddStep(-2 , 0.2); // 10
    if (min <  2  && max >  2  ) sueBinning.AddStep( 2 , 0.1); // 10
    if (min <  4  && max >  4  ) sueBinning.AddStep( 4 , 0.2); // 20
    if (min < 10  && max >  10 ) sueBinning.AddStep( 10, 0.5); // 12
    if (min < 25  && max >  25 ) sueBinning.AddStep( 25, 1.0); // 15
    if (min < 50  && max >  50 ) sueBinning.AddStep( 50, 2.5); // 10
    if (min < 100 && max >  100) sueBinning.AddStep(100, 5.0); // 10
    if (min < 200 && max >  200) sueBinning.AddStep(200,10.0); // 10
    
    if      ( max <=-200 ) sueBinning.AddStep(max,20.0); 
    else if ( max <=-100 ) sueBinning.AddStep(max,10.0);
    else if ( max <= -50 ) sueBinning.AddStep(max, 5.0);
    else if ( max <= -25 ) sueBinning.AddStep(max, 2.5); 
    else if ( max <= -10 ) sueBinning.AddStep(max, 1.0); 
    else if ( max <= -4  ) sueBinning.AddStep(max, 0.5); 
    else if ( max <= -2  ) sueBinning.AddStep(max, 0.2); 
    else if ( max <=  2  ) sueBinning.AddStep(max, 0.1); 
    else if ( max <=  4  ) sueBinning.AddStep(max, 0.2); 
    else if ( max <=  10 ) sueBinning.AddStep(max, 0.5); 
    else if ( max <=  25 ) sueBinning.AddStep(max, 1.0); 
    else if ( max <=  50 ) sueBinning.AddStep(max, 2.5); 
    else if ( max <=  100) sueBinning.AddStep(max, 5.0); 
    else if ( max <=  200) sueBinning.AddStep(max,10.0); 
    else                   sueBinning.AddStep(max,20.0); 
    
    TArrayD sueBinsArray;
    sueBinning.CreateBinEdges(sueBinsArray);
    GetHistogramRanges()->SetHistoPtSumSubArr(sueBinsArray);
  }
  
  // Exoticity
  if ( GetHistogramRanges()->GetHistoExoticityArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoExoticityMin();
    Float_t max = GetHistogramRanges()->GetHistoExoticityMax();
    
    if ( min < -1 ) min = -1;
    if ( max >  1 ) max =  1;

    TCustomBinning fBinning;
    fBinning.SetMinimum(min);
    
    if ( min < 0    && max > 0    ) fBinning.AddStep(0.000,0.2000); // 5
    if ( min < 0.5  && max > 0.5  ) fBinning.AddStep(0.500,0.1000); // 5
    if ( min < 0.7  && max > 0.7  ) fBinning.AddStep(0.700,0.0500); // 4
    if ( min < 0.8  && max > 0.8  ) fBinning.AddStep(0.800,0.0200); // 5
    if ( min < 0.85 && max > 0.85 ) fBinning.AddStep(0.850,0.0100); // 5
    if ( min < 0.9  && max > 0.9  ) fBinning.AddStep(0.900,0.0050); // 20 
    if ( min < 1.0  && max > 1.0  ) fBinning.AddStep(1.002,0.0020); // 51 
    
    if      ( max <= 0    ) fBinning.AddStep(max,0.2000); // 5
    else if ( max <= 0.5  ) fBinning.AddStep(max,0.1000); // 5
    else if ( max <= 0.7  ) fBinning.AddStep(max,0.0500); // 4
    else if ( max <= 0.8  ) fBinning.AddStep(max,0.0200); // 5
    else if ( max <= 0.85 ) fBinning.AddStep(max,0.0100); // 5
    else if ( max <= 0.9  ) fBinning.AddStep(max,0.0050); // 20 
    else if ( max <= 1.0  ) fBinning.AddStep(max,0.0020); // 51 

    TArrayD fBinsArray;
    fBinning.CreateBinEdges(fBinsArray);
    GetHistogramRanges()->SetHistoExoticityArr(fBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoMassArr().GetSize() == 0 )
  {
    Float_t min = GetHistogramRanges()->GetHistoMassMin();
    Float_t max = GetHistogramRanges()->GetHistoMassMax();
    Float_t binWidth = (max - min) / GetHistogramRanges()->GetHistoMassBins();
    
    TCustomBinning massBinning;
    massBinning.SetMinimum(min);
    massBinning.AddStep(max, binWidth);  
    
    TArrayD massBinsArray;
    massBinning.CreateBinEdges(massBinsArray);
    GetHistogramRanges()->SetHistoMassArr(massBinsArray);
  }
  
  if ( GetHistogramRanges()->GetHistoCentralityArr().GetSize() == 0 )
    AliAnaCaloTrackCorrBaseClass::SetEventCentralityBins();
}

//_________________________________________________________
/// Check if there is any track attached to this cluster.
/// \param cluster: pointer to calorimeter cluster.
/// \param event: AliVEvent pointer. Needed to get the tracks or the magnetic field.
/// \return kTRUE if cluster is matched by a track.
Bool_t AliAnaCaloTrackCorrBaseClass::IsTrackMatched(AliVCluster * cluster, AliVEvent* event) 
{ 
  Bool_t bRes = kFALSE, bEoP = kFALSE;
  return GetCaloPID()->IsTrackMatched(cluster, fCaloUtils, event, bEoP, bRes); 
} 

//_________________________________________________________
/// Check if there is any track attached to this cluster.
/// \param cluster: pointer to calorimeter cluster.
/// \param event: AliVEvent pointer. Needed to get the tracks or the magnetic field.
/// \param bEoP: If rejection is due to E over P cut, set it true, else false
/// \param bRes: If rejection is due to residual eta-phi cut, set it true, else false
/// \return kTRUE if cluster is matched by a track.
Bool_t AliAnaCaloTrackCorrBaseClass::IsTrackMatched(AliVCluster * cluster, AliVEvent* event, 
                                                    Bool_t & bEoP, Bool_t & bRes) 
{ 
  return GetCaloPID()->IsTrackMatched(cluster, fCaloUtils, event, bEoP, bRes); 
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

  printf("pT/E range          = [%2.2f,%2.2f]\n", fMinPt,fMaxPt) ;
  printf("Pair time cut       = %2.2f",fPairTimeCut);
  printf("Max Photon pT       =     %3.2f\n", fMaxPt) ;
  printf("Calo %s (%d) settings:\n",fCalorimeterString.Data(),fCalorimeter);
  printf("\t nSM %d, nRCU %d, First SM %d, Last SM %d; first TRD SM %d\n",
          fNModules,fNRCU,fFirstModule,fLastModule,fTRDSMCovered);
  printf("\t nMax cols %d, nMax Rows %d; full SM nMax Cols %d, nMax Rows %d; Rows Full SM: Min %d, Max %d\n",
          fNMaxCols,fNMaxRows,fNMaxColsFull,fNMaxRowsFull,fNMaxRowsFullMin,fNMaxRowsFullMax); 
    
  //printf("Check PID           =     %d\n",    fCheckCaloPID) ;
  printf("Recalculate PID     =     %d\n",    fRecalculateCaloPID) ;
  printf("Check Fiducial cut  =     %d\n",    fCheckFidCut) ;
  printf("Check Real Calo Acc =     %d\n",    fCheckRealCaloAcc) ;
  printf("Check MC labels     =     %d\n",    fDataMC);
  printf("Make plots?         =     %d\n",    fMakePlots);
  printf("Debug Level         =     %d\n",    fDebug);
  
  printf("Do mix event        =     %d\n",    fDoOwnMix);
  printf("\t N: z vert %d, RP %d, Centr %d, Mix event Pool size %d, Use track mult %d\n",
         fNZvertBin,fNrpBin,fNCentrBin,fNmaxMixEv,fUseTrackMultBins);
  printf("Fill histo: pile-up %d, high mult %d, embed %d, generated particles %d",
         fFillPileUpHistograms,fFillHighMultHistograms,fFillEmbedHistograms,fFillGenPartHisto);
  printf("Select embedded clusters/tracks %d\n",fSelectEmbededSignal);

  printf("Pt histograms bin array non constant? %d \n",fHistoPtBinNonConstantInArray);
} 

//_______________________________________________________________
/// Set the calorimeter for the analysis. A string.
//_______________________________________________________________
void AliAnaCaloTrackCorrBaseClass::SetCalorimeter(TString calo)
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


