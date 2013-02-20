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
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors (CTS).
// Not all MC particles/tracks/clusters are kept, some kinematical/fiducial restrictions are done.
// Mother class of : AliCaloTrackESDReader: Fills ESD data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackMCReader : Fills Kinematics data in 3 TObjArrays (PHOS, EMCAL, CTS)
//                 : AliCaloTrackAODReader: Fills AOD data in 3 TObjArrays (PHOS, EMCAL, CTS) 
//-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TFile.h>
#include <TGeoManager.h>

// ---- ANALYSIS system ----
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliMixedEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliTriggerAnalysis.h"
#include "AliESDVZERO.h"
#include "AliVCaloCells.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// ---- Detectors ----
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

// ---- CaloTrackCorr ---
#include "AliCalorimeterUtils.h"
#include "AliCaloTrackReader.h"

ClassImp(AliCaloTrackReader)
  
  
//________________________________________
AliCaloTrackReader::AliCaloTrackReader() : 
TObject(),                   fEventNumber(-1), //fCurrentFileName(""),
fDataType(0),                fDebug(0), 
fFiducialCut(0x0),           fCheckFidCut(kFALSE), 
fComparePtHardAndJetPt(0),   fPtHardAndJetPtFactor(0),
fComparePtHardAndClusterPt(0),fPtHardAndClusterPtFactor(0),
fCTSPtMin(0),                fEMCALPtMin(0),                  fPHOSPtMin(0), 
fCTSPtMax(0),                fEMCALPtMax(0),                  fPHOSPtMax(0),
fUseEMCALTimeCut(1),         fUseParamTimeCut(0),             fUseTrackTimeCut(0),
fEMCALTimeCutMin(-10000),    fEMCALTimeCutMax(10000),
fEMCALParamTimeCutMin(),     fEMCALParamTimeCutMax(),
fTrackTimeCutMin(-10000),    fTrackTimeCutMax(10000),
fUseTrackDCACut(0),
fAODBranchList(0x0),
fCTSTracks(0x0),             fEMCALClusters(0x0),             fPHOSClusters(0x0),
fEMCALCells(0x0),            fPHOSCells(0x0),
fInputEvent(0x0),            fOutputEvent(0x0),fMC(0x0),
fFillCTS(0),                 fFillEMCAL(0),                   fFillPHOS(0),
fFillEMCALCells(0),          fFillPHOSCells(0), 
fRecalculateClusters(kFALSE),fSelectEmbeddedClusters(kFALSE),
fTrackStatus(0),             fTrackFilterMask(0),             
fESDtrackCuts(0),            fESDtrackComplementaryCuts(0),   fConstrainTrack(kFALSE),
fSelectHybridTracks(0),      fSelectSPDHitTracks(kFALSE),
fTrackMult(0),               fTrackMultEtaCut(0.9),
fReadStack(kFALSE),          fReadAODMCParticles(kFALSE), 
fDeltaAODFileName(""),       fFiredTriggerClassName(""),      
fEventTriggerMask(0),        fMixEventTriggerMask(0),         fEventTriggerAtSE(0), 
fAnaLED(kFALSE),
fTaskName(""),               fCaloUtils(0x0), 
fMixedEvent(NULL),           fNMixedEvent(0),                 fVertex(NULL), 
fListMixedTracksEvents(),    fListMixedCaloEvents(),
fLastMixedTracksEvent(-1),   fLastMixedCaloEvent(-1),
fWriteOutputDeltaAOD(kFALSE),fOldAOD(kFALSE),                 fCaloFilterPatch(kFALSE),
fEMCALClustersListName(""),  fZvtxCut(0.),                    
fAcceptFastCluster(kFALSE),  fRemoveLEDEvents(kTRUE), 
fDoEventSelection(kFALSE),   fDoV0ANDEventSelection(kFALSE),
fDoVertexBCEventSelection(kFALSE),
fDoRejectNoTrackEvents(kFALSE),
fUseEventsWithPrimaryVertex(kFALSE),
fTriggerAnalysis (0x0),      fTimeStampEventSelect(0),
fTimeStampEventFracMin(0),   fTimeStampEventFracMax(0),
fTimeStampRunMin(0),         fTimeStampRunMax(0),
fNPileUpClusters(-1),        fNNonPileUpClusters(-1),         fNPileUpClustersCut(3),
fVertexBC(-200),             fRecalculateVertexBC(0),
fCentralityClass(""),            fCentralityOpt(0),
fEventPlaneMethod(""),       fImportGeometryFromFile(kFALSE), fImportGeometryFilePath("")
{
  //Ctor

  //Initialize parameters
  InitParameters();
}

//_______________________________________
AliCaloTrackReader::~AliCaloTrackReader() 
{
  //Dtor
  
  delete fFiducialCut ;
	
  if(fAODBranchList)
  {
    fAODBranchList->Delete();
    delete fAODBranchList ;
  }  
  
  if(fCTSTracks)
  {
    if(fDataType!=kMC)fCTSTracks->Clear() ; 
    else              fCTSTracks->Delete() ; 
    delete fCTSTracks ;
  }
  
  if(fEMCALClusters)
  {
    if(fDataType!=kMC)fEMCALClusters->Clear("C") ; 
    else              fEMCALClusters->Delete() ; 
    delete fEMCALClusters ;
  }
  
  if(fPHOSClusters)
  {
    if(fDataType!=kMC)fPHOSClusters->Clear("C") ; 
    else              fPHOSClusters->Delete() ; 
    delete fPHOSClusters ;
  }
  
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++) 
    {
      delete [] fVertex[i] ;
      
    }
    delete [] fVertex ;
	}
  
  delete fESDtrackCuts;
  delete fESDtrackComplementaryCuts;
  delete fTriggerAnalysis;
  
  //  Pointers not owned, done by the analysis frame
  //  if(fInputEvent)  delete fInputEvent ;
  //  if(fOutputEvent) delete fOutputEvent ;
  //  if(fMC)          delete fMC ;  
  //  Pointer not owned, deleted by maker
  //  if (fCaloUtils) delete fCaloUtils ;
  
}

//________________________________________________________________________
Bool_t  AliCaloTrackReader::AcceptDCA(const Float_t pt, const Float_t dca)
{
  // Accept track if DCA is smaller than function
  
  Float_t cut = fTrackDCACut[0]+fTrackDCACut[1]/TMath::Power(pt,fTrackDCACut[2]);
  
  if(TMath::Abs(dca) < cut)
    return kTRUE;
  else
    return kFALSE;
  
}

//________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndJetPt()
{
  // Check the event, if the requested ptHard is much smaller than the jet pT, then there is a problem.
  // Only for PYTHIA.
  
  //printf("AliCaloTrackReader::ComparePtHardAndJetPt() - GenHeaderName : %s\n",GetGenEventHeader()->ClassName());
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    TParticle * jet =  0;
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Int_t nTriggerJets =  pygeh->NTriggerJets();
    Float_t ptHard = pygeh->GetPtHard();
    
    if(fDebug > 1) 
      printf("AliCaloTrackReader::ComparePtHardAndJetPt() - Njets: %d, pT Hard %f\n",nTriggerJets, ptHard);
    
    Float_t tmpjet[]={0,0,0,0};
    for(Int_t ijet = 0; ijet< nTriggerJets; ijet++)
    {
      pygeh->TriggerJet(ijet, tmpjet);
      jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
      
      if(fDebug > 1) 
        printf("AliCaloTrackReader::ComparePtHardAndJetPt() - jet %d; pycell jet pT %f\n",ijet, jet->Pt());
      
      //Compare jet pT and pt Hard
      if(jet->Pt() > fPtHardAndJetPtFactor * ptHard)
      {
        printf("AliCaloTrackReader::ComparePtHardAndJetPt() - Reject jet event with : pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n",
                ptHard, jet->Pt(), fPtHardAndJetPtFactor);
        return kFALSE;
      }
    }
    
    if(jet) delete jet; 
  }
  
  return kTRUE ;
  
}

//____________________________________________________________________
Bool_t AliCaloTrackReader::ComparePtHardAndClusterPt()
{
  // Check the event, if the requested ptHard is smaller than the calorimeter cluster E, then there is a problem.
  // Only for PYTHIA.
  
  if(!strcmp(GetGenEventHeader()->ClassName(), "AliGenPythiaEventHeader"))
  {
    AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) GetGenEventHeader();
    Float_t ptHard = pygeh->GetPtHard();
    
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++) 
    {
      AliVCluster * clus = fInputEvent->GetCaloCluster(iclus) ; 
      Float_t ecluster = clus->E();
      
      if(ecluster > fPtHardAndClusterPtFactor*ptHard) 
      {
        printf("AliCaloTrackReader::ComparePtHardAndClusterPt() - Reject : ecluster %2.2f, calo %d, factor %2.2f, ptHard %f\n",ecluster,clus->GetType(),fPtHardAndClusterPtFactor,ptHard);

        return kFALSE;
      }
    }
    
  }
  
  return kTRUE ;
  
}

//____________________________________________
AliStack* AliCaloTrackReader::GetStack() const 
{
  //Return pointer to stack
  if(fMC)
    return fMC->Stack();
  else
  {
    if(fDebug > 1) printf("AliCaloTrackReader::GetStack() - Stack is not available\n"); 
    return 0x0 ;
  }
}

//__________________________________________________
TString AliCaloTrackReader::GetFiredTriggerClasses() 
{ 
  // List of triggered classes in a TString

  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (GetInputEvent());
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (GetInputEvent());
    
  if     (esdevent) return esdevent->GetFiredTriggerClasses();
  else if(aodevent) return aodevent->GetFiredTriggerClasses();
  else              return ""; // Mixed Event, MC event, does not have this trigger info

}

//______________________________________________
AliHeader* AliCaloTrackReader::GetHeader() const 
{
  //Return pointer to header
  if(fMC)
  {
    return fMC->Header();
  }
  else
  {
    printf("AliCaloTrackReader::Header is not available\n"); 
    return 0x0 ;
  }
}

//______________________________________________________________
AliGenEventHeader* AliCaloTrackReader::GetGenEventHeader() const 
{
  //Return pointer to Generated event header
  if     (ReadStack() && fMC)
  {
    return fMC->GenEventHeader();
  }
  else if(ReadAODMCParticles() && GetAODMCHeader())
  {
    //printf("AliCaloTrackReader::GetGenEventHeader() - N headers %d\n",GetAODMCHeader()->GetNCocktailHeaders());
    if( GetAODMCHeader()->GetNCocktailHeaders() > 0)
      return GetAODMCHeader()->GetCocktailHeader(0) ;
    else 
      return 0x0;
  }
  else 
  {
    //printf("AliCaloTrackReader::GetGenEventHeader() - MC header not available! \n");
    return 0;
  }
}

//____________________________________________________________________
TClonesArray* AliCaloTrackReader::GetAODMCParticles() const 
{
  //Return list of particles in AOD. Do it for the corresponding input event.
  
  TClonesArray * rv = NULL ; 
  if(fDataType == kAOD)
  {
    //Normal input AOD
    AliAODEvent * evt = dynamic_cast<AliAODEvent*> (fInputEvent) ;
    if(evt)
      rv = (TClonesArray*)evt->FindListObject("mcparticles");
    else  
      printf("AliCaloTrackReader::GetAODMCParticles() - Null AOD event \n"); 
  } 
  else 
  {
    printf("AliCaloTrackReader::GetAODMCParticles() - Input are not AODs\n"); 
  }
  
  return rv ; 
}

//________________________________________________________
AliAODMCHeader* AliCaloTrackReader::GetAODMCHeader() const 
{
  //Return MC header in AOD. Do it for the corresponding input event.
  
  AliAODMCHeader *mch = NULL;
  
  if(fDataType == kAOD)
  {
    AliAODEvent * aod = dynamic_cast<AliAODEvent*> (fInputEvent);
    if(aod) mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
  }
  else 
  {
    printf("AliCaloTrackReader::GetAODMCHeader() - Input are not AODs\n");
  }
  
  return mch;
}

//___________________________________________________________
Int_t AliCaloTrackReader::GetVertexBC(const AliVVertex * vtx)
{
  // Get the vertex BC
    
  Int_t vertexBC=vtx->GetBC();
  if(!fRecalculateVertexBC) return vertexBC;
  
  // In old AODs BC not stored, recalculate it
  // loop over the global track and select those which have small DCA to primary vertex (e.g. primary).
  // If at least one of these primaries has valid BC != 0, then this vertex is a pile-up candidate.
  // Execute after CTS
  Double_t bz  = fInputEvent->GetMagneticField();
  Bool_t   bc0 = kFALSE;
  Int_t    ntr = GetCTSTracks()->GetEntriesFast();
  //printf("N Tracks %d\n",ntr);
  
  for(Int_t i = 0 ; i < ntr ; i++)
  {
    AliVTrack * track =  (AliVTrack*) (GetCTSTracks()->At(i));
    
    //Check if has TOF info, if not skip
    ULong_t status  = track->GetStatus();
    Bool_t  okTOF   = (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ;
    vertexBC        = track->GetTOFBunchCrossing(bz);
    Float_t pt      = track->Pt();
    
    if(!okTOF) continue;
    
    // Get DCA x, y
    Double_t dca[2]   = {1e6,1e6};
    Double_t covar[3] = {1e6,1e6,1e6};
    track->PropagateToDCA(vtx,bz,100.,dca,covar);
    
    if(AcceptDCA(pt,dca[0]))
    {
      if     (vertexBC !=0 && fVertexBC != AliVTrack::kTOFBCNA) return vertexBC;
      else if(vertexBC == 0)                                    bc0 = kTRUE;
    }
  }
  
  if( bc0 ) vertexBC = 0 ;
  else      vertexBC = AliVTrack::kTOFBCNA ;
  
  return vertexBC;
  
}

//_____________________________
void AliCaloTrackReader::Init()
{
  //Init reader. Method to be called in AliAnaPartCorrMaker

  //printf(" AliCaloTrackReader::Init() %p \n",gGeoManager);

  if(fReadStack && fReadAODMCParticles)
  {
    printf("AliCaloTrackReader::Init() - Cannot access stack and mcparticles at the same time, change them \n");
    fReadStack          = kFALSE;
    fReadAODMCParticles = kFALSE;
  }
  
  // Init geometry, I do not like much to do it like this ...
  if(fImportGeometryFromFile && !gGeoManager) 
  {
    if(fImportGeometryFilePath=="") // If not specified, set a default location
    fImportGeometryFilePath = "$ALICE_ROOT/OADB/EMCAL/geometry_2011.root"; // "$ALICE_ROOT/EVE/alice-data/default_geo.root"

    printf("AliCaloTrackReader::Init() - Import %s\n",fImportGeometryFilePath.Data());
    TGeoManager::Import(fImportGeometryFilePath) ; // default need file "geometry.root" in local dir!!!!
  }

  if(!fESDtrackCuts)
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //initialize with TPC only tracks

}

//_______________________________________
void AliCaloTrackReader::InitParameters()
{
  //Initialize the parameters of the analysis.
  fDataType   = kESD ;
  fCTSPtMin   = 0.1 ;
  fEMCALPtMin = 0.1 ;
  fPHOSPtMin  = 0.1 ;
  fCTSPtMax   = 1000. ;
  fEMCALPtMax = 1000. ;
  fPHOSPtMax  = 1000. ;
  
  //Track DCA cuts
  // dca_xy cut = 0.0105+0.0350/TMath::Power(pt,1.1);
  fTrackDCACut[0] = 0.0105; 
  fTrackDCACut[1] = 0.0350; 
  fTrackDCACut[2] = 1.1;
  
  //Do not filter the detectors input by default.
  fFillEMCAL      = kFALSE;
  fFillPHOS       = kFALSE;
  fFillCTS        = kFALSE;
  fFillEMCALCells = kFALSE;
  fFillPHOSCells  = kFALSE;
  
  fReadStack             = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fReadAODMCParticles    = kFALSE; // Check in the constructor of the other readers if it was set or in the configuration file
  fDeltaAODFileName      = "deltaAODPartCorr.root";
  fFiredTriggerClassName = "";
  fEventTriggerMask      = AliVEvent::kAny;
  fMixEventTriggerMask   = AliVEvent::kAnyINT;
  fEventTriggerAtSE      = kTRUE; // Use only events that pass event selection at SE base class
  
  fAcceptFastCluster = kTRUE;
  fAnaLED            = kFALSE;
  
  //We want tracks fitted in the detectors:
  //fTrackStatus=AliESDtrack::kTPCrefit;
  //fTrackStatus|=AliESDtrack::kITSrefit; 
  fTrackStatus     = 0;
  fTrackFilterMask = 128; //For AODs, but what is the difference between fTrackStatus and fTrackFilterMask?
  
  fESDtrackCuts = 0;
  fESDtrackComplementaryCuts = 0;
  
  fConstrainTrack = kFALSE ; // constrain tracks to vertex
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0; 
  fV0Mul[0] = 0;   fV0Mul[1] = 0; 
  
  fZvtxCut   = 10.;
  
  fNMixedEvent = 1;
  
  fPtHardAndJetPtFactor     = 7.;
  fPtHardAndClusterPtFactor = 1.;

  //Centrality
  fCentralityClass  = "V0M";
  fCentralityOpt    = 10;
  fCentralityBin[0] = fCentralityBin[1]=-1;
  
  fEventPlaneMethod = "V0";

  // Allocate memory (not sure this is the right place)
  fCTSTracks       = new TObjArray();
  fEMCALClusters   = new TObjArray();
  fPHOSClusters    = new TObjArray(); 
  fTriggerAnalysis = new AliTriggerAnalysis;
  fAODBranchList   = new TList ;

  fImportGeometryFromFile = kFALSE;
  
  fPileUpParamSPD[0] = 3   ; fPileUpParamSPD[1] = 0.8 ;
  fPileUpParamSPD[2] = 3.0 ; fPileUpParamSPD[3] = 2.0 ; fPileUpParamSPD[4] = 5.0;
  
  // Parametrized time cut (LHC11d)
  fEMCALParamTimeCutMin[0] =-5; fEMCALParamTimeCutMin[1] =-1 ; fEMCALParamTimeCutMin[2] = 3.5 ; fEMCALParamTimeCutMin[3] = 1.  ;   
  fEMCALParamTimeCutMax[0] = 5; fEMCALParamTimeCutMax[1] = 50; fEMCALParamTimeCutMax[2] = 0.45; fEMCALParamTimeCutMax[3] = 1.25;   

  // Parametrized time cut (LHC11c)
  //fEMCALParamTimeCutMin[0] =-5;   fEMCALParamTimeCutMin[1] =-1 ; fEMCALParamTimeCutMin[2] = 1.87; fEMCALParamTimeCutMin[3] = 0.4;   
  //fEMCALParamTimeCutMax[0] = 3.5; fEMCALParamTimeCutMax[1] = 50; fEMCALParamTimeCutMax[2] = 0.15; fEMCALParamTimeCutMax[3] = 1.6;
  
  fTimeStampRunMin = -1;
  fTimeStampRunMax = 1e12;
  fTimeStampEventFracMin = -1;
  fTimeStampEventFracMax = 2;

  for(Int_t i = 0; i < 19; i++)
  {
    fEMCalBCEvent   [i] = 0;
    fEMCalBCEventCut[i] = 0;
    fTrackBCEvent   [i] = 0;
    fTrackBCEventCut[i] = 0;
  }
  
}

//___________________________________________________________
Bool_t AliCaloTrackReader::IsInTimeWindow(const Double_t tof, const Float_t energy) const
{
  // Cluster time selection window
  
  // Parametrized cut depending on E
  if(fUseParamTimeCut)
  {
    Float_t minCut= fEMCALParamTimeCutMin[0]+fEMCALParamTimeCutMin[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMin[2])/fEMCALParamTimeCutMin[3]);
    Float_t maxCut= fEMCALParamTimeCutMax[0]+fEMCALParamTimeCutMax[1]*TMath::Exp(-(energy-fEMCALParamTimeCutMax[2])/fEMCALParamTimeCutMax[3]);
    //printf("tof %f, minCut %f, maxCut %f\n",tof,minCut,maxCut);
    if( tof < minCut || tof > maxCut )  return kFALSE ;
  }
  
  //In any case, the time should to be larger than the fixed window ...
  if( tof < fEMCALTimeCutMin  || tof > fEMCALTimeCutMax )  return kFALSE ;
  
  return kTRUE ;
}

//________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPD() const
{
  // Check if event is from pile-up determined by SPD
  // Default values: (3, 0.8, 3., 2., 5.)
  return fInputEvent->IsPileupFromSPD((Int_t) fPileUpParamSPD[0] , fPileUpParamSPD[1] , 
                                              fPileUpParamSPD[2] , fPileUpParamSPD[3] , fPileUpParamSPD[4] ); 
  //printf("Param : %d, %2.2f, %2.2f, %2.2f, %2.2f\n",(Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);

}

//__________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCal() const
{
  // Check if event is from pile-up determined by EMCal
  if(fNPileUpClusters > fNPileUpClustersCut) return kTRUE ;
  else                                       return kFALSE;
}

//________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndEMCal() const
{
  // Check if event is from pile-up determined by SPD and EMCal
  if( IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//_______________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDOrEMCal() const
{
  // Check if event is from pile-up determined by SPD or EMCal
  if( IsPileUpFromSPD() || IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromSPDAndNotEMCal() const
{
  // Check if event is from pile-up determined by SPD and not by EMCal
  if( IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                          return kFALSE;
}

//___________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromEMCalAndNotSPD() const
{
  // Check if event is from pile-up determined by EMCal, not by SPD
  if( !IsPileUpFromSPD() && IsPileUpFromEMCal()) return kTRUE ;
  else                                           return kFALSE;
}

//______________________________________________________________
Bool_t AliCaloTrackReader::IsPileUpFromNotSPDAndNotEMCal() const
{
  // Check if event not from pile-up determined neither by SPD nor by EMCal
  if( !IsPileUpFromSPD() && !IsPileUpFromEMCal()) return kTRUE ;
  else                                            return kFALSE;
}

//________________________________________________________
void AliCaloTrackReader::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n",    GetName(), GetTitle() ) ;
  printf("Task name      : %s\n",          fTaskName.Data()) ;
  printf("Data type      : %d\n",          fDataType) ;
  printf("CTS Min pT     : %2.1f GeV/c\n", fCTSPtMin) ;
  printf("EMCAL Min pT   : %2.1f GeV/c\n", fEMCALPtMin) ;
  printf("PHOS Min pT    : %2.1f GeV/c\n", fPHOSPtMin) ;
  printf("CTS Max pT     : %2.1f GeV/c\n", fCTSPtMax) ;
  printf("EMCAL Max pT   : %2.1f GeV/c\n", fEMCALPtMax) ;
  printf("PHOS Max pT    : %2.1f GeV/c\n", fPHOSPtMax) ;
  printf("EMCAL Time Cut: %3.1f < TOF  < %3.1f\n", fEMCALTimeCutMin, fEMCALTimeCutMax);
  printf("Use CTS         =     %d\n",     fFillCTS) ;
  printf("Use EMCAL       =     %d\n",     fFillEMCAL) ;
  printf("Use PHOS        =     %d\n",     fFillPHOS) ;
  printf("Use EMCAL Cells =     %d\n",     fFillEMCALCells) ;
  printf("Use PHOS  Cells =     %d\n",     fFillPHOSCells) ;
  printf("Track status    =     %d\n", (Int_t) fTrackStatus) ;
  printf("AODs Track filter mask  =  %d or hybrid %d, SPD hit %d\n", (Int_t) fTrackFilterMask,fSelectHybridTracks,fSelectSPDHitTracks) ;
  printf("Track Mult Eta Cut =  %d\n", (Int_t) fTrackMultEtaCut) ;
  printf("Write delta AOD =     %d\n",     fWriteOutputDeltaAOD) ;
  printf("Recalculate Clusters = %d\n",    fRecalculateClusters) ;
  
  printf("Use Triggers selected in SE base class %d; If not what trigger Mask? %d; Trigger max for mixed %d \n", 
         fEventTriggerAtSE, fEventTriggerMask,fMixEventTriggerMask);
  
  if(fComparePtHardAndClusterPt)
	  printf("Compare jet pt and pt hard to accept event, factor = %2.2f",fPtHardAndJetPtFactor);
  
  if(fComparePtHardAndClusterPt)
	  printf("Compare cluster pt and pt hard to accept event, factor = %2.2f",fPtHardAndClusterPtFactor);
  
  printf("Read Kine from, stack? %d, AOD ? %d \n", fReadStack, fReadAODMCParticles) ;
  printf("Delta AOD File Name =     %s\n", fDeltaAODFileName.Data()) ;
  printf("Centrality: Class %s, Option %d, Bin [%d,%d] \n", fCentralityClass.Data(),fCentralityOpt,fCentralityBin[0], fCentralityBin[1]) ;
  
  printf("    \n") ;
  
} 

//_________________________________________________________________________
Bool_t AliCaloTrackReader::FillInputEvent(const Int_t iEntry, 
                                          const char * /*currentFileName*/) 
{
  //Fill the event counter and input lists that are needed, called by the analysis maker.
  
  fEventNumber = iEntry;
  //fCurrentFileName = TString(currentFileName);
  if(!fInputEvent)
  {
	  if(fDebug >= 0) printf("AliCaloTrackReader::FillInputEvent() - Input event not available, skip event analysis\n");
	  return kFALSE;
  }
  
  //Select events only fired by a certain trigger configuration if it is provided
  Int_t eventType = 0;
  if(fInputEvent->GetHeader())
	  eventType = ((AliVHeader*)fInputEvent->GetHeader())->GetEventType();
  
  if (GetFiredTriggerClasses().Contains("FAST")  && !GetFiredTriggerClasses().Contains("ALL") && !fAcceptFastCluster) 
  {
    if(fDebug > 0)  printf("AliCaloTrackReader::FillInputEvent - Do not count events from fast cluster, trigger name %s\n",fFiredTriggerClassName.Data());
    return kFALSE;
  }
  
  //-------------------------------------------------------------------------------------
  // Reject event if large clusters with large energy
  // Use only for LHC11a data for the moment, and if input is clusterizer V1 or V1+unfolding
  // If clusterzer NxN or V2 it does not help
  //-------------------------------------------------------------------------------------
  Int_t run = fInputEvent->GetRunNumber();
  if( fRemoveLEDEvents && run > 146857  && run < 146861 )
  {
    //printf("Event %d\n",GetEventNumber());
    
    // Count number of cells with energy larger than 0.1 in SM3, cut on this number
    Int_t ncellsSM3 = 0;
    for(Int_t icell = 0; icell < fInputEvent->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      Int_t absID = fInputEvent->GetEMCALCells()->GetCellNumber(icell);
      Int_t sm    = GetCaloUtils()->GetEMCALGeometry()->GetSuperModuleNumber(absID);
      if(fInputEvent->GetEMCALCells()->GetAmplitude(icell) > 0.1 && sm==3) ncellsSM3++;
    }
    
    Int_t ncellcut = 21;
    if(fFiredTriggerClassName.Contains("EMC")) ncellcut = 35;
    
    if(ncellsSM3 >= ncellcut)
    {
      if(fDebug > 0) printf(" AliCaloTrackReader::FillInputEvent() - reject event with ncells in SM3 %d\n",ncellsSM3);
      return kFALSE;
    }
  }// Remove LED events
  
  // Reject pure LED events?
  if( fFiredTriggerClassName  !="" && !fAnaLED)
  {
    //printf("Event type %d\n",eventType);
    if(eventType!=7)
      return kFALSE; //Only physics event, do not use for simulated events!!!
    
    if(fDebug > 0) 
      printf("AliCaloTrackReader::FillInputEvent() - FiredTriggerClass <%s>, selected class <%s>, compare name %d\n",
             GetFiredTriggerClasses().Data(),fFiredTriggerClassName.Data(), GetFiredTriggerClasses().Contains(fFiredTriggerClassName));
    
    if( !GetFiredTriggerClasses().Contains(fFiredTriggerClassName) ) return kFALSE;
    else if(fDebug > 0) printf("AliCaloTrackReader::FillInputEvent() - Accepted triggered event\n");
  }
  else if(fAnaLED)
  {
    //	  kStartOfRun =       1,    // START_OF_RUN
    //	  kEndOfRun =         2,    // END_OF_RUN
    //	  kStartOfRunFiles =  3,    // START_OF_RUN_FILES
    //	  kEndOfRunFiles =    4,    // END_OF_RUN_FILES
    //	  kStartOfBurst =     5,    // START_OF_BURST
    //	  kEndOfBurst =       6,    // END_OF_BURST
    //	  kPhysicsEvent =     7,    // PHYSICS_EVENT
    //	  kCalibrationEvent = 8,    // CALIBRATION_EVENT
    //	  kFormatError =      9,    // EVENT_FORMAT_ERROR
    //	  kStartOfData =      10,   // START_OF_DATA
    //	  kEndOfData =        11,   // END_OF_DATA
    //	  kSystemSoftwareTriggerEvent   = 12, // SYSTEM_SOFTWARE_TRIGGER_EVENT
    //	  kDetectorSoftwareTriggerEvent = 13  // DETECTOR_SOFTWARE_TRIGGER_EVENT
    
	  if(eventType!=7 && fDebug > 1 )printf("AliCaloTrackReader::FillInputEvent() - DO LED, Event Type <%d>, 8 Calibration \n",  eventType);
	  if(eventType!=8)return kFALSE;
  }
  
  //In case of analysis of events with jets, skip those with jet pt > 5 pt hard	
  if(fComparePtHardAndJetPt) 
  {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
  }
  
  if(fComparePtHardAndClusterPt) 
  {
    if(!ComparePtHardAndClusterPt()) return kFALSE ;
  }
  
  //Fill Vertex array
  FillVertexArray();
  //Reject events with Z vertex too large, only for SE analysis, if not, cut on the analysis code
  if(!GetMixedEvent() && TMath::Abs(fVertex[0][2]) > fZvtxCut) return kFALSE;  
  
  //------------------------------------------------------
  //Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------
  if(fDataType==kESD && fTimeStampEventSelect)
  {
    AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
    if(esd)
    {
      Int_t timeStamp = esd->GetTimeStamp();
      Float_t timeStampFrac = 1.*(timeStamp-fTimeStampRunMin) / (fTimeStampRunMax-fTimeStampRunMin);
      
      //printf("stamp0 %d, max0 %d, frac %f\n", timeStamp-fTimeStampRunMin,fTimeStampRunMax-fTimeStampRunMin, timeStampFrac);
      
      if(timeStampFrac < fTimeStampEventFracMin || timeStampFrac > fTimeStampEventFracMax) return kFALSE;
    }
    //printf("\t accept time stamp\n");
  }

  
  //------------------------------------------------------
  //Event rejection depending on vertex, pileup, v0and
  //------------------------------------------------------

  if(fUseEventsWithPrimaryVertex)
  {
    if( !CheckForPrimaryVertex() )              return kFALSE;
    if( TMath::Abs(fVertex[0][0] ) < 1.e-6 && 
        TMath::Abs(fVertex[0][1] ) < 1.e-6 && 
        TMath::Abs(fVertex[0][2] ) < 1.e-6    ) return kFALSE;
  }
  
  //printf("Reader : IsPileUp %d, Multi %d\n",IsPileUpFromSPD(),fInputEvent->IsPileupFromSPDInMultBins());
  
  if(fDoEventSelection)
  {
    if(!fCaloFilterPatch)
    {
      // Do not analyze events with pileup
      Bool_t bPileup = IsPileUpFromSPD();
      //IsPileupFromSPDInMultBins() // method to try
      //printf("pile-up %d, %d, %2.2f, %2.2f, %2.2f, %2.2f\n",bPileup, (Int_t) fPileUpParamSPD[0], fPileUpParamSPD[1], fPileUpParamSPD[2], fPileUpParamSPD[3], fPileUpParamSPD[4]);
      if(bPileup) return kFALSE;
      
      if(fDoV0ANDEventSelection)
      {
        Bool_t bV0AND = kTRUE; 
        AliESDEvent* esd = dynamic_cast<AliESDEvent*> (fInputEvent);
        if(esd) 
          bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0AND);
        //else bV0AND = //FIXME FOR AODs
        if(!bV0AND) return kFALSE;
      }
    }//CaloFilter patch
    else
    { 
      if(fInputEvent->GetNumberOfCaloClusters() > 0) 
      {
        AliVCluster * calo = fInputEvent->GetCaloCluster(0);
        if(calo->GetNLabels() == 4)
        {
          Int_t * selection = calo->GetLabels();
          Bool_t bPileup = selection[0];
          if(bPileup) return kFALSE;
          
          Bool_t bGoodV = selection[1]; 
          if(fUseEventsWithPrimaryVertex && !bGoodV) return kFALSE;
          
          if(fDoV0ANDEventSelection)
          {
            Bool_t bV0AND = selection[2]; 
            if(!bV0AND) return kFALSE;
          }
          
          fTrackMult = selection[3];
          if(fTrackMult == 0) return kFALSE;
        } 
        else 
        {
          //First filtered AODs, track multiplicity stored there.  
          fTrackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
          if(fTrackMult == 0) return kFALSE;          
        }
      }//at least one cluster
      else 
      {
        //printf("AliCaloTrackReader::FillInputEvent() - No clusters in event\n");
        //Remove events with  vertex (0,0,0), bad vertex reconstruction
        if(fUseEventsWithPrimaryVertex && TMath::Abs(fVertex[0][0]) < 1.e-6 && TMath::Abs(fVertex[0][1]) < 1.e-6 && TMath::Abs(fVertex[0][2]) < 1.e-6) return kFALSE;
        
        //First filtered AODs, track multiplicity stored there.  
        fTrackMult = (Int_t) ((AliAODHeader*)fInputEvent->GetHeader())->GetCentrality();
        if(fTrackMult == 0) return kFALSE;
      }// no cluster
    }// CaloFileter patch
  }// Event selection/AliceSoft/AliRoot/trunk/PWG/CaloTrackCorrBase/AliCaloTrackReader.h
  //------------------------------------------------------
  
  //Check if there is a centrality value, PbPb analysis, and if a centrality bin selection is requested
  //If we need a centrality bin, we select only those events in the corresponding bin.
  if(GetCentrality() && fCentralityBin[0]>=0 && fCentralityBin[1]>=0 && fCentralityOpt==100)
  {
    Int_t cen = GetEventCentrality();
    if(cen > fCentralityBin[1] || cen < fCentralityBin[0]) return kFALSE; //reject events out of bin.
  }
    
  //Fill the arrays with cluster/tracks/cells data
  
  if(!fEventTriggerAtSE)
  {
    // In case of mixing analysis, accept MB events, not only Trigger
    // Track and cluster arrays filled for MB in order to create the pool in the corresponding analysis
    // via de method in the base class FillMixedEventPool()
    
    AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(manager->GetInputEventHandler());
    
    if(!inputHandler) return kFALSE ;  // to content coverity
    
    UInt_t isTrigger = inputHandler->IsEventSelected() & fEventTriggerMask;
    UInt_t isMB      = inputHandler->IsEventSelected() & fMixEventTriggerMask;
    
    if(!isTrigger && !isMB) return kFALSE;
    
    //printf("Selected triggered event : %s\n",GetFiredTriggerClasses().Data());
  }
  
  // Get the main vertex BC, in case not available
  // it is calculated in FillCTS checking the BC of tracks
  // with DCA small (if cut applied, if open)
  fVertexBC=fInputEvent->GetPrimaryVertex()->GetBC();
  
  if(fFillCTS)
  {
    FillInputCTS();
    //Accept events with at least one track
    if(fTrackMult == 0 && fDoRejectNoTrackEvents) return kFALSE ;
  }
  
  if(fDoVertexBCEventSelection)
  {
    if(fVertexBC!=0 && fVertexBC!=AliVTrack::kTOFBCNA) return kFALSE ;
  }
    
  if(fFillEMCALCells)
    FillInputEMCALCells();
  
  if(fFillPHOSCells)
    FillInputPHOSCells();
  
  if(fFillEMCAL)
    FillInputEMCAL();
  
  if(fFillPHOS)
    FillInputPHOS();
  
  FillInputVZERO();

  
  return kTRUE ;
}

//___________________________________
void AliCaloTrackReader::ResetLists() 
{
  //  Reset lists, called by the analysis maker 
  
  if(fCTSTracks)       fCTSTracks     -> Clear();
  if(fEMCALClusters)   fEMCALClusters -> Clear("C");
  if(fPHOSClusters)    fPHOSClusters  -> Clear("C");
  
  fV0ADC[0] = 0;   fV0ADC[1] = 0; 
  fV0Mul[0] = 0;   fV0Mul[1] = 0; 
  
}

//____________________________________________________________
void AliCaloTrackReader::SetInputEvent(AliVEvent* const input)  
{
  fInputEvent  = input;
  fMixedEvent = dynamic_cast<AliMixedEvent*>(GetInputEvent()) ; 
  if (fMixedEvent) 
    fNMixedEvent = fMixedEvent->GetNumberOfEvents() ; 
  
  //Delete previous vertex
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++) 
    {
      delete [] fVertex[i] ; 
    }
    delete [] fVertex ;
  }
  
  fVertex = new Double_t*[fNMixedEvent] ; 
  for (Int_t i = 0; i < fNMixedEvent; i++) 
  {
    fVertex[i] = new Double_t[3] ; 
    fVertex[i][0] = 0.0 ; 
    fVertex[i][1] = 0.0 ; 
    fVertex[i][2] = 0.0 ; 
  }
}

//__________________________________________________
Int_t AliCaloTrackReader::GetEventCentrality() const 
{
  //Return current event centrality
  
  if(GetCentrality())
  {
    if     (fCentralityOpt==100) return (Int_t) GetCentrality()->GetCentralityPercentile(fCentralityClass); // 100 bins max
    else if(fCentralityOpt==10)  return GetCentrality()->GetCentralityClass10(fCentralityClass);// 10 bins max
    else if(fCentralityOpt==20)  return GetCentrality()->GetCentralityClass5(fCentralityClass); // 20 bins max
    else 
    {
      printf("AliCaloTrackReader::GetEventCentrality() - Unknown centrality option %d, use 10, 20 or 100\n",fCentralityOpt);
      return -1;
    } 
  }
  else return -1;
  
}

//_____________________________________________________
Double_t AliCaloTrackReader::GetEventPlaneAngle() const 
{
  //Return current event centrality
  
  if(GetEventPlane())
  {
    Float_t ep =  GetEventPlane()->GetEventplane(GetEventPlaneMethod(), GetInputEvent());
    
    if(GetEventPlaneMethod()=="Q" && (ep < 0 || ep > TMath::Pi())) 
    {
      if(fDebug > 0 ) printf("AliCaloTrackReader::GetEventPlaneAngle() -  Bad EP for <Q> method : %f\n",ep);
      return -1000;
    }
    else if(GetEventPlaneMethod().Contains("V0")  ) 
    {
      if((ep > TMath::Pi()/2 || ep < -TMath::Pi()/2))
      {
        if(fDebug > 0 ) printf("AliCaloTrackReader::GetEventPlaneAngle() -  Bad EP for <%s> method : %f\n",GetEventPlaneMethod().Data(), ep);
        return -1000;
      }
      
      ep+=TMath::Pi()/2; // put same range as for <Q> method
      
    }
  
    //printf("AliCaloTrackReader::GetEventPlaneAngle() = %f\n",ep);
    if(fDebug > 0 )
    {
      if     (ep > TMath::Pi()) printf("AliCaloTrackReader::GetEventPlaneAngle() - Too large angle = %f\n",ep);
      else if(ep < 0          ) printf("AliCaloTrackReader::GetEventPlaneAngle() - Negative angle = %f\n" ,ep);
    }
    
    return ep;
  }
  else
  {
    if(fDataType!=kMC && fDebug > 0) printf("AliCaloTrackReader::GetEventPlaneAngle() -  No EP pointer\n");
    return -1000;
  } 
  
}

//__________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3]) const 
{
  //Return vertex position to be used for single event analysis
  vertex[0]=fVertex[0][0];  
  vertex[1]=fVertex[0][1];  
  vertex[2]=fVertex[0][2];
}

//____________________________________________________________
void AliCaloTrackReader::GetVertex(Double_t vertex[3], 
                                   const Int_t evtIndex) const 
{
  //Return vertex position for mixed event, recover the vertex in a particular event.
  
  vertex[0]=fVertex[evtIndex][0];  vertex[1]=fVertex[evtIndex][1];  vertex[2]=fVertex[evtIndex][2];
  
}

//________________________________________
void AliCaloTrackReader::FillVertexArray() 
{
  
  //Fill data member with vertex
  //In case of Mixed event, multiple vertices
  
  //Delete previous vertex
  if(fVertex)
  {
    for (Int_t i = 0; i < fNMixedEvent; i++) 
    {
      delete [] fVertex[i] ; 
    }
    delete [] fVertex ;  
  }
  
  fVertex = new Double_t*[fNMixedEvent] ; 
  for (Int_t i = 0; i < fNMixedEvent; i++) 
  {
    fVertex[i] = new Double_t[3] ; 
    fVertex[i][0] = 0.0 ; 
    fVertex[i][1] = 0.0 ; 
    fVertex[i][2] = 0.0 ; 
  }          
  
  if (!fMixedEvent) 
  { //Single event analysis
    if(fDataType!=kMC)
    {
      
      if(fInputEvent->GetPrimaryVertex())
      {
        fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex[0]); 
      }
      else 
      {
        printf("AliCaloTrackReader::FillVertexArray() - NULL primary vertex\n");
        fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
      }//Primary vertex pointer do not exist
      
    } else
    {//MC read event 
      fVertex[0][0]=0.;   fVertex[0][1]=0.;   fVertex[0][2]=0.;
    }
    
    if(fDebug > 1)
      printf("AliCaloTrackReader::FillVertexArray() - Single Event Vertex : %f,%f,%f\n",fVertex[0][0],fVertex[0][1],fVertex[0][2]);
    
  } else 
  { // MultiEvent analysis
    for (Int_t iev = 0; iev < fNMixedEvent; iev++) 
    {
      if (fMixedEvent->GetVertexOfEvent(iev))
        fMixedEvent->GetVertexOfEvent(iev)->GetXYZ(fVertex[iev]);
      else
      { // no vertex found !!!!
        AliWarning("No vertex found");
      }
      
      if(fDebug > 1)
        printf("AliCaloTrackReader::FillVertexArray() - Multi Event %d Vertex : %f,%f,%f\n",iev,fVertex[iev][0],fVertex[iev][1],fVertex[iev][2]);
      
    }
  }
  
}

//_____________________________________
void AliCaloTrackReader::FillInputCTS() 
{
  //Return array with Central Tracking System (CTS) tracks
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS()\n");
  
  Double_t pTrack[3] = {0,0,0};
  
  Int_t nTracks = fInputEvent->GetNumberOfTracks() ;
  fTrackMult    = 0;
  Int_t nstatus = 0;
  Double_t bz   = GetInputEvent()->GetMagneticField();

  for(Int_t i = 0; i < 19; i++)
  {
    fTrackBCEvent   [i] = 0;
    fTrackBCEventCut[i] = 0;
  }
  
  Bool_t   bc0  = kFALSE;
  if(fRecalculateVertexBC) fVertexBC=AliVTrack::kTOFBCNA;
  
  for (Int_t itrack =  0; itrack <  nTracks; itrack++)
  {////////////// track loop
    AliVTrack * track = (AliVTrack*)fInputEvent->GetTrack(itrack) ; // retrieve track from esd
    
    //Select tracks under certain conditions, TPCrefit, ITSrefit ... check the set bits
    ULong_t status = track->GetStatus();

    if (fTrackStatus && !((status & fTrackStatus) == fTrackStatus))
      continue ;
    
    nstatus++;
    
    Float_t dcaTPC =-999;
    
    if     (fDataType==kESD)
    {
      AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
      
      if(esdTrack)
      {
        if(fESDtrackCuts->AcceptTrack(esdTrack))
        {
          track->GetPxPyPz(pTrack) ;

          if(fConstrainTrack)
          {
            if(esdTrack->GetConstrainedParam())
            {
              const AliExternalTrackParam* constrainParam = esdTrack->GetConstrainedParam();
              esdTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
              esdTrack->GetConstrainedPxPyPz(pTrack);
            }
            else continue;
            
          } // use constrained tracks
          
          if(fSelectSPDHitTracks)
          {//Not much sense to use with TPC only or Hybrid tracks
            if(!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1)) continue ;
          }
        }
        // Complementary track to global : Hybrids (make sure that the previous selection is for Global)
        else  if(fESDtrackComplementaryCuts && fESDtrackComplementaryCuts->AcceptTrack(esdTrack))
        {
          // constrain the track
          if(esdTrack->GetConstrainedParam())
          {
            esdTrack->Set(esdTrack->GetConstrainedParam()->GetX(),esdTrack->GetConstrainedParam()->GetAlpha(),esdTrack->GetConstrainedParam()->GetParameter(),esdTrack->GetConstrainedParam()->GetCovariance());
          
            track->GetPxPyPz(pTrack) ;

          }
          else continue;
        }
        else continue;
      }
    } // ESD
    else if(fDataType==kAOD)
    {
      AliAODTrack *aodtrack = dynamic_cast <AliAODTrack*>(track);
      
      if(aodtrack)
      {
       if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS():AOD track type: %d (primary %d), hybrid? %d \n",
                              aodtrack->GetType(),AliAODTrack::kPrimary,
                              aodtrack->IsHybridGlobalConstrainedGlobal());
        
        
        if (fSelectHybridTracks)
        {
          if (!aodtrack->IsHybridGlobalConstrainedGlobal())       continue ;
        }
        else 
        {
          if ( aodtrack->TestFilterBit(fTrackFilterMask)==kFALSE) continue ;
        }
        
        if(fSelectSPDHitTracks)
        {//Not much sense to use with TPC only or Hybrid tracks
          if(!aodtrack->HasPointOnITSLayer(0) && !aodtrack->HasPointOnITSLayer(1)) continue ;
        }
        
        if (aodtrack->GetType()!= AliAODTrack::kPrimary)          continue ;
        
        if (fDebug > 2 ) printf("AliCaloTrackReader::FillInputCTS(): \t accepted track! \n");
        
        //In case of AODs, TPC tracks cannot be propagated back to primary vertex,
        // info stored here
        dcaTPC = aodtrack->DCA();
        
        track->GetPxPyPz(pTrack) ;
        
      } // aod track exists
      else continue ;
      
    } // AOD
    
    TLorentzVector momentum(pTrack[0],pTrack[1],pTrack[2],0);
    
    Bool_t okTOF  = ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) ;
    Double_t tof  = -1000;
    Int_t trackBC = -1000 ;
    
    if(okTOF)
    {
      trackBC = track->GetTOFBunchCrossing(bz);
      SetTrackEventBC(trackBC+9);

      tof = track->GetTOFsignal()*1e-3;
    }
                
    if(fUseTrackDCACut)
    {
      //normal way to get the dca, cut on dca_xy
      if(dcaTPC==-999)
      {
        Double_t dca[2]   = {1e6,1e6};
        Double_t covar[3] = {1e6,1e6,1e6};
        Bool_t okDCA = track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),bz,100.,dca,covar);
        if( okDCA) okDCA = AcceptDCA(momentum.Pt(),dca[0]);
        if(!okDCA)
        {
          //printf("AliCaloTrackReader::FillInputCTS() - Reject track pt %2.2f, dca_xy %2.4f, BC %d\n",momentum.Pt(),dca[0],trackBC);
          continue ;
        }
      }
    }// DCA cuts
    
    if(okTOF)
    {
      //SetTrackEventBCcut(bc);
      SetTrackEventBCcut(trackBC+9);
      
      //After selecting tracks with small DCA, pointing to vertex, set vertex BC depeding on tracks BC
      if(fRecalculateVertexBC)
      {
        if     (trackBC !=0 && trackBC != AliVTrack::kTOFBCNA) fVertexBC = trackBC;
        else if(trackBC == 0)                                  bc0 = kTRUE;
      }

      //In any case, the time should to be larger than the fixed window ...
      if( fUseTrackTimeCut && (trackBC!=0 || tof < fTrackTimeCutMin  || tof > fTrackTimeCutMax) )
      {
        //printf("Remove track time %f and bc = %d\n",tof,trackBC);
        continue ;
      }
      //else printf("Accept track time %f and bc = %d\n",tof,trackBC);
      
    }
    
    //Count the tracks in eta < 0.9
    //printf("Eta %f cut  %f\n",TMath::Abs(track->Eta()),fTrackMultEtaCut);
    if(TMath::Abs(track->Eta())< fTrackMultEtaCut) fTrackMult++;
    
    if(fCTSPtMin > momentum.Pt() || fCTSPtMax < momentum.Pt()) continue ;
    
    if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"CTS")) continue;
    
    if(fDebug > 2 && momentum.Pt() > 0.1)
      printf("AliCaloTrackReader::FillInputCTS() - Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
             momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
    
    if (fMixedEvent)  track->SetID(itrack);

    fCTSTracks->Add(track);
    
  }// track loop
	
  if(fVertexBC ==0 || fVertexBC == AliVTrack::kTOFBCNA)
  {
    if( bc0 ) fVertexBC = 0 ;
    else      fVertexBC = AliVTrack::kTOFBCNA ;
  }
  
  if(fDebug > 1)
    printf("AliCaloTrackReader::FillInputCTS()   - aod entries %d, input tracks %d, pass status %d, multipliticy %d\n", fCTSTracks->GetEntriesFast(), nTracks, nstatus, fTrackMult);//fCTSTracksNormalInputEntries);
  
}

//__________________________________________________________________
void AliCaloTrackReader::FillInputEMCALAlgorithm(AliVCluster * clus, 
                                                 const Int_t iclus) 
{
  //Fill the EMCAL data in the array, do it
    
  Int_t vindex = 0 ;  
  if (fMixedEvent) 
    vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
  
  //Reject clusters with bad channels, close to borders and exotic;
  if(!GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,GetCaloUtils()->GetEMCALGeometry(),GetEMCALCells(),fInputEvent->GetBunchCrossNumber())) return;
  
  //Mask all cells in collumns facing ALICE thick material if requested
  if(GetCaloUtils()->GetNMaskCellColumns())
  {
    Int_t absId   = -1;
    Int_t iSupMod = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Bool_t shared = kFALSE;
    GetCaloUtils()->GetEMCALRecoUtils()->GetMaxEnergyCell(GetCaloUtils()->GetEMCALGeometry(), GetEMCALCells(),clus,absId,iSupMod,ieta,iphi,shared);
    if(GetCaloUtils()->MaskFrameCluster(iSupMod, ieta)) return;
  }
  
  if(fSelectEmbeddedClusters)
  {
    if(clus->GetNLabels()==0 || clus->GetLabel() < 0) return;
    //else printf("Embedded cluster,  %d, n label %d label %d  \n",iclus,clus->GetNLabels(),clus->GetLabel());
  }
  
  //Float_t pos[3];
  //clus->GetPosition(pos);
  //printf("Before Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
  
  if(fRecalculateClusters)
  {
    //Recalibrate the cluster energy 
    if(GetCaloUtils()->IsRecalibrationOn())
    {
      Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, GetEMCALCells());
      
      clus->SetE(energy);
      //printf("Recalibrated Energy %f\n",clus->E());  
      
      GetCaloUtils()->RecalculateClusterShowerShapeParameters(GetEMCALCells(),clus);
      GetCaloUtils()->RecalculateClusterPID(clus);
      
    } // recalculate E
    
    //Recalculate distance to bad channels, if new list of bad channels provided
    GetCaloUtils()->RecalculateClusterDistanceToBadChannel(GetEMCALCells(),clus);
    
    //Recalculate cluster position
    if(GetCaloUtils()->IsRecalculationOfClusterPositionOn())
    {
      GetCaloUtils()->RecalculateClusterPosition(GetEMCALCells(),clus); 
      //clus->GetPosition(pos);
      //printf("After  Corrections: e %f, x %f, y %f, z %f\n",clus->E(),pos[0],pos[1],pos[2]);
    }
    
    // Recalculate TOF
    if(GetCaloUtils()->GetEMCALRecoUtils()->IsTimeRecalibrationOn()) 
    {
      Double_t tof      = clus->GetTOF();
      Float_t  frac     =-1;
      Int_t    absIdMax = GetCaloUtils()->GetMaxEnergyCell(fEMCALCells, clus,frac);
      
      if(fDataType==AliCaloTrackReader::kESD)
      { 
        tof = fEMCALCells->GetCellTime(absIdMax);
      }
      
      GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
      
      clus->SetTOF(tof);
      
    }// Time recalibration    
  }
  
  //Correct non linearity
  if(GetCaloUtils()->IsCorrectionOfClusterEnergyOn())
  {
    GetCaloUtils()->CorrectClusterEnergy(clus) ;
    //printf("Linearity Corrected Energy %f\n",clus->E());  
    
    //In case of MC analysis, to match resolution/calibration in real data
    Float_t rdmEnergy = GetCaloUtils()->GetEMCALRecoUtils()->SmearClusterEnergy(clus);
    // printf("\t Energy %f, smeared %f\n", clus->E(),rdmEnergy);
    clus->SetE(rdmEnergy);
  }
  
  Double_t tof = clus->GetTOF()*1e9;

  Int_t bc = TMath::Nint(tof/50) + 9;
  //printf("tof %2.2f, bc+5=%d\n",tof,bc);
  
  SetEMCalEventBC(bc);
  
  if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) return ;

  TLorentzVector momentum ;
  
  clus->GetMomentum(momentum, fVertex[vindex]);
  
  if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) return ;
  
  SetEMCalEventBCcut(bc);

  if(!IsInTimeWindow(tof,clus->E()))
  {
    fNPileUpClusters++ ;
    if(fUseEMCALTimeCut) return ;
  }
  else
    fNNonPileUpClusters++;
  
  if(fDebug > 2 && momentum.E() > 0.1) 
    printf("AliCaloTrackReader::FillInputEMCAL() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
           momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
  
  if (fMixedEvent) 
    clus->SetID(iclus) ; 
  
  //Correct MC label for AODs
  
  fEMCALClusters->Add(clus);
  
}

//_______________________________________
void AliCaloTrackReader::FillInputEMCAL() 
{
  //Return array with EMCAL clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputEMCAL()\n");
  
  // First recalibrate cells, time or energy
  //  if(GetCaloUtils()->IsRecalibrationOn())
  //    GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCells(GetCaloUtils()->GetEMCALGeometry(), 
  //                                                          GetEMCALCells(), 
  //                                                          fInputEvent->GetBunchCrossNumber());
  
  fNPileUpClusters    = 0; // Init counter
  fNNonPileUpClusters = 0; // Init counter
  for(Int_t i = 0; i < 19; i++)
  {
    fEMCalBCEvent   [i] = 0;
    fEMCalBCEventCut[i] = 0;
  }
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  if(fEMCALClustersListName=="")
  {
    Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++) 
    {
      AliVCluster * clus = 0;
      if ( (clus = fInputEvent->GetCaloCluster(iclus)) ) 
      {
        if (IsEMCALCluster(clus))
        {          
          FillInputEMCALAlgorithm(clus, iclus);
        }//EMCAL cluster
      }// cluster exists
    }// cluster loop
    
    //Recalculate track matching
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent);
    
  }//Get the clusters from the input event
  else
  {
    TClonesArray * clusterList = 0x0; 
    
    if      (fInputEvent->FindListObject(fEMCALClustersListName))
    {
      clusterList = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(fEMCALClustersListName));
    }
    else if(fOutputEvent)
    {
      clusterList = dynamic_cast<TClonesArray*> (fOutputEvent->FindListObject(fEMCALClustersListName));
    }
    
    if(!clusterList)
    {
        printf("AliCaloTrackReader::FillInputEMCAL() - Wrong name of list with clusters?  <%s>\n",fEMCALClustersListName.Data());
        return;
    }
    
    Int_t nclusters = clusterList->GetEntriesFast();
    for (Int_t iclus =  0; iclus <  nclusters; iclus++) 
    {
      AliVCluster * clus = dynamic_cast<AliVCluster*> (clusterList->At(iclus));
      //printf("E %f\n",clus->E());
      if (clus) FillInputEMCALAlgorithm(clus, iclus);
      else printf("AliCaloTrackReader::FillInputEMCAL() - Null cluster in list!\n");
    }// cluster loop
    
    // Recalculate the pile-up time, in case long time clusters removed during clusterization
    //printf("Input event INIT : Pile-up clusters %d, NO pile-up %d\n",fNPileUpClusters,fNNonPileUpClusters);

    fNPileUpClusters    = 0; // Init counter
    fNNonPileUpClusters = 0; // Init counter
    for(Int_t i = 0; i < 19; i++)
    {
      fEMCalBCEvent   [i] = 0;
      fEMCalBCEventCut[i] = 0;
    }
    
    for (Int_t iclus =  0; iclus < fInputEvent->GetNumberOfCaloClusters(); iclus++)
    {
      AliVCluster * clus = 0;
      
      if ( (clus = fInputEvent->GetCaloCluster(iclus)) )
      {
        if (IsEMCALCluster(clus))
        {
          
          Float_t  frac     =-1;
          Int_t    absIdMax = GetCaloUtils()->GetMaxEnergyCell(fEMCALCells, clus,frac);
          Double_t tof = clus->GetTOF();
          GetCaloUtils()->GetEMCALRecoUtils()->RecalibrateCellTime(absIdMax,fInputEvent->GetBunchCrossNumber(),tof);
          tof*=1e9;

          //printf("Input event cluster : AbsIdMax %d, E %2.2f, time %2.2f \n", absIdMax,clus->E(),tof);

          //Reject clusters with bad channels, close to borders and exotic;
          if(!GetCaloUtils()->GetEMCALRecoUtils()->IsGoodCluster(clus,GetCaloUtils()->GetEMCALGeometry(),GetEMCALCells(),fInputEvent->GetBunchCrossNumber()))  continue;

          Int_t bc = TMath::Nint(tof/50) + 9;
          SetEMCalEventBC(bc);
          
          if(fEMCALPtMin > clus->E() || fEMCALPtMax < clus->E()) continue ;

          TLorentzVector momentum ;
          
          clus->GetMomentum(momentum, fVertex[0]);
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) return ;

          SetEMCalEventBCcut(bc);

          if(!IsInTimeWindow(tof,clus->E()))
            fNPileUpClusters++ ;
          else
            fNNonPileUpClusters++;
          
        }
      }
    }
    
    //printf("Input event : Pile-up clusters %d, NO pile-up %d\n",fNPileUpClusters,fNNonPileUpClusters);
    
    // Recalculate track matching, not necessary if already done in the reclusterization task.
    // in case it was not done ...
    GetCaloUtils()->RecalculateClusterTrackMatching(fInputEvent,clusterList);
    
  }
    
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputEMCAL() - aod entries %d, n pile-up clusters %d, n non pile-up %d \n",  fEMCALClusters->GetEntriesFast(),fNPileUpClusters,fNNonPileUpClusters);
  
}

//______________________________________
void AliCaloTrackReader::FillInputPHOS() 
{
  //Return array with PHOS clusters in aod format
  
  if(fDebug > 2 ) printf("AliCaloTrackReader::FillInputPHOS()\n");
  
  //Loop to select clusters in fiducial cut and fill container with aodClusters
  Int_t nclusters = fInputEvent->GetNumberOfCaloClusters();
  for (Int_t iclus = 0; iclus < nclusters; iclus++) 
  {
    AliVCluster * clus = 0;
    if ( (clus = fInputEvent->GetCaloCluster(iclus)) ) 
    {
      if (IsPHOSCluster(clus))
      {
        //Check if the cluster contains any bad channel and if close to calorimeter borders
        Int_t vindex = 0 ;  
        if (fMixedEvent) 
          vindex = fMixedEvent->EventIndexForCaloCluster(iclus);
        if( GetCaloUtils()->ClusterContainsBadChannel("PHOS",clus->GetCellsAbsId(), clus->GetNCells())) 
          continue;
        if(!GetCaloUtils()->CheckCellFiducialRegion(clus, fInputEvent->GetPHOSCells(), fInputEvent, vindex)) 
          continue;
        
        if(fRecalculateClusters)
        {
          //Recalibrate the cluster energy 
          if(GetCaloUtils()->IsRecalibrationOn()) 
          {
            Float_t energy = GetCaloUtils()->RecalibrateClusterEnergy(clus, (AliAODCaloCells*)GetPHOSCells());
            clus->SetE(energy);
          }
        }
        
        TLorentzVector momentum ;
        
        clus->GetMomentum(momentum, fVertex[vindex]);      
        
        if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(momentum,"PHOS")) continue;
        
        if(fPHOSPtMin > momentum.E() || fPHOSPtMax < momentum.E())          continue;
        
        if(fDebug > 2 && momentum.E() > 0.1) 
          printf("AliCaloTrackReader::FillInputPHOS() - Selected clusters E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
                 momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());        
        
        
        if (fMixedEvent) 
        {
          clus->SetID(iclus) ; 
        }              
        
        fPHOSClusters->Add(clus);	
        
      }//PHOS cluster
    }//cluster exists
  }//esd cluster loop
  
  if(fDebug > 1) printf("AliCaloTrackReader::FillInputPHOS()  - aod entries %d\n",  fPHOSClusters->GetEntriesFast());
  
}

//____________________________________________
void AliCaloTrackReader::FillInputEMCALCells() 
{
  //Return array with EMCAL cells in aod format
  
  fEMCALCells = fInputEvent->GetEMCALCells(); 
  
}

//___________________________________________
void AliCaloTrackReader::FillInputPHOSCells() 
{
  //Return array with PHOS cells in aod format
  
  fPHOSCells = fInputEvent->GetPHOSCells(); 
  
}

//_______________________________________
void AliCaloTrackReader::FillInputVZERO()
{
  //Fill VZERO information in data member, add all the channels information.
  AliVVZERO* v0 = fInputEvent->GetVZEROData();
  //printf("Init V0: ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  
  if (v0) 
  {
    AliESDVZERO* esdV0 = dynamic_cast<AliESDVZERO*> (v0);
    for (Int_t i = 0; i < 32; i++)
    {
      if(esdV0)
      {//Only available in ESDs
        fV0ADC[0] += (Int_t)esdV0->GetAdcV0C(i);
        fV0ADC[1] += (Int_t)esdV0->GetAdcV0A(i);
      }
      
      fV0Mul[0] += (Int_t)v0->GetMultiplicityV0C(i);
      fV0Mul[1] += (Int_t)v0->GetMultiplicityV0A(i);
    }
    if(fDebug > 0)
      printf("V0: ADC (%d,%d), Multiplicity (%d,%d) \n",fV0ADC[0],fV0ADC[1],fV0Mul[0],fV0Mul[1]);
  }
  else
  {
    if(fDebug > 0)
      printf("Cannot retrieve V0 ESD! Run w/ null V0 charges\n ");
  }
}


//___________________________________________________________________
Bool_t AliCaloTrackReader::IsEMCALCluster(AliVCluster* cluster) const 
{
  // Check if it is a cluster from EMCAL. For old AODs cluster type has
  // different number and need to patch here
  
  if(fDataType==kAOD && fOldAOD)
  {
    if (cluster->GetType() == 2) return kTRUE;
    else                         return kFALSE;
  }
  else 
  {
    return cluster->IsEMCAL();
  }
  
}

//___________________________________________________________________
Bool_t AliCaloTrackReader::IsPHOSCluster(AliVCluster * cluster) const 
{
  //Check if it is a cluster from PHOS.For old AODs cluster type has
  // different number and need to patch here
  
  if(fDataType==kAOD && fOldAOD)
  {
    Int_t type = cluster->GetType();
    if (type == 0 || type == 1) return kTRUE;
    else                        return kFALSE;
  }
  else 
  {
    return cluster->IsPHOS();
  }
  
}

//________________________________________________
Bool_t AliCaloTrackReader::CheckForPrimaryVertex()
{
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //Only for ESDs ...
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (fInputEvent);
  if(!event) return kTRUE;
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0) 
  {
    return kTRUE;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1) 
  {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0) 
    {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1)
    {
      //      cout<<"bad vertex type::"<< event->GetPrimaryVertex()->GetName() << endl;
      return kFALSE;
    }
  }
  
  return kFALSE;  
  
}

//____________________________________________________________
void  AliCaloTrackReader::SetTrackCuts(AliESDtrackCuts * cuts) 
{ 
  // Set Track cuts
  
  if(fESDtrackCuts) delete fESDtrackCuts ;
  
  fESDtrackCuts = cuts ; 
  
}		  

//_________________________________________________________________________
void  AliCaloTrackReader::SetTrackComplementaryCuts(AliESDtrackCuts * cuts)
{
  // Set Track cuts for complementary tracks (hybrids)
  
  if(fESDtrackComplementaryCuts) delete fESDtrackComplementaryCuts ;
  
  fESDtrackComplementaryCuts = cuts ;
  
}


