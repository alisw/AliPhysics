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
// This analysis provides a new list of clusters to be used in other analysis
//
// Author: Gustavo Conesa Balbastre,
//         Adapted from analysis class from Deepa Thomas
//
//
//_________________________________________________________________________

// --- Root ---
#include "TString.h"
#include "TRefArray.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"

// --- AliRoot Analysis Steering
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliGeomManager.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"

// --- EMCAL
#include "AliEMCALRecParam.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALRecoUtils.h"

#include "AliAnalysisTaskEMCALClusterize.h"

ClassImp(AliAnalysisTaskEMCALClusterize)

//________________________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize(const char *name) 
  : AliAnalysisTaskSE(name)
  , fGeom(0), fGeomName("EMCAL_FIRSTYEARV1"),  fGeomMatrixSet(kFALSE), fLoadGeomMatrices(kFALSE)
  , fCalibData(0),       fPedestalData(0),     fOCDBpath("raw://")
  , fDigitsArr(0),       fClusterArr(0),       fCaloClusterArr(0)
  , fRecParam(0),        fClusterizer(0),      fUnfolder(0),           fJustUnfold(kFALSE) 
  , fOutputAODBranch(0), fOutputAODBranchName("newEMCALClusters")
  , fFillAODFile(kTRUE), fFillAODHeader(0),    fFillAODCaloCells(0)
  , fRun(-1),            fRecoUtils(0),        fConfigName("")
  , fCellLabels(),       fCellSecondLabels(),  fMaxEvent(1000000000)
  
  {
  //ctor
  for(Int_t i = 0; i < 10;    i++)  fGeomMatrix[i] =  0;
    for(Int_t j = 0; j < 24*48*11; j++)  {
      fCellLabels[j]       = -1;
      fCellSecondLabels[j] = -1;
    }  
  fDigitsArr       = new TClonesArray("AliEMCALDigit",200);
  fClusterArr      = new TObjArray(100);
  fCaloClusterArr  = new TObjArray(100);
  fRecParam        = new AliEMCALRecParam;
  fBranchNames     = "ESD:AliESDHeader.,EMCALCells.";
  fRecoUtils       = new AliEMCALRecoUtils();
}

//________________________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEMCALClusterize")
  , fGeom(0), fGeomName("EMCAL_FIRSTYEARV1"),   fGeomMatrixSet(kFALSE), fLoadGeomMatrices(kFALSE)
  , fCalibData(0),        fPedestalData(0),     fOCDBpath("raw://")
  , fDigitsArr(0),        fClusterArr(0),       fCaloClusterArr(0)
  , fRecParam(0),         fClusterizer(0),      fUnfolder(0),           fJustUnfold(kFALSE)
  , fOutputAODBranch(0),  fOutputAODBranchName("newEMCALClusters")
  , fFillAODFile(kFALSE), fFillAODHeader(0),    fFillAODCaloCells(0)
  , fRun(-1),             fRecoUtils(0),        fConfigName("") 
  , fCellLabels(),        fCellSecondLabels(),  fMaxEvent(1000000000)
{
  // Constructor
  for(Int_t i = 0; i < 10;    i++)  fGeomMatrix[i] =  0;
  for(Int_t j = 0; j < 24*48*11; j++)  {
    fCellLabels[j]       = -1;
    fCellSecondLabels[j] = -1;
  }
  fDigitsArr       = new TClonesArray("AliEMCALDigit",200);
  fClusterArr      = new TObjArray(100);
  fCaloClusterArr  = new TObjArray(100);
  fRecParam        = new AliEMCALRecParam;
  fBranchNames     = "ESD:AliESDHeader.,EMCALCells.";
  fRecoUtils       = new AliEMCALRecoUtils();
}


//________________________________________________________________________
AliAnalysisTaskEMCALClusterize::~AliAnalysisTaskEMCALClusterize()
{
  //dtor 
  
  if (fDigitsArr){
    fDigitsArr->Clear("C");
    delete fDigitsArr; 
  }
  
  if (fClusterArr){
    fClusterArr->Delete();
    delete fClusterArr;
  }
  
  if (fCaloClusterArr){
    fCaloClusterArr->Delete();
    delete fCaloClusterArr; 
  }

  if(fClusterizer) delete fClusterizer;
  if(fUnfolder)    delete fUnfolder;   
  if(fRecoUtils)   delete fRecoUtils;

}

//-------------------------------------------------------------------
void AliAnalysisTaskEMCALClusterize::Init()
{
  //Init analysis with configuration macro if available
  
  if(gROOT->LoadMacro(fConfigName) >=0){
    printf("Configure analysis with %s\n",fConfigName.Data());
    AliAnalysisTaskEMCALClusterize *clus = (AliAnalysisTaskEMCALClusterize*)gInterpreter->ProcessLine("ConfigEMCALClusterize()");
    fGeomName         = clus->fGeomName; 
    fLoadGeomMatrices = clus->fLoadGeomMatrices;
    fOCDBpath         = clus->fOCDBpath;
    fRecParam         = clus->fRecParam;
    fJustUnfold       = clus->fJustUnfold;
    fFillAODFile      = clus->fFillAODFile;
    fRecoUtils        = clus->fRecoUtils; 
    fConfigName       = clus->fConfigName;
    fOutputAODBranchName = clus->fOutputAODBranchName;
    for(Int_t i = 0; i < 10; i++) fGeomMatrix[i] = clus->fGeomMatrix[i] ;

  }

}  

//-------------------------------------------------------------------
void AliAnalysisTaskEMCALClusterize::UserCreateOutputObjects()
{
  // Init geometry, create list of output clusters

  fGeom =  AliEMCALGeometry::GetInstance(fGeomName) ;	
  if(fOutputAODBranchName.Length()!=0){
    fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);
    fOutputAODBranch->SetName(fOutputAODBranchName);
    AddAODBranch("TClonesArray", &fOutputAODBranch);
  }
  else {
    AliFatal("fOutputAODBranchName not set\n");
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODCaloCells() {
  // Put calo cells in standard branch  
  AliVEvent * event = InputEvent();
  AliVCaloCells &eventEMcells = *(event->GetEMCALCells());
  Int_t nEMcell = eventEMcells.GetNumberOfCells() ;
  
  AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
  aodEMcells.CreateContainer(nEMcell);
  aodEMcells.SetType(AliVCaloCells::kEMCALCell);
  Double_t calibFactor = 1.;   
  for (Int_t iCell = 0; iCell < nEMcell; iCell++) { 
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fGeom->GetCellIndex(eventEMcells.GetCellNumber(iCell),imod,iTower,iIphi,iIeta); 
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
    
    if(fRecoUtils->IsRecalibrationOn()){ 
      calibFactor = fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
    }
    
    if(!fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi)){ //Channel is not declared as bad
      aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),eventEMcells.GetAmplitude(iCell)*calibFactor);
      //printf("GOOD channel\n");
    }
    else {
      aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),0);
      //printf("BAD channel\n");
    }
  }
  aodEMcells.Sort();
  
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODHeader() {
  //Put event header information in standard AOD branch
  
  AliVEvent* event      = InputEvent();
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);

  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  AliAODHeader* header = AODEvent()->GetHeader();
  header->SetRunNumber(event->GetRunNumber());
  
  if(esdevent){
    TTree* tree = fInputHandler->GetTree();
    if (tree) {
      TFile* file = tree->GetCurrentFile();
      if (file) header->SetESDFileName(file->GetName());
    }
  }
  else if (aodevent) header->SetESDFileName(aodevent->GetHeader()->GetESDFileName());
  
  header->SetBunchCrossNumber(event->GetBunchCrossNumber());
  header->SetOrbitNumber(event->GetOrbitNumber());
  header->SetPeriodNumber(event->GetPeriodNumber());
  header->SetEventType(event->GetEventType());
  
  //Centrality
  if(event->GetCentrality()){
    header->SetCentrality(new AliCentrality(*(event->GetCentrality())));
  }
  else{
    header->SetCentrality(0);
  }
  
  //Trigger  
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection
  if      (esdevent) header->SetFiredTriggerClasses(esdevent->GetFiredTriggerClasses());
  else if (aodevent) header->SetFiredTriggerClasses(aodevent->GetFiredTriggerClasses());
  header->SetTriggerMask(event->GetTriggerMask()); 
  header->SetTriggerCluster(event->GetTriggerCluster());
  if(esdevent){
    header->SetL0TriggerInputs(esdevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(esdevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(esdevent->GetHeader()->GetL2TriggerInputs());    
  }
  else if (aodevent){
    header->SetL0TriggerInputs(aodevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(aodevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(aodevent->GetHeader()->GetL2TriggerInputs());    
  }
  
  header->SetMagneticField(event->GetMagneticField());
  //header->SetMuonMagFieldScale(esdevent->GetCurrentDip()/6000.); 
  
  header->SetZDCN1Energy(event->GetZDCN1Energy());
  header->SetZDCP1Energy(event->GetZDCP1Energy());
  header->SetZDCN2Energy(event->GetZDCN2Energy());
  header->SetZDCP2Energy(event->GetZDCP2Energy());
  header->SetZDCEMEnergy(event->GetZDCEMEnergy(0),event->GetZDCEMEnergy(1));
  
  Float_t diamxy[2]={event->GetDiamondX(),event->GetDiamondY()};
  Float_t diamcov[3];
  event->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  if      (esdevent) header->SetDiamondZ(esdevent->GetDiamondZ(),esdevent->GetSigma2DiamondZ());
  else if (aodevent) header->SetDiamondZ(aodevent->GetDiamondZ(),aodevent->GetSigma2DiamondZ());
  
  //
  //
  Int_t nVertices = 1 ;/* = prim. vtx*/;
  Int_t nCaloClus = event->GetNumberOfCaloClusters();
  
  AODEvent()->ResetStd(0, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  Int_t jVertices=0;
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  event->GetPrimaryVertex()->GetXYZ(pos);
  Float_t chi = 0;
  if      (esdevent){
    esdevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = esdevent->GetPrimaryVertex()->GetChi2toNDF();
  }
  else if (aodevent){
    aodevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = aodevent->GetPrimaryVertex()->GetChi2perNDF();//Different from ESD?
  }
  
  AliAODVertex * primary = new(vertices[jVertices++])
  AliAODVertex(pos, covVtx, chi, NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(event->GetPrimaryVertex()->GetName());
  primary->SetTitle(event->GetPrimaryVertex()->GetTitle());
  
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  Int_t eventN = Entry();
  if(aodIH) 
    eventN = aodIH->GetReadEntry();
  //printf("Clusterizer --- Event %d-- \n",eventN);

  if (eventN > fMaxEvent) return;

  //Remove the contents of output list set in the previous event 
  fOutputAODBranch->Clear("C");
  
//  if(fOutputAODBranchName.Length()!=0){
//    //Remove the contents of output list set in the previous event 
//    fOutputAODBranch->Clear("C");
//  }
//  else{ 
//    //Reset and put new clusters in standard branch, also cells
//    AODEvent()->ResetStd(0, 1, 0, 0, 0, AODEvent()->GetNumberOfCaloClusters(), 0, 0);
//    fOutputAODBranch = AODEvent()->GetCaloClusters();// Put new clusters in standard branch
//  }
  
  //Magic line to write events to AOD file
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(fFillAODFile);
  LoadBranches();
  
  AccessOCDB();
  
  //Get the event
  AliVEvent   * event    = 0;
  AliESDEvent * esdevent = 0;
    
  //Fill output event with headerø

  //Check if input event are embedded events
  //If so, take output event
  if (aodIH && aodIH->GetMergeEvents()) {
    //printf("AliAnalysisTaskEMCALClusterize::UserExec() - Use embedded events\øn");
    event  = AODEvent();
    
    if(!aodIH->GetMergeEMCALCells()) 
      AliFatal("Events merged but not EMCAL cells, check analysis settings!");
    
//    printf("InputEvent  N Clusters %d, N Cells %d\n",InputEvent()->GetNumberOfCaloClusters(),
//           InputEvent()->GetEMCALCells()->GetNumberOfCells());
//    printf("MergedEvent  N Clusters %d, N Cells %d\n",aodIH->GetEventToMerge()->GetNumberOfCaloClusters(),
//           aodIH->GetEventToMerge()->GetEMCALCells()->GetNumberOfCells());
//    printf("OutputEvent N Clusters %d, N Cells %d\n", AODEvent()->GetNumberOfCaloClusters(),
//           AODEvent()->GetEMCALCells()->GetNumberOfCells());
    
  }
  else {
    event =  InputEvent();
    esdevent = dynamic_cast<AliESDEvent*>(InputEvent());
    if(fFillAODCaloCells) FillAODCaloCells();   
    if(fFillAODHeader)    FillAODHeader();
  }
  
  if (!event) {
    Error("UserExec","Event not available");
    return;
  }

  //-------------------------------------------------------------------------------------
  //Set the geometry matrix, for the first event, skip the rest
  //-------------------------------------------------------------------------------------
  if(!fGeomMatrixSet){
    if(fLoadGeomMatrices){
      for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
        if(fGeomMatrix[mod]){
          if(DebugLevel() > 1) 
            fGeomMatrix[mod]->Print();
          fGeom->SetMisalMatrix(fGeomMatrix[mod],mod) ;  
        }
        fGeomMatrixSet=kTRUE;
      }//SM loop
    }//Load matrices
    else if(!gGeoManager){
      Info("UserExec","Get geo matrices from data");
      //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
      if(!strcmp(event->GetName(),"AliAODEvent")) {
        if(DebugLevel() > 1) 
          Warning("UserExec","Use ideal geometry, values geometry matrix not kept in AODs.");
      }//AOD
      else {	
        if(DebugLevel() > 1) 
          Info("UserExec","AliAnalysisTaskEMCALClusterize Load Misaligned matrices.");
        AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event) ;
        if(!esd) {
          Error("UserExec","This event does not contain ESDs?");
          return;
        }
        for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
          if(DebugLevel() > 1) 
            esd->GetEMCALMatrix(mod)->Print();
          if(esd->GetEMCALMatrix(mod)) fGeom->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
        } 
        fGeomMatrixSet=kTRUE;
      }//ESD
    }//Load matrices from Data 
    
    //Recover time dependent corrections, put then in recalibration histograms. Do it once
    fRecoUtils->SetTimeDependentCorrections(InputEvent()->GetRunNumber());

  }//first event

  //Get the list of cells needed for unfolding and reclustering
  AliVCaloCells *cells= event->GetEMCALCells();
  Int_t nClustersOrg = 0;
  //-------------------------------------------
  //---------Unfolding clusters----------------
  //-------------------------------------------
  if (fJustUnfold) {
    
    //Fill the array with the EMCAL clusters, copy them
    for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++)
    {
      AliVCluster *clus = event->GetCaloCluster(i);
      if(clus->IsEMCAL()){        
        //printf("Org Cluster %d, E %f\n",i, clus->E());
        
        //recalibrate/remove bad channels/etc if requested
        if(fRecoUtils->ClusterContainsBadChannel(fGeom,clus->GetCellsAbsId(), clus->GetNCells())){
          //printf("Remove cluster\n");
          continue;
        } 
        
        if(fRecoUtils->IsRecalibrationOn()){
          //printf("Energy before %f ",clus->E());
          fRecoUtils->RecalibrateClusterEnergy(fGeom, clus, cells);
          //printf("; after %f\n",clus->E());
        }
        //Cast to ESD or AOD, needed to create the cluster array
        AliESDCaloCluster * esdCluster = dynamic_cast<AliESDCaloCluster*> (clus);
        AliAODCaloCluster * aodCluster = dynamic_cast<AliAODCaloCluster*> (clus);
        if     (esdCluster){
          fCaloClusterArr->Add( new AliESDCaloCluster(*esdCluster) );   
        }//ESD
        else if(aodCluster){
          fCaloClusterArr->Add( new AliAODCaloCluster(*aodCluster) );   
        }//AOD
        else 
          Warning("UserExec()"," - Wrong CaloCluster type?");
        nClustersOrg++;
      }
    }
    
    //Do the unfolding
    fUnfolder->UnfoldClusters(fCaloClusterArr, cells);
    
    //CLEAN-UP
    fUnfolder->Clear();
    
    
    //printf("nClustersOrg %d\n",nClustersOrg);
  }
  //-------------------------------------------
  //---------- Recluster cells ----------------
  //-------------------------------------------
  
  else{
    //-------------------------------------------------------------------------------------
    //Transform CaloCells into Digits
    //-------------------------------------------------------------------------------------
    
    //In case of MC, first loop on the clusters and fill MC label to array
    //.....................................................................
    
//    for(Int_t j = 0; j < 24*48*11; j++)  {
//      if(fCellLabels[j]!=-1) printf("Not well initialized cell %d, label %d\n",j,fCellLabels[j]) ;
//    }
    
    Int_t nClusters = event->GetNumberOfCaloClusters();
    if(aodIH) 
      nClusters = aodIH->GetEventToMerge()->GetNumberOfCaloClusters(); //Get clusters directly from embedded signal
    for (Int_t i = 0; i < nClusters; i++)
    {
      AliVCluster *clus = 0;
      if(aodIH) clus = aodIH->GetEventToMerge()->GetCaloCluster(i); //Get clusters directly from embedded signal
      else      clus = event->GetCaloCluster(i);

      if(!clus) {
        printf("AliEMCALReclusterize::UserExec() - No Clusters\n");
        return;
      }
      
      if(clus->IsEMCAL()){   
        //printf("Cluster Signal %d, energy %f\n",clus->GetID(),clus->E());
        Int_t label = clus->GetLabel();
        Int_t label2 = -1 ;
        if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;
        UShort_t * index    = clus->GetCellsAbsId() ;
        for(Int_t icell=0; icell < clus->GetNCells(); icell++ ){
          fCellLabels[index[icell]]=label;
          fCellSecondLabels[index[icell]]=label2;
          //printf("1) absID %d, label %d\n",index[icell], fCellLabels[index[icell]]);
        }
        nClustersOrg++;
      }
    } 

    // Create digits 
    //.................
    Int_t idigit =  0;
    Int_t id     = -1;
    Float_t amp  = -1; 
    Float_t time = -1; 
    
    Double_t cellAmplitude = 0;
    Double_t cellTime      = 0;
    Short_t  cellNumber    = 0;
        
    TTree *digitsTree = new TTree("digitstree","digitstree");
    digitsTree->Branch("EMCAL","TClonesArray", &fDigitsArr, 32000);
    
    
    for (Int_t icell = 0; icell < cells->GetNumberOfCells(); icell++)
    {
      if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
        break;
      
      time = cellTime;
      amp  = cellAmplitude;
      id   = cellNumber;
      
      Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
      fGeom->GetCellIndex(id,imod,iTower,iIphi,iIeta); 
      fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
      
      //Do not include bad channels found in analysis?
      if( fRecoUtils->IsBadChannelsRemovalSwitchedOn() && 
          fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi)){
          fCellLabels[id]=-1; //reset the entry in the array for next event
        //printf("Remove channel %d\n",id);
        continue;
      }
             
      //Recalibrate?
      if(fRecoUtils->IsRecalibrationOn()){ 
        //printf("CalibFactor %f times %f for id %d\n",fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi),amp,id);
        amp *=fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
      }
           
      //Create the digit, put a fake primary deposited energy to trick the clusterizer when checking the most likely primary
      new((*fDigitsArr)[idigit]) AliEMCALDigit( fCellLabels[id], fCellLabels[id],id, amp, time,AliEMCALDigit::kHG,idigit, 0, 0, 1); 
      //if(fCellLabels[id]>=0)printf("2) Digit cell %d, label %d\n",id,fCellLabels[id]) ;
      //else                  printf("2) Digit cell %d, no label, amp %f \n",id,amp) ;
      fCellLabels[id]=-1; //reset the entry in the array for next event
      
      //AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->New(idigit);
      //digit->SetId(id);
      //digit->SetAmplitude(amp);
      //digit->SetTime(time);
      //digit->SetTimeR(time);
      //digit->SetIndexInList(idigit);
      //digit->SetType(AliEMCALDigit::kHG);
      
      //printf("Digit: Id %d, amp %f, time %e, index %d\n",id, amp,time,idigit);
      idigit++;
    }
    
    //Fill the tree with digits
    digitsTree->Fill();
    
    //-------------------------------------------------------------------------------------
    //Do the clusterization
    //-------------------------------------------------------------------------------------        
    TTree *clustersTree = new TTree("clustertree","clustertree");

    fClusterizer->SetInput(digitsTree);
    fClusterizer->SetOutput(clustersTree);
    fClusterizer->Digits2Clusters("");

    //-------------------------------------------------------------------------------------
    //Transform the recpoints into AliVClusters
    //-------------------------------------------------------------------------------------
    
    clustersTree->SetBranchStatus("*",0); //disable all branches
    clustersTree->SetBranchStatus("EMCALECARP",1); //Enable only the branch we need
    
    TBranch *branch = clustersTree->GetBranch("EMCALECARP");
    branch->SetAddress(&fClusterArr);
    branch->GetEntry(0);

    RecPoints2Clusters(fDigitsArr, fClusterArr, fCaloClusterArr);
    
    //---CLEAN UP-----
    fClusterizer->Clear();
    fDigitsArr  ->Clear("C");
    fClusterArr ->Delete(); // Do not Clear(), it leaks, why?

    clustersTree->Delete("all");
    digitsTree  ->Delete("all");
  }
  
  //Recalculate track-matching for the new clusters, only with ESDs
  if(esdevent)fRecoUtils->FindMatches(esdevent,fCaloClusterArr);

  
  //-------------------------------------------------------------------------------------
  //Put the new clusters in the AOD list
  //-------------------------------------------------------------------------------------
  
  Int_t kNumberOfCaloClusters   = fCaloClusterArr->GetEntries();
  //printf("New clusters %d, Org clusters %d\n",kNumberOfCaloClusters, nClustersOrg);
  for(Int_t i = 0; i < kNumberOfCaloClusters; i++){
    AliAODCaloCluster *newCluster = (AliAODCaloCluster *) fCaloClusterArr->At(i);
    //if(Entry()==0) Info("UserExec","newCluster E %f\n", newCluster->E());
    
    //Add matched track, if any, only with ESDs
    if(esdevent){
      Int_t trackIndex = fRecoUtils->GetMatchedTrackIndex(i);
      if(trackIndex >= 0){
        newCluster->AddTrackMatched(event->GetTrack(trackIndex));
        if(DebugLevel() > 1) 
          Info("UserExec","Matched Track index %d to new cluster %d \n",trackIndex,i);
      }
    }
    
    //In case of new bad channels, recalculate distance to bad channels
    if(fRecoUtils->IsBadChannelsRemovalSwitchedOn()){
      //printf("DistToBad before %f ",newCluster->GetDistanceToBadChannel());
      fRecoUtils->RecalculateClusterDistanceToBadChannel(fGeom, event->GetEMCALCells(), newCluster);
      //printf("; after %f \n",newCluster->GetDistanceToBadChannel());
    }
    
//    if(newCluster->GetNLabels()>0){
//      printf("3) MC: N labels %d, label %d ; ", newCluster->GetNLabels(), newCluster->GetLabel() );
//      UShort_t * newindex    = newCluster->GetCellsAbsId() ;
//      for(Int_t icell=0; icell < newCluster->GetNCells(); icell++ ){
//       printf("\t absID %d\n",newindex[icell]);
//     }
//    }
    
    //printf("Cluster %d, energy %f\n",newCluster->GetID(),newCluster->E());
    
    newCluster->SetID(i);
    new((*fOutputAODBranch)[i])  AliAODCaloCluster(*newCluster);
  }
  
  //if(fOutputAODBranchName.Length()!=0)
  fOutputAODBranch->Expand(kNumberOfCaloClusters); // resize TObjArray to 'remove' slots
  
  //---CLEAN UP-----
  fCaloClusterArr->Delete(); // Do not Clear(), it leaks, why?
}      

//_____________________________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::AccessOCDB()
{
  //Access to OCDB stuff

  AliVEvent * event = InputEvent();
  if (!event)
  {
    Warning("AccessOCDB","Event not available!!!");
    return kFALSE;
  }

  if (event->GetRunNumber()==fRun)
    return kTRUE;
  fRun = event->GetRunNumber();

  if(DebugLevel() > 1 )
    Info("AccessODCD()"," Begin");

  fGeom = AliEMCALGeometry::GetInstance(fGeomName);
  
  
  AliCDBManager *cdb = AliCDBManager::Instance();
    

  if (fOCDBpath.Length()){
    cdb->SetDefaultStorage(fOCDBpath.Data());
    Info("AccessOCDB","Default storage %s",fOCDBpath.Data());
  }
  
  cdb->SetRun(event->GetRunNumber());

  //
  // EMCAL from RAW OCDB
  if (fOCDBpath.Contains("alien:"))
  {
    cdb->SetSpecificStorage("EMCAL/Calib/Data","alien://Folder=/alice/data/2010/OCDB");
    cdb->SetSpecificStorage("EMCAL/Calib/Pedestals","alien://Folder=/alice/data/2010/OCDB");
  }

  TString path = cdb->GetDefaultStorage()->GetBaseFolder();
  
  // init parameters:
  
  //Get calibration parameters	
  if(!fCalibData)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
      AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
    
    if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
  }
  
  if(!fCalibData)
    AliFatal("Calibration parameters not found in CDB!");
    
  //Get calibration parameters	
  if(!fPedestalData)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
      AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals");
    
    if (entry) fPedestalData =  (AliCaloCalibPedestal*) entry->GetObject();
  }
    
  if(!fPedestalData)
    AliFatal("Dead map not found in CDB!");

  InitClusterization();
  
  return kTRUE;
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::InitClusterization()
{
  //Select clusterization/unfolding algorithm and set all the needed parameters
  
  if (fJustUnfold){
    // init the unfolding afterburner 
    delete fUnfolder;
    fUnfolder =  new AliEMCALAfterBurnerUF(fRecParam->GetW0(),fRecParam->GetLocMaxCut());
    return;
 }

  //First init the clusterizer
  delete fClusterizer;
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1 (fGeom, fCalibData, fPedestalData);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) 
    fClusterizer = new AliEMCALClusterizerNxN(fGeom, fCalibData, fPedestalData);
  else if(fRecParam->GetClusterizerFlag() > AliEMCALRecParam::kClusterizerNxN) {
   AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(fGeom, fCalibData, fPedestalData);
   clusterizer->SetNRowDiff(2);
   clusterizer->SetNColDiff(2);
   fClusterizer = clusterizer;
  } else{
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
  }
    
  //Now set the parameters
  fClusterizer->SetECAClusteringThreshold( fRecParam->GetClusteringThreshold() );
  fClusterizer->SetECALogWeight          ( fRecParam->GetW0()                  );
  fClusterizer->SetMinECut               ( fRecParam->GetMinECut()             );    
  fClusterizer->SetUnfolding             ( fRecParam->GetUnfold()              );
  fClusterizer->SetECALocalMaxCut        ( fRecParam->GetLocMaxCut()           );
  fClusterizer->SetTimeCut               ( fRecParam->GetTimeCut()             );
  fClusterizer->SetTimeMin               ( fRecParam->GetTimeMin()             );
  fClusterizer->SetTimeMax               ( fRecParam->GetTimeMax()             );
  fClusterizer->SetInputCalibrated       ( kTRUE                               );
  fClusterizer->SetJustClusters          ( kTRUE                               );  
  //In case of unfolding after clusterization is requested, set the corresponding parameters
  if(fRecParam->GetUnfold()){
    Int_t i=0;
    for (i = 0; i < 8; i++) {
      fClusterizer->SetSSPars(i, fRecParam->GetSSPars(i));
    }//end of loop over parameters
    for (i = 0; i < 3; i++) {
      fClusterizer->SetPar5  (i, fRecParam->GetPar5(i));
      fClusterizer->SetPar6  (i, fRecParam->GetPar6(i));
    }//end of loop over parameters
    
    fClusterizer->InitClusterUnfolding();
  }// to unfold
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::RecPoints2Clusters(TClonesArray *digitsArr, TObjArray *recPoints, TObjArray *clusArray)
{
  // Restore clusters from recPoints
  // Cluster energy, global position, cells and their amplitude fractions are restored
  
  for(Int_t i = 0; i < recPoints->GetEntriesFast(); i++)
  {
    AliEMCALRecPoint *recPoint = (AliEMCALRecPoint*) recPoints->At(i);
    
    const Int_t ncells = recPoint->GetMultiplicity();
    Int_t ncells_true = 0;
    
    // cells and their amplitude fractions
    UShort_t   absIds[ncells];  
    Double32_t ratios[ncells];
    
    for (Int_t c = 0; c < ncells; c++) {
      AliEMCALDigit *digit = (AliEMCALDigit*) digitsArr->At(recPoint->GetDigitsList()[c]);
      
      absIds[ncells_true] = digit->GetId();
      ratios[ncells_true] = recPoint->GetEnergiesList()[c]/digit->GetAmplitude();
      
      if (ratios[ncells_true] > 0.001) ncells_true++;
    }
    
    if (ncells_true < 1) {
      Warning("AliAnalysisTaskEMCALClusterize::RecPoints2Clusters", "skipping cluster with no cells");
      continue;
    }
    
    TVector3 gpos;
    Float_t g[3];
    
    // calculate new cluster position
    recPoint->EvalGlobalPosition(fRecParam->GetW0(), digitsArr);
    recPoint->GetGlobalPosition(gpos);
    gpos.GetXYZ(g);
    
    // create a new cluster
    AliAODCaloCluster *clus = new AliAODCaloCluster();
    clus->SetType(AliVCluster::kEMCALClusterv1);
    clus->SetE(recPoint->GetEnergy());
    clus->SetPosition(g);
    clus->SetNCells(ncells_true);
    clus->SetCellsAbsId(absIds);
    clus->SetCellsAmplitudeFraction(ratios);
    clus->SetDispersion(recPoint->GetDispersion());
    clus->SetChi2(-1); //not yet implemented
    clus->SetTOF(recPoint->GetTime()) ; //time-of-flight
    clus->SetNExMax(recPoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recPoint->GetElipsAxis(elipAxis);
    clus->SetM02(elipAxis[0]*elipAxis[0]) ;
    clus->SetM20(elipAxis[1]*elipAxis[1]) ;
    clus->SetDistanceToBadChannel(recPoint->GetDistanceToBadTower()); 
    
    //MC
    Int_t  parentMult  = 0;
    Int_t *parentList = recPoint->GetParents(parentMult);
    clus->SetLabel(parentList, parentMult); 
    
    //Write the second major contributor to each MC cluster.
    Int_t iNewLabel ;
    for ( Int_t iLoopCell = 0 ; iLoopCell < clus->GetNCells() ; iLoopCell++ ){
      
      Int_t idCell = clus->GetCellAbsId(iLoopCell) ;
      iNewLabel = 1 ; //iNewLabel makes sure we  don't write twice the same label.
      for ( UInt_t iLoopLabels = 0 ; iLoopLabels < clus->GetNLabels() ; iLoopLabels++ )
      {
	if ( fCellSecondLabels[idCell] == -1 )  iNewLabel = 0;  // -1 is never a good second label.
	if ( fCellSecondLabels[idCell] == clus->GetLabelAt(iLoopLabels) )  iNewLabel = 0;
      }
      if (iNewLabel == 1) 
      {
	Int_t * newLabelArray = (Int_t*)malloc((clus->GetNLabels()+1)*sizeof(Int_t)) ;
	for ( UInt_t iLoopNewLabels = 0 ; iLoopNewLabels < clus->GetNLabels() ; iLoopNewLabels++ ){
	  newLabelArray[iLoopNewLabels] = clus->GetLabelAt(iLoopNewLabels) ;
	}
	newLabelArray[clus->GetNLabels()] = fCellSecondLabels[idCell] ;
	clus->SetLabel(newLabelArray,clus->GetNLabels()+1) ;
      }
      
    }
    
    clusArray->Add(clus);
  } // recPoints loop
}
