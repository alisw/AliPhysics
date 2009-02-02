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

/* $Id$ */

//_________________________________________________________________________
//--
//-- Yves Schutz (SUBATECH) 
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
// 
//-- Aleksei Pavlinov : added staf for EMCAL jet trigger 9Apr 25, 2008)
//                    : fgDigitsArr should read just once at event

// --- ROOT system ---
#include <TList.h>
#include <TClonesArray.h>
#include <TH2.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALReconstructor.h"

#include "AliCodeTimer.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDtrack.h"
#include "AliEMCALLoader.h"
#include "AliEMCALRawUtils.h"
#include "AliEMCALDigit.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALPID.h"
#include "AliEMCALTrigger.h"
#include "AliRawReader.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCAL.h"
#include "AliEMCALHistoUtilities.h"
#include "AliESDVZERO.h"

#include "AliRunLoader.h"
#include "AliRun.h"

ClassImp(AliEMCALReconstructor) 

const AliEMCALRecParam* AliEMCALReconstructor::fgkRecParam = 0;  // EMCAL rec. parameters
AliEMCALRawUtils* AliEMCALReconstructor::fgRawUtils = 0;   // EMCAL raw utilities class
AliEMCALClusterizer* AliEMCALReconstructor::fgClusterizer = 0;   // EMCAL clusterizer class
TClonesArray*     AliEMCALReconstructor::fgDigitsArr = 0;  // shoud read just once at event
//____________________________________________________________________________
AliEMCALReconstructor::AliEMCALReconstructor() 
  : fDebug(kFALSE), fList(0), fGeom(0) 
{
  // ctor

  fgRawUtils = new AliEMCALRawUtils;
  fgClusterizer = new AliEMCALClusterizerv1;

  //To make sure we match with the geometry in a simulation file,
  //let's try to get it first.  If not, take the default geometry
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  } else {
    AliInfo(Form("Using default geometry in reconstruction"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

} 

//____________________________________________________________________________
AliEMCALReconstructor::AliEMCALReconstructor(const AliEMCALReconstructor & rec)
  : AliReconstructor(rec),
    fDebug(rec.fDebug),
    fList(rec.fList),
    fGeom(rec.fGeom)
{
  //copy ctor
}

//____________________________________________________________________________
AliEMCALReconstructor::~AliEMCALReconstructor()
{
  // dtor
  delete fGeom;
  AliCodeTimer::Instance()->Print();
} 

//____________________________________________________________________________
void AliEMCALReconstructor::Init()
{
  // Trigger hists - Oct 24, 2007
  fList = AliEMCALHistoUtilities::GetTriggersListOfHists(kTRUE);
}

//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESD object to retrieve the tracks reconstructed by 
  // the global tracking.
  // Works on the current event.

  AliCodeTimerAuto("")

  ReadDigitsArrayFromTree(digitsTree);
  fgClusterizer->InitParameters();
  fgClusterizer->SetOutput(clustersTree);

  if(fgDigitsArr && fgDigitsArr->GetEntries()) {

    fgClusterizer->SetInput(digitsTree);
    
    if(Debug())
      fgClusterizer->Digits2Clusters("deb all") ;
    else
      fgClusterizer->Digits2Clusters("");
    
    fgClusterizer->Clear();

  }

}

//____________________________________________________________________________
void AliEMCALReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const

{
  // Conversion from raw data to
  // EMCAL digits.
  // Works on a single-event basis

  rawReader->Reset() ; 

  TClonesArray *digitsArr = new TClonesArray("AliEMCALDigit",200);
  Int_t bufsize = 32000;
  digitsTree->Branch("EMCAL", &digitsArr, bufsize);

  //must be done here because, in constructor, option is not yet known
  fgRawUtils->SetOption(GetOption());

  fgRawUtils->SetRawFormatHighLowGainFactor(GetRecParam()->GetHighLowGainFactor());
  fgRawUtils->SetRawFormatOrder(GetRecParam()->GetOrderParameter());
  fgRawUtils->SetRawFormatTau(GetRecParam()->GetTau());
  fgRawUtils->SetNoiseThreshold(GetRecParam()->GetNoiseThreshold());
  fgRawUtils->SetNPedSamples(GetRecParam()->GetNPedSamples());

  fgRawUtils->Raw2Digits(rawReader,digitsArr);

  digitsTree->Fill();
  digitsArr->Delete();
  delete digitsArr;

}


//____________________________________________________________________________
void AliEMCALReconstructor::FillESD(TTree* digitsTree, TTree* clustersTree, 
				    AliESDEvent* esd) const
{
  // Called by AliReconstruct after Reconstruct() and global tracking and vertexing 
  // and V0 
  // Works on the current event
  //  printf(" ## AliEMCALReconstructor::FillESD() is started ### \n ");
  //return;

  //######################################################
  //#########Calculate trigger and set trigger info###########
  //######################################################
 
  AliEMCALTrigger tr;
  //   tr.SetPatchSize(1);  // create 4x4 patches
  tr.SetSimulation(kFALSE); // Reconstruction mode
  tr.SetDigitsList(fgDigitsArr);
  // Get VZERO total multiplicity for jet trigger simulation 
  // The simulation of jey trigger will be incorrect if no VZERO data 
  // at ESD
  AliESDVZERO* vZero = esd->GetVZEROData();
  if(vZero) {
    tr.SetVZER0Multiplicity(vZero->GetMTotV0A() + vZero->GetMTotV0C());
  }
  //
  tr.Trigger();

  Float_t maxAmp2x2  = tr.Get2x2MaxAmplitude();
  Float_t maxAmpnxn  = tr.GetnxnMaxAmplitude();
  Float_t ampOutOfPatch2x2  = tr.Get2x2AmpOutOfPatch() ;
  Float_t ampOutOfPatchnxn  = tr.GetnxnAmpOutOfPatch() ;

  Int_t iSM2x2      = tr.Get2x2SuperModule();
  Int_t iSMnxn      = tr.GetnxnSuperModule();
  Int_t iModulePhi2x2 = tr.Get2x2ModulePhi();
  Int_t iModulePhinxn = tr.GetnxnModulePhi();
  Int_t iModuleEta2x2 = tr.Get2x2ModuleEta();
  Int_t iModuleEtanxn = tr.GetnxnModuleEta();

  AliDebug(2, Form("Trigger 2x2 max amp %f, out amp %f, SM %d, iphi %d ieta %d",  maxAmp2x2, ampOutOfPatch2x2, iSM2x2,iModulePhi2x2, iModuleEta2x2));
  AliDebug(2, Form("Trigger 4x4 max amp %f , out amp %f, SM %d, iphi %d, ieta %d",  maxAmpnxn, ampOutOfPatchnxn, iSMnxn,iModulePhinxn, iModuleEtanxn));

  TVector3    pos2x2(-1,-1,-1);
  TVector3    posnxn(-1,-1,-1);

  Int_t iAbsId2x2 = fGeom->GetAbsCellIdFromCellIndexes( iSM2x2, iModulePhi2x2, iModuleEta2x2) ; // should be changed to Module
  Int_t iAbsIdnxn = fGeom->GetAbsCellIdFromCellIndexes( iSMnxn, iModulePhinxn, iModuleEtanxn) ;
  fGeom->GetGlobal(iAbsId2x2, pos2x2);
  fGeom->GetGlobal(iAbsIdnxn, posnxn);
  //printf(" iAbsId2x2 %i iAbsIdnxn %i \n", iAbsId2x2, iAbsIdnxn);
  
  TArrayF triggerPosition(6);
  triggerPosition[0] = pos2x2(0) ;   
  triggerPosition[1] = pos2x2(1) ;   
  triggerPosition[2] = pos2x2(2) ;  
  triggerPosition[3] = posnxn(0) ;   
  triggerPosition[4] = posnxn(1) ;   
  triggerPosition[5] = posnxn(2) ;
  //printf(" triggerPosition ");
  //for(int i=0; i<6; i++) printf(" %i %f : ", i, triggerPosition[i]);

  TArrayF triggerAmplitudes(4);
  triggerAmplitudes[0] = maxAmp2x2 ;   
  triggerAmplitudes[1] = ampOutOfPatch2x2 ;    
  triggerAmplitudes[2] = maxAmpnxn ;   
  triggerAmplitudes[3] = ampOutOfPatchnxn ;   
  //printf("\n triggerAmplitudes ");
  //for(int i=0; i<4; i++) printf(" %i %f : ", i, triggerAmplitudes[i]);
  //printf("\n");
  //tr.Print("");
  //
  // Trigger jet staff
  //
  if(tr.GetNJetThreshold()>0) {
    // Jet phi/eta
    Int_t n0 = triggerPosition.GetSize();
    const TH2F *hpatch = tr.GetJetMatrixE();
    triggerPosition.Set(n0 + 2);
    for(Int_t i=0; i<2; i++) triggerPosition[n0+i] = hpatch->GetMean(i+1);   
    // Add jet ampitudes
    n0 = triggerAmplitudes.GetSize();
    triggerAmplitudes.Set(n0 + tr.GetNJetThreshold());
    Double_t *ampJet = tr.GetL1JetThresholds();
    for(Int_t i=0; i<tr.GetNJetThreshold(); i++){
      triggerAmplitudes[n0 + i] = Float_t(ampJet[i]);
    }
  }
  esd->AddEMCALTriggerPosition(triggerPosition);
  esd->AddEMCALTriggerAmplitudes(triggerAmplitudes);
  // Fill trigger hists
  AliEMCALHistoUtilities::FillTriggersListOfHists(fList,&triggerPosition,&triggerAmplitudes);

  //########################################
  //##############Fill CaloCells###############
  //########################################

  TClonesArray *digits = new TClonesArray("AliEMCALDigit",1000);
  TBranch *branchdig = digitsTree->GetBranch("EMCAL");
  if (!branchdig) { 
    AliError("can't get the branch with the PHOS digits !");
    return;
  }
  branchdig->SetAddress(&digits);
  digitsTree->GetEvent(0);
  Int_t nDigits = digits->GetEntries(), idignew = 0 ;
  AliDebug(1,Form("%d digits",nDigits));

  AliESDCaloCells &emcCells = *(esd->GetEMCALCells());
  emcCells.CreateContainer(nDigits);
  emcCells.SetType(AliESDCaloCells::kEMCALCell);
  for (Int_t idig = 0 ; idig < nDigits ; idig++) {
    const AliEMCALDigit * dig = (const AliEMCALDigit*)digits->At(idig);
    if(dig->GetAmp() > 0 ){
      emcCells.SetCell(idignew,dig->GetId(),dig->GetAmp(), dig->GetTime());   
      idignew++;
    }
  }
  emcCells.SetNumberOfCells(idignew);
  emcCells.Sort();

  //------------------------------------------------------------
  //-----------------CLUSTERS-----------------------------
  //------------------------------------------------------------
  TObjArray *clusters = new TObjArray(100);
  TBranch *branch = clustersTree->GetBranch("EMCALECARP");
  branch->SetAddress(&clusters);
  clustersTree->GetEvent(0);

  Int_t nClusters = clusters->GetEntries(),  nClustersNew=0;
  AliDebug(1,Form("%d clusters",nClusters));
  esd->SetFirstEMCALCluster(esd->GetNumberOfCaloClusters()); // Put after Phos clusters 


  //######################################################
  //#######################TRACK MATCHING###############
  //######################################################
  //Fill list of integers, each one is index of track to which the cluster belongs.

  // step 1 - initialize array of matched track indexes
  Int_t *matchedTrack = new Int_t[nClusters];
  for (Int_t iclus = 0; iclus < nClusters; iclus++)
    matchedTrack[iclus] = -1;  // neg. index --> no matched track
  
  // step 2, change the flag for all matched clusters found in tracks
  Int_t iemcalMatch = -1;
  Int_t endtpc = esd->GetNumberOfTracks();
  for (Int_t itrack = 0; itrack < endtpc; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack) ; // retrieve track
    iemcalMatch = track->GetEMCALcluster();
    if(iemcalMatch >= 0) matchedTrack[iemcalMatch] = itrack;
  } 
  
  //########################################
  //##############Fill CaloClusters#############
  //########################################
  esd->SetNumberOfEMCALClusters(nClusters);
  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) {
    const AliEMCALRecPoint * clust = (const AliEMCALRecPoint*)clusters->At(iClust);
    //if(clust->GetClusterType()== AliESDCaloCluster::kEMCALClusterv1) nRP++; else nPC++;
    if (Debug()) clust->Print();
    // Get information from EMCAL reconstruction points
    Float_t xyz[3];
    TVector3 gpos;
    clust->GetGlobalPosition(gpos);
    for (Int_t ixyz=0; ixyz<3; ixyz++)
      xyz[ixyz] = gpos[ixyz];
    Float_t elipAxis[2];
    clust->GetElipsAxis(elipAxis);
       //Create digits lists
    Int_t cellMult = clust->GetMultiplicity();
    //TArrayS digiList(digitMult);
    Float_t *amplFloat = clust->GetEnergiesList();
    Int_t   *digitInts = clust->GetAbsId();
    TArrayS absIdList(cellMult);
    TArrayD fracList(cellMult);

    Int_t newCellMult = 0;
    for (Int_t iCell=0; iCell<cellMult; iCell++) {
      if (amplFloat[iCell] > 0) {
      absIdList[newCellMult] = (UShort_t)(digitInts[iCell]);
      //Uncomment when unfolding is done
      //if(emcCells.GetCellAmplitude(digitInts[iCell])>0)
      //fracList[newCellMult] = amplFloat[iCell]/(emcCells.GetCellAmplitude(digitInts[iCell])*calibration);//get cell calibration value 
      //else
      fracList[newCellMult] = 0; 
      newCellMult++;
      }
    }

    absIdList.Set(newCellMult);
    fracList.Set(newCellMult);
    
    if(newCellMult > 0) { // accept cluster if it has some digit
      nClustersNew++;
      //Primaries
      Int_t  parentMult  = 0;
      Int_t *parentList =  clust->GetParents(parentMult);
      // fills the ESDCaloCluster
      AliESDCaloCluster * ec = new AliESDCaloCluster() ;
      ec->SetClusterType(AliESDCaloCluster::kEMCALClusterv1);
      ec->SetPosition(xyz);
      ec->SetE(clust->GetEnergy());
      ec->SetNCells(newCellMult);
      //Change type of list from short to ushort
      UShort_t *newAbsIdList  = new UShort_t[newCellMult];
      Double_t *newFracList  = new Double_t[newCellMult];
      for(Int_t i = 0; i < newCellMult ; i++) {
        newAbsIdList[i]=absIdList[i];
        newFracList[i]=fracList[i];
      }
      ec->SetCellsAbsId(newAbsIdList);
      ec->SetCellsAmplitudeFraction(newFracList);
      ec->SetClusterDisp(clust->GetDispersion());
      ec->SetClusterChi2(-1); //not yet implemented
      ec->SetM02(elipAxis[0]*elipAxis[0]) ;
      ec->SetM20(elipAxis[1]*elipAxis[1]) ;
      ec->SetTOF(clust->GetTime()) ; //time-of-fligh
      ec->SetNExMax(clust->GetNExMax());          //number of local maxima
      TArrayI arrayTrackMatched(1);// Only one track, temporal solution.
      arrayTrackMatched[0]= matchedTrack[iClust];
      ec->AddTracksMatched(arrayTrackMatched);

      TArrayI arrayParents(parentMult,parentList);
      ec->AddLabels(arrayParents);

      // add the cluster to the esd object
      esd->AddCaloCluster(ec);
      delete ec;
      delete [] newAbsIdList ;
      delete [] newFracList ;
   }
 } // cycle on clusters

 delete [] matchedTrack;

 esd->SetNumberOfEMCALClusters(nClustersNew);
 //if(nClustersNew != nClusters)
 //printf(" ##### nClusters %i -> new %i ##### \n", nClusters, nClustersNew );

 //Fill ESDCaloCluster with PID weights
  AliEMCALPID *pid = new AliEMCALPID;
  //pid->SetPrintInfo(kTRUE);
  pid->SetReconstructor(kTRUE);
  pid->RunPID(esd);
  delete pid;
  
  delete digits;
  delete clusters;
  
  // printf(" ## AliEMCALReconstructor::FillESD() is ended : ncl %i -> %i ### \n ",nClusters, nClustersNew); 
}

//__________________________________________________________________________
void AliEMCALReconstructor::ReadDigitsArrayFromTree(TTree *digitsTree) const
{
  // See AliEMCALClusterizer::SetInput(TTree *digitsTree);
  if(fgDigitsArr) {
    // Clear previous digits 
    fgDigitsArr->Delete();
    delete fgDigitsArr;
  }
  // Read the digits from the input tree
  TBranch *branch = digitsTree->GetBranch("EMCAL");
  if (!branch) { 
    AliError("can't get the branch with the EMCAL digits !");
    return;
  }
  fgDigitsArr = new TClonesArray("AliEMCALDigit",100);
  branch->SetAddress(&fgDigitsArr);
  branch->GetEntry(0);
}

