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
// --- ROOT system ---
#include "TGeoManager.h"
#include "TGeoMatrix.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliAltroMapping.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSTracker.h"
#include "AliRawReader.h"
#include "AliPHOSCalibData.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliPHOSTrigger.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSRawFitterv0.h"
#include "AliPHOSRawFitterv1.h"
#include "AliPHOSRawFitterv2.h"
#include "AliPHOSRawFitterv3.h"
#include "AliPHOSRawFitterv4.h"
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSCpvRawDigiProducer.h"
#include "AliPHOSPulseGenerator.h"
#include "AliPHOSTriggerRawDigit.h"
#include "AliPHOSTriggerRawDigiProducer.h"
#include "AliPHOSTriggerParameters.h"
#include <fstream>

ClassImp(AliPHOSReconstructor)

Bool_t AliPHOSReconstructor::fgDebug = kFALSE ; 
TClonesArray*     AliPHOSReconstructor::fgDigitsArray = 0;   // Array of PHOS digits
TObjArray*        AliPHOSReconstructor::fgEMCRecPoints = 0;   // Array of EMC rec.points
TObjArray*        AliPHOSReconstructor::fgCPVRecPoints = 0;   // Array of CPV rec.points
AliPHOSCalibData * AliPHOSReconstructor::fgCalibData  = 0 ;
TClonesArray*     AliPHOSReconstructor::fgTriggerDigits = 0;   // Array of PHOS trigger digits

//____________________________________________________________________________
AliPHOSReconstructor::AliPHOSReconstructor() :
  fGeom(NULL),fClusterizer(NULL),fTSM(NULL),fPID(NULL),fTmpDigLG(NULL)
{
  // ctor
  fGeom          = AliPHOSGeometry::GetInstance("IHEP","");
  fClusterizer   = new AliPHOSClusterizerv1      (fGeom);
  fTSM           = new AliPHOSTrackSegmentMakerv1(fGeom);
  fPID           = new AliPHOSPIDv1              (fGeom);
  fTmpDigLG      = new TClonesArray("AliPHOSDigit",100);
  fgDigitsArray  = new TClonesArray("AliPHOSDigit",100);
  fgEMCRecPoints = new TObjArray(100) ;
  fgCPVRecPoints = new TObjArray(100) ;
  if (!fgCalibData)
    fgCalibData = new AliPHOSCalibData(-1); //use AliCDBManager's run number
  
  fgTriggerDigits  = new TClonesArray("AliPHOSTriggerRawDigit",100);
  
  AliInfo(Form("PHOS bad channel map contains %d bad channel(s).\n",
               fgCalibData->GetNumOfEmcBadChannels()));
 
}

//____________________________________________________________________________
AliPHOSReconstructor::~AliPHOSReconstructor()
{
  // dtor
  //  delete fGeom; // RS: fGeom is a static object owned by AliPHOSGeometry::fGeom, don't delete
  delete fClusterizer;
  delete fTSM;
  delete fPID;
  delete fTmpDigLG;
  delete fgDigitsArray;
  delete fgEMCRecPoints;
  delete fgCPVRecPoints;
  delete fgTriggerDigits;
} 

//____________________________________________________________________________
void AliPHOSReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  // 'single-event' local reco method called by AliReconstruction; 
  // Only the clusterization is performed,; the rest of the reconstruction is done in FillESD because the track
  // segment maker needs access to the AliESDEvent object to retrieve the tracks reconstructed by 
  // the global tracking.

  fClusterizer->InitParameters();
  fClusterizer->SetInput(digitsTree);
  fClusterizer->SetOutput(clustersTree);
  if ( Debug() ) 
    fClusterizer->Digits2Clusters("deb all") ; 
  else 
    fClusterizer->Digits2Clusters("") ;
}

//____________________________________________________________________________
void AliPHOSReconstructor::FillESD(TTree* digitsTree, TTree* clustersTree, 
				   AliESDEvent* esd) const
{
  // This method produces PHOS rec-particles,
  // then it creates AliESDtracks out of them and
  // write tracks to the ESD

  // do current event; the loop over events is done by AliReconstruction::Run()
  fTSM->SetESD(esd) ; 
  fTSM->SetInput(clustersTree);
  if ( Debug() ) 
    fTSM->Clusters2TrackSegments("deb all") ;
  else 
    fTSM->Clusters2TrackSegments("") ;
  
  fPID->SetInput(clustersTree, fTSM->GetTrackSegments()) ; 
  fPID->SetESD(esd) ; 
  if ( Debug() ) 
    fPID->TrackSegments2RecParticles("deb all") ;
  else 
    fPID->TrackSegments2RecParticles("") ;

  TClonesArray *recParticles  = fPID->GetRecParticles();
  Int_t nOfRecParticles = recParticles->GetEntriesFast();
  
  AliDebug(2,Form("%d rec. particles, option %s",nOfRecParticles,GetOption()));
  
  // Read digits array

  TBranch *branch = digitsTree->GetBranch("PHOS");
  if (!branch) { 
    AliError("can't get the branch with the PHOS digits !");
    return;
  }
  branch->SetAddress(&fgDigitsArray);
  branch->GetEntry(0);

  // Get the clusters array

  TBranch *emcbranch = clustersTree->GetBranch("PHOSEmcRP");
  if (!emcbranch) { 
    AliError("can't get the branch with the PHOS EMC clusters !");
    return;
  }

  emcbranch->SetAddress(&fgEMCRecPoints);
  emcbranch->GetEntry(0);
  
  TBranch *cpvbranch = clustersTree->GetBranch("PHOSCpvRP");
  if (cpvbranch) { 
    cpvbranch->SetAddress(&fgCPVRecPoints);
    cpvbranch->GetEntry(0);
  }
  
  // Trigger
  
  TBranch *tbranch = digitsTree->GetBranch("TPHOS");
  if (tbranch) { 
    
    tbranch->SetAddress(&fgTriggerDigits);
    tbranch->GetEntry(0);
    
    AliESDCaloTrigger* trgESD = esd->GetCaloTrigger("PHOS");
    
    if (trgESD) {
      trgESD->Allocate(fgTriggerDigits->GetEntriesFast());
      
      for (Int_t i = 0; i < fgTriggerDigits->GetEntriesFast(); i++) {
	AliPHOSTriggerRawDigit* tdig = (AliPHOSTriggerRawDigit*)fgTriggerDigits->At(i);
	
	Int_t mod,modX,modZ;
	tdig->GetModXZ(mod,modX,modZ);
	
	const Int_t relId[4] = {5-mod,0,modX+1,modZ+1};
	Int_t absId;
	
	fGeom->RelToAbsNumbering(relId,absId);
	trgESD->Add(mod,absId,tdig->GetAmp(),0.,(Int_t*)NULL,0,tdig->GetL1Threshold(),tdig->GetType());
      }
    }  
  }
  
//   //#########Calculate trigger and set trigger info###########

//   AliPHOSTrigger tr ;
//   //   tr.SetPatchSize(1);//create 4x4 patches
//   tr.SetSimulation(kFALSE);
//   tr.Trigger(fgDigitsArray);
  
//   Float_t maxAmp2x2  = tr.Get2x2MaxAmplitude();
//   Float_t maxAmpnxn  = tr.GetnxnMaxAmplitude();
//   Float_t ampOutOfPatch2x2  = tr.Get2x2AmpOutOfPatch() ;
//   Float_t ampOutOfPatchnxn  = tr.GetnxnAmpOutOfPatch() ;

//   Int_t iSM2x2      = tr.Get2x2SuperModule();
//   Int_t iSMnxn      = tr.GetnxnSuperModule();
//   Int_t iCrystalPhi2x2 = tr.Get2x2CrystalPhi();
//   Int_t iCrystalPhinxn = tr.GetnxnCrystalPhi();
//   Int_t iCrystalEta2x2 = tr.Get2x2CrystalEta();
//   Int_t iCrystalEtanxn = tr.GetnxnCrystalEta();

//   AliDebug(2, Form("Trigger 2x2 max amp %f, out amp %f, SM %d, iphi %d ieta %d",  
// 		   maxAmp2x2, ampOutOfPatch2x2, iSM2x2,iCrystalPhi2x2, iCrystalEta2x2));
//   AliDebug(2, Form("Trigger 4x4 max amp %f , out amp %f, SM %d, iphi %d, ieta %d",
// 		   maxAmpnxn, ampOutOfPatchnxn, iSMnxn,iCrystalPhinxn, iCrystalEtanxn));

//   // Attention! PHOS modules in order to calculate AbsId need to be 1-5 not 0-4 as returns trigger.
//   Int_t iRelId2x2 []= {iSM2x2+1,0,iCrystalPhi2x2,iCrystalEta2x2};
//   Int_t iAbsId2x2 =-1;
//   Int_t iRelIdnxn []= {iSMnxn+1,0,iCrystalPhinxn,iCrystalEtanxn};
//   Int_t iAbsIdnxn =-1;
//   TVector3    pos2x2(-1,-1,-1);
//   TVector3    posnxn(-1,-1,-1);
//   fGeom->RelToAbsNumbering(iRelId2x2, iAbsId2x2);
//   fGeom->RelToAbsNumbering(iRelIdnxn, iAbsIdnxn);
//   fGeom->RelPosInAlice(iAbsId2x2, pos2x2);
//   fGeom->RelPosInAlice(iAbsIdnxn, posnxn);

//   TArrayF triggerPosition(6);
//   triggerPosition[0] = pos2x2(0) ;   
//   triggerPosition[1] = pos2x2(1) ;   
//   triggerPosition[2] = pos2x2(2) ;  
//   triggerPosition[3] = posnxn(0) ;   
//   triggerPosition[4] = posnxn(1) ;   
//   triggerPosition[5] = posnxn(2) ;  

//   TArrayF triggerAmplitudes(4);
//   triggerAmplitudes[0] = maxAmp2x2 ;   
//   triggerAmplitudes[1] = ampOutOfPatch2x2 ;    
//   triggerAmplitudes[2] = maxAmpnxn ;   
//   triggerAmplitudes[3] = ampOutOfPatchnxn ;   

//   //esd->SetPHOSTriggerCells(triggerPosition);
//   esd->AddPHOSTriggerPosition(triggerPosition);
//   esd->AddPHOSTriggerAmplitudes(triggerAmplitudes);
  

  //########################################
  //############# Fill CaloCells ###########
  //########################################
  Int_t nDigits = fgDigitsArray->GetEntries();
  Int_t idignew = 0 ;
  AliDebug(1,Form("%d digits",nDigits));

  const Int_t knEMC = fGeom->GetNModules()*fGeom->GetNPhi()*fGeom->GetNZ();
  AliESDCaloCells &phsCells = *(esd->GetPHOSCells());
  phsCells.CreateContainer(nDigits);
  phsCells.SetType(AliESDCaloCells::kPHOSCell);

  // Add to CaloCells only EMC digits with non-zero energy 
  for (Int_t idig = 0 ; idig < nDigits ; idig++) {
    const AliPHOSDigit * dig = (const AliPHOSDigit*)fgDigitsArray->At(idig);
    if(dig->GetId() <= knEMC && 
       Calibrate(dig->GetEnergy(),dig->GetId()) > GetRecoParam()->GetEMCMinE() ){
      Int_t primary = dig->GetPrimary(1) ;
      phsCells.SetCell(idignew,dig->GetId(), Calibrate(dig->GetEnergy(),dig->GetId()),
                                             CalibrateT(dig->GetTime(),dig->GetId(),dig->IsLG()),
                                             primary,0.,!dig->IsLG()) ;
					     
      idignew++;
    }
    //Fill CPV digits
    if(dig->GetId() > knEMC && 
       Calibrate(dig->GetEnergy(),dig->GetId()) > GetRecoParam()->GetCPVMinE() ){
      Int_t primary = dig->GetPrimary(1) ;
      //NB! CPV cell ID can not fit Short_t <32000
      //So, for CPV we subtract EMC part and set negative ID
      phsCells.SetCell(idignew,-(dig->GetId()-56*64*5), Calibrate(dig->GetEnergy(),dig->GetId()),
                                             0.,
                                             primary,0.,0) ;					     
      idignew++;
    }    
  }
  phsCells.SetNumberOfCells(idignew);
  phsCells.Sort();

  //########################################
  //############## Fill CaloClusters #######
  //########################################

  for (Int_t recpart = 0 ; recpart < nOfRecParticles ; recpart++) {
    AliPHOSRecParticle  *rp    = static_cast<AliPHOSRecParticle*>(recParticles->At(recpart));
    if (Debug()) 
      rp->Print();
    // Get track segment and EMC rec.point associated with this rec.particle
    AliPHOSTrackSegment *ts    = static_cast<AliPHOSTrackSegment *>(fTSM->GetTrackSegments()
								    ->At(rp->GetPHOSTSIndex()));

    AliPHOSEmcRecPoint  *emcRP = static_cast<AliPHOSEmcRecPoint *>(fgEMCRecPoints->At(ts->GetEmcIndex()));
    AliESDCaloCluster   *ec    = new AliESDCaloCluster() ; 
    
    Float_t xyz[3];
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = rp->GetPos()[ixyz];
    
    AliDebug(2,Form("Global position xyz=(%f,%f,%f)",xyz[0],xyz[1],xyz[2]));
   
    // Create cell lists

    Int_t     cellMult   = emcRP->GetDigitsMultiplicity();
    Int_t    *digitsList = emcRP->GetDigitsList();
    Float_t  *rpElist    = emcRP->GetEnergiesList() ;
    UShort_t *absIdList  = new UShort_t[cellMult];
    Double_t *fracList   = new Double_t[cellMult];

    for (Int_t iCell=0; iCell<cellMult; iCell++) {
      AliPHOSDigit *digit = static_cast<AliPHOSDigit *>(fgDigitsArray->At(digitsList[iCell]));
      absIdList[iCell] = (UShort_t)(digit->GetId());
      if (digit->GetEnergy() > 0)
 	fracList[iCell] = rpElist[iCell]/(Calibrate(digit->GetEnergy(),digit->GetId()));
      else
 	fracList[iCell] = 0;
    }

    //Primaries
    Int_t  primMult  = 0;
    Int_t *primList =  emcRP->GetPrimaries(primMult);

    Float_t energy=0.;
    if (GetRecoParam()->EMCEcore2ESD())
      energy = emcRP->GetCoreEnergy();
    else
      energy = rp->Energy();
    //Apply nonlinearity correction
    if(GetRecoParam()->GetEMCEnergyCorrectionOn())
      energy=CorrectNonlinearity(energy) ;

    // fills the ESDCaloCluster
    ec->SetType(AliVCluster::kPHOSNeutral);
    ec->SetPosition(xyz);                       //rec.point position in MARS
    ec->SetE(energy);                           //total or core particle energy
    ec->SetDispersion(emcRP->GetDispersion());  //cluster dispersion
    ec->SetPID(rp->GetPID()) ;            //array of particle identification
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima
    ec->SetEmcCpvDistance(ts->GetCpvDistance("r")); //Only radius, what about separate x,z????
    ec->SetTrackDistance(ts->GetCpvDistance("x"),ts->GetCpvDistance("z")); 
    ec->SetChi2(-1);                     //not yet implemented
    ec->SetTOF(emcRP->GetTime());               //Time of flight - already calibrated in EMCRecPoint

    //Cells contributing to clusters
    ec->SetNCells(cellMult);
    ec->SetCellsAbsId(absIdList);
    ec->SetCellsAmplitudeFraction(fracList);

    //Distance to the nearest bad crystal
    ec->SetDistanceToBadChannel(emcRP->GetDistanceToBadCrystal()); 
  
    //Array of MC indeces
    TArrayI arrayPrim(primMult,primList);
    ec->AddLabels(arrayPrim);
    
    //Matched ESD track
    if(ts->GetTrackIndex()>=0){ //there is a match
      TArrayI arrayTrackMatched(1);
      arrayTrackMatched[0]= ts->GetTrackIndex();
      ec->AddTracksMatched(arrayTrackMatched);
    
      Int_t index = esd->AddCaloCluster(ec);

      //Set pointer to this cluster in ESD track
      Int_t nt=esd->GetNumberOfTracks();
      for (Int_t itr=0; itr<nt; itr++) {
        AliESDtrack *esdTrack=esd->GetTrack(itr);
        if(!esdTrack->IsPHOS())
          continue ;
        if(esdTrack->GetPHOScluster()==-recpart){ //we store negative cluster number
          esdTrack->SetPHOScluster(index) ;
//no garatie that only one track matched this cluster
//      break ;
        }
      }
    }
    else{ //Story empty list
      TArrayI arrayTrackMatched(0);
      ec->AddTracksMatched(arrayTrackMatched);      
      esd->AddCaloCluster(ec);
    }
 
    delete ec;   
    delete [] fracList;
    delete [] absIdList;
  }

  //Fill CPV clusters
  Int_t nOfCPVclu=fgCPVRecPoints->GetEntriesFast() ;
  Int_t nOfTS = fTSM->GetTrackSegments()->GetEntriesFast() ;
  for (Int_t recpoint = 0 ; recpoint < nOfCPVclu ; recpoint++) {
    AliPHOSCpvRecPoint  *cpvRP = static_cast<AliPHOSCpvRecPoint *>(fgCPVRecPoints->At(recpoint));
    AliESDCaloCluster   *ec    = new AliESDCaloCluster() ; 
    
    Float_t xyz[3];
    TVector3 pos ; 
    fGeom->GetGlobalPHOS(cpvRP, pos) ; 
    for (Int_t ixyz=0; ixyz<3; ixyz++) 
      xyz[ixyz] = pos[ixyz];
       
    // Create cell lists
    Int_t     cellMult   = cpvRP->GetDigitsMultiplicity();
    Int_t    *digitsList = cpvRP->GetDigitsList();
    Float_t  *rpElist    = cpvRP->GetEnergiesList() ;
    UShort_t *absIdList  = new UShort_t[cellMult];
    Double_t *fracList   = new Double_t[cellMult];

    for (Int_t iCell=0; iCell<cellMult; iCell++) {
      AliPHOSDigit *digit = static_cast<AliPHOSDigit *>(fgDigitsArray->At(digitsList[iCell]));
      absIdList[iCell] = (UShort_t)(digit->GetId());
      if (digit->GetEnergy() > 0)
 	fracList[iCell] = rpElist[iCell]/(Calibrate(digit->GetEnergy(),digit->GetId()));
      else
 	fracList[iCell] = 0;
    }

    //Primaries
    Int_t  primMult  = 0;
    Int_t *primList =  cpvRP->GetPrimaries(primMult);

    Float_t energy = cpvRP->GetEnergy();
    // fills the ESDCaloCluster
    ec->SetType(AliVCluster::kPHOSCharged);
    ec->SetPosition(xyz);                       //rec.point position in MARS
    ec->SetE(energy);                           //total or core particle energy
    ec->SetDispersion(cpvRP->GetDispersion());  //cluster dispersion
    ec->SetM02(cpvRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(cpvRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(cpvRP->GetNExMax());          //number of local maxima
    ec->SetChi2(-1);                     //not yet implemented

    //Cells contributing to clusters
    ec->SetNCells(cellMult);
    ec->SetCellsAbsId(absIdList);
    ec->SetCellsAmplitudeFraction(fracList);

    //Distance to the nearest bad crystal
    ec->SetDistanceToBadChannel(cpvRP->GetDistanceToBadCrystal()); 
  
    //Array of MC indeces
    TArrayI arrayPrim(primMult,primList);
    ec->AddLabels(arrayPrim);
    
    //Matched CPV-EMC clusters
    //CPV TrackSerments go after EMC-related TSs
    if(recpoint+nOfRecParticles<nOfTS){
      AliPHOSTrackSegment *ts    = static_cast<AliPHOSTrackSegment *>(fTSM->GetTrackSegments()
								    ->At(recpoint+nOfRecParticles));  
      if(ts){
	if(ts->GetEmcIndex()>=0){
          TArrayI arrayTrackMatched(1);
          arrayTrackMatched[0]= ts->GetEmcIndex();
          ec->AddTracksMatched(arrayTrackMatched);
          ec->SetTrackDistance(ts->GetCpvDistance("x"),ts->GetCpvDistance("z")); 
	}
	else{//No match
          TArrayI arrayTrackMatched(0);
          ec->AddTracksMatched(arrayTrackMatched);
          ec->SetTrackDistance(999.,999.); 
	}
      }
    }
    Int_t index = esd->AddCaloCluster(ec);

    //Set pointer to this cluster in ESD track
//    Int_t nt=esd->GetNumberOfTracks();
//    for (Int_t itr=0; itr<nt; itr++) {
//      AliESDtrack *esdTrack=esd->GetTrack(itr);
//      if(!esdTrack->IsPHOS())
//        continue ;
//      if(esdTrack->GetPHOScluster()==-recpart){ //we store negative cluster number
//        esdTrack->SetPHOScluster(index) ;
//no garatie that only one track matched this cluster
//      break ;
//      }
//    }
 
    delete ec;   
    delete [] fracList;
    delete [] absIdList;
  }
  
  
  fgDigitsArray ->Clear("C");
  fgEMCRecPoints->Clear("C");
  recParticles  ->Clear();

  //Store PHOS misalignment matrixes
  FillMisalMatrixes(esd) ;

}

//____________________________________________________________________________
AliTracker* AliPHOSReconstructor::CreateTracker() const
{
  // creates the PHOS tracker
  return new AliPHOSTracker();
}

//____________________________________________________________________________
void  AliPHOSReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  // Converts raw data to
  // PHOS digits
  // Works on a single-event basis
  rawReader->Reset() ; 

  // Create a new array of PHOS and CPV digits and fill it in PHOS and CPV raw data decoders
  TClonesArray *digits = new TClonesArray("AliPHOSDigit",1);
  digits->SetName("DIGITS");
  Int_t bufsize = 32000;
  digitsTree->Branch("PHOS", &digits, bufsize);

  ConvertDigitsEMC(rawReader,digits);
  ConvertDigitsCPV(rawReader,digits);

  AliDebug(2,Form("Number of created digits = %d",digits->GetEntriesFast()));

  // Create a new array of PHOS trigger digits and fill it from raw data
  TClonesArray *tdigits = new TClonesArray("AliPHOSTriggerRawDigit",1);
  tdigits->SetName("TDIGITS");
  digitsTree->Branch("TPHOS", &tdigits, bufsize);  

  rawReader->Reset();
  AliPHOSTriggerRawDigiProducer tdp(rawReader);
  
  AliPHOSTriggerParameters* parameters = (AliPHOSTriggerParameters*)AliPHOSRecoParam::GetTriggerParameters();
  
  tdp.SetTriggerParameters(parameters);
  tdp.ProcessEvent(tdigits);
  
  if (AliLog::GetGlobalDebugLevel() == 1) {
    Int_t modMax=-111;
    Int_t colMax=-111;
    Int_t rowMax=-111;
    Float_t eMax=-333;
    //!!!for debug!!!
    
    Int_t relId[4];
    for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
      AliPHOSDigit* digit = (AliPHOSDigit*)digits->At(iDigit);
      if(digit->GetEnergy()>eMax) {
	fGeom->AbsToRelNumbering(digit->GetId(),relId);
	eMax=digit->GetEnergy();
	modMax=relId[0];
	rowMax=relId[2];
	colMax=relId[3];
      }
    }
    
    AliDebug(1,Form("Digit with max. energy:  modMax %d colMax %d rowMax %d  eMax %f\n\n",
		    modMax,colMax,rowMax,eMax));
  }
  
  digitsTree->Fill();
  
  digits->Delete();
  delete digits;
  
  tdigits->Delete();
  delete tdigits;
}
//==================================================================================
void  AliPHOSReconstructor::ConvertDigitsEMC(AliRawReader* rawReader, TClonesArray* digits) const
{
  // Converts CPV raw data to PHOS EMC digits
  // Works on a single-event basis
  AliPHOSRawFitterv0 * fitter ;

  const TObjArray* maps = AliPHOSRecoParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  AliAltroMapping *mapping[20];
  for(Int_t i = 0; i < 20; i++) {
    mapping[i] = (AliAltroMapping*)maps->At(i);
  }

  if      (strcmp(GetRecoParam()->EMCFitterVersion(),"v0")==0) 
    fitter=new AliPHOSRawFitterv0();
  else if (strcmp(GetRecoParam()->EMCFitterVersion(),"v1")==0) 
    fitter=new AliPHOSRawFitterv1();
  else if (strcmp(GetRecoParam()->EMCFitterVersion(),"v2")==0) 
    fitter=new AliPHOSRawFitterv2();
  else if (strcmp(GetRecoParam()->EMCFitterVersion(),"v3")==0) 
    fitter=new AliPHOSRawFitterv3();
  else
    fitter=new AliPHOSRawFitterv4();

  fitter->SubtractPedestals(GetRecoParam()->EMCSubtractPedestals());
  fitter->SetAmpOffset     (GetRecoParam()->GetGlobalAltroOffset());
  fitter->SetAmpThreshold  (GetRecoParam()->GetGlobalAltroThreshold());

  AliPHOSRawDigiProducer rdp(rawReader,mapping);

  rdp.SetEmcMinAmp(GetRecoParam()->GetEMCRawDigitThreshold()); // in ADC
  rdp.SetCpvMinAmp(GetRecoParam()->GetCPVMinE());
  rdp.SetSampleQualityCut(GetRecoParam()->GetEMCSampleQualityCut());
  rdp.SetSubtractL1phase(GetRecoParam()->GetSubtractL1phase());
  rdp.MakeDigits(digits,fTmpDigLG,fitter);

  delete fitter ;
}
//==================================================================================
void  AliPHOSReconstructor::ConvertDigitsCPV(AliRawReader* rawReader, TClonesArray* digits) const
{
  // Converts CPV raw data to PHOS CPV digits
  // Works on a single-event basis
  AliPHOSCpvRawDigiProducer rdp(rawReader);
  rdp.SetCalibData(fgCalibData);
  rdp.SetCpvMinAmp(GetRecoParam()->GetCPVMinE());
  rdp.SetTurbo(kTRUE);
  rdp.MakeDigits(digits);
    
}
//==================================================================================
Float_t AliPHOSReconstructor::Calibrate(Float_t amp, Int_t absId)const{
  // Calibrate EMC digit, i.e. multiply its Amp by a factor read from CDB

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]){ //CPV
    Float_t calibration = fgCalibData->GetADCchannelCpv(module,row,column);//corrected by sevdokim
    return amp*calibration ;
  }
  else{ //EMC
    Float_t calibration = fgCalibData->GetADCchannelEmc(module,column,row);
    return amp*calibration ;
  }
}
//==================================================================================
Float_t AliPHOSReconstructor::CalibrateT(Float_t time, Int_t absId,Bool_t isLG)const{
  // Calibrate EMC digit, i.e. multiply its Amp by a factor read from CDB

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]){ //CPV
    return 0. ;
  }
  else{ //EMC
    if(isLG)
      time += fgCalibData->GetLGTimeShiftEmc(module,column,row);
    else
      time += fgCalibData->GetTimeShiftEmc(module,column,row);
    return time ;
  }
}
//==================================================================================
void AliPHOSReconstructor::FillMisalMatrixes(AliESDEvent* esd)const{
  //Store PHOS matrixes in ESD Header

  //Check, if matrixes was already stored
  for(Int_t mod=0 ;mod<5; mod++){
    if(esd->GetPHOSMatrix(mod)!=0)
      return ;
  }

  //Create and store matrixes
  if(!gGeoManager){
    AliError("Can not store misal. matrixes: no gGeoManager! \n") ;
    return ;
  }
  //Note, that owner of copied marixes will be header
  char path[255] ;
  TGeoHMatrix * m ;
  for(Int_t mod=0; mod<5; mod++){
    snprintf(path,255,"/ALIC_1/PHOS_%d",mod+1) ; //In Geometry modules numbered 1,2,.,5 without CPV
    if (gGeoManager->CheckPath(path)){
      gGeoManager->cd(path) ;
      m = gGeoManager->GetCurrentMatrix() ;
      esd->SetPHOSMatrix(new TGeoHMatrix(*m),mod) ;
    }
    else{
      snprintf(path,255,"/ALIC_1/PHOC_%d",mod+1) ; //In Geometry modules numbered 1,2,3 with CPV
      if (gGeoManager->CheckPath(path)){
        gGeoManager->cd(path) ;
        m = gGeoManager->GetCurrentMatrix() ;
        esd->SetPHOSMatrix(new TGeoHMatrix(*m),mod) ;
      }
      else{
        snprintf(path,255,"/ALIC_1/PHOH_%d",mod+1) ; //In Geometry modules numbered 1,2,.,5
        if (gGeoManager->CheckPath(path)){
          gGeoManager->cd(path) ;
          m = gGeoManager->GetCurrentMatrix() ;
          esd->SetPHOSMatrix(new TGeoHMatrix(*m),mod) ;
        }
        else{
          esd->SetPHOSMatrix(NULL,mod) ;
        }
      }
    }
  }

}
//==================================================================================
Float_t AliPHOSReconstructor::CorrectNonlinearity(Float_t en){

  //For backward compatibility, if no RecoParameters found
  if(!GetRecoParam()){
    return 0.0241+1.0504*en+0.000249*en*en ;
  }

  if(strcmp(GetRecoParam()->GetNonlinearityCorrectionVersion(),"NoCorrection")==0){
    return en ;
  }
  if(strcmp(GetRecoParam()->GetNonlinearityCorrectionVersion(),"Gustavo2005")==0){
    const Float_t *par=GetRecoParam()->GetNonlinearityParams() ;
    return par[0]+par[1]*en + par[2]*en*en ;
  }
  if(strcmp(GetRecoParam()->GetNonlinearityCorrectionVersion(),"Henrik2010")==0){
     const Float_t *par=GetRecoParam()->GetNonlinearityParams() ;
     return en*(par[0]+par[1]*TMath::Exp(-en*par[2]))*(1.+par[3]*TMath::Exp(-en*par[4]))*(1.+par[6]/(en*en+par[5])) ;
  }
  //For backward compatibility
  if(strcmp(GetRecoParam()->GetNonlinearityCorrectionVersion(),"")==0){
    return 0.0241+1.0504*en+0.000249*en*en ;
  }
  return en ;
}

void AliPHOSReconstructor::readTRUParameters(AliPHOSTriggerParameters* parameters) const
{
  //Read trigger parameters.

  TString path(gSystem->Getenv("ALICE_ROOT"));
  path += "/PHOS/macros/Trigger/OCDB/";
  
  for (Int_t mod = 2; mod < 5; ++mod) { // module
    for (Int_t tru = 0; tru < 4; tru++) { // tru row
      for (Int_t branch = 0; branch < 2; branch++) { // branch
	
	// Open the Appropriate pedestal file
	TString fileName = path;
	fileName += "pedestal_m";
	fileName = fileName += mod;
	fileName+="_r";
	fileName+=tru;
	fileName+="_b";
	fileName+=branch;
	fileName+=".dat";
	std::ifstream instream;
	instream.open(fileName.Data());
	
	// Read pedestals from file
	if( ! instream.is_open() )
	  Printf("E-TRUPedestals: could not open %s", fileName.Data());
	else
	  {
	    Int_t ped[112];
	    
	    char ch_s[36];
	    char *ch_s_p = ch_s;
	    //Int_t nlines = 0;
	    
	    Int_t t_ped_0 =0;
	    Int_t t_ped_1 =0;
	    Int_t t_ped_2 =0;
	    
	    for(Int_t n=0; n<112; n++)
	      {
		instream.getline(ch_s_p,36);
		if (ch_s_p[23]<='9' && ch_s_p[23]>='0')
		  {
		    t_ped_0 = ch_s_p[23]-'0';
		  }
		else if (ch_s_p[23]>='A' && ch_s_p[23]<='Z')
		  {
		    t_ped_0 = ch_s_p[23]-'A'+10;
		    
		  }
		  
		if (ch_s_p[22]<='9' && ch_s_p[22]>='0')
		  {
		    t_ped_1 = ch_s_p[22]-'0';
		  }
		else if (ch_s_p[22]<='Z' && ch_s_p[22]>='A')
		  {
		    t_ped_1 = ch_s_p[22]-'A'+10;
		  }
		
		if (ch_s_p[21]<='9' && ch_s_p[21]>='0')
		  {
		    t_ped_2 = ch_s_p[21]-'0';
		  }
		else if (ch_s_p[21]<='Z' && ch_s_p[21]>='A')
		  {
		    t_ped_2 = ch_s_p[21]-'A'+10;
		  }
		
		ped[n]=t_ped_2*256+t_ped_1*16+t_ped_0;
		
		
	      }
	    for (Int_t xrow = 0; xrow < 8; xrow++){
	      for (Int_t zcol = 0; zcol < 14; zcol++){
		Int_t pedestal = ped[zcol*8+xrow];
		
		if( pedestal < 612 && pedestal > 412 ) // resonable
		  parameters->SetTRUPedestal(pedestal, mod, tru, branch, xrow, zcol);
		else // unresonable
		  continue;
	      }
	    }
	  } // Ends read of pedestals from branch from file.
	instream.close();
      }// end branch
    }// end tru
    
  }// end for mod
}




