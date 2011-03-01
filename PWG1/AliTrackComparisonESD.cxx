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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ANALYSIS task to perrorm TPC calibration                                  //

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliTrackComparisonESD.h"
#include "AliTrackComparison.h"
#include "TChain.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackPointArray.h"
#include "AliESDtrackCuts.h"
#include "AliTracker.h"
#include "AliESDCaloCluster.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliEMCALGeometry.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "AliHMPIDParam.h"
//#include <TGeoHMatrix>
#include "AliGeomManager.h"
//#include "AliCDBManager.h"
//#include "AliGRPManager.h"

ClassImp(AliTrackComparisonESD)

//________________________________________________________________________
AliTrackComparisonESD::AliTrackComparisonESD()
  :AliAnalysisTask(),
   fESD(0),
   fESDCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
   fESDfriend(0),
   fCurrentRun(-1),
   fDebugOutputPath(""),
   //   fOcdbPath("local:///lustre/alice/alien/alice/data/2010/OCDB"),
   fOutput(0),
   fEMCAL(0),
   fHMPID(0),
   fTOF(0),
   fGeom(0),
   fCutX(10),
   fCutY(10),
   fCutZ(10)
{
  //
  // default constructor
  // 
  for(Int_t i=0; i<4; i++) fTransMatrix[i]=0;
  
}

//________________________________________________________________________
AliTrackComparisonESD::AliTrackComparisonESD(const char *name) 
  :AliAnalysisTask(name,""),
   fESD(0),
   fESDCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010()),
   fESDfriend(0),
   fCurrentRun(-1),
   fDebugOutputPath(""),
   //   fOcdbPath("local:///lustre/alice/alien/alice/data/2010/OCDB"),
   fOutput(0),
   fEMCAL(0),
   fHMPID(0),
   fTOF(0),
   fGeom(0),
   fCutX(10),
   fCutY(10),
   fCutZ(10)
{
  //
  // Constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TObjArray::Class());

  for(Int_t i=0; i<4; i++) fTransMatrix[i]=0;
}

//________________________________________________________________________
AliTrackComparisonESD::~AliTrackComparisonESD() {
  //
  // destructor
  //
  printf("AliTrackComparisonESD::~AliTrackComparisonESD");
  if(fOutput) delete fOutput; fOutput=0;
  if(fEMCAL)  delete fEMCAL;  fEMCAL=0;
  if(fHMPID)  delete fHMPID;  fHMPID=0;
  if(fTOF)    delete fTOF;    fTOF=0;
  for(Int_t i=0; i<4; i++)
    {
      if(fTransMatrix[i]) {delete fTransMatrix[i];fTransMatrix[i]=0;}
    }
}

//________________________________________________________________________
void AliTrackComparisonESD::ConnectInputData(Option_t *) {
  //
  //
  //
  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    //Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      //Printf("ERROR: Could not get ESDInputHandler");
    } 
    else {
      //esdH->SetReadFriends(kTRUE);
      //esdH->SetActiveBranches("ESDfriend");
      fESD = esdH->GetEvent();
      //Printf("*** CONNECTED NEW EVENT ****");
    }
  }

}

//________________________________________________________________________
void AliTrackComparisonESD::CreateOutputObjects() {
  //
  //
  //
  //OpenFile(0, "RECREATE");
  TFile *ftmp = OpenFile(0);
  if(!ftmp)AliError(Form("File %s not found!",ftmp->GetName()));
  fOutput=new TObjArray(0);
  fOutput->SetOwner(kTRUE);

  fEMCAL = new AliTrackComparison("EMCAL","EMCAL");
  fEMCAL->SetLayerID(20);
  fEMCAL->SetFillAll(kFALSE);
  fEMCAL->Init();
  fHMPID = new AliTrackComparison("HMPID","HMPID");
  fHMPID->SetRangeDY(-5,5);
  fHMPID->SetRangeDZ(-5,5);
  fHMPID->SetLayerID(18);
  fHMPID->SetFillAll(kFALSE);
  fHMPID->Init();
  fTOF   = new AliTrackComparison("TOF","TOF");
  fTOF->SetRangeDY(-5,5);
  fTOF->SetRange1Pt(-1,1);
  fTOF->SetLayerID(15);
  fTOF->SetFillAll(kFALSE);
  fTOF->Init();

  fOutput->Add(fEMCAL);
  fOutput->Add(fHMPID);
  fOutput->Add(fTOF);

  Double_t rotationMatrix[4][9] = {{-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996},
				   {-0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996},
				   {-0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990},
				   {-0.345864,  0.938278,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990}};

  Double_t translationMatrix[4][3] = {{0.351659, 447.576446,  176.269742},
				      {1.062577, 446.893974, -173.728870},
				      {-154.213287, 419.306156,  176.753692},
				      {-153.018950, 418.623681, -173.243605}};
  for(Int_t imx=0; imx<4; imx++)
    {
      fTransMatrix[imx] = new TGeoHMatrix();
      fTransMatrix[imx]->SetRotation(rotationMatrix[imx]);
      fTransMatrix[imx]->SetTranslation(translationMatrix[imx]);
      fTransMatrix[imx]->Print();
    }


  PostData(0,fOutput);
}

//________________________________________________________________________
Bool_t AliTrackComparisonESD::SetupEvent() {
  //
  // Setup Event
  //
  // check if something to be done

  if(!fESD)
    return kFALSE;

  if (fCurrentRun == fESD->GetRunNumber())
    return kTRUE;
  else
    fCurrentRun = fESD->GetRunNumber();

//   // OCDB
//   printf("setting run to %d\n",fCurrentRun);
//   AliCDBManager::Instance()->SetDefaultStorage(fOcdbPath.Data());
//   AliCDBManager::Instance()->SetRun(fCurrentRun); 

//   // magnetic field
//   if ( !TGeoGlobalMagField::Instance()->GetField() ) {
//     printf("Loading field map...\n");
//     AliGRPManager grpMan;
//     if( !grpMan.ReadGRPEntry() ) { 
//       printf("Cannot get GRP entry\n"); 
//       return kFALSE;
//     }
//     if( !grpMan.SetMagField() ) { 
//       printf("Problem with magnetic field setup\n"); 
//       return kFALSE;
//     }
//   }

//   // geometry
//   printf("Loading geometry...\n");
//   AliGeomManager::LoadGeometry();
//   if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD TOF PHOS EMCAL HMPID") ) {
//     //printf("Problem with align objects\n");
//   }
  fGeom =  AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1");
  printf("%s\n",fGeom->GetName());
  if(!fGeom)
    {
      printf("EMCAL geometry not loaded!\n");
      return kFALSE;
    }
  return kTRUE;
}

//________________________________________________________________________
void AliTrackComparisonESD::Exec(Option_t *) {
  //
  // Exec function
  // Loop over tracks and call  Process function
  if (!fESD) {
    //Printf("ERROR: fESD not available");
    return;
  }

  if(!SetupEvent()) return;


  fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
    //Printf("ERROR: fESDfriend not available");
    return;
  }


  if ( fESDfriend->GetNumberOfTracks() <=0 ) {
    //Printf("ERROR: fESDfriend Tracks not available");
    return;
  }


  //Get primary vertex
  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  if(!vertex) AliError("No primary vertex found!\n");
  Double_t vPos[3];
  vertex->GetXYZ(vPos);
  if(TMath::Abs(vPos[2])>7) return;

  
  //Get EMCAL clusters and cells
  TRefArray *clusters = new TRefArray();
  Int_t nclusters = fESD->GetEMCALClusters(clusters);
  AliESDCaloCells *cells = fESD->GetEMCALCells();
  RecalClusterPos(clusters,cells);
//   Float_t pos[3];
//   for(Int_t icl=0; icl<nclusters; icl++)
//     {
//       AliESDCaloCluster *cluster = (AliESDCaloCluster*) clusters->At(icl);
//       cluster->GetPosition(pos);
//       printf("cluster %d pass00 pos: (%5.3f,%5.3f,%5.3f,%5.3f)\n",icl,pos[0],pos[1],pos[2],TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]));
//     }



  //Loop over tracks
  Int_t ntracks = fESD->GetNumberOfTracks();
  for(Int_t itr=0; itr<ntracks; itr++)
    {
      AliESDtrack *track = fESD->GetTrack(itr);
      if(!track || !fESDCuts->AcceptTrack(track)) continue;

      AliESDfriendTrack *friendTrack = fESDfriend->GetTrack(itr);
      if(!friendTrack) continue;
      //printf(" --- %d < %d || %p | %p -- %p \n", itr, fESDfriend->GetNumberOfTracks(), track, fESDfriend, friendTrack);
      ProcessTOF(track,friendTrack,vPos);
      if(nclusters>0)ProcessEMCAL(track,friendTrack,clusters,vPos);
      ProcessHMPID(track,friendTrack,vPos);
    }//End of track loop

  delete clusters;
  PostData(0,fOutput);
}

//________________________________________________________________________
void AliTrackComparisonESD::ProcessTOF(AliESDtrack *track, AliESDfriendTrack *friendTrack, Double_t *vPos){
  //
  // Process TPC-TOF extrapolation
  //

  //printf("Enter function!\n");
  if (track->GetTOFsignal()<=0)  return;
  if (!friendTrack) return;
  if (!friendTrack->GetTPCOut()) return;

  AliExternalTrackParam *pTPC = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
  if(!pTPC) return;

  const AliTrackPointArray *points=friendTrack->GetTrackPointArray();
  if (!points) return;
  Int_t npoints = points->GetNPoints();
  if(npoints>1000) return; //the default value is more than 30000, why not -1???
  AliTrackPoint point;
  //
  Int_t counter=0;
  for (Int_t ipoint=0;ipoint<npoints;ipoint++){
    //cout<<"(npoints,ipoint) = ("<<npoints<<","<<ipoint<<")"<<endl;
    if(!points->GetPoint(point,ipoint)) continue;
    //Float_t xyz[3];
    //point.GetXYZ(xyz);
    //Float_t r=10;
    //Float_t r=TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);

    //printf("fVolumeID %d | LayerID %d\n",point.GetVolumeID(),AliGeomManager::VolUIDToLayer(point.GetVolumeID()));

    //if (r<350) continue;
    //if (r>400) continue;
    //cout<<"r="<<r<<endl;
    if(AliGeomManager::VolUIDToLayer(point.GetVolumeID())==AliGeomManager::kTOF)
      {
	counter++;
	fTOF->AddTracks(pTPC,&point,track->GetMass(),track->P(),vPos);
      }
  }
  //Printf("# of track points in TOF: %d!\n",counter);
}

//________________________________________________________________________
void AliTrackComparisonESD::ProcessEMCAL(AliESDtrack *track, AliESDfriendTrack *friendTrack, TRefArray *clusters, Double_t *vPos){
  if(clusters->GetEntriesFast()==0) return;

  Double_t rEMCal = 438;

  AliExternalTrackParam *pTPC = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
  if(!pTPC) return;

  Double_t trPos[3];
  Float_t clPos[3];
  AliExternalTrackParam *pTest = new AliExternalTrackParam(*pTPC);
  if(!AliTracker::PropagateTrackToBxByBz(pTest, rEMCal , track->GetMass(), 1 , kFALSE,0.99,-1)) return;
  if(!pTest->GetXYZ(trPos)) return;

  AliExternalTrackParam *p0=0;
  AliTrackPoint *p1=new AliTrackPoint();
  Int_t nclusters = clusters->GetEntries();
  for(Int_t icl=0; icl<nclusters; icl++)
    {
      AliESDCaloCluster *cluster = (AliESDCaloCluster*) clusters->At(icl);
      if(!cluster) continue;
      cluster->GetPosition(clPos);
      if( TMath::Abs(clPos[0]-trPos[0])>fCutX || TMath::Abs(clPos[1]-trPos[1])>fCutY || TMath::Abs(clPos[2]-trPos[2]>fCutZ) ) continue;
      
      p0 = pTPC;
      //printf("cluster pos: (%5.3f,%5.3f,%5.3f,%5.3f)\n",clPos[0],clPos[1],clPos[2],TMath::Sqrt(clPos[0]*clPos[0]+clPos[1]*clPos[1]));
      p1->SetXYZ(clPos[0],clPos[1],clPos[2],0);
      //printf("Found EMCAL point!\n");
      fEMCAL->AddTracks(p0,p1,track->GetMass(),cluster->E(),vPos);
    }

  delete pTest;
  delete p1;
}

//________________________________________________________________________
void AliTrackComparisonESD::ProcessHMPID(AliESDtrack *track, AliESDfriendTrack *friendTrack, Double_t *vPos){
  //
  // Process TPC-TOF extrapolation
  //
  if (track->GetHMPIDsignal()<=0)  return;

  AliExternalTrackParam *pTPC = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
  if(!pTPC) return;

  Int_t q, nph, ch;
  Float_t x, y;
  track->GetHMPIDmip(x,y,q,nph);
  Double_t pHmp[3]={0}, pHmp3=0;
  if (track->GetOuterHmpPxPyPz(pHmp)) 
    pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);

  ch = track->GetHMPIDcluIdx()/1000000;

  AliHMPIDParam *pParam = AliHMPIDParam::Instance(); 
  TVector3 vG = pParam->Lors2Mars(ch,x,y);

  AliTrackPoint *p1 = new AliTrackPoint();
  p1->SetXYZ(vG.X(),vG.Y(),vG.Z());
  //printf("Found HMPID point!\n");
  fHMPID->AddTracks(pTPC,p1,track->GetMass(),pHmp3,vPos);
  delete p1;
}


//________________________________________________________________________
void AliTrackComparisonESD::RecalClusterPos(TRefArray *clusters, AliESDCaloCells *cells){
  Double_t iLocal[3], iGlobal[3];
  Float_t cPos[3];
  //Float_t pos[3];
  Int_t nclusters = clusters->GetEntries();
  for(Int_t icl=0; icl<nclusters; icl++)
    {
      AliESDCaloCluster *cluster = (AliESDCaloCluster*) clusters->At(icl);
      UShort_t *absId = cluster->GetCellsAbsId();
      Int_t nCells = cluster->GetNCells();
      for(Int_t i=0;i<3;i++)cPos[i]=0;
      Double_t wTot=0;
      for(Int_t iCell=0; iCell<nCells; iCell++)
	{
	  Double_t cellEnergy = cells->GetCellAmplitude(absId[iCell]);
	  Double_t dist = 1.31*(log(cluster->E())+4.82+0.5);
	  fGeom->RelPosCellInSModule(absId[iCell],dist,iLocal[0],iLocal[1],iLocal[2]);
	  //fGeom->GetGlobal(iLocal,iGlobal,fGeom->GetSuperModuleNumber(absId[iCell]));
	  Int_t sm = fGeom->GetSuperModuleNumber(absId[iCell]);
	  //matrix[sm]->Print();
	  //cout<<"sm = "<<sm<<endl;
	  fTransMatrix[sm]->LocalToMaster(iLocal,iGlobal);

	  Double_t w = TMath::Max( 0., 4.5 + TMath::Log( cellEnergy / cluster->E() ));
	  if(w>0.0)
	    {
	      wTot += w;
	      for(Int_t i=0; i<3; i++ )
		  cPos[i] += (w*iGlobal[i]);
	    }
	}//End of cell loop   
      if(wTot>0)
	{
	  for(int i=0; i<3; i++ ) 
	    {
	      cPos[i] /= wTot;
	      cluster->SetPosition(cPos);
	    }
	}
      //cluster->GetPosition(pos);
      //printf("cluster %d pass10 pos: (%5.3f,%5.3f,%5.3f)\n",icl,pos[0],pos[1],pos[2]);
    }//End of cluster loop
}

//________________________________________________________________________
void AliTrackComparisonESD::Terminate(Option_t */*option*/) {
  //
  // Terminate
  //
  AliInfo("AliTrackComparisonESD::Terminate()\n");
  
}

//________________________________________________________________________
void AliTrackComparisonESD::FinishTaskOutput(){
  //
  // According description in AliAnalisysTask this method is call 
  // on the slaves before sending data
  //
  Terminate("slave");
  if(!fDebugOutputPath.IsNull()) { 
    RegisterDebugOutput();
  }
}

//________________________________________________________________________
Long64_t AliTrackComparisonESD::Merge(TCollection *li) {
  if(li) return 1;
  else return 1;
}

//________________________________________________________________________
void AliTrackComparisonESD::Analyze() {
  //
  // Analyze the content of the task
  //

}

//________________________________________________________________________
void AliTrackComparisonESD::RegisterDebugOutput(){
  //
  //
  //
}
