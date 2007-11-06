#include "AliHMPIDTracker.h"     //class header
#include "AliHMPIDCluster.h"     //GetTrackPoint(),PropagateBack() 
#include "AliHMPIDParam.h"       //GetTrackPoint(),PropagateBack()
#include "AliHMPIDRecon.h"       //Recon()
#include "AliHMPIDReconHTA.h"    //ReconHTA()
#include <AliESDEvent.h>         //PropagateBack(),Recon()  
#include <AliESDtrack.h>         //Intersect()  
#include <AliRun.h>              //GetTrackPoint(),PropagateBack()  
#include <AliTrackPointArray.h>  //GetTrackPoint()
#include <AliAlignObj.h>         //GetTrackPoint()
#include <AliCDBManager.h>       //PropageteBack()
#include <AliCDBEntry.h>         //PropageteBack()
//.
// HMPID base class fo tracking
//.
//.
//.
ClassImp(AliHMPIDTracker)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDTracker::AliHMPIDTracker():AliTracker()
{
// ctor. Create TObjArray of TClonesArray of AliHMPIDCluster  
// 
//  
  fClu=new TObjArray(AliHMPIDParam::kMaxCh+1);  fClu->SetOwner(kTRUE);
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++) fClu->AddAt(new TClonesArray("AliHMPIDCluster"),i);
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
Bool_t AliHMPIDTracker::GetTrackPoint(Int_t idx, AliTrackPoint& point) const
{
// Interface callback methode invoked from AliReconstruction::WriteAlignmentData() to get position of MIP cluster in MARS associated to a current track.
// MIP cluster is reffered by index which is stored in AliESDtrack  ???????
// Arguments: idx- cluster index which is stored by HMPID in AliESDtrack
//            point- reference to the object where to store the point     
//   Returns: status of operation  if FALSE then AliReconstruction::WriteAlignmentData() do not store this point to array of points for current track. 
  if(idx<0) return kFALSE; //no MIP cluster assigned to this track in PropagateBack()
  Int_t iCham=idx/1000000; Int_t iClu=idx%1000000;
  point.SetVolumeID(AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,iCham-1));//layer and chamber number
  TClonesArray *pArr=(TClonesArray*)(*fClu)[iCham];
  AliHMPIDCluster *pClu=(AliHMPIDCluster*)pArr->UncheckedAt(iClu);//get pointer to cluster
  Double_t mars[3];
  AliHMPIDParam::Instance()->Lors2Mars(iCham,pClu->X(),pClu->Y(),mars);
  point.SetXYZ(mars[0],mars[1],mars[2]);
  return kTRUE;
}//GetTrackPoint()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::IntTrkCha(AliESDtrack *pTrk,Float_t &xPc,Float_t &yPc)
{
// Static method to find intersection in between given track and HMPID chambers
// Arguments: pTrk- ESD track; xPc,yPc- track intersection with PC in LORS [cm]
//   Returns: intersected chamber ID or -1
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  Float_t xRa=0,yRa=0,theta=0,phi=0;                                                            //track intersection at PC and angles at RAD, LORS  
  for(Int_t i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){                              //chambers loop
    Double_t p1[3],n1[3]; pParam->Norm(i,n1); pParam->Point(i,p1,AliHMPIDParam::kRad);          //point & norm  for middle of radiator plane
    Double_t p2[3],n2[3]; pParam->Norm(i,n2); pParam->Point(i,p2,AliHMPIDParam::kPc);           //point & norm  for entrance to PC plane
    if(pTrk->Intersect(p1,n1,-GetBz())==kFALSE) continue;                                       //try to intersect track with the middle of radiator
    if(pTrk->Intersect(p2,n2,-GetBz())==kFALSE) continue;                                       //try to intersect track with PC
    pParam->Mars2LorsVec(i,n1,theta,phi);                                                       //track angles at RAD
    pParam->Mars2Lors   (i,p1,xRa,yRa);                                                         //TRKxRAD position
    pParam->Mars2Lors   (i,p2,xPc,yPc);                                                         //TRKxPC position
    if(AliHMPIDParam::IsInside(xPc,yPc,pParam->DistCut())==kFALSE) continue;                    //not in active area  
    pTrk->SetHMPIDtrk      (xRa,yRa,theta,phi);                                                 //store track intersection info
    pTrk->SetHMPIDcluIdx   (i,0);
    return i;
  }                                                                                             //chambers loop
  return -1; //no intersection with HMPID chambers
}//IntTrkCha()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::LoadClusters(TTree *pCluTree)
{
// Interface callback methode invoked from AliReconstruction::RunTracking() to load HMPID clusters before PropagateBack() gets control. Done once per event.
// Arguments: pCluTree- pointer to clusters tree got by AliHMPIDLoader::LoadRecPoints("read") then AliHMPIDLoader::TreeR()
//   Returns: error code (currently ignored in AliReconstruction::RunTraking())    
  for(int i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++) pCluTree->SetBranchAddress(Form("HMPID%d",i),&((*fClu)[i]));
  pCluTree->GetEntry(0);
  return 0;  
}//LoadClusters()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::PropagateBack(AliESDEvent *pEsd)
{
// Interface pure virtual in AliTracker. Invoked from AliReconstruction::RunTracking() after invocation of AliTracker::LoadClusters() once per event
// Agruments: pEsd - pointer to ESD
//   Returns: error code    
  AliCDBEntry *pNmeanEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Nmean"); //contains TObjArray of 21 TF1
  AliCDBEntry *pQthreEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Qthre"); //contains TObjArray of 7 TF1
  if(!pNmeanEnt) AliFatal("No Nmean C6F14 ");
  if(!pQthreEnt) AliFatal("No Qthre");
    
  return Recon(pEsd,fClu,(TObjArray*)pNmeanEnt->GetObject(),(TObjArray*)pQthreEnt->GetObject());  
}//PropagateBack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::Recon(AliESDEvent *pEsd,TObjArray *pClus,TObjArray *pNmean, TObjArray *pQthre)
{
// Static method to reconstruct Theta Ckov for all valid tracks of a given event.
// Arguments: pEsd- pointer ESD; pClu- pointer to clusters for all chambers; pNmean - pointer to all function Nmean=f(time)
//   Returns: error code, 0 if no errors   
  AliHMPIDRecon recon;                                                                       //instance of reconstruction class, nothing important in ctor
  Float_t xPc,yPc;
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){                                       //ESD tracks loop
    AliESDtrack *pTrk = pEsd->GetTrack(iTrk);                                                    //get reconstructed track    
    Int_t cham=IntTrkCha(pTrk,xPc,yPc);                                                          //get chamber intersected by this track 
    if(cham<0) continue;                                                                         //no intersection at all, go after next track
    Double_t nmean=((TF1*)pNmean->At(3*cham))->Eval(pEsd->GetTimeStamp());                       //C6F14 Nmean for this chamber
    Double_t qthre=((TF1*)pQthre->At(cham))  ->Eval(pEsd->GetTimeStamp());                       //Qthre for this chamber
    recon.SetImpPC(xPc,yPc);                                                                     //store track impact to PC
    recon.CkovAngle(pTrk,(TClonesArray *)pClus->At(cham),nmean,qthre);                           //search for Cerenkov angle of this track
//    Printf("AliHMPIDTracker::Recon: nmean %f, qthre %f",nmean,qthre);
  }                                                                                              //ESD tracks loop
  return 0; // error code: 0=no error;
}//Recon()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::ReconHiddenTrk(Int_t iCh,AliESDtrack *pTrk,TClonesArray *pCluLst,TObjArray *pNmean,TObjArray *pQthre)
{
// Static method to reconstruct Theta Ckov for all valid tracks of a given event.
// Arguments: pEsd- pointer ESD; pClu- pointer to clusters for all chambers; pNmean - pointer to all function Nmean=f(time), pQthre - pointer to all function Qthre=f(time)
//   Returns: error code, 0 if no errors
  AliHMPIDReconHTA reconHTA;                                                                          //instance of reconstruction class, nothing important in ctor
  Double_t nmean=((TF1*)pNmean->At(3*iCh))->Eval(0);                                            //C6F14 Nmean for this chamber
  Double_t qthre=((TF1*)pQthre->At(iCh))  ->Eval(0);                                            //C6F14 Nmean for this chamber
  if(pCluLst->GetEntriesFast()<4) return 1;                                                     //min 4 clusters (3 + 1 mip) to find a ring! 
  if(reconHTA.CkovHiddenTrk(pTrk,pCluLst,nmean,qthre)) return 0;                                   //search for track parameters and Cerenkov angle of this track
  else return 1;                                                                                // error code: 0=no error,1=fit not performed;
}//Recon()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
