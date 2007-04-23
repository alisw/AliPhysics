#include "AliHMPIDTracker.h"     //class header
#include "AliHMPIDCluster.h"     //GetTrackPoint(),PropagateBack() 
#include "AliHMPIDParam.h"       //GetTrackPoint(),PropagateBack()
#include "AliHMPIDRecon.h"       //PropagateBack()
#include <AliESD.h>              //PropagateBack()  
#include <AliRun.h>              //GetTrackPoint(),PropagateBack()  
#include <AliTrackPointArray.h>  //GetTrackPoint()
#include <AliAlignObj.h>         //GetTrackPoint()
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
//.
// HMPID base class fo tracking
//.
//.
//.
ClassImp(AliHMPIDTracker)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDTracker::AliHMPIDTracker():AliTracker()
{
// ctor
  fClu=new TObjArray(AliHMPIDDigit::kMaxCh+1);
  fClu->SetOwner(kTRUE);
  for(int i=AliHMPIDDigit::kMinCh;i<=AliHMPIDDigit::kMaxCh;i++) fClu->AddAt(new TClonesArray("AliHMPIDCluster"),i);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
Bool_t AliHMPIDTracker::GetTrackPoint(Int_t idx, AliTrackPoint& point) const
{
// Interface callback methode invoked from AliReconstruction::WriteAlignmentData() to get position of MIP cluster in MARS associated to a current track.
// MIP cluster is reffered by index which is stored in AliESDtrack  ???????
// Arguments: idx- cluster index which is stored by HMPID in AliESDtrack
//            point- reference to the object where to store the point     
//   Returns: status of operation  if FALSE then AliReconstruction::WriteAlignmentData() do not store this point to array of points for current track. 
  if(idx<0) return kFALSE; //no MIP cluster assigned to this track in PropagateBack()
  Int_t iCham=idx/1000000;
  Int_t iClu=idx%1000000;
  point.SetVolumeID(AliAlignObj::LayerToVolUID(AliAlignObj::kHMPID,iCham-1));//layer and chamber number
  TClonesArray *pArr=(TClonesArray*)(*fClu)[iCham];
  AliHMPIDCluster *pClu=(AliHMPIDCluster*)pArr->UncheckedAt(iClu);//get pointer to cluster
  Double_t mars[3];
  AliHMPIDParam::Instance()->Lors2Mars(iCham,pClu->X(),pClu->Y(),mars);
  point.SetXYZ(mars[0],mars[1],mars[2]);
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::LoadClusters(TTree *pCluTree)
{
// Interface callback methode invoked from AliReconstruction::RunTracking() to load HMPID clusters before PropagateBack() gets control 
// Arguments: pCluTree- pointer to clusters tree got by AliHMPIDLoader::LoadRecPoints("read") then AliHMPIDLoader::TreeR()
//   Returns: error code (currently ignored in AliReconstruction::RunTraking())    
  AliDebug(1,"Start.");
  for(int i=AliHMPIDDigit::kMinCh;i<=AliHMPIDDigit::kMaxCh;i++) pCluTree->SetBranchAddress(Form("HMPID%d",i),&((*fClu)[i]));
  pCluTree->GetEntry(0);
  AliDebug(1,"Stop."); 
  return 0;  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::PropagateBack(AliESD *pEsd)
{
// This method defined as pure virtual in AliTracker. It is invoked from AliReconstruction::RunTracking() after invocation of AliTracker::LoadClusters()
// Agruments: pEsd - pointer to ESD
//   Returns: error code    
  AliCDBEntry *pNmeanEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Nmean"); //contains TObjArray of 21 TF1
  AliCDBEntry *pQthreEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Qthre"); //contains TObjArray of 7 TF1
  if(!pNmeanEnt) AliFatal("No Nmean C6F14 ");
  if(!pQthreEnt) AliFatal("No Qthre");
    
  return Recon(pEsd,fClu,(TObjArray*)pNmeanEnt->GetObject());  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::Recon(AliESD *pEsd,TObjArray *pCluAll,TObjArray *pNmean)
{
// Interface callback method invoked by AliRecontruction::RunTracking() during tracking after TOF. It's done just once per event
// Arguments: pEsd - pointer to Event Summary Data class instance which contains a list of tracks
//   Returns: error code, 0 if no errors   
  Int_t iNtracks=pEsd->GetNumberOfTracks();  AliDebugClass(1,Form("Start with %i tracks",iNtracks));
  
  AliHMPIDRecon recon;                                                                       //instance of reconstruction class, nothing important in ctor
  Double_t xPc,yPc;
  for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){                                                        //ESD tracks loop
    AliESDtrack *pTrk = pEsd->GetTrack(iTrk);                                                    //get next reconstructed track    
    Int_t cham=IntTrkCha(pTrk,xPc,yPc);                                                          //get chamber intersected by this track 
    if(cham<0) continue;                                                                         //no intersection at all, go after next track
    Double_t nmean=((TF1*)pNmean->At(3*cham))->Eval(pEsd->GetTimeStamp());                       //C6F14 Nmean for this chamber
    recon.SetImpPC(xPc,yPc);                                                                     //store track impact to PC
    recon.CkovAngle(pTrk,(TClonesArray *)pCluAll->At(cham),nmean);                               //search for Cerenkov angle of this track
  }                                                                                              //ESD tracks loop
  AliDebugClass(1,"Stop pattern recognition");
  return 0; // error code: 0=no error;
}//PropagateBack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::IntTrkCha(AliESDtrack *pTrk,Double_t &xToPc,Double_t &yToPc)
{
// Static method to find intersection in between given track and HMPID chambers
// Arguments: pTrk    - ESD track 
//            xPc,yPc - track intersection with PC, LORS
//   Returns: intersected chamber ID or -1
  
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  Float_t xRa=0,yRa=0,xPc=0,yPc=0,theta=0,phi=0;                                                //track intersection point and angles, LORS  
  for(Int_t i=AliHMPIDDigit::kMinCh;i<=AliHMPIDDigit::kMaxCh;i++){                              //chambers loop
    Double_t p1[3],n1[3]; pParam->Norm(i,n1); pParam->Lors2Mars(i,0,0,p1,AliHMPIDParam::kRad);  //point & norm  for RAD
    Double_t p2[3],n2[3]; pParam->Norm(i,n2); pParam->Lors2Mars(i,0,0,p2,AliHMPIDParam::kPc);   //point & norm  for PC
      
    if(pTrk->Intersect(p1,n1,-GetBz())==kFALSE) continue;                                       //try to intersect track with the middle of radiator
    if(pTrk->Intersect(p2,n2,-GetBz())==kFALSE) continue;                                       //try to intersect track with PC
      
    pParam->Mars2LorsVec(i,n1,theta,phi);                                                       //track angles at RAD
    pParam->Mars2Lors   (i,p1,xRa,yRa);                                                         //TRKxRAD position
    pParam->Mars2Lors   (i,p2,xPc,yPc);                                                         //TRKxPC position
    xToPc=(Double_t)xPc;yToPc=(Double_t)yPc;                                                    //conversion float->double only
    if(AliHMPIDDigit::IsInside(xPc,yPc,pParam->DistCut())==kFALSE) continue;             //not in active area  
    pTrk->SetHMPIDtrk      (xRa,yRa,theta,phi);                                                 //store track intersection info
    pTrk->SetHMPIDcluIdx   (i,0);
    return i;
  }                                                                                             //chambers loop
  return -1; //no intersection with HMPID chambers
}
