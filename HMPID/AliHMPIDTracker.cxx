#include "AliHMPIDTracker.h"     //class header
#include "AliHMPIDtrack.h"       //class header
#include "AliHMPIDCluster.h"     //GetTrackPoint(),PropagateBack() 
#include "AliHMPIDParam.h"       //GetTrackPoint(),PropagateBack()
#include "AliHMPIDPid.h"         //Recon(),reconHTA()
#include "AliHMPIDRecon.h"       //Recon()
#include "AliHMPIDRecoParamV1.h"   //Recon()
#include "AliHMPIDReconstructor.h"//Recon()
#include "AliHMPIDReconHTA.h"    //ReconHTA()
#include <AliLog.h>              //Recon()  
#include <AliESDEvent.h>         //PropagateBack(),Recon()  
#include <AliESDtrack.h>         //Intersect() 
#include <AliTracker.h> 
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
AliHMPIDTracker::AliHMPIDTracker():
  AliTracker(),
  fClu(new TObjArray(AliHMPIDParam::kMaxCh+1))  
{
// ctor. Create TObjArray of TClonesArray of AliHMPIDCluster  
// 
//  
  fClu->SetOwner(kTRUE);
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
  iClu = iClu%1000; //GetHMPIDcluIdx -> 1e+6*ch + 1e+3*clusize + cluIdx;
  point.SetVolumeID(AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,iCham));//layer and chamber number
  TClonesArray *pArr=(TClonesArray*)(*fClu)[iCham];
  AliHMPIDCluster *pClu=(AliHMPIDCluster*)pArr->UncheckedAt(iClu);//get pointer to cluster
  Float_t xyz[3];
  pClu->GetGlobalXYZ(xyz);
  Float_t cov[6];
  pClu->GetGlobalCov(cov);
  point.SetXYZ(xyz,cov);
  return kTRUE;
}//GetTrackPoint()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::IntTrkCha(AliESDtrack *pTrk,Float_t &xPc,Float_t &yPc,Float_t &xRa,Float_t &yRa,Float_t &theta,Float_t &phi)
{
// Static method to find intersection in between given track and HMPID chambers
// Arguments: pTrk- ESD track; xPc,yPc- track intersection with PC in LORS [cm]
//   Returns: intersected chamber ID or -1
  AliHMPIDtrack *hmpTrk = new AliHMPIDtrack(*pTrk);                                             //create a hmpid track to be used for propagation and matching 
  for(Int_t i=AliHMPIDParam::kMinCh;i<=AliHMPIDParam::kMaxCh;i++){                              //chambers loop
    Int_t chInt = IntTrkCha(i,hmpTrk,xPc,yPc,xRa,yRa,theta,phi);
    if(chInt>=0) {delete hmpTrk;hmpTrk=0x0;return chInt;}
  }                                                                                             //chambers loop
  delete hmpTrk; hmpTrk=0x0;
  return -1;                                                                                    //no intersection with HMPID chambers
}//IntTrkCha()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::IntTrkCha(Int_t ch,AliHMPIDtrack *pTrk,Float_t &xPc,Float_t &yPc,Float_t &xRa,Float_t &yRa,Float_t &theta,Float_t &phi)
{
// Static method to find intersection in between given track and HMPID chambers
// Arguments: pTrk- HMPID track; xPc,yPc- track intersection with PC in LORS [cm]
//   Returns: intersected chamber ID or -1
    AliHMPIDParam *pParam=AliHMPIDParam::Instance();
    Double_t p1[3],n1[3]; 
    pParam->Norm(ch,n1); 
    pParam->Point(ch,p1,AliHMPIDParam::kRad);                                                    //point & norm  for middle of radiator plane
    Double_t p2[3],n2[3]; 
    pParam->Norm(ch,n2);
    pParam->Point(ch,p2,AliHMPIDParam::kPc);                                                     //point & norm  for entrance to PC plane
    if(pTrk->Intersect(p1,n1)==kFALSE) return -1;                                                //try to intersect track with the middle of radiator
    if(pTrk->Intersect(p2,n2)==kFALSE) return -1;   
    pParam->Mars2LorsVec(ch,n1,theta,phi);                                                       //track angles at RAD
    pParam->Mars2Lors   (ch,p1,xRa,yRa);                                                         //TRKxRAD position
    pParam->Mars2Lors   (ch,p2,xPc,yPc);                                                         //TRKxPC position
    if(AliHMPIDParam::IsInside(xPc,yPc,pParam->DistCut())==kTRUE) return ch;                     //return intersected chamber  
  return -1;                                                                                     //no intersection with HMPID chambers
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
  AliCDBEntry *pNmeanEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Nmean"); //contains TObjArray of 42 TF1 + 1 EPhotMean
  AliCDBEntry *pQthreEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Qthre"); //contains TObjArray of 42 (7ch * 6sec) TF1
  if(!pNmeanEnt) AliError("No Nmean C6F14 ");
  if(!pQthreEnt) AliError("No Qthre");
    
  return Recon(pEsd,fClu,(TObjArray*)pNmeanEnt->GetObject(),(TObjArray*)pQthreEnt->GetObject());  
}//PropagateBack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::Recon(AliESDEvent *pEsd,TObjArray *pClus,TObjArray *pNmean, TObjArray *pQthre)
{
// Static method to reconstruct Theta Ckov for all valid tracks of a given event.
// Arguments: pEsd- pointer ESD; pClu- pointer to clusters for all chambers; pNmean - pointer to all function Nmean=f(time)
//   Returns: error code, 0 if no errors   
  
  AliHMPIDRecon recon;                                                                           //instance of reconstruction class, nothing important in ctor
  AliHMPIDParam *pParam = AliHMPIDParam::Instance();                                             //Instance of AliHMPIDParam
  Float_t xPc,yPc,xRa,yRa,theta,phi;
  Double_t cluLORS[2]={0};
//  Double_t cluMARS[3]={0},trkMARS[3]={0};
//  Double_t bestcluMARS[3]={0,0,0};
//  Double_t radClu,radInitTrk;   
  Int_t nMipClusTot=0;
//  Double_t d3d=0;
  Double_t qthre = 0;   Double_t nmean=0; Int_t hvsec=0;
  Int_t nClusCh[AliHMPIDParam::kMaxCh+1];

  UInt_t tsmin = (UInt_t)((TF1*)pQthre->At(0))->GetXmin();                                        //
  UInt_t tsmax = (UInt_t)((TF1*)pQthre->At(0))->GetXmax();                                        //
  UInt_t ts = pEsd->GetTimeStamp();
  
  if(ts<tsmin || ts>tsmax) {
    AliWarning(Form(" in HMPID time stamp out of range!!! Please check!!! ts = %i",ts));
    return 1;
  }
   
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){                                        //loop on the ESD tracks in the event
           
//    Double_t bestChi2=99999;chi2=99999;                                                          //init. track matching params
    Double_t dmin=999999,bz=0,distCut=1,distParams[5]={1};

    Bool_t isOkDcut=kFALSE;
    Bool_t isOkQcut=kFALSE;
    Bool_t isMatched=kFALSE;
    
    AliHMPIDCluster *bestHmpCluster=0x0;                                                          //the best matching cluster
    AliESDtrack *pTrk = pEsd->GetTrack(iTrk);                                                     //get reconstructed track   
    
    if(!pTrk->IsOn(AliESDtrack::kTPCout)) continue;
 
    if(pTrk->IsOn(AliESDtrack::kTPCrefit)) continue;
 
    AliHMPIDtrack *hmpTrk = new AliHMPIDtrack(*pTrk);                                             //create a hmpid track to be used for propagation and matching 
    bz=AliTracker::GetBz();  
//initial flags for HMPID ESD infos    
    pTrk->SetHMPIDtrk(0,0,0,0);                                                                //no intersection found
    pTrk->SetHMPIDmip(0,0,0,0);                                                                //store mip info in any case 
    pTrk->SetHMPIDcluIdx(99,99999);                                                            //chamber not found, mip not yet considered
    pTrk->SetHMPIDsignal(AliHMPIDRecon::kNotPerformed);                                        //ring reconstruction not yet performed
    
    Int_t ipCh=IntTrkCha(pTrk,xPc,yPc,xRa,yRa,theta,phi);                                        //find the intersected chamber for this track 
    if(ipCh<0) {delete hmpTrk;hmpTrk=0x0;continue;}                                                                         //no intersection at all, go after next track

    pTrk->SetHMPIDtrk(xPc,yPc,theta,phi);                                                        //store initial infos
    pTrk->SetHMPIDcluIdx(ipCh,9999);                                                             //set chamber, index of cluster + cluster size
    
// track intersects the chamber ipCh: find the MIP          
    
    TClonesArray *pMipCluLst=(TClonesArray *)pClus->At(ipCh);                                   //get the list of clusters
    nMipClusTot = pMipCluLst->GetEntries();                                                     //total number of clusters in the given chamber
    nClusCh[ipCh] = nMipClusTot;
    
    if(nMipClusTot==0) {delete hmpTrk;hmpTrk=0x0;continue;}                                                                         
    
    Int_t index=-1;                                                                             //index of the "best" matching cluster
    
    for (Int_t iClu=0; iClu<nMipClusTot;iClu++) {                                               //clusters loop
      
      AliHMPIDCluster *pClu=(AliHMPIDCluster*)pMipCluLst->UncheckedAt(iClu);                    //get the cluster
// evaluate qThre
      if(pQthre->GetEntriesFast()==pParam->kMaxCh+1) {
        qthre=((TF1*)pQthre->At(pClu->Ch()))->Eval(ts);                                         //
      } else {                                                                                  // in the past just 1 qthre
        hvsec = pParam->InHVSector(pClu->Y());                                                  //  per chamber
        if(hvsec>=0){
          qthre=((TF1*)pQthre->At(6*ipCh+hvsec))->Eval(ts);                                     //
        }

       }                                                                                            //
//
      if(pClu->Q()<qthre) continue;                                                                      //charge compartible with MIP clusters      
      isOkQcut = kTRUE;
//
      cluLORS[0]=pClu->X(); cluLORS[1]=pClu->Y();                                            //get the LORS coordinates of the cluster
      Double_t dist = TMath::Sqrt((xPc-cluLORS[0])*(xPc-cluLORS[0])+(yPc-cluLORS[1])*(yPc-cluLORS[1]));
      
      if(dist<dmin) {
        dmin = dist;
        index=iClu;
        bestHmpCluster=pClu;
      }
    } // clusters loop

    if(!isOkQcut) {
      pTrk->SetHMPIDsignal(pParam->kMipQdcCut);
      delete hmpTrk;hmpTrk=0x0; continue;                                                                     
    }

    Double_t radius = (pParam->Lors2Mars(ipCh,pParam->SizeAllX()/2,pParam->SizeAllY()/2)).Mag(); 
    
    if(!AliTracker::PropagateTrackToBxByBz(hmpTrk,radius,pTrk->GetMass(),1,kFALSE)) {delete hmpTrk;hmpTrk=0x0;continue;}
              
    if(!hmpTrk->PropagateTo(bestHmpCluster)) {delete hmpTrk;hmpTrk=0x0;continue;}

    Int_t cluSiz = bestHmpCluster->Size();
    pTrk->SetHMPIDmip(bestHmpCluster->X(),bestHmpCluster->Y(),(Int_t)bestHmpCluster->Q(),0);  //store mip info in any case 
    pTrk->SetHMPIDcluIdx(ipCh,index+1000*cluSiz);                                             //set chamber, index of cluster + cluster size

    if(AliHMPIDReconstructor::GetRecoParam())                                                 //retrieve distance cut
    {
      if(AliHMPIDReconstructor::GetRecoParam()->IsFixedDistCut()==kTRUE)                      //distance cut is fixed number
      { 
        distCut=AliHMPIDReconstructor::GetRecoParam()->GetHmpTrackMatchingDist();
      }
      else 
      {
        for(Int_t ipar=0;ipar<5;ipar++) distParams[ipar]=AliHMPIDReconstructor::GetRecoParam()->GetHmpTrackMatchingDistParam(ipar);      //prevision: distance cut is function of momentum
        distCut=distParams[0]+distParams[1]*TMath::Power(distParams[2]*pTrk->GetP(),distParams[3]); //prevision: change functional form to be more precise
      }
    }
    else {distCut=pParam->DistCut();}
      
    if(dmin < distCut) {
      isOkDcut = kTRUE;
    }   
    
    if(!isOkDcut) {
      pTrk->SetHMPIDsignal(pParam->kMipDistCut);                                                //closest cluster with enough charge is still too far from intersection
    }
    
    if(isOkQcut*isOkDcut) isMatched = kTRUE;                                                    // MIP-Track matched !!    
    
    if(!isMatched) {delete hmpTrk;hmpTrk=0x0;continue;}                                           // If matched continue...

    Bool_t isOk = hmpTrk->Update(bestHmpCluster,0.1,0);
    if(!isOk) {delete hmpTrk;hmpTrk=0x0;continue;}
    pTrk->SetOuterHmpParam(hmpTrk,AliESDtrack::kHMPIDout);                 

    FillResiduals(hmpTrk,bestHmpCluster,kFALSE);
 
    
    //evaluate nMean
    if(pNmean->GetEntries()==21) {                                                              //for backward compatibility
      nmean=((TF1*)pNmean->At(3*ipCh))->Eval(ts);                                               //C6F14 Nmean for this chamber
    } else {
      Int_t iRad     = pParam->Radiator(yRa);                                                   //evaluate the radiator involved
      if(iRad < 0) {
	nmean = -1;
      } else {
      Double_t tLow  = ((TF1*)pNmean->At(6*ipCh+2*iRad  ))->Eval(ts);                           //C6F14 low  temp for this chamber
      Double_t tHigh = ((TF1*)pNmean->At(6*ipCh+2*iRad+1))->Eval(ts);                           //C6F14 high temp for this chamber
      Double_t tExp  = pParam->FindTemp(tLow,tHigh,yRa);                                        //estimated temp for that chamber at that y
      nmean = pParam->NIdxRad(AliHMPIDParam::Instance()->GetEPhotMean(),tExp);                  //mean ref idx @ a given temp
      }
      if(nmean < 0){                                                                            //track didn' t pass through the radiator
         pTrk->SetHMPIDsignal(AliHMPIDRecon::kNoRad);                                           //set the appropriate flag
         pTrk->SetHMPIDcluIdx(ipCh,index+1000*cluSiz);                                          //set index of cluster
         delete hmpTrk;hmpTrk=0x0; 
         continue;
      }
    }
    //
    recon.SetImpPC(xPc,yPc);                                                                     //store track impact to PC
    recon.CkovAngle(pTrk,(TClonesArray *)pClus->At(ipCh),index,nmean,xRa,yRa);                   //search for Cerenkov angle of this track

    if(pTrk->GetHMPIDsignal()<0) {delete hmpTrk;hmpTrk=0x0;continue;}
        
    AliHMPIDPid pID;
    Double_t prob[5];
    pID.FindPid(pTrk,5,prob);
    pTrk->SetHMPIDpid(prob);
//      Printf(" Prob e- %6.2f mu %6.2f pi %6.2f k %6.2f p %6.2f",prob[0]*100,prob[1]*100,prob[2]*100,prob[3]*100,prob[4]*100);
    delete hmpTrk;hmpTrk=0x0;
  }//iTrk

  return 0; // error code: 0=no error;
}//Recon()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDTracker::ReconHiddenTrk(AliESDEvent *pEsd,TObjArray *pClus,TObjArray *pNmean, TObjArray *pQthre)
{
// Static method to reconstruct Theta Ckov for all valid tracks of a given event.
// Arguments: pEsd- pointer ESD; pClu- pointer to clusters for all chambers; pNmean - pointer to all function Nmean=f(time), pQthre - pointer to all function Qthre=f(time)
//   Returns: error code, 0 if no errors
  
  AliHMPIDReconHTA reconHTA;                                                                     //instance of reconstruction class, nothing important in ctor
  
  AliHMPIDParam *pParam = AliHMPIDParam::Instance();                                             //Instance of AliHMPIDParam
  
  UInt_t tsmin = (UInt_t)((TF1*)pQthre->At(0))->GetXmin();                                        //
  UInt_t tsmax = (UInt_t)((TF1*)pQthre->At(0))->GetXmax();                                        //
  UInt_t ts = pEsd->GetTimeStamp();

  if(ts<tsmin || ts>tsmax) {
    AliWarning(Form(" in HMPID time stamp out of range!!! Please check!!! ts = %i",ts));
    return 1;
  }
   
  for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){                                        //loop on the ESD tracks in the event
    
    AliESDtrack *pTrk = pEsd->GetTrack(iTrk);                                                     //here it is simulated or just empty track
    Int_t ipCh;
    ipCh = pTrk->GetHMPIDcluIdx();ipCh/=1000000;
    if(ipCh<0) continue;

    TClonesArray *pMipCluLst=(TClonesArray *)pClus->At(ipCh);                                   //get the list of clusters
    Int_t nMipClusTot = pMipCluLst->GetEntries();                                               //total number of clusters in the given chamber
    
    Double_t qMip=-1;
    Int_t chMip=-1;    
    Double_t xMip = 0;
    Double_t yMip = 0;
    Int_t indMip  = 0;
    Int_t cluMipSiz = 0;

    for (Int_t iClu=0; iClu<nMipClusTot;iClu++) {                                               //clusters loop
      
      AliHMPIDCluster *pClu=(AliHMPIDCluster*)pMipCluLst->UncheckedAt(iClu);                    //get the cluster
      Double_t qClus = pClu->Q();
      if(qClus>qMip) {
        qMip  = qClus;
        chMip = pClu->Ch();
        xMip = pClu->X();
        yMip = pClu->Y();
        indMip = iClu;
        cluMipSiz = pClu->Size();
      }
    }//clus loop

    if(chMip<0) return 1;    
    
    Int_t hvsec;
    Double_t qthre=0;
// evaluate qThre
    if(pQthre->GetEntriesFast()==pParam->kMaxCh+1) {                                              // just for backward compatibility
      qthre=((TF1*)pQthre->At(chMip))->Eval(ts);                                                  //
    } else {                                                                                      // in the past just 1 qthre
      hvsec = pParam->InHVSector(yMip);                                                           //  per chamber
      if(hvsec>=0) qthre=((TF1*)pQthre->At(6*chMip+hvsec))->Eval(ts);                             //
    }
//
    if(qMip<qthre) {
      pTrk->SetHMPIDmip(xMip,yMip,(Int_t)qMip,0);                                                 //store mip info in any case 
      pTrk->SetHMPIDcluIdx(chMip,indMip+1000*cluMipSiz);                                                          
      pTrk->SetHMPIDsignal(pParam->kMipQdcCut);
      return 1;                                                                                   //charge compatible with MIP clusters      
    }
      
    pTrk->SetHMPIDmip(xMip,yMip,(Int_t)qMip,0);                                                   //store mip info in any case 
    pTrk->SetHMPIDcluIdx(chMip,indMip+1000*cluMipSiz);                                                          

    Double_t yRa = yMip;                                                                        //just an approx...
    Double_t nmean;
    //evaluate nMean
    if(pNmean->GetEntries()==21) {                                                              //for backward compatibility
      nmean=((TF1*)pNmean->At(3*chMip))->Eval(ts);                                              //C6F14 Nmean for this chamber
    } else {
      Int_t iRad     = pParam->Radiator(yRa);                                                   //evaluate the radiator involved
      if(iRad < 0) {
	nmean = -1;
      } else {
      Double_t tLow  = ((TF1*)pNmean->At(6*chMip+2*iRad  ))->Eval(ts);                          //C6F14 low  temp for this chamber
      Double_t tHigh = ((TF1*)pNmean->At(6*chMip+2*iRad+1))->Eval(ts);                          //C6F14 high temp for this chamber
      Double_t tExp  = pParam->FindTemp(tLow,tHigh,yRa);                                        //estimated temp for that chamber at that y
      nmean = pParam->NIdxRad(AliHMPIDParam::Instance()->GetEPhotMean(),tExp);                  //mean ref idx @ a given temp
      }
      if(nmean < 0){                                                                            //track didn' t pass through the radiator
         pTrk->SetHMPIDsignal(AliHMPIDRecon::kNoRad);                                           //set the appropriate flag
         return 1;
      }
    }
    //
    if(!reconHTA.CkovHiddenTrk(pTrk,(TClonesArray *)pClus->At(ipCh),indMip,nmean)) {                 //search for track parameters and Cerenkov angle of this track
      AliHMPIDPid pID;
      Double_t prob[5];
      pID.FindPid(pTrk,5,prob);
      pTrk->SetHMPIDpid(prob);
    }
//      Printf(" Prob e- %6.2f mu %6.2f pi %6.2f k %6.2f p %6.2f",prob[0]*100,prob[1]*100,prob[2]*100,prob[3]*100,prob[4]*100);
  }//iTrk

  return 0; // error code: 0=no error;

}//Recon()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDTracker::FillClusterArray(TObjArray* array) const {
  
 // Publishes all pointers to clusters known to the tracker into the
  // passed object array.
  // The ownership is not transfered - the caller is not expected to delete
  // the clusters
 
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){    
    TClonesArray *pCluArr=(TClonesArray*)(*fClu)[iCh];
    for (Int_t iClu=0; iClu<pCluArr->GetEntriesFast();iClu++){
      AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluArr->UncheckedAt(iClu);    
      array->AddLast(pClu);
    }//cluster loop in iCh
    pCluArr->Delete();
  }//Ch loop
    
  return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
