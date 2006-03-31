#include "AliRICHTracker.h" //class header
#include "AliRICH.h"
#include "AliRICHRecon.h"
#include <AliESD.h>
#include <TVector3.h>
#include <TTree.h>          //EsdPrint() 
#include <TFile.h>          //EsdPrint()    
#include "AliRICHHelix.h"
#include <AliMagF.h>
#include <AliStack.h>
#include <TParticle.h>
#include <TMath.h>
#include <AliRun.h>
#include <TNtupleD.h>            //RecWithStack();
#include <AliTrackPointArray.h>  //GetTrackPoint()
#include <AliAlignObj.h>         //GetTrackPoint()
ClassImp(AliRICHTracker)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliRICHTracker::AliRICHTracker():AliTracker()
{
// AliRICHTracker is created from AliReconstraction::Run() which invokes AliReconstraction::CreateTrackers() 
// which in turn invokes AliRICHReconstructor::CreateTracker(). 
// Note that this is done just once per session before AliReconstruction::Run() goes to events loop.
  AliRICHParam::Instance()->CdbRead(0,0); 
  for(Int_t i=0;i<5;i++)fErrPar[i]=0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
Bool_t AliRICHTracker::GetTrackPoint(Int_t idx, AliTrackPoint& point) const
{
// Interface callback methode invoked from AliReconstruction::WriteAlignmentData() to get position of MIP cluster in MARS associated to a current track.
// MIP cluster is reffered by index which is stored in AliESDtrack  ???????
// Arguments: idx- cluster index which is stored by RICH in AliESDtrack
//            point- reference to the object where to store the point     
//   Returns: status of operation  if FALSE then AliReconstruction::WriteAlignmentData() do not store this point to array of points for current track. 
  if(idx<0) return kFALSE; //no MIP cluster assigned to this track in PropagateBack()
  Int_t iCham=idx/1000000;
  Int_t iClu=idx%1000000;
  point.SetVolumeID(AliAlignObj::LayerToVolUID(AliAlignObj::kRICH,iCham-1));//layer and chamber number
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));  
  AliRICHCluster *pClu=(AliRICHCluster*)pRich->Clus(iCham)->UncheckedAt(iClu);//get pointer to cluster
  TVector3 mars=AliRICHParam::Instance()->Lors2Mars(iCham,pClu->X(),pClu->Y());
  point.SetXYZ(mars.X(),mars.Y(),mars.Z());
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHTracker::LoadClusters(TTree *pCluTree)
{
// Interface callback methode invoked from AliReconstruction::RunTracking() to load RICH clusters for RICH
// Arguments: pCluTree- pointer to clusters tree got by AliRICHLoader::LoadRecPoints("read") then AliRICHLoader::TreeR()
//   Returns: error code (currently ignored in AliReconstruction::RunTraking())    
  AliDebug(1,"Start.");  pCluTree->GetEntry(0);  AliDebug(1,"Stop."); return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHTracker::PropagateBack(AliESD *pESD)
{
// Interface callback methode invoked by AliRecontruction::RunTracking() during tracking after TOF. It's done just once per event
// Arguments: pESD - pointer to Event Summary Data class instance which contains a list of tracks
//   Returns: error code, 0 if no errors   
  Int_t iNtracks=pESD->GetNumberOfTracks();
  AliDebug(1,Form("Start with %i tracks",iNtracks));
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));  
  AliRICHRecon recon;
   
  for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrk);// get next reconstructed track    
    AliRICHHelix helix(pTrack->X3(),pTrack->P3(),(Int_t)pTrack->GetSign(),-GetBz()/10); //construct helix out of track running parameters  
     //Printf(" magnetic field %f charged %f\n",GetBz(),pTrack->GetSign()); helix.Print("Track");
    Int_t iChamber=helix.RichIntersect(AliRICHParam::Instance());    
    if(!iChamber) continue;                                         //no intersection with chambers, ignore this track go after the next one
  
    //find MIP cluster candidate (closest to track intersection point cluster with large enough QDC)
    Double_t    dR=9999,   dX=9999,   dY=9999; //distance between track-PC intersection point and current cluster
    Double_t dRmip=9999,dXmip=9999,dYmip=9999; //distance between track-PC intersection point and nearest cluster
    Int_t   iMipId=-1;                         //index of this nearest cluster
    for(Int_t iClu=0;iClu<pRich->Clus(iChamber)->GetEntries();iClu++){//clusters loop for intersected chamber
      AliRICHCluster *pClu=(AliRICHCluster*)pRich->Clus(iChamber)->UncheckedAt(iClu);//get pointer to current cluster
      if(pClu->Q()<AliRICHParam::QthMIP()) continue; //to low QDC, go after another one
      pClu->DistXY(helix.PosPc(),dX,dY); dR=TMath::Sqrt(dX*dX+dY*dY); //get distance for current cluster
      if(dR<dRmip){iMipId=iClu; dRmip=dR; dXmip=dX; dYmip=dY; }       //current cluster is closer, overwrite data for min cluster
    }//clusters loop for intersected chamber

    pTrack->SetRICHthetaPhi(helix.Ploc().Theta(),helix.Ploc().Phi()); //store track impact angles with respect to RICH planes
    pTrack->SetRICHdxdy(dXmip,dYmip);                                 //distance between track-PC intersection and closest cluster with Qdc>100
    
    if(iMipId==-1)                        {pTrack->SetRICHsignal(kMipQdcCut);  continue;} //no cluster with enough QDC found
    if(dRmip>AliRICHParam::DmatchMIP())   {pTrack->SetRICHsignal(kMipDistCut); continue;} //closest cluster with enough carge is still too far 
  
    pTrack->SetRICHcluster(iMipId+1000000*iChamber);                                //set mip cluster index
    pTrack->SetRICHsignal(recon.ThetaCerenkov(&helix,pRich->Clus(iChamber),iMipId));//search for mean Cerenkov angle for this track
    pTrack->SetRICHnclusters(iMipId);                                               //on return iMipId is number of photon clusters accepted in reconstruction

    AliDebug(1,Form("Ch=%i PC Intersection=(%5.2f,%5.2f) cm MIP cluster dist=(%5.2f,%5.2f)=%5.2f cm ThetaCkov=%f",
                     iChamber,helix.PosPc().X(),helix.PosPc().Y(),            dXmip,dYmip,dRmip,     pTrack->GetRICHsignal()));
    
//here comes PID calculations    
    if(pTrack->GetRICHsignal()>0) {
      AliDebug(1,Form("Start to assign the probabilities"));
      Double_t sigmaPID[AliPID::kSPECIES];
      Double_t richPID[AliPID::kSPECIES];
      for (Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) {
        sigmaPID[iPart] = 0;
        fErrPar[iPart] = 0;
        for(Int_t iphot=0;iphot<pRich->Clus(iChamber)->GetEntries();iphot++) {
          recon.SetPhotonIndex(iphot);
          if(recon.GetPhotonFlag() == 2) {
            Double_t theta_g=recon.GetTrackTheta();
            Double_t phi_g=(recon.GetPhiPoint()-recon.GetTrackPhi());
            Double_t sigma2 = AliRICHParam::Instance()->SigmaSinglePhoton(iPart,pTrack->GetP(),theta_g,phi_g).Mag2(); 
            if(sigma2>0) sigmaPID[iPart] += 1/sigma2;
          }
        }
        if (sigmaPID[iPart]>0)
          sigmaPID[iPart] *= (Double_t)(iMipId-recon.GetPhotBKG())/(Double_t)(iMipId); // n total phots, m are background...the sigma are scaled..
          if(sigmaPID[iPart]>0) sigmaPID[iPart] = 1/TMath::Sqrt(sigmaPID[iPart])*0.001; // sigma from parametrization are in mrad...
          else                  sigmaPID[iPart] = 0;
          fErrPar[iPart]=sigmaPID[iPart];
        AliDebug(1,Form("sigma for %s is %f rad",AliPID::ParticleName(iPart),sigmaPID[iPart]));
      }
      CalcProb(pTrack->GetRICHsignal(),pTrack->GetP(),sigmaPID,richPID);
      pTrack->SetRICHpid(richPID);         
      AliDebug(1,Form("PROBABILITIES ---> %f - %f - %f - %f - %f",richPID[0],richPID[1],richPID[2],richPID[3],richPID[4]));
    }    
  }//ESD tracks loop
  AliDebug(1,"Stop pattern recognition");
  return 0; // error code: 0=no error;
}//PropagateBack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHTracker::RecWithStack(TNtupleD *hn)
{
// Reconstruction for particles from STACK. This methode is to be used for RICH standalone when no other detectors are switched on, so normal tracking is not available.
// Arguments: hn- output ntuple where to store all variables
//   Returns: none       
  AliDebug(1,"Start.");  
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  
//  pRich->GetLoader()->GetRunLoader()->LoadHeader();
  if(!pRich->GetLoader()->GetRunLoader()->TreeK()) pRich->GetLoader()->GetRunLoader()->LoadKinematics();
  AliStack *pStack =   pRich->GetLoader()->GetRunLoader()->Stack();
  if(!pStack) {AliDebug(1,Form("No STACK found in AliRoot"));return;}
  Int_t iNtracks=pStack->GetNtrack();
  AliDebug(1,Form(" Start reconstruction with %i track(s) from Stack",iNtracks));
  
  Double_t hnvec[20];
  
  Double_t b=GetFieldMap()->SolenoidField()/10;// magnetic field in Tesla
  AliDebug(1,Form("Start with simulated %i tracks in %f Tesla field",iNtracks,b));
  TVector3 x0(0,0,0); TVector3 p0(0,0,0);//tmp storage for AliRICHHelix
  

  if(pRich->GetLoader()->LoadRecPoints()) {AliDebug(1,Form("No clusters found in RICH"));return;}
  pRich->GetLoader()->TreeR()->GetEntry(0);

  AliRICHRecon recon;
  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//stack particles loop
    TParticle *pParticle = pStack->Particle(iTrackN);
    if(!pParticle) {AliDebug(1,Form("Not a valid TParticle pointer. Track skipped"));continue;}
    AliDebug(1,Form(" PDG code : %i",pParticle->GetPdgCode()));
//
// problem of PDG code of some extra particles to be solved!!!!!!!!!
//
// found problem! Look in TRD directory : codes from Fluka are :
//
//    if ((pdg_code == 10010020) ||
//        (pdg_code == 10010030) ||
//        (pdg_code == 50000050) ||
//        (pdg_code == 50000051) ||
//        (pdg_code == 10020040)) {
//
    if(pParticle->GetPdgCode()>=50000050||pParticle->GetPdgCode()==0||pParticle->GetPdgCode()>10000) {AliDebug(1,Form("A photon as track... Track skipped"));continue;}
    if(!pParticle->GetPDG()) continue;
//
// to be updated for us!!
//
    AliDebug(1,Form("Track %i is a %s with charge %i and momentum %f",
            iTrackN,pParticle->GetPDG()->GetName(),Int_t(pParticle->GetPDG()->Charge()),pParticle->P()));
//    if(pParticle->GetMother(0)!=-1) continue; //consider only primaries
    if(pParticle->GetPDG()->Charge()==0||TMath::Abs(Int_t(pParticle->GetPDG()->Charge()))!=3) continue; //to avoid photons from stack...
    hnvec[0]=pParticle->P();
    hnvec[1]=pParticle->GetPDG()->Charge();
    
    p0.SetMagThetaPhi(pParticle->P(),pParticle->Theta(),pParticle->Phi());   x0.SetXYZ(pParticle->Vx(),pParticle->Vy(),pParticle->Vz());
    AliRICHHelix helix(x0,p0,TMath::Sign(1,(Int_t)pParticle->GetPDG()->Charge()),b);   
    Int_t iChamber=helix.RichIntersect(AliRICHParam::Instance());        
    if(!iChamber) continue;// no intersection with RICH found
    
    hnvec[2]=helix.Ploc().Theta();
    hnvec[3]=helix.Ploc().Phi();
    hnvec[4]=helix.PosPc().X();
    hnvec[5]=helix.PosPc().Y();
    
    Double_t dX,dY,dR,dRmip=9999;      //min distance between clusters and track position on PC 
    Int_t iMipId=-1;       //index of that min distance cluster 
    for(Int_t iClu=0;iClu<pRich->Clus(iChamber)->GetEntries();iClu++){//clusters loop for intersected chamber
      AliRICHCluster *pClu=(AliRICHCluster*)pRich->Clus(iChamber)->UncheckedAt(iClu);//get pointer to current cluster
      pClu->DistXY(helix.PosPc(),dX,dY); dR=TMath::Sqrt(dX*dX+dY*dY);//ditance between current cluster and helix intersection with PC
      if(dR<dRmip){dRmip=dR; hnvec[6]=pClu->X();hnvec[7]=pClu->Y();hnvec[8]=pClu->Q();
                                                  iMipId=1000000*iChamber+iClu;}//find cluster nearest to the track 
      
    }//clusters loop for intersected chamber
    
    
    hnvec[9]=recon.ThetaCerenkov(&helix,pRich->Clus(iChamber),iMipId); //search for mean Cerenkov angle for this track
    hnvec[10]=iMipId;//on return from ThetaCerenkov() contains number of photon candidates accepted
    hnvec[11]=(Double_t)iMipId;
    hnvec[12]=(Double_t)iChamber;
    hnvec[13]=(Double_t)pParticle->GetPdgCode();
    if(hn) hn->Fill(hnvec);
    AliDebug(1,Form("FINAL Theta Cerenkov=%f",hnvec[9]));
  }//stack particles loop
  
  pRich->GetLoader()->UnloadRecPoints();
  AliDebug(1,"Stop.");  
}//RecWithStack
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHTracker::CalcProb(Double_t thetaCer,Double_t pmod, Double_t *sigmaPID, Double_t *richPID)
{
// Calculates probability to be a electron-muon-pion-kaon-proton
// from the given Cerenkov angle and momentum assuming no initial particle composition
// (i.e. apriory probability to be the particle of the given sort is the same for all sorts)
  Double_t height[AliPID::kSPECIES];Double_t totalHeight=0;
  Double_t thetaTh[AliPID::kSPECIES];
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++){
    height[iPart]=0;
    Double_t mass = AliRICHParam::fgMass[iPart];
    Double_t refIndex=AliRICHParam::Instance()->IdxC6F14(AliRICHParam::EckovMean());
    Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(refIndex*pmod);
    thetaTh[iPart]=0;
    if(cosThetaTh>=1) continue;
    thetaTh[iPart] = TMath::ACos(cosThetaTh);
//    Double_t sinThetaThNorm = TMath::Sin(thetaTh)/TMath::Sqrt(1-1/(refIndex*refIndex));
//    Double_t sigmaThetaTh = (0.014*(1/sinThetaThNorm-1) + 0.0043)*1.25;
//    height[iPart] = TMath::Gaus(thetaCer,thetaTh,sigmaThetaTh);
    if(sigmaPID[iPart]>0) height[iPart] = TMath::Gaus(thetaCer,thetaTh[iPart],sigmaPID[iPart],kTRUE);
    else                  height[iPart] = 0;
    totalHeight +=height[iPart];
    AliDebugClass(1,Form(" Particle %s with mass %f with height %f and thetaTH %f",AliPID::ParticleName(iPart),mass,height[iPart],thetaTh[iPart]));
    AliDebugClass(1,Form(" partial height %15.14f total height %15.14f",height[iPart],totalHeight));
  }
  if(totalHeight<1e-5) {for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)richPID[iPart]=1.0/AliPID::kSPECIES;return;}
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) richPID[iPart] = height[iPart]/totalHeight;
  Int_t iPartNear = TMath::LocMax(AliPID::kSPECIES,richPID);
  if(TMath::Abs(thetaCer-thetaTh[iPartNear])>5.*sigmaPID[iPartNear]) for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)richPID[iPart]=1.0/AliPID::kSPECIES;
  //last line is to check if the nearest thetacerenkov to the teorethical one is within 5 sigma, otherwise no response (equal prob to every particle

}//CalcProb
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHTracker::EsdPrint()
{
// Reads a set of ESD files and print out some information
// Arguments: probCut - cut on probability 
//   Returns: none
    
  TFile *pFile=TFile::Open("AliESDs.root","read");              if(!pFile) {Printf("ERROR: AliESDs.root does not exist!");return;} 
  TTree *pTr=(TTree*)pFile->Get("esdTree");                     if(!pTr)   {Printf("ERROR: AliESDs.root, no ESD tree inside!");return;} 
  AliESD *pEsd=new AliESD;  pTr->SetBranchAddress("ESD", &pEsd); 
  
  Int_t iNevt=pTr->GetEntries();  Printf("This ESD contains %i events",iNevt);
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//ESD events loop
    pTr->GetEvent(iEvt);
    Int_t iNtracks=pEsd->GetNumberOfTracks(); Printf("ESD contains %i tracks created in Bz=%.2f Tesla",iNtracks,pEsd->GetMagneticField()/10.);
    for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){//ESD tracks loop
      AliESDtrack *pTrack = pEsd->GetTrack(iTrk);// get next reconstructed track
      Double_t dx,dy;         pTrack->GetRICHdxdy(dx,dy);
      Double_t theta,phi;     pTrack->GetRICHthetaPhi(theta,phi);
      Printf("Track %2i Q=%4.1f P=%.3f GeV RICH: ChamberCluster %7i Track-Mip=(%7.2f,%7.2f)=%5.2f cm ThetaCer %7.1f rad",iTrk,pTrack->GetSign(),pTrack->GetP(),
                   pTrack->GetRICHcluster(),dx,dy,TMath::Sqrt(dx*dx+dy*dy),pTrack->GetRICHsignal());
    }//ESD tracks loop
  }//ESD events loop
  delete pEsd;  pFile->Close();//close AliESDs.root
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHTracker::MatrixPrint(Double_t probCut)
{
// Reads a set of 3  ESD files from current directory and prints out the matrix of probabilities to pion kaon or proton completely blindly withou nay assumption on the contents of files.
// Normally it implies that those 3 ESDs contain only particles of the same sort namly pions, kaons and protons in that order.  
// Arguments: probCut - cut on probability 
//   Returns: none
  for(Int_t iFile=0;iFile<3;iFile++){
    TFile *pFile=TFile::Open(Form("Esd%1i.root",iFile+1),"read"); if(!pFile) {Printf("ERROR: Esd%1i.root does not exist!",iFile+1);return;} 
    TTree *pTr=(TTree*)pFile->Get("esdTree");                     if(!pTr)   {Printf("ERROR: Esd%1i.root, no ESD tree inside!",iFile+1);return;} 
    AliESD *pEsd=new AliESD;  pTr->SetBranchAddress("ESD", &pEsd); 
    Int_t iProtCnt=0,iKaonCnt=0,iPionCnt=0,iUnreconCnt=0,iTrkCnt=0; //counters
    
    for(Int_t iEvt=0;iEvt<pTr->GetEntries();iEvt++){//ESD events loop
      pTr->GetEvent(iEvt);
      iTrkCnt+=pEsd->GetNumberOfTracks();
      for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//ESD tracks loop
        AliESDtrack *pTrack = pEsd->GetTrack(iTrk);// get next reconstructed track
        Double_t dx,dy;         pTrack->GetRICHdxdy(dx,dy);
        Double_t theta,phi;     pTrack->GetRICHthetaPhi(theta,phi);
        Double_t prob[5];       pTrack->GetRICHpid(prob);
        if(pTrack->GetRICHsignal()>0){
          if(prob[4]>probCut)                         iProtCnt++; 
          if(prob[3]>probCut)                         iKaonCnt++;
          if((prob[0]+prob[1]+prob[2])>probCut)       iPionCnt++;
        } else
          iUnreconCnt++;       
      }//ESD tracks loop
      
    }//ESD events loop
    Printf("Bz=%5.2f Events=%i Total tracks=%i No recognized tracks=%i Pion=%i Kaon=%i Proton=%i ProbCut=%.2f",
        0.1*pEsd->GetMagneticField(),pTr->GetEntries(),iTrkCnt,iUnreconCnt,iPionCnt,iKaonCnt,iProtCnt,probCut);
    delete pEsd;  pFile->Close();//close AliESDs.root
  }//files loop
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
