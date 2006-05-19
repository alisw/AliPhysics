#include "AliRICHTracker.h" //class header
#include "AliRICH.h"
#include "AliRICHRecon.h"
#include "AliRICHParam.h"
#include <AliESD.h>
#include <TGeoManager.h>     //EsdQA()
#include <TVector3.h>
#include <TTree.h>          //EsdQA() 
#include <TFile.h>          //EsdQA()  
#include <TProfile.h>   //EsdQA() 
#include "AliRICHHelix.h"
#include <AliMagF.h>
#include <AliStack.h>
#include <TParticle.h>
#include <TMath.h>
#include <AliRun.h>
#include <TNtupleD.h>            //RecWithStack();
#include <AliTrackPointArray.h>  //GetTrackPoint()
#include <AliAlignObj.h>         //GetTrackPoint()
#include <TH1F.h>                //EsdQA()  
#include <TH2F.h>                //EsdQA()  
#include <TCanvas.h>             //EsdQA()  
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
  AliPID pid; // needed to retrive all the PID info
   
  for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrk);// get next reconstructed track    
    Double_t mom[3], pos[3];
    pTrack->GetPxPyPz(mom); TVector3 mom3(mom[0],mom[1],mom[2]);
    pTrack->GetXYZ(pos); TVector3 pos3(pos[0],pos[1],pos[2]);
    AliRICHHelix helix(pos3,mom3,(Int_t)pTrack->GetSign(),-0.1*GetBz()); //construct helix out of track running parameters  
     //Printf(" magnetic field %f charged %f\n",GetBz(),pTrack->GetSign()); helix.Print("Track");
    Int_t iChamber=helix.RichIntersect(AliRICHParam::Instance());    
    if(!iChamber) continue;                                         //no intersection with chambers, ignore this track go after the next one
  
    //find MIP cluster candidate (closest to track intersection point cluster with large enough QDC)
    Double_t    dR=9999,   dX=9999,   dY=9999;                     //distance between track-PC intersection point and current cluster
    Double_t mipDr=9999,mipDx=9999,mipDy=9999,mipX=9999,mipY=9999; //nearest cluster parameters
    Int_t   iMipId=-1;                         //index of this nearest cluster
    for(Int_t iClu=0;iClu<pRich->Clus(iChamber)->GetEntries();iClu++){//clusters loop for intersected chamber
      AliRICHCluster *pClu=(AliRICHCluster*)pRich->Clus(iChamber)->UncheckedAt(iClu);//get pointer to current cluster
      if(pClu->Q()<AliRICHParam::QthMIP()) continue; //to low QDC, go after another one
      pClu->DistXY(helix.PosPc(),dX,dY); dR=TMath::Sqrt(dX*dX+dY*dY); //get distance for current cluster
      if(dR<mipDr){iMipId=iClu; mipDr=dR; mipDx=dX; mipDy=dY; mipX=pClu->X(); mipY=pClu->Y();} //current cluster is closer, overwrite data for min cluster
    }//clusters loop for intersected chamber

    pTrack->SetRICHthetaPhi(helix.Ploc().Theta(),helix.Ploc().Phi()); //store track impact angles with respect to RICH planes
    pTrack->SetRICHdxdy(mipDx,mipDy);                                 //distance between track-PC intersection and closest cluster with Qdc>100
    pTrack->SetRICHmipXY(mipX,mipY);                                  //position of that closest cluster with Qdc>100
    
    if(iMipId==-1)                        {pTrack->SetRICHsignal(kMipQdcCut);  continue;} //no cluster with enough QDC found
    if(mipDr>AliRICHParam::DmatchMIP())   {pTrack->SetRICHsignal(kMipDistCut); continue;} //closest cluster with enough carge is still too far 
  
    pTrack->SetRICHcluster(iMipId+1000000*iChamber);                                //set mip cluster index
    pTrack->SetRICHsignal(recon.ThetaCerenkov(&helix,pRich->Clus(iChamber),iMipId));//search for mean Cerenkov angle for this track
    pTrack->SetRICHchi2(recon.GetRingSigma2());
    pTrack->SetRICHnclusters(iMipId);                                               //on return iMipId is number of photon clusters accepted in reconstruction

    AliDebug(1,Form("Ch=%i PC Intersection=(%5.2f,%5.2f) cm MIP cluster dist=(%5.2f,%5.2f)=%5.2f cm ThetaCkov=%f",
                     iChamber,helix.PosPc().X(),helix.PosPc().Y(),            mipDx,mipDy,mipDr,     pTrack->GetRICHsignal()));
    
//here comes PID calculations    
//    CalcProb(pTrack->GetRICHsignal(),pTrack->GetP(),sigmaPID,richPID);
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
void AliRICHTracker::EsdQA(Bool_t isPrint)
{
// Reads a set of ESD files and print out or plot some information
// Arguments: isPrint is a flag to choose between printing (isPrint = kTRUE) and plotting (isPrint = kFALSE)
//   Returns: none
    
  TFile *pFile=TFile::Open("AliESDs.root","read");              if(!pFile) {Printf("ERROR: AliESDs.root does not exist!");return;} 
  TTree *pTr=(TTree*)pFile->Get("esdTree");                     if(!pTr)   {Printf("ERROR: AliESDs.root, no ESD tree inside!");return;} 
  AliESD *pEsd=new AliESD;  pTr->SetBranchAddress("ESD", &pEsd); 
  
  Int_t iNevt=pTr->GetEntries();  Printf("This ESD contains %i events",iNevt);
   //AliRICHParam *pPar;
  
  TH1D *pDx,*pDy,*pProbEl,*pProbMu,*pProbPi,*pProbKa,*pProbPr; TH2F *pThP; TProfile *pChiTh;
  if(!isPrint){
    TH1::AddDirectory(kFALSE);    
    pDx    =new TH1D("dX"  ,"distance between Xmip and Xtrack;cm",300,-1.5,1.5);
    pDy    =new TH1D("dY"  ,"distance between Ymip and Ytrack;cm",300,-1.5,1.5);
    pThP   =new TH2F("tvsp","#theta_{Ckov} radian;P GeV"              ,65 ,-0.5,6.0,75,0,0.75); pThP->SetStats(0);
    pProbPi=new TH1D("RichProbPion","HMPID PID probability for  e #mu #pi"   ,101,0.05,1.05); pProbPi->SetLineColor(kRed);
    pProbEl=new TH1D("RichProbEle" ,""   ,101,0.05,1.05);                                     pProbEl->SetLineColor(kGreen); 
    pProbMu=new TH1D("RichProbMuon" ,""  ,101,0.05,1.05);                                     pProbMu->SetLineColor(kBlue);
    pProbKa=new TH1D("RichProbKaon","HMPID PID probability for K"     ,101,0.05,1.05);
    pProbPr=new TH1D("RichProbProton","HMPID PID probability for p"   ,101,0.05,1.05);
    pChiTh =new TProfile("RichChiTh","#chi^{2};#theta_{C}"            ,80 ,0,0.8 , -2,2);
    
    //if(!gGeoManager)TGeoManager::Import("geometry.root");
    //pPar = AliRICHParam::Instance();
  }
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//ESD events loop
    pTr->GetEvent(iEvt);
    Int_t iNtracks=pEsd->GetNumberOfTracks(); Printf("ESD contains %i tracks created in Bz=%.2f Tesla",iNtracks,pEsd->GetMagneticField()/10.);
    for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){//ESD tracks loop
      AliESDtrack *pTrack = pEsd->GetTrack(iTrk);// get next reconstructed track
      Float_t dx,dy;         pTrack->GetRICHdxdy(dx,dy);
      Float_t theta,phi;     pTrack->GetRICHthetaPhi(theta,phi);
      TString comment;
      if(isPrint){
        if(pTrack->GetRICHsignal()>0)                         comment="OK";
        else if(pTrack->GetRICHsignal()==kMipQdcCut)          comment="no enough QDC";
        else if(pTrack->GetRICHsignal()==kMipDistCut)         comment="nearest cluster is too far";
        else if(pTrack->GetRICHsignal()==-1)                  comment="no intersection";
        Double_t pid[5]; pTrack->GetRICHpid(pid);           
        Printf("Tr %2i Q=%4.1f P=%.3f GeV Tr-Mip=%5.2f cm Th=%7.3f rad Prob : e=%.4f mu=%.4f pi=%.4f K=%.4f p=%.4f %s" ,
                     iTrk,pTrack->GetSign(),pTrack->GetP(),TMath::Sqrt(dx*dx+dy*dy),pTrack->GetRICHsignal(), 
                     pid[0],pid[1],pid[2],pid[3],pid[4], comment.Data());
      }else{//collect hists
    
        pDx->Fill(dx);  pDy->Fill(dy);
        pThP->Fill(pTrack->GetP(),pTrack->GetRICHsignal());
        pChiTh->Fill(pTrack->GetRICHsignal(),pTrack->GetRICHchi2());
        Double_t prob[5]; pTrack->GetRICHpid(prob);
        pProbEl->Fill(prob[0]); pProbMu->Fill(prob[1]); pProbPi->Fill(prob[2]);  pProbKa->Fill(prob[3]);  pProbPr->Fill(prob[4]); 
      }//if plot
    }//ESD tracks loop
  }//ESD events loop
  delete pEsd;  pFile->Close();//close AliESDs.root
  if(!isPrint){
    TCanvas *pC=new TCanvas("c","Quality",1200,1500); pC->Divide(2,2);
    TF1 *pPion = new TF1("RICHtheor","acos(sqrt(x*x+[0]*[0])/(x*[1]))",1.2,6);  pPion->SetLineWidth(1);pPion->SetParameter(1,1.292);
    AliPID ppp;                   pPion->SetParameter(0,AliPID::ParticleMass(AliPID::kPion))  ; pPion->SetLineColor(kRed);
    TF1 *pKaon = (TF1*)pPion->Clone();  pKaon->SetParameter(0,AliPID::ParticleMass(AliPID::kKaon))  ; pKaon->SetLineColor(kGreen);
    TF1 *pProt = (TF1*)pPion->Clone();  pProt->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)); pProt->SetLineColor(kBlue);
    
    
    pC->cd(1);        pDx->Draw();
    pC->cd(2);        pDy->Draw();
    pC->cd(3);        pThP->Draw();       pPion->Draw("same"); pKaon->Draw("same"); pProt->Draw("same");
    pC->cd(4);        pChiTh->Draw();
    
    TCanvas *pC2=new TCanvas("c2","Quality 2",1200,1500); pC2->Divide(2,2);
    pC2->cd(1);        pProbPi->Draw(); pProbMu->Draw("same"); pProbEl->Draw("same");
    pC2->cd(2);        pProbKa->Draw();
    pC2->cd(3);        pProbPr->Draw();
    
    
  }
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
        Float_t dx,dy;         pTrack->GetRICHdxdy(dx,dy);
        Float_t theta,phi;     pTrack->GetRICHthetaPhi(theta,phi);
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
