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
  Int_t iNtracks=pESD->GetNumberOfTracks();  AliDebug(1,Form("Start with %i tracks",iNtracks));
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));  
  AliRICHRecon recon;
  Int_t nphots =0;
   
  for(Int_t iTrk=0;iTrk<iNtracks;iTrk++){//ESD tracks loop
    AliESDtrack *pTrk = pESD->GetTrack(iTrk);// get next reconstructed track    
    Double_t mom[3], pos[3];
    pTrk->GetPxPyPz(mom); TVector3 mom3(mom[0],mom[1],mom[2]);
    pTrk->GetXYZ(pos);    TVector3 pos3(pos[0],pos[1],pos[2]);
    AliRICHHelix helix(pos3,mom3,(Int_t)pTrk->GetSign(),-0.1*GetBz()); //construct helix out of track running parameters  
     //Printf(" magnetic field %f charged %f\n",GetBz(),pTrack->GetSign()); helix.Print("Track");
    Int_t iChamber=helix.RichIntersect(AliRICHParam::Instance());    
    if(!iChamber) continue;                                         //no intersection with chambers, ignore this track go after the next one
  
    //find MIP cluster candidate (closest to track intersection point cluster with large enough QDC)
    Double_t    dR=9999,   dX=9999,   dY=9999;                                             //distance between track-PC intersection point and current cluster
    Double_t mipDr=9999,mipDx=9999,mipDy=9999,mipX=9999,mipY=9999; Int_t mipQ=0;           //nearest cluster parameters
    Int_t   iMipId=-1;                                                                     //index of this nearest cluster
    for(Int_t iClu=0;iClu<pRich->Clus(iChamber)->GetEntries();iClu++){                     //clusters loop for intersected chamber
      AliRICHCluster *pClu=(AliRICHCluster*)pRich->Clus(iChamber)->UncheckedAt(iClu);      //get pointer to current cluster
      if(pClu->Q()<AliRICHParam::QthMIP()) continue;                                       //to low QDC, go after another one
      pClu->DistXY(helix.PosPc(),dX,dY); dR=TMath::Sqrt(dX*dX+dY*dY);                      //get distance for current cluster
      if(dR<mipDr){iMipId=iClu; mipDr=dR; mipDx=dX; mipDy=dY; mipX=pClu->X(); mipY=pClu->Y(); mipQ=pClu->Q();} //current cluster is closer, overwrite data for min cluster
    }//clusters loop for intersected chamber

    pTrk->SetRICHthetaPhi(helix.Ploc().Theta(),helix.Ploc().Phi()); //store track impact angles with respect to RICH planes
    pTrk->SetRICHdxdy(mipDx,mipDy);                                 //distance between track-PC intersection and closest cluster with Qdc>100
    pTrk->SetRICHmipXY(mipX,mipY);                                  //position of that closest cluster with Qdc>100
    pTrk->SetRICHnclusters(1000000*mipQ);                           //charge of that closest cluster with Qdc>100 
    
    if(iMipId==-1)                        {pTrk->SetRICHsignal(kMipQdcCut);  continue;} //no cluster with enough QDC found
    if(mipDr>AliRICHParam::DmatchMIP())   {pTrk->SetRICHsignal(kMipDistCut); continue;} //closest cluster with enough carge is still too far 
  
    pTrk->SetRICHcluster(iMipId+1000000*iChamber);                                //set mip cluster index
    pTrk->SetRICHsignal(recon.ThetaCerenkov(&helix,pRich->Clus(iChamber),nphots));//search for mean Cerenkov angle for this track
    pTrk->SetRICHnclusters(1000000*mipQ+nphots);                                  //on return nphots is number of photon clusters accepted in reconstruction
    pTrk->SetRICHchi2(recon.GetRingSigma2());

    AliDebug(1,Form("Ch=%i PC Intersection=(%5.2f,%5.2f) cm MIP cluster dist=(%5.2f,%5.2f)=%5.2f cm ThetaCkov=%f",
                     iChamber,helix.PosPc().X(),helix.PosPc().Y(),            mipDx,mipDy,mipDr,     pTrk->GetRICHsignal()));
    
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
// Reads ESD file and print out or plot some information for QA
// Arguments: isPrint is a flag to choose between printing (isPrint = kTRUE) and plotting (isPrint = kFALSE)
//   Returns: none
    
  TFile *pFile=TFile::Open("AliESDs.root","read");              if(!pFile) {Printf("ERROR: AliESDs.root does not exist!");return;} 
  TTree *pTr=(TTree*)pFile->Get("esdTree");                     if(!pTr)   {Printf("ERROR: AliESDs.root, no ESD tree inside!");return;} 
  AliESD *pEsd=new AliESD;  pTr->SetBranchAddress("ESD", &pEsd); 
  
  TH1D     *pProbEl=0,*pProbMu=0,*pProbPi=0,*pProbKa=0,*pProbPr=0,*pMom=0,*pMipQ=0; 
  TH2F     *pThP=0,*pDxDy=0; 
  TProfile *pChiTh=0;
  if(!isPrint){
    TH1::AddDirectory(kFALSE);    
    pProbEl=new TH1D("RiProbE" ,"Prob e"       ,101  ,0   ,1.05); pProbEl->SetLineColor(kGreen); 
    pProbPi=new TH1D("RiProbPi","Prob #pi"     ,101  ,0   ,1.05); pProbPi->SetLineColor(kRed);
    pProbMu=new TH1D("RiProbMu","Prob #mu"     ,101  ,0   ,1.05); pProbMu->SetLineColor(kBlue);
    pProbKa=new TH1D("RiProbK" ,"Prob K"       ,101  ,0   ,1.05);
    pProbPr=new TH1D("RiProbP" ,"Prob p"       ,101  ,0   ,1.05);
    pMom   =new TH1D("pMom"    ,"Track P, GeV" ,200  ,0   ,20  );
    pMipQ  =new TH1D("RiMipQ"  ,"Mip Q, ADC"   ,2000 ,0   ,4000);
    pThP   =new TH2F("RiThP"   ,"#theta_{Ckov} radian;P GeV"              ,65 ,-0.5,6.0,75,0,0.75); pThP->SetStats(0);
    pDxDy  =new TH2F("RiDxDy"  ,"distance between mip and track;cm",300,-2.5,2.5, 300,-2.5,2.5);
    pChiTh =new TProfile("RiChiTh","#chi^{2};#theta_{C}"            ,80 ,0,0.8 , -2,2);
  }
  
  Int_t iEvtCnt=0,iTrkCnt=0,iGoodCnt=0;   Float_t bz=0;
  for(Int_t iEvt=0;iEvt<pTr->GetEntries();iEvt++){//ESD events loop
    pTr->GetEvent(iEvt); iEvtCnt++; if(isPrint) Printf("");
    bz=pEsd->GetMagneticField()/10.;
    for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){//ESD tracks loop
      AliESDtrack *pTrk=pEsd->GetTrack(iTrk); iTrkCnt++;         //get next reconstructed track and increment total tracks counter
      
      Float_t mom           =pTrk->GetP();                       //track momentum
      Double_t sign         =pTrk->GetSign();                    //track sign 
      Float_t ckov          =pTrk->GetRICHsignal();              //Theta ckov for this track, rad
      Float_t chi2          =pTrk->GetRICHchi2();                //Theta ckov error for this track, rad^2 
      Int_t   qdc           =pTrk->GetRICHnclusters()/1000000;   //Mip candidate charge, qdc 
      Int_t   nphot         =pTrk->GetRICHnclusters()%1000000;   //number of photon candidates
      Float_t dx,dy;         pTrk->GetRICHdxdy(dx,dy);           //distance between mip position and track instersection
      Float_t theta,phi;     pTrk->GetRICHthetaPhi(theta,phi);   //track inclination angles in LORS
      Double_t pid[5];       pTrk->GetRICHpid(pid);              //pid vector
      
      if(ckov>0) iGoodCnt++;
      if(isPrint){
        TString comment;
        if(ckov>0)                         comment="OK";
        else if(ckov==kMipQdcCut)          comment="small QDC";
        else if(ckov==kMipDistCut)         comment="mip too far";
        else if(ckov==-1)                  comment="no intersection";
        Printf("Tr=%2i Q=%4.1f P=%.3f R=%4.2f Th=%6.3f MipQ= %4i Nph=%2i" " rad Prob : e=%.4f mu=%.4f pi=%.4f K=%.4f p=%.4f %s" ,
                     iTrk,sign,mom,TMath::Sqrt(dx*dx+dy*dy),ckov,qdc,nphot,           pid[0],pid[1],pid[2],pid[3],pid[4], comment.Data());
      }else{//collect hists
                                 pMom->Fill(mom); 
        pMipQ->Fill(qdc);
                                 pDxDy->Fill(dx,dy);
        pThP->Fill(mom,ckov); 
                                 pChiTh->Fill(ckov,chi2);
        pProbEl->Fill(pid[0]); 
                                 pProbMu->Fill(pid[1]); 
        pProbPi->Fill(pid[2]);  
                                 pProbKa->Fill(pid[3]);  
        pProbPr->Fill(pid[4]); 
      }//if plot
    }//ESD tracks loop
  }//ESD events loop
  delete pEsd;  pFile->Close();//close AliESDs.root
  
  TString summary=Form("Events: %i Tracks %i Good RICH: %i Mag Fld %.2f Tesla",iEvtCnt,iTrkCnt,iGoodCnt,bz);
  if(isPrint){
    Printf(summary.Data());
  }else{  
    TCanvas *pC=new TCanvas("c",summary.Data()); pC->Divide(2,2);
    TF1 *pPion = new TF1("RITheor","acos(sqrt(x*x+[0]*[0])/(x*[1]))",1.2,6); pPion->SetLineWidth(1);
                                                                    pPion->SetParameter(1,1.292);                                  //ref idx
    AliPID ppp;                        pPion->SetLineColor(kRed);   pPion->SetParameter(0,AliPID::ParticleMass(AliPID::kPion));    //mass
    TF1 *pKaon = (TF1*)pPion->Clone(); pKaon->SetLineColor(kGreen); pKaon->SetParameter(0,AliPID::ParticleMass(AliPID::kKaon)); 
    TF1 *pProt = (TF1*)pPion->Clone(); pProt->SetLineColor(kBlue);  pProt->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)); 
    
    pC->cd(1);        pDxDy->Draw();                                                                     //distance between mip and track intersection 
    pC->cd(2);        pMipQ->Draw();
    pC->cd(3);        pThP->Draw();       pPion->Draw("same"); pKaon->Draw("same"); pProt->Draw("same"); //Theta Ckov versus p + theoretical curves
    pC->cd(4);        pChiTh->Draw();                                                                    //Theta Ckov error versus theta Ckov
    
    TCanvas *pC2=new TCanvas("c2",summary.Data()); pC2->Divide(2,2);
    pC2->cd(1);        pProbPi->Draw(); pProbMu->Draw("same"); pProbEl->Draw("same");
    pC2->cd(2);        pProbKa->Draw();
    pC2->cd(3);        pProbPr->Draw();
    pC2->cd(4);        pMom->Draw();
  }
}//EsdQA()
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
