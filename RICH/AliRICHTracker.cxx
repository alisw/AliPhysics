#include "AliRICHTracker.h"
#include <AliESD.h>
#include <TVector3.h>
#include <TTree.h>
#include "AliRICH.h"
#include "AliRICHHelix.h"
#include <AliMagF.h>
#include "AliRICHRecon.h"
#include <AliStack.h>
#include <TParticle.h>
#include <TMath.h>
#include <AliRun.h>
ClassImp(AliRICHTracker)
//__________________________________________________________________________________________________
Int_t AliRICHTracker::PropagateBack(AliESD *pESD)
{
// Interface callback methode invoked by AliRecontruction during tracking after TOF
// It steers to different way to provide the final reconstructed information sutable for analisys:
// 1. AliESD  - reconstructed tracks are used     
// 2. RICH private ntuple for debug- stack particles used instead of reconstructed tracks     
  AliDebug(1,"Start pattern recognition");
  if(pESD->GetNumberOfTracks())
    RecWithESD(pESD);
  else
    RecWithStack(0);
  AliDebug(1,"Stop pattern recognition");
  return 0; // error code: 0=no error;
}//PropagateBack()
//__________________________________________________________________________________________________
void AliRICHTracker::RecWithESD(AliESD *pESD)
{
//recontruction from ESD- primary way to reconstruct particle ID signal from tracks provided by core detectors

  Int_t iNtracks=pESD->GetNumberOfTracks();
  Double_t b=GetFieldMap()->SolenoidField()/10;// magnetic field in Tesla
  AliDebug(1,Form("Start with %i tracks in %f Tesla field",iNtracks,b));
  
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  
  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
//  if((pTrack->GetStatus()&AliESDtrack::kTOFout)==0) continue; //ignore tracks not recontructed by TOF
//    pTrack->GetXYZ(xb); 
//    pTrack->GetPxPyPz(pb); 
    Int_t status=pTrack->GetStatus()&AliESDtrack::kTOFout;//get running track parameters
    Int_t charge = (Int_t)(-TMath::Sign(1.,pTrack->GetSign()*b));
    AliDebug(1,Form("Track %i pmod=%f charge=%i stat=%i",iTrackN,pTrack->GetP(),charge,status));
    AliRICHHelix helix(pTrack->X3(),pTrack->P3(),charge,b);   
    Int_t iChamber=helix.RichIntersect(pRich->P());        
    AliDebug(1,Form("intersection with %i chamber found",iChamber));
    if(!iChamber) continue;//intersection with no chamber found
//find MIP cluster candidate (cluster which is closest to track intersection point)    
    Double_t distMip=9999,distX=0,distY=0; //min distance between clusters and track position on PC 
    Int_t iMipId=0; //index of that min distance cluster
    Double_t chargeMip=0; //charge of the MIP
    for(Int_t iClusN=0;iClusN<pRich->Clusters(iChamber)->GetEntries();iClusN++){//clusters loop for intersected chamber
      AliRICHCluster *pClus=(AliRICHCluster*)pRich->Clusters(iChamber)->UncheckedAt(iClusN);//get pointer to current cluster
      Double_t distCurrent=pClus->DistTo(helix.PosPc());//distance between current cluster and helix intersection with PC
      if(distCurrent<distMip){
        distMip=distCurrent;
        iMipId=iClusN;
        distX=pClus->DistX(helix.PosPc());
        distY=pClus->DistY(helix.PosPc());
        chargeMip=pClus->Q();
      }//find cluster nearest to the track       
      AliDebug(1,Form("Ploc (%f,%f,%f) dist= %f",helix.Ploc().Mag(),helix.Ploc().Theta()*TMath::RadToDeg(),
                                       helix.Ploc().Phi()*TMath::RadToDeg(),pClus->DistTo(helix.PosPc())));
    }//clusters loop for intersected chamber
    
    AliDebug(1,Form("Min distance cluster: %i dist is %f",iMipId,distMip));
//
// HERE CUTS ON GOLD RINGS....
//
    if(distMip>1||chargeMip<100) {
      //track not accepted for pattern recognition
      pTrack->SetRICHsignal(-999.); //to be improved by flags...
      continue;
    }
//
    AliRICHRecon recon(&helix,pRich->Clusters(iChamber),iMipId); //actual job is done there

    Double_t thetaCerenkov=recon.ThetaCerenkov(); //search for mean Cerenkov angle for this track
    
    pTrack->SetRICHcluster(iMipId+100000*iChamber);
    pTrack->SetRICHdxdy(distX,distY);
    pTrack->SetRICHthetaPhi(helix.Ploc().Theta(),helix.Ploc().Phi());
    pTrack->SetRICHsignal(thetaCerenkov);
    pTrack->SetRICHnclusters(recon.GetHoughPhotons());
    
    AliDebug(1,Form("FINAL Theta Cerenkov=%f",pTrack->GetRICHsignal()));
//
    if(pTrack->GetRICHsignal()>0) {
      AliDebug(1,Form("Start to assign the probabilities"));
      Double_t sigmaPID[AliPID::kSPECIES];
      Double_t richPID[AliPID::kSPECIES];
      for (Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) {
        sigmaPID[iPart] = 0;
        for(Int_t iphot=0;iphot<pRich->Clusters(iChamber)->GetEntries();iphot++) {
          recon.SetPhotonIndex(iphot);
          if(recon.GetPhotonFlag() == 2) {
            Double_t sigma = AliRICHParam::SigmaSinglePhoton(iPart,pTrack->GetP(),recon.GetTrackTheta(),recon.GetPhiPoint()-recon.GetTrackPhi()).Mag();
            sigmaPID[iPart] += 1/(sigma*sigma);
          }
        }
        sigmaPID[iPart] = 1/TMath::Sqrt(sigmaPID[iPart])*0.001;
        AliDebug(1,Form("sigma for %s is %f rad",AliPID::ParticleName(iPart),sigmaPID[iPart]));
      }
      CalcProb(thetaCerenkov,pTrack->GetP(),sigmaPID,richPID);
      pTrack->SetRICHpid(richPID);         
      AliDebug(1,Form("PROBABILITIES ---> %f - %f - %f - %f - %f",richPID[0],richPID[1],richPID[2],richPID[3],richPID[4]));
    }
    
    
  }//ESD tracks loop
  AliDebug(1,"Stop.");  
} //RecWithESD
//__________________________________________________________________________________________________
void AliRICHTracker::RecWithStack(TNtupleD *hn)
{
//Reconstruction for particles from STACK. This methode is to be used for RICH standalone when no other detectors are switched on,
//so normal tracking is not available   
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
//
// to be updated for us!!
//
    AliDebug(1,Form("Track %i is a %s with charge %i and momentum %f",
            iTrackN,pParticle->GetPDG()->GetName(),Int_t(pParticle->GetPDG()->Charge()),pParticle->P()));
//    if(pParticle->GetMother(0)!=-1) continue; //consider only primaries
    if(pParticle->GetPDG()->Charge()==0||TMath::Abs(Int_t(pParticle->GetPDG()->Charge()))!=3) continue; //to avoid photons from stack...
    hnvec[0]=pParticle->P();
    hnvec[1]=pParticle->GetPDG()->Charge();
    p0.SetMagThetaPhi(pParticle->P(),pParticle->Theta(),pParticle->Phi());
    x0.SetXYZ(pParticle->Vx(),pParticle->Vy(),pParticle->Vz());
    AliRICHHelix helix(x0,p0,TMath::Sign(1,(Int_t)pParticle->GetPDG()->Charge()),b);   
    Int_t iChamber=helix.RichIntersect(pRich->P());        
    hnvec[2]=helix.Ploc().Theta();
    hnvec[3]=helix.Ploc().Phi();
    AliDebug(1,Form("intersection with %i chamber found",iChamber));
    if(!iChamber) continue;// no intersection with RICH found
    hnvec[4]=helix.PosPc().X();
    hnvec[5]=helix.PosPc().Y();
    Double_t distMip=9999;   //min distance between clusters and track position on PC 
    Double_t mipX=-1;      //min distance between clusters and track position on PC 
    Double_t mipY=-1;      //min distance between clusters and track position on PC 
    Double_t chargeMip=-1; // charge MIP to find
    Int_t iMipId=-1;       //index of that min distance cluster 
    for(Int_t iClusN=0;iClusN<pRich->Clusters(iChamber)->GetEntries();iClusN++){//clusters loop for intersected chamber
      AliRICHCluster *pClus=(AliRICHCluster*)pRich->Clusters(iChamber)->UncheckedAt(iClusN);//get pointer to current cluster
      Double_t distCurrent=pClus->DistTo(helix.PosPc());//ditance between current cluster and helix intersection with PC
      if(distCurrent<distMip){distMip=distCurrent;mipX=pClus->X();
                                                  mipY=pClus->Y();
                                                  chargeMip=pClus->Q();iMipId=1000000*iChamber+iClusN;}//find cluster nearest to the track 
      
      AliDebug(1,Form("Ploc (%f,%f,%f) dist= %f",helix.Ploc().Mag(),helix.Ploc().Theta()*TMath::RadToDeg(),
                                                                    helix.Ploc().Phi()*TMath::RadToDeg(),pClus->DistTo(helix.PosPc())));
    }//clusters loop for intersected chamber
    
    AliDebug(1,Form("Min distance cluster: %i dist is %f",iMipId,distMip));
    hnvec[6]=mipX;hnvec[7]=mipY;
    hnvec[8]=chargeMip;
    AliRICHRecon recon(&helix,pRich->Clusters(iChamber),iMipId);
    Double_t thetaCerenkov=recon.ThetaCerenkov(); //search for mean Cerenkov angle for this track
    hnvec[9]=thetaCerenkov;
    hnvec[10]=recon.GetHoughPhotons();
    hnvec[11]=(Double_t)iMipId;
    hnvec[12]=(Double_t)iChamber;
    hnvec[13]=(Double_t)pParticle->GetPdgCode();
    if(hn) hn->Fill(hnvec);
    AliDebug(1,Form("FINAL Theta Cerenkov=%f",thetaCerenkov));
  }//stack particles loop
  
  pRich->GetLoader()->UnloadRecPoints();
  AliDebug(1,"Stop.");  
}//RecWithStack
//__________________________________________________________________________________________________
Int_t AliRICHTracker::LoadClusters(TTree *pTree)
{
// Load clusters for RICH
  AliDebug(1,"Start.");  pTree->GetEntry(0);  AliDebug(1,"Stop."); return 0;
}
//__________________________________________________________________________________________________
void AliRICHTracker::CalcProb(Double_t thetaCer,Double_t pmod, Double_t *sigmaPID, Double_t *richPID)
{
// Calculates probability to be a electron-muon-pion-kaon-proton
// from the given Cerenkov angle and momentum assuming no initial particle composition
// (i.e. apriory probability to be the particle of the given sort is the same for all sorts)
  Double_t height[AliPID::kSPECIES];Double_t totalHeight=0;
  Double_t thetaTh[AliPID::kSPECIES];
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++){
    height[iPart]=0;
    Double_t mass = AliPID::ParticleMass(iPart);
    Double_t refIndex=AliRICHParam::RefIdxC6F14(AliRICHParam::MeanCkovEnergy());
    Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(refIndex*pmod);
    thetaTh[iPart]=0;
    if(cosThetaTh>=1) continue;
    thetaTh[iPart] = TMath::ACos(cosThetaTh);
//    Double_t sinThetaThNorm = TMath::Sin(thetaTh)/TMath::Sqrt(1-1/(refIndex*refIndex));
//    Double_t sigmaThetaTh = (0.014*(1/sinThetaThNorm-1) + 0.0043)*1.25;
//    height[iPart] = TMath::Gaus(thetaCer,thetaTh,sigmaThetaTh);
    height[iPart] = TMath::Gaus(thetaCer,thetaTh[iPart],sigmaPID[iPart],kTRUE);
    totalHeight +=height[iPart];
    AliDebug(1,Form(" Particle %s with mass %f with height %f and thetaTH %f",AliPID::ParticleName(iPart),mass,height[iPart],thetaTh[iPart]));
    AliDebug(1,Form(" partial height %15.14f total height %15.14f",height[iPart],totalHeight));
  }
  if(totalHeight<1e-5) {for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)richPID[iPart]=1.0/AliPID::kSPECIES;return;}
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) richPID[iPart] = height[iPart]/totalHeight;
  Int_t iPartNear = TMath::LocMax(AliPID::kSPECIES,richPID);
  if(TMath::Abs(thetaCer-thetaTh[iPartNear])/sigmaPID[iPartNear]>3) for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++)richPID[iPart]=1.0/AliPID::kSPECIES;
  //last line is to check if the nearest thetacerenkov to the teorethical one is within 3 sigma, otherwise no response (equal prob to every particle

}//CalcProb
//__________________________________________________________________________________________________
