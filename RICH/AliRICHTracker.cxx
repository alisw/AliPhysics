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
ClassImp(AliRICHTracker)
//__________________________________________________________________________________________________
Int_t AliRICHTracker::PropagateBack(AliESD *pESD)
{
  //invoked by AliRecontruction for RICH
  //if ESD doesn't contain tracks, try to reconstruct with particle from STACK 
  //(this case is just to forsee a standalone RICH simulation
  TNtupleD *hn=0;
  AliDebug(1,"Start pattern recognition");
  if(pESD->GetNumberOfTracks()) RecWithESD(pESD);
  else
    RecWithStack(hn);
  AliDebug(1,"Stop pattern recognition");

  return 0; // error code: 0=no error;
} //pure virtual from AliTracker
//__________________________________________________________________________________________________
void AliRICHTracker::RecWithESD(AliESD *pESD)
{
  //recontruction from ESD
  //
  Int_t iNtracks=pESD->GetNumberOfTracks();
  Double_t b=GetFieldMap()->SolenoidField()/10;// magnetic field in Tesla
  AliDebug(1,Form("Start with %i tracks in %f Tesla field",iNtracks,b));
  Double_t xb[3],pb[3];//tmp storage for track parameters
  TVector3 x0(0,0,0); TVector3 p0(0,0,0);//tmp storage for AliRICHHelix
  
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  
  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
//  if((pTrack->GetStatus()&AliESDtrack::kTOFout)==0) continue; //ignore tracks not recontructed by TOF
    pTrack->GetXYZ(xb); 
    pTrack->GetPxPyPz(pb); 
    Int_t status=pTrack->GetStatus()&AliESDtrack::kTOFout;//get running track parameters
    AliDebug(1,Form("Track %i pmod=%f mass=%f stat=%i",iTrackN,pTrack->GetP(),pTrack->GetMass(),status));
    x0.SetXYZ(xb[0],xb[1],xb[2]); p0.SetXYZ(xb[0],xb[1],xb[2]);
    AliRICHHelix helix(x0,p0,pTrack->GetSign(),b);   
    Int_t iChamber=helix.RichIntersect(pRich->P());        
    AliDebug(1,Form("intersection with %i chamber found",iChamber));
    if(!iChamber) continue;//intersection with no chamber found
    
    Double_t distMip=9999; //min distance between clusters and track position on PC 
    Int_t iMipId=0; //index of that min distance cluster 
    for(Int_t iClusN=0;iClusN<pRich->Clusters(iChamber)->GetEntries();iClusN++){//clusters loop for intersected chamber
      AliRICHcluster *pClus=(AliRICHcluster*)pRich->Clusters(iChamber)->UncheckedAt(iClusN);//get pointer to current cluster
      Double_t distCurrent=pClus->DistTo(helix.PosPc());//ditance between current cluster and helix intersection with PC
      if(distCurrent<distMip){distMip=distCurrent;iMipId=iClusN;}//find cluster nearest to the track 
      
      AliDebug(1,Form("Ploc (%f,%f,%f) dist= %f",helix.Ploc().Mag(),helix.Ploc().Theta()*TMath::RadToDeg(),
                                                                    helix.Ploc().Phi()*TMath::RadToDeg(),pClus->DistTo(helix.PosPc())));
    }////clusters loop for intersected chamber
    
    AliDebug(1,Form("Min distance cluster: %i dist is %f",iMipId,distMip));
    
    AliRICHRecon recon(&helix,pRich->Clusters(iChamber),iMipId);
    Double_t thetaCerenkov=recon.ThetaCerenkov(); //search for mean Cerenkov angle for this track
    AliDebug(1,Form("FINAL Theta Cerenkov=%f",thetaCerenkov));
    pTrack->SetRICHsignal(thetaCerenkov);
        
//    Double_t richPID[5]={0.2,0.2,0.2,0.2,0.2}; //start with equal probs for (e,mu,pi,k,p)
//    CalcProb(thetaCerenkov,p0.Mag(),richPID);
//    pTrack->SetRICHpid(richPID);         
    
  }//ESD tracks loop
  AliDebug(1,"Stop.");  
} //RecWithESD
//__________________________________________________________________________________________________
void AliRICHTracker::RecWithStack(TNtupleD *hn)
{
  // reconstruction for particles from STACK
  //
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  
//  pRich->GetLoader()->GetRunLoader()->LoadHeader();
  pRich->GetLoader()->GetRunLoader()->LoadKinematics();
  AliStack *pStack =   pRich->GetLoader()->GetRunLoader()->Stack();
  Int_t iNtracks=pStack->GetNtrack();
  AliDebug(1,Form(" Start reconstruction with %i track(s) from Stack",iNtracks));
  
  Double_t hnvec[12];
  
  Double_t b=GetFieldMap()->SolenoidField()/10;// magnetic field in Tesla
  AliDebug(1,Form("Start with simulated %i tracks in %f Tesla field",iNtracks,b));
  TVector3 x0(0,0,0); TVector3 p0(0,0,0);//tmp storage for AliRICHHelix
  

  if(pRich->GetLoader()->LoadRecPoints()) {AliDebug(1,Form("No clusters found in RICH"));return;}
  pRich->GetLoader()->TreeR()->GetEntry(0);

  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    TParticle *pParticle = pStack->Particle(iTrackN);
    AliDebug(1,Form("Track %i is a %s with charge %i and momentum %f",
            iTrackN,pParticle->GetPDG()->GetName(),Int_t(pParticle->GetPDG()->Charge()),pParticle->P()));
//    if(pParticle->GetMother(0)!=-1) continue; //consider only primaries
    if(pParticle->GetPDG()->Charge()==0||TMath::Abs(Int_t(pParticle->GetPDG()->Charge()))!=3) continue; //to avoid photons from stack...
    hnvec[0]=pParticle->P();
    hnvec[1]=pParticle->GetPDG()->Charge();
    hnvec[2]=pParticle->Theta();
    hnvec[3]=pParticle->Phi();
    p0.SetMagThetaPhi(pParticle->P(),pParticle->Theta(),pParticle->Phi());
    x0.SetXYZ(pParticle->Vx(),pParticle->Vy(),pParticle->Vz());
    AliRICHHelix helix(x0,p0,TMath::Sign(1,(Int_t)pParticle->GetPDG()->Charge()),b);   
    Int_t iChamber=helix.RichIntersect(pRich->P());        
    AliDebug(1,Form("intersection with %i chamber found",iChamber));
    if(!iChamber) continue;// no intersection with RICH found
    hnvec[4]=helix.PosPc().X();
    hnvec[5]=helix.PosPc().Y();
    Double_t distMip=9999;   //min distance between clusters and track position on PC 
    Double_t mipX=kBad;      //min distance between clusters and track position on PC 
    Double_t mipY=kBad;      //min distance between clusters and track position on PC 
    Double_t chargeMip=kBad; // charge MIP to find
    Int_t iMipId=kBad;       //index of that min distance cluster 
    for(Int_t iClusN=0;iClusN<pRich->Clusters(iChamber)->GetEntries();iClusN++){//clusters loop for intersected chamber
      AliRICHcluster *pClus=(AliRICHcluster*)pRich->Clusters(iChamber)->UncheckedAt(iClusN);//get pointer to current cluster
      Double_t distCurrent=pClus->DistTo(helix.PosPc());//ditance between current cluster and helix intersection with PC
      if(distCurrent<distMip){distMip=distCurrent;mipX=pClus->X();
                                                  mipY=pClus->Y();
                                                  chargeMip=pClus->Q();iMipId=1000000*iChamber+iClusN;}//find cluster nearest to the track 
      
      AliDebug(1,Form("Ploc (%f,%f,%f) dist= %f",helix.Ploc().Mag(),helix.Ploc().Theta()*TMath::RadToDeg(),
                                                                    helix.Ploc().Phi()*TMath::RadToDeg(),pClus->DistTo(helix.PosPc())));
    }////clusters loop for intersected chamber
    
    AliDebug(1,Form("Min distance cluster: %i dist is %f",iMipId,distMip));
    hnvec[6]=mipX;hnvec[7]=mipY;
    hnvec[8]=chargeMip;
    AliRICHRecon recon(&helix,pRich->Clusters(iChamber),iMipId);
    Double_t thetaCerenkov=recon.ThetaCerenkov(); //search for mean Cerenkov angle for this track
    hnvec[9]=thetaCerenkov;
    hnvec[10]=recon.GetHoughPhotons();
    hnvec[11]=(Double_t)iMipId;
    if(hn) hn->Fill(hnvec);
    AliDebug(1,Form("FINAL Theta Cerenkov=%f",thetaCerenkov));
//    pTrack->SetRICHsignal(thetaCerenkov);

//    Double_t richPID[5]={0.2,0.2,0.2,0.2,0.2}; //start with equal probs for (e,mu,pi,k,p)
//    CalcProb(thetaCerenkov,p0.Mag(),richPID);
    
  }//ESD tracks loop
  
 pRich->GetLoader()->UnloadRecPoints();

  AliDebug(1,"Stop.");  
} //RecWithStack

Int_t AliRICHTracker::LoadClusters(TTree *pTree)
{
// Load clusters for RICH
  AliDebug(1,"Start.");  pTree->GetEntry(0);  AliDebug(1,"Stop."); return 0;
}

//__________________________________________________________________________________________________
void AliRICHTracker::CalcProb(Double_t thetaCer,Double_t pmod, Double_t *richPID)
{
// 
  Double_t height[AliPID::kSPECIES];Double_t totalHeight=0;
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++){
    Double_t mass = AliPID::ParticleMass(iPart);
    Double_t refIndex=AliRICHParam::IndOfRefC6F14(6.755);
    Double_t cosThetaTh = TMath::Sqrt(mass*mass+pmod*pmod)/(refIndex*pmod);
    if(cosThetaTh>=1) {break;}
    Double_t thetaTh = TMath::ACos(cosThetaTh);
    Double_t sinThetaThNorm = TMath::Sin(thetaTh)/TMath::Sqrt(1-1/(refIndex*refIndex));
    Double_t sigmaThetaTh = (0.014*(1/sinThetaThNorm-1) + 0.0043)*1.25;
    height[iPart] = TMath::Gaus(thetaCer,thetaTh,sigmaThetaTh);
    totalHeight +=height[iPart];
  }
  for(Int_t iPart=0;iPart<AliPID::kSPECIES;iPart++) richPID[iPart] = height[iPart]/totalHeight;    
}//CalcProb
