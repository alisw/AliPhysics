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

ClassImp(AliRICHTracker)


Int_t AliRICHTracker::PropagateBack(AliESD *pESD)
{
  Int_t iNtracks=pESD->GetNumberOfTracks();
  Double_t b=GetFieldMap()->SolenoidField()/10;// magnetic field in Tesla
  AliDebug(1,Form("Start with %i tracks in %f Tesla field",iNtracks,b));
  Double_t xb[3],pb[3];//tmp storage for track parameters
  TVector3 x0(0,0,0); TVector3 p0(0,0,0);//tmp storage for AliRICHHelix
  
  AliRICH *pRich=((AliRICH*)gAlice->GetDetector("RICH"));
  
//  pRich->GetLoader()->GetRunLoader()->LoadHeader();pRich->GetLoader()->GetRunLoader()->LoadKinematics();
//  AliStack *pStack =   pRich->GetLoader()->GetRunLoader()->Stack();
//  TParticle *pParticle = pStack->Particle(0);
//  p0.SetMagThetaPhi(pParticle->P(),pParticle->Theta(),pParticle->Phi());

  for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
    AliESDtrack *pTrack = pESD->GetTrack(iTrackN);// get next reconstructed track
//  if((pTrack->GetStatus()&AliESDtrack::kTOFout)==0) continue; //ignore tracks not recontructed by TOF
    pTrack->GetXYZ(xb); 
      pTrack->GetPxPyPz(pb); 
          Int_t status=pTrack->GetStatus()&AliESDtrack::kTOFout;//get running track parameters
    AliDebug(1,Form("Track %i pmod=%f mass=%f stat=%i",iTrackN,pTrack->GetP(),pTrack->GetMass(),status));
//      x.Print();p.Print();
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
        
    
  }//ESD tracks loop
  AliDebug(1,"Stop.");  
  return 0;
} //pure virtual from AliTracker

Int_t AliRICHTracker::LoadClusters(TTree *pTree)
{
// Load clusters for RICH
  AliDebug(1,"Start.");  pTree->GetEntry(0);  AliDebug(1,"Stop."); return 0;
}
