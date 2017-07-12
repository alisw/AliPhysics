#include "AliAODConversionMother.h"
#include "AliKFConversionMother.h"
#include "AliKFParticle.h"

// Author D. Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAODConversionMother)

AliAODConversionMother::AliAODConversionMother():
AliAODConversionParticle(),
    fMCLabel(-1),
    fChi2(-1),
    fOpeningAngle(-1),
    fAlpha(-1),
    fWeight(1),
    fdcaBetweenPhotons(1),
    fdcaZPrimVtx(100),
    fdcaRPrimVtx(100),
    fQuality(0),
    fTrueMeson(0)
{
	fLabel[0] = -1;
	fLabel[1] = -1;
	fLabel[2] = 0;
   
   fProductionVtx[0]=0;
   fProductionVtx[1]=0;
   fProductionVtx[2]=0;
  
}

AliAODConversionMother::AliAODConversionMother(const AliKFConversionMother *kf):
AliAODConversionParticle(),
fMCLabel(kf->GetMCLabel()),
fChi2(kf->GetChi2()),
fOpeningAngle(kf->GetOpeningAngle()),
fAlpha(kf->GetAlpha()),
fWeight(1),
fdcaBetweenPhotons(100),
fdcaZPrimVtx(100),
fdcaRPrimVtx(100),
fQuality(0),
fTrueMeson(0)
{
    // Set 4momentu
    SetPxPyPzE(kf->GetPx(),kf->GetPy(),kf->GetPz(),kf->GetE());

   fProductionVtx[0]=0;
   fProductionVtx[1]=0;
   fProductionVtx[2]=0;
  
     //Set Decay Photon Labels
    fLabel[0]=kf->GetGammaLabel(0);
    fLabel[1]=kf->GetGammaLabel(1);
    fLabel[2]=0;
}

AliAODConversionMother::AliAODConversionMother(const AliAODConversionPhoton *y1, const AliAODConversionPhoton *y2):
AliAODConversionParticle(),
fMCLabel(-1),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1),
fWeight(1),
fdcaBetweenPhotons(1),
fdcaZPrimVtx(100),
fdcaRPrimVtx(100),
fQuality(0),
fTrueMeson(0)
{
    // Set 4momentum
    SetPxPyPzE(y1->Px()+y2->Px(),y1->Py()+y2->Py(),y1->Pz()+y2->Pz(),y1->E()+y2->E());

    // Calculate Opening Angle
    TVector3 v1(y1->Px(),y1->Py(),y1->Pz());
    TVector3 v2(y2->Px(),y2->Py(),y2->Pz());
    fOpeningAngle=v1.Angle(v2);
    fdcaBetweenPhotons = CalculateDistanceBetweenPhotons(y1,y2,fProductionVtx);
    DetermineMesonQuality(y1,y2);
    // Calculate Alpha
    if((y1->E()+y2->E()) != 0){
		fAlpha=(y1->E()-y2->E())/(y1->E()+y2->E());
    }

    // Set Chi2 to the mean chi2 of gammas
 //   fChi2=0.5*(y1->GetChi2perNDF()+y2->GetChi2perNDF());

    //Set Decay Photon Labels
    fLabel[0]=-1;
    fLabel[1]=-1;
    fLabel[2]=0;
}

AliAODConversionMother::AliAODConversionMother(const AliAODConversionMother *meson, const AliAODConversionPhoton *gamma):
AliAODConversionParticle(),
fMCLabel(-1),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1),
fWeight(1),
fdcaBetweenPhotons(1),
fdcaZPrimVtx(100),
fdcaRPrimVtx(100),
fQuality(0),
fTrueMeson(0)
{
    // Set 4momentum
    SetPxPyPzE(meson->Px()+gamma->Px(),meson->Py()+gamma->Py(),meson->Pz()+gamma->Pz(),meson->E()+gamma->E());

    // Calculate Opening Angle
    TVector3 v1(meson->Px(),meson->Py(),meson->Pz());
    TVector3 v2(gamma->Px(),gamma->Py(),gamma->Pz());
    fOpeningAngle=v1.Angle(v2);
     
	fProductionVtx[0]=0;
	fProductionVtx[1]=0;
	fProductionVtx[2]=0;

    // Calculate Alpha
    if((meson->E()+gamma->E()) != 0){
		fAlpha=(meson->E()-gamma->E())/(meson->E()+gamma->E());
    }

    // Set Chi2 to the mean chi2 of gammas
	// fChi2=0.5*(y1->GetChi2perNDF()+y2->GetChi2perNDF());

    //Set Decay Photon Labels
    fLabel[0]=-1;
    fLabel[1]=-1;
    fLabel[2]=0;
}


AliAODConversionMother::AliAODConversionMother(const AliAODConversionMother *meson1, const AliAODConversionMother *meson2):
AliAODConversionParticle(),
fMCLabel(-1),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1),
fWeight(1),
fdcaBetweenPhotons(1),
fdcaZPrimVtx(100),
fdcaRPrimVtx(100),
fQuality(0),
fTrueMeson(0)
{
	// Set 4momentum
	SetPxPyPzE(meson1->Px()+meson2->Px(),meson1->Py()+meson2->Py(),meson1->Pz()+meson2->Pz(),meson1->E()+meson2->E());

	// Calculate Opening Angle
	TVector3 v1(meson1->Px(),meson1->Py(),meson1->Pz());
	TVector3 v2(meson2->Px(),meson2->Py(),meson2->Pz());
	fOpeningAngle=v1.Angle(v2);

	fProductionVtx[0]=0;
	fProductionVtx[1]=0;
	fProductionVtx[2]=0;

	// Calculate Alpha
	if((meson1->E()+meson2->E()) != 0){
		fAlpha=(meson1->E()-meson2->E())/(meson1->E()+meson2->E());
	}

	// Set Chi2 to the mean chi2 of gammas
	// fChi2=0.5*(y1->GetChi2perNDF()+y2->GetChi2perNDF());

	//Set Decay Photon Labels
	fLabel[0]=-1;
	fLabel[1]=-1;
	fLabel[2]=0;
}


AliAODConversionMother::AliAODConversionMother(const AliV0ParticleStrange *v0, const AliAODConversionPhoton *gamma):
AliAODConversionParticle(),
fMCLabel(-1),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1),
fWeight(1),
fdcaBetweenPhotons(1),
fdcaZPrimVtx(100),
fdcaRPrimVtx(100),
fQuality(0),
fTrueMeson(0)
{
    // Set 4momentum
    SetPxPyPzE(v0->Px()+gamma->Px(),v0->Py()+gamma->Py(),v0->Pz()+gamma->Pz(),v0->E()+gamma->E());

    // Calculate Opening Angle
    TVector3 v1(v0->Px(),v0->Py(),v0->Pz());
    TVector3 v2(gamma->Px(),gamma->Py(),gamma->Pz());
    fOpeningAngle=v1.Angle(v2);
     
	fProductionVtx[0]=0;
	fProductionVtx[1]=0;
	fProductionVtx[2]=0;

    // Calculate Alpha
    if((v0->E()+gamma->E()) != 0){
		fAlpha=(v0->E()-gamma->E())/(v0->E()+gamma->E());
    }

    // Set Chi2 to the mean chi2 of gammas
	// fChi2=0.5*(y1->GetChi2perNDF()+y2->GetChi2perNDF());

    //Set Decay Photon Labels
    fLabel[0]=-1;
    fLabel[1]=-1;
    fLabel[2]=0;
}


AliAODConversionMother::~AliAODConversionMother() {
    // empty standard destructor
}

TParticle *AliAODConversionMother::GetMCParticle(AliMCEvent *mcEvent){
    if(!mcEvent){AliError("MCEvent not defined");return 0x0;}

    if(fMCLabel>-1){
        return mcEvent->Particle(fMCLabel);
    }
    return 0x0;
}

Bool_t AliAODConversionMother::IsTrueMeson(AliMCEvent *mcEvent,Int_t pdgcode){
    TParticle *part=GetMCParticle(mcEvent);

    if(part){
	// Check if it is a true photon
	if(part->GetPdgCode()==pdgcode){
	    return kTRUE;
	}
    }
    return kFALSE;
}

Float_t AliAODConversionMother::CalculateDistanceBetweenPhotons(const AliAODConversionPhoton* y1, const AliAODConversionPhoton* y2 , Double_t prodPoint[3]){

   TVector3 a(y1->GetConversionX(),y1->GetConversionY(),y1->GetConversionZ());
   TVector3 b(y1->GetPx(),y1->GetPy(),y1->GetPz());
   TVector3 c(y2->GetConversionX(),y2->GetConversionY(),y2->GetConversionZ());
   TVector3 d(y2->GetPx(),y2->GetPy(),y2->GetPz());
   
   TVector3 n = b.Cross(d);
   TVector3 nn = n.Unit();
   
   Double_t dist = 0;
   if (n.Mag() == 0){
      TVector3 e = a-c;
      if (d.Mag() != 0){
         dist = TMath::Abs((e.Cross(d)).Mag())/TMath::Abs(d.Mag());
      }
      prodPoint[0] = 0;
      prodPoint[1] = 0;
      prodPoint[2] = 0;
   } else {   
      dist = TMath::Abs(n.Dot(c-a))/TMath::Abs(n.Mag());
      Double_t lambda = (b.Dot(d) * (a-c).Dot(d) - d.Dot(d) * (a-c).Dot(b))/(b.Dot(b) * d.Dot(d) - TMath::Power(b.Dot(d),2));
      Double_t mu = ((a-c).Dot(d) * b.Dot(b) - (a-c).Dot(b) * b.Dot(d) )/(b.Dot(b) * d.Dot(d) - TMath::Power(b.Dot(d),2));
      
      TVector3 S1 = a + lambda* b;
      TVector3 S2 = c + mu* d;
      TVector3 Prod = S1 + 0.5*dist*(S2-S1).Unit();
      prodPoint[0] = Prod(0);
      prodPoint[1] = Prod(1);
      prodPoint[2] = Prod(2);
      
   }
   if (dist > 1000) dist = 999.;
   return dist;
}

///________________________________________________________________________
void AliAODConversionMother::CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex){

   Double_t primCo[3] = {primVertex->GetX(),primVertex->GetY(),primVertex->GetZ()};

   Double_t absoluteP = TMath::Sqrt(TMath::Power(Px(),2) + TMath::Power(Py(),2) + TMath::Power(Pz(),2));   
   Double_t p[3] = {Px()/absoluteP,Py()/absoluteP,Pz()/absoluteP};
   Double_t CP[3];
   
   CP[0] =  fProductionVtx[0] - primCo[0];
   CP[1] =  fProductionVtx[1] - primCo[1];
   CP[2] =  fProductionVtx[2] - primCo[2];
   
   Double_t Lambda = - (CP[0]*p[0]+CP[1]*p[1]+CP[2]*p[2])/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
   
   Double_t S[3];
   S[0] = fProductionVtx[0] + p[0]*Lambda;
   S[1] = fProductionVtx[1] + p[1]*Lambda;
   S[2] = fProductionVtx[2] + p[2]*Lambda;
   
   fdcaRPrimVtx = TMath::Sqrt( TMath::Power(primCo[0]-S[0],2) + TMath::Power(primCo[1]-S[1],2));
   fdcaZPrimVtx = primCo[2]-S[2];
   
   
//    cout << "DCA z: " << dca[1] << "\t DCA r: " << dca[0] << "\t DCA 3d: " << TMath::Sqrt(dca[1]*dca[1] + dca[0]*dca[0]) << endl;
   
   
   return;
}

///________________________________________________________________________
void AliAODConversionMother::DetermineMesonQuality(const AliAODConversionPhoton* y1, const AliAODConversionPhoton* y2){
   UChar_t photonQA1 = y1->GetPhotonQuality();
   UChar_t photonQA2 = y2->GetPhotonQuality();
   
   if (photonQA1 == 0 || photonQA2 == 0){
      fQuality = 0;
      return;
   } 
   if (photonQA1 == 1 && photonQA2 == 1){
      fQuality = 1;
      return;
   }               
   if (photonQA1 == 2 && photonQA2 == 2){
      fQuality = 4;
      return;
   }               
   if (photonQA1 == 3 && photonQA2 == 3){
      fQuality = 6;
      return;
   }         
   if (photonQA1 == 1){
       if (photonQA2 == 2){
        fQuality = 2;  
        return;
       }
       if (photonQA2 == 3){
        fQuality = 3;  
        return;
       }
   }
   if (photonQA2 == 1){
       if (photonQA1 == 2){
        fQuality = 2;  
        return;
       }
       if (photonQA1 == 3){
        fQuality = 3;  
        return;
       }
   }
   if ((photonQA1 == 2 && photonQA2 == 3)|| (photonQA1 == 3 && photonQA2 == 2)){
        fQuality = 5;  
        return;
   }
   
}

