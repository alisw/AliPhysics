#include "AliAODConversionMother.h"
#include "AliKFConversionMother.h"
#include "AliKFParticle.h"


// Author D. Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAODConversionMother)

AliAODConversionMother::AliAODConversionMother() :
AliAODConversionParticle(),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1)
{
	fLabel[0] = -1;
	fLabel[1] = -1;
}

AliAODConversionMother::AliAODConversionMother(AliKFConversionMother *kf):
AliAODConversionParticle(),
fChi2(kf->GetChi2()),
fOpeningAngle(kf->GetOpeningAngle()),
fAlpha(kf->GetAlpha())
{
    // Set 4momentu
    SetPxPyPzE(kf->GetPx(),kf->GetPy(),kf->GetPz(),kf->GetE());

     //Set Decay Photon Labels
    fLabel[0]=kf->GetGammaLabel(0);
    fLabel[1]=kf->GetGammaLabel(1);
}

AliAODConversionMother::AliAODConversionMother(AliAODConversionPhoton *y1,AliAODConversionPhoton *y2):
AliAODConversionParticle(),
fChi2(-1),
fOpeningAngle(-1),
fAlpha(-1)
{
    // Set 4momentum
    SetPxPyPzE(y1->Px()+y2->Px(),y1->Py()+y2->Py(),y1->Pz()+y2->Pz(),y1->E()+y2->E());

    // Calculate Opening Angle
    TVector3 v1(y1->Px(),y1->Py(),y1->Pz());
    TVector3 v2(y2->Px(),y2->Py(),y2->Pz());
    fOpeningAngle=v1.Angle(v2);

    // Calculate Alpha
    if((y1->E()+y2->E()) != 0){
	fAlpha=TMath::Abs((y1->E()-y2->E())/(y1->E()+y2->E()));
    }

    // Set Chi2 to the mean chi2 of gammas
    fChi2=0.5*(y1->GetChi2perNDF()+y2->GetChi2perNDF());

    //Set Decay Photon Labels
    fLabel[0]=-1;
    fLabel[1]=-1;

}

AliAODConversionMother::~AliAODConversionMother() {
// empty standard destructor

}

