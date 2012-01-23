#include "AliKFConversionMother.h"
#include "AliKFParticle.h"
#include "TMath.h"
#include "TVector3.h"

using namespace std;

ClassImp(AliKFConversionMother)

AliKFConversionMother::AliKFConversionMother() :
AliKFParticle(),
fOpeningAngle(-1),
fAlpha(-1)

{
  //Default constructor
  fLabel[0] = -1;
  fLabel[1] = -1;
}

/*AliKFConversionMother::AliKFConversionMother(const AliKFParticle& d1, const AliKFParticle& d2) :
AliKFParticle(d1,d2),
fOpeningAngle(-1),
fAlpha(-1)

{
  //Default constructor
  fLabel[0] = -1;
  fLabel[1] = -1;

	// Calculate Opening Angle
	TVector3 v1(d1.GetPx(),d1.GetPy(),d1.GetPz());
	TVector3 v2(d2.GetPx(),d2.GetPy(),d2.GetPz());
      	fOpeningAngle=v1.Angle(v2);
	// Calculate Alpha
	if((d1.GetE()+d2.GetE()) != 0){
	    fAlpha=TMath::Abs((d1.GetE()-d2.GetE())/(d1.GetE()+d2.GetE()));
	}
}*/

AliKFConversionMother::AliKFConversionMother(const AliKFConversionPhoton& d1, const AliKFConversionPhoton& d2) :
AliKFParticle(d1,d2),
fOpeningAngle(-1),
fAlpha(-1)

{
  //Default constructor
    fLabel[0] = -1;
    fLabel[1] = -1;

	// Calculate Opening Angle
	TVector3 v1(d1.GetPx(),d1.GetPy(),d1.GetPz());
	TVector3 v2(d2.GetPx(),d2.GetPy(),d2.GetPz());
      	fOpeningAngle=v1.Angle(v2);
	// Calculate Alpha
	if((d1.GetE()+d2.GetE()) != 0){
	    fAlpha=TMath::Abs((d1.GetE()-d2.GetE())/(d1.GetE()+d2.GetE()));
	}
}

AliKFConversionMother::AliKFConversionMother(const AliKFConversionMother & original) :
AliKFParticle(original),
fOpeningAngle(original.fOpeningAngle),
fAlpha(original.fAlpha)
{
  //Copy constructor
  fLabel[0] = original.fLabel[0];
  fLabel[1] = original.fLabel[1];
}


AliKFConversionMother & AliKFConversionMother::operator = (const AliKFConversionMother & /*source*/)
{
  // assignment operator
  return *this;
}

Double_t AliKFConversionMother::GetRapidity()
{
    Double_t rapidity;
    if(GetE() - GetPz() <= 0 || GetE() + GetPz() <= 0){
	cout << "Error: |Pz| > E !!!! " << endl;
	rapidity=8.;
    } else {
	rapidity = 0.5*(TMath::Log((GetE() +GetPz()) / (GetE()-GetPz())));
    }

return rapidity;
}

Double_t AliKFConversionMother::Phi() const
{
    Double_t phi = AliKFParticle::GetPhi();
    if (phi < 0.) phi += 2. * TMath::Pi();
    return phi;
}
