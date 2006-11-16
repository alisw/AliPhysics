// #include "EVE/Alieve/EventAlieve.h"
// #include "Reve/RGTopFrame.h"
// #include "Reve/Cascade.h"

// #include "AliESD.h"
// #include "AliESDtrack.h"
// #include "AliESDcascade.h"
// #include "AliESDVertex.h"

// using namespace Reve;
// using namespace Alieve;



Reve::Cascade* esd_make_cas(Reve::TrackRnrStyle* rnrStyle, AliESDVertex* primVtx, 
			    AliESDcascade* cas, AliESDtrack* neg, AliESDtrack* pos,
			    AliESDtrack* bach,Int_t i) {

  Reve::Cascade* myCas = new Reve::Cascade(rnrStyle);
  myCas->SetESDIndex(i);

  static Double_t vx,vy,vz, px,py,pz;
  cas->GetBPxPyPz(px,py,pz);
  myCas->SetBachP(px,py,pz);
  cas->GetNPxPyPz(px,py,pz);
  myCas->SetNegP(px,py,pz);
  cas->GetPPxPyPz(px,py,pz);
  myCas->SetPosP(px,py,pz);

  cas->GetXYZ(vx,vy,vz); // v0 decay vertex
  myCas->SetV0vtx(vx,vy,vz);
  cas->GetXYZcascade(vx,vy,vz); // cascade decay vertex
  myCas->SetCascadeVtx(vx,vy,vz);

  Double_t primx = primVtx->GetXv(),
    primy = primVtx->GetYv(),
    primz = primVtx->GetZv();

  myCas->SetCasCosPointingAngle( cas->GetCascadeCosineOfPointingAngle(primx,primy,primz) );
  myCas->SetDecayLength(primx, primy, primz);
  myCas->SetDCA_v0_Bach(cas->GetDcaXiDaughters());

  Float_t p = neg->GetP(), mc = neg->GetMass();
  Float_t betaNeg = p/TMath::Sqrt(p*p + mc*mc);
  p = pos->GetP(); mc = pos->GetMass();
  Float_t betaPos = p/TMath::Sqrt(p*p + mc*mc);
  p = bach->GetP(); mc = bach->GetMass();
  Float_t betaBach = p/TMath::Sqrt(p*p + mc*mc);
  if (bach->GetSign()<0) betaBach = -betaBach; // sign is stored is this parameter

  myCas->SetBeta(betaNeg, betaPos, betaBach);

  return myCas;
}




Reve::CascadeList* esd_Cascade(Double_t min_pt=0.1, Double_t max_pt=100)
{

  AliESD* esd = Alieve::Event::AssertESD();
  AliESDVertex* primVertex =(AliESDVertex*) esd->GetVertex();

  Reve::CascadeList* cont = new Reve::CascadeList("ESD cascade"); 
  cont->SetMainColor(Color_t(3)); // green
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  gReve->AddRenderElement(cont);

  Int_t count = 0;
  //for (Int_t n=0; count<3; n++) {
  for (Int_t n=0; n<esd->GetNumberOfCascades(); n++) {

    AliESDcascade *cas = esd->GetCascade(n);
    Int_t negInd = cas->GetNindex();
    Int_t posInd = cas->GetPindex();
    Int_t bachInd = cas->GetBindex();
    AliESDtrack* negTr = esd->GetTrack(negInd);
    AliESDtrack* posTr = esd->GetTrack(posInd);
    AliESDtrack* bachTr = esd->GetTrack(bachInd);

    if (cas) {
      Reve::Cascade* myCas = esd_make_cas(rnrStyle, primVertex, cas,
					  negTr, posTr, bachTr, n);
      if (myCas) {
	gReve->AddRenderElement(cont, myCas);
	count++;
      }
    }
  }

  cont->SetTitle("CascadeList");
  cont->UpdateItems();

  cont->MakeCascades();
  gReve->Redraw3D();

  return cont;
}
