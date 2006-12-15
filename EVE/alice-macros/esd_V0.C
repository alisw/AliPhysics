

// #include "EVE/Alieve/EventAlieve.h"
// #include "Reve/RGTopFrame.h"
// #include "Reve/V0.h"

// #include "AliESD.h"
// #include "AliESDtrack.h"
// #include "AliESDv0.h"
// #include "AliESDVertex.h"

// using namespace Reve;
// using namespace Alieve;


Reve::V0* esd_make_v0(Reve::TrackRnrStyle* rnrStyle, AliESDVertex* primVtx,
		      AliESDtrack* neg, AliESDtrack* pos, AliESDv0* v0, Int_t i)
{
  if (! v0->GetOnFlyStatus())
  { // v0 on fly do not have the momentum vector filled...
    Reve::RecTrack  rcPos;
    Reve::RecTrack  rcNeg;
    Reve::RecV0 rcV0;

    Double_t p[3];
    v0->GetNPxPyPz(p[0], p[1], p[2]);
    rcV0.P_pos.Set(p);
    v0->GetPPxPyPz(p[0], p[1], p[2]);
    rcV0.P_neg.Set(p);

    v0->GetPxPyPz(p[0], p[1], p[2]);
    Double_t v[3];
    v0->GetXYZ(v[0], v[1], v[2]);


    //   printf(" %f %f %f / %f %f %f    %i\n",p[0], p[1], p[2],
    // 	 v[0], v[1], v[2], v0->GetOnFlyStatus());

    rcV0.V_neg.Set(v); //original track vertices at dca not stored 
    rcV0.V_pos.Set(v);
    rcV0.V_ca.Set(v);

    rcV0.d_label[0] = v0->GetNindex();
    rcV0.d_label[1] = v0->GetPindex();

    Double_t ep = neg->GetP(), mc = neg->GetMass();
    rcNeg.beta = ep/TMath::Sqrt(ep*ep + mc*mc);
    ep = pos->GetP(); mc = pos->GetMass();
    rcPos.beta = ep/TMath::Sqrt(ep*ep + mc*mc);


    Reve::V0* myV0 = new Reve::V0(&rcNeg, &rcPos, &rcV0, rnrStyle);
    char ch[50];
    //   sprintf(ch,"ESDv0%i",i); 
    //   myV0->SetName(ch);
    //   myV0->SetTitle(ch);
    myV0->SetESDIndex(i);
    myV0->SetDaughterDCA(v0->GetDcaV0Daughters());

    Double_t primx = primVtx->GetXv(),
      primy = primVtx->GetYv(),
      primz = primVtx->GetZv();
    myV0->SetCosPointingAngle(v0->GetV0CosineOfPointingAngle(primx,primy,primz));

    myV0->SetDecayLength(primVtx->GetXv(), primVtx->GetYv(), primVtx->GetZv());

    return myV0;
  } else {
    return 0;
  }

}


Reve::V0List* esd_V0(Double_t min_pt=0.1, Double_t max_pt=100)
{

  AliESD* esd = Alieve::Event::AssertESD();
  AliESDVertex* primVertex =(AliESDVertex*) esd->GetVertex();

  Reve::V0List* cont = new Reve::V0List("ESD v0"); 
  cont->SetMainColor(Color_t(3)); // green
  Reve::TrackRnrStyle* rnrStyle = cont->GetRnrStyle();
  rnrStyle->SetMagField( esd->GetMagneticField() );

  gReve->AddRenderElement(cont);

  Int_t count = 0;
  //for (Int_t n=0; count<3; n++) {
  for (Int_t n=0; n<esd->GetNumberOfV0s(); n++) {

    AliESDv0 *v0 = esd->GetV0(n);
    if (v0->GetOnFlyStatus()) continue;

    Int_t negInd = v0->GetNindex();
    Int_t posInd = v0->GetPindex();
    AliESDtrack* negTr = esd->GetTrack(negInd);
    AliESDtrack* posTr = esd->GetTrack(posInd);
    
    Reve::V0* myV0 = esd_make_v0(rnrStyle, primVertex, negTr,posTr, v0, n);
    if (myV0) {
      gReve->AddRenderElement(cont, myV0);
      count++;
    }
  }

  cont->SetTitle("testV0List ");
  cont->UpdateItems();

  cont->MakeV0s();
  cont->MakeMarkers();
  gReve->Redraw3D();

  return cont;
}
















