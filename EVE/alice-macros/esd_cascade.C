// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class AliEveCascade;
class AliEveCascadeList;

AliEveCascade* esd_make_cas(TEveTrackPropagator* rnrStyle, AliESDVertex* primVtx,
			    AliESDcascade* cas, AliESDtrack* neg, AliESDtrack* pos,
			    AliESDtrack* bach,Int_t i)
{
  AliEveCascade* myCas = new AliEveCascade(rnrStyle);
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


AliEveCascadeList* esd_cascade(Double_t min_pt=0.1, Double_t max_pt=100)
{
  printf("THIS SCRIPT DOES NOT WORK.\n"
	 "AliEveCascade classes have been temporarily removed.\n"
	 "They need to be cleaned up.\n");
  return;

  AliESDEvent* esd = AliEveEventManager::AssertESD();
  AliESDVertex* primVertex =(AliESDVertex*) esd->GetVertex();

  CascadeList* cont = new CascadeList("ESD cascade");
  cont->SetMainColor(3); // green
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );

  gEve->AddElement(cont);

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
      AliEveCascade* myCas = esd_make_cas(rnrStyle, primVertex, cas,
					  negTr, posTr, bachTr, n);
      if (myCas) {
	gEve->AddElement(myCas, cont);
	count++;
      }
    }
  }

  cont->SetTitle("CascadeList");

  cont->MakeCascades();
  gEve->Redraw3D();

  return cont;
}
