// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void esd_v0_init_rectrack(TEveRecTrack& rt, AliExternalTrackParam* tp)
{
  Double_t      pbuf[3], vbuf[3];

  rt.fSign = tp->GetSign();
  tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
  tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
  // Double_t ep = at->GetP(), mc = at->GetMass();
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveV0* esd_make_v0(TEveTrackPropagator* rnrStyle, AliESDVertex* primVtx,
		      AliESDtrack* neg, AliESDtrack* pos, AliESDv0* v0, Int_t i)
{
  TEveRecTrack  rcPos;
  TEveRecTrack  rcNeg;
  TEveRecV0     rcV0;

  Double_t p[3];
  v0->GetNPxPyPz(p[0], p[1], p[2]);
  rcV0.fPPos.Set(p);
  v0->GetPPxPyPz(p[0], p[1], p[2]);
  rcV0.fPNeg.Set(p);

  v0->GetPxPyPz(p[0], p[1], p[2]);

  Double_t v[3];

  v0->GetXYZ(v[0], v[1], v[2]);
  rcV0.fVCa.Set(v);

  v0->GetParamN()->GetXYZ(v);  rcV0.fVNeg.Set(v);
  v0->GetParamP()->GetXYZ(v);  rcV0.fVPos.Set(v);

  rcV0.fV0Birth.Set(primVtx->GetXv(), primVtx->GetYv(), primVtx->GetZv());

  // Simulation data not directly available in AliESDv0
  //rcV0.fDLabel[0] = v0->GetNindex();
  //rcV0.fDLabel[1] = v0->GetPindex();

  esd_v0_init_rectrack(rcNeg, v0->GetParamN());
  rcNeg.fIndex = v0->GetNindex();
  esd_v0_init_rectrack(rcPos, v0->GetParamP());
  rcPos.fIndex = v0->GetPindex();

  AliEveV0* myV0 = new AliEveV0(&rcNeg, &rcPos, &rcV0, rnrStyle);
  myV0->SetElementName(Form("ESDv0 %d", i));
  myV0->SetElementTitle(Form("OnFly: %d\nDCA %f",
                             v0->GetOnFlyStatus(),
                             v0->GetDcaV0Daughters()));
  myV0->SetESDIndex(i);
  myV0->SetOnFlyStatus(v0->GetOnFlyStatus());
  myV0->SetDaughterDCA(v0->GetDcaV0Daughters());

  return myV0;
}


AliEveV0List* esd_V0(Bool_t onFly=kFALSE)
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliESDVertex* primVertex = (AliESDVertex*) esd->GetPrimaryVertex();

  AliEveV0List* cont = new AliEveV0List("ESD v0");
  cont->SetMainColor(3); // green
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );

  gEve->AddElement(cont);

  Int_t count = 0;
  for (Int_t n=0; n<esd->GetNumberOfV0s(); ++n)
  {
    AliESDv0 *v0 = esd->GetV0(n);

    if (v0->GetOnFlyStatus() != onFly) continue;

    Int_t negInd = v0->GetNindex();
    Int_t posInd = v0->GetPindex();
    AliESDtrack* negTr = esd->GetTrack(negInd);
    AliESDtrack* posTr = esd->GetTrack(posInd);

    AliEveV0* myV0 = esd_make_v0(rnrStyle, primVertex, negTr,posTr, v0, n);
    if (myV0)
    {
      gEve->AddElement(myV0, cont);
      ++count;
    }
  }

  cont->SetTitle("test");

  cont->MakeV0s();
  gEve->Redraw3D();

  return cont;
}
















