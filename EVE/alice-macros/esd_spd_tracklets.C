// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// To use when per-line id is supported

class TrackletId : public TObject
{
public:
  // label, phi, theta
  // virtual void Print(const Option_t* opt="") {}
};

TEveStraightLineSet* esd_spd_tracklets(Float_t rad=8)
{
  AliESDEvent         * esd = AliEveEventManager::AssertESD();
  AliESDVertex   * pv  = esd->GetPrimaryVertexSPD();
  AliMultiplicity* mul = esd->GetMultiplicity();

  Double_t pvx[3], pve[3];
  pv->GetXYZ(pvx);
  pv->GetSigmaXYZ(pve);

  TEveStraightLineSet* ls = new TEveStraightLineSet("SPD tracklets");

  for (Int_t i=0; i<mul->GetNumberOfTracklets(); ++i)
  {
    using namespace TMath;
    Float_t dr[3];
    Float_t phi = mul->GetPhi(i);
    dr[0] = rad*Cos(phi);
    dr[1] = rad*Sin(phi);
    dr[2] = rad/Tan(mul->GetTheta(i));
    ls->AddLine(pvx[0], pvx[1], pvx[2],
		pvx[0]+dr[0], pvx[1]+dr[1], pvx[2]+dr[2]);
  }

  ls->SetElementTitle(Form("N=%d", mul->GetNumberOfTracklets()));
  gEve->AddElement(ls);
  gEve->Redraw3D();

  return ls;
}
