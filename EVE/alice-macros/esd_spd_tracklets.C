// To use when per-line id is supported
class TrackletId : public TObject
{
public:
  // label, phi, theta
  // virtual void Print(const Option_t* opt="") {}
};

Reve::StraightLineSet* esd_spd_tracklets(Float_t rad=8)
{
  AliESD         * esd = Alieve::Event::AssertESD();
  AliESDVertex   * pv  = esd->GetPrimaryVertex();
  AliMultiplicity* mul = esd->GetMultiplicity();

  Double_t pvx[3], pve[3];
  pv->GetXYZ(pvx);
  pv->GetSigmaXYZ(pve);

  Reve::StraightLineSet* ls = new Reve::StraightLineSet();

  for (Int_t i=0; i<mul->GetNumberOfTracklets(); ++i)
  {
    using namespace TMath;
    Float_t dr[3];
    dr[0] = rad*Cos(mul->GetPhi(i));
    dr[1] = rad*Sin(mul->GetPhi(i));
    dr[2] = rad/Tan(mul->GetTheta(i));
    ls->AddLine(pvx[0], pvx[1], pvx[2],
		pvx[0]+dr[0], pvx[1]+dr[1], pvx[2]+dr[2]);
  }

  gReve->AddRenderElement(ls);
  gReve->Redraw3D();

  return ls;
}
