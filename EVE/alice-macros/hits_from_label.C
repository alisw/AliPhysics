// $Id$

void hits_from_label(Int_t label=0)
{
  Reve::PointSet* h;
  Reve::LoadMacro("its_hits.C");
  h = its_hits("fX:fY:fZ", Form("ITS.fTrack==%d", label));
  h->SetMarkerSize(1);

  Reve::LoadMacro("tpc_hits.C");
  h = tpc_hits("TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
	       Form("TPC2.fArray.fTrackID==%d", label));
  h->SetMarkerSize(1);

  gReve->Redraw3D();
}
