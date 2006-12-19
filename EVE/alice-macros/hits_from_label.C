// $Id$

void hits_from_label(Int_t label=0)
{
  Reve::PointSet* h;
  Reve::LoadMacro("its_hits.C");
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  h = its_hits("fX:fY:fZ", Form("ITS.fTrack==%d", label));
  char form[1000];
  sprintf(form,"ITS.fTrack==%d", label);
  h = its_hits("fX:fY:fZ", form);
  h->SetMarkerSize(1);

  Reve::LoadMacro("tpc_hits.C");
  sprintf(form,"TPC2.fArray.fTrackID==%d", label);
  h = tpc_hits("TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",form);
  //PH  h = tpc_hits("TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
  //PH	       Form("TPC2.fArray.fTrackID==%d", label));
  h->SetMarkerSize(1);

  Reve::LoadMacro("trd_hits.C");
  sprintf(form,"TRD.fTrack==%d", label);
  h = trd_hits("fX:fY:fZ", form);
  h->SetMarkerSize(1);

  gReve->Redraw3D();
}
