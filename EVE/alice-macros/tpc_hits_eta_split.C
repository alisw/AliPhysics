// $Id$

void tpc_hits_eta_split(const char *varexp    = "TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
			const char *selection = "TPC2.fArray.fR>80",
			Option_t *option      = "goff")
{
  // Extracts 'major' TPC hits (not the compressed ones).
  // This gives ~2.5% of all hits.

  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  rl->LoadHits("TPC");

  TTree* ht = rl->GetTreeH("TPC", false);
  ht->SetEstimate(400*ht->GetEntries());
  ht->Draw(varexp, selection, option);

  gReve->DisableRedraw();

  Reve::PointSetArray* l = new Reve::PointSetArray
    ("TPC hits - Eta Slices", "");
  l->SetMainColor((Color_t)3);
  TGListTreeItem *ti = gReve->AddRenderElement(l);

  l->InitBins(ti, "Eta", 20, -2, 2);

  Double_t *vr = ht->GetV1(), *vphi = ht->GetV2(), *vz = ht->GetV3();
  Long64_t nr  = ht->GetSelectedRows();
  while(nr-- > 0) {
    using namespace TMath;
    Double_t ctg = *vz / *vr;
    Double_t eta = -Log(Hypot(1,ctg)-Abs(ctg)); if(ctg < 0) eta = -eta;
    Double_t cos_theta = *vz / Hypot(*vr, *vz);
    Double_t eta1       = -0.5*Log( (1.0-cos_theta)/(1.0+cos_theta) );
    if(Abs(eta1 - eta) > 0.01) printf("etadiff %lf %lf\n", eta1, eta);
    l->Fill(eta, *vr * Cos(*vphi), *vr * Sin(*vphi), *vz);
    ++vr; ++vphi; ++vz;
  }

  l->CloseBins(20, 1); // Full circle, size 1.

  gReve->DrawRenderElement(l);
  gReve->EnableRedraw();
}
