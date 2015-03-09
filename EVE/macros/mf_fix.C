// Fix magnetic field ... so that we can display
// old data from alien.

void mf_fix(Double_t fSol=-1., Double_t fDip=1.0)
{
  AliMagF* fld = new AliMagF("Map", "Map", fSol,fDip, AliMagF::k5kG, AliMagF::kNoBeamField, -1.0);
  TGeoGlobalMagField::Instance()->SetField(fld);
  TGeoGlobalMagField::Instance()->Lock();
}
