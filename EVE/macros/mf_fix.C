// Fix magnetic field ... so that we can display
// old data from alien.

void mf_fix(Double_t f=-1.)
{
  AliMagF* fld = new AliMagF("Map", "Map", f, 1.0, AliMagF::k5kG, AliMagF::kNoBeamField, -1.0);
  TGeoGlobalMagField::Instance()->SetField(fld);
  TGeoGlobalMagField::Instance()->Lock();
}
