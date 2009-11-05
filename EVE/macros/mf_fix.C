// Fix magnetic field ... so that we can display
// old data from alien.

void mf_fix()
{
  AliMagF* fld = new AliMagF("Map", "Map", -1.0, 1.0, AliMagF::k5kG, AliMagF::kNoBeamField, -1.0);
  TGeoGlobalMagField::Instance()->SetField(fld);
  TGeoGlobalMagField::Instance()->Lock();
}
