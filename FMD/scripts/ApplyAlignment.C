//____________________________________________________________________
//
// $Id$
//
// Read in the geometry, and get alignment data from CDB, and  apply
// that to the geometry. 
//
/** Apply alignment to a geometry 
    @ingroup FMD_simple_script
 */
void
ApplyAlignment()
{
  gSystem->Load("libFMDutil");
  TGeoManager::Import("geometry.root");

  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://cdb");
  AliCDBEntry*   align = cdb->Get("FMD/Align/Data");
  if (align) {
    Info("ApplyAlignment","Got alignment data from CDB");
    TClonesArray* array = dynamic_cast<TClonesArray*>(align->GetObject());
    if (!array) {
      Warning("ApplyAlignement", "Invalid align data from CDB");
    }
    else {
      Int_t nAlign = array->GetEntries();
      for (Int_t i = 0; i < nAlign; i++) {
	AliAlignObjParams* a = static_cast<AliAlignObjParams*>(array->At(i));
	if (!a->ApplyToGeometry()) {
	  Warning("ApplyAlignement", "Failed to apply alignment to %s", 
		  a->GetVolPath());
	}
      }
    }
  }
  TCanvas* c = new TCanvas("Geometry", "Geometry");
  c->SetFillColor(1);
  gGeoManager->GetTopVolume()->Draw();
}
//____________________________________________________________________
//
// EOF
// 
