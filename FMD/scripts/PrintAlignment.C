//____________________________________________________________________
//
// $Id$
//
// Read in the geometry, and get alignment data from CDB, and  apply
// that to the geometry. 
//
/** Print alignment to a geometry 
    @ingroup FMD_simple_script
 */
#include <iomanip>
void
PrintAlignment()
{
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBEntry*   align = cdb->Get("FMD/Align/Data");
  if (!align) {
    Error("PrintAlignment","didn't alignment data from CDB");
    return;
  }
  
  TClonesArray* array = dynamic_cast<TClonesArray*>(align->GetObject());
  if (!array) {
    Warning("PrintAlignement", "Invalid align data from CDB");
    return;
  }
  Int_t nAlign = array->GetEntries();
  for (Int_t i = 0; i < nAlign; i++) {
    AliAlignObjParams* a = static_cast<AliAlignObjParams*>(array->At(i));
    Double_t ang[3];
    Double_t trans[3];
    a->GetAngles(ang);
    a->GetTranslation(trans);
    std::cout << a->GetVolPath() << "\n" 
	      << "  translation: "
	      << "(" << std::setw(12) << trans[0] 
	      << "," << std::setw(12) << trans[1] 
	      << "," << std::setw(12) << trans[2] << ")\n"
	      << "  rotation:    "
	      << "(" << std::setw(12) << ang[0] 
	      << "," << std::setw(12) << ang[1] 
	      << "," << std::setw(12) << ang[2]  << ")" << std::endl;
    // a->Print();
  }
}
//____________________________________________________________________
//
// EOF
// 
