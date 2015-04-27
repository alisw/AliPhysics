
const char** PPSys() { return x;}
UShort_t*    PPSNN() { return x;}
const char*  PPTrg() { return x;}
const char*  PASys() { return x;}
UShort_t*    PASNN() { return x;}
const char*  PATrg() { return x;}
const char*  APTrg() { const char* x[] = { "CENTV0X", 0 };            return x;}


void DrawAll(UShort_t which=1)
{
  const char* fwd = 0;
  if (gSystem->Getenv("FWD"))
    fwd = gSystem->Getenv("FWD");
  else 
    fwd = gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  gROOT->SetMacroPath(Form("%s/dndeta:%s", fwd, gROOT->GetMacroPath()));

  if (!gROOT->GetClass("Drawer"))  gROOT->LoadMacro("Drawer.C+g");

  switch(which) {
  case 1: {
    const char* s[] = { "pp", 0 };                 
    UShort_t    e[] = { 900, 2760, 7000, 8000, 0 };
    const char* t[] = { "INEL", "NSD", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  case 2: {
    const char* s[] = { "pPb", "Pbp", 0 };         
    UShort_t    e[] = { 5023, 0 };                 
    const char* t[] = { "CENTZNX", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  case 3: {
    const char* s[] = { "pPb", "Pbp", 0 };         
    UShort_t    e[] = { 5023, 0 };                 
    const char* t[] = { "CENTV0X", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  case 4: {
    const char* s[] = { "Pbp", 0 };         
    UShort_t    e[] = { 5023, 0 };                 
    const char* t[] = { "CENTV0C", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  case 5: {
    const char* s[] = { "pPb", "Pbp", 0 };         
    UShort_t    e[] = { 5023, 0 };                 
    const char* t[] = { "CENTV0M", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  case 6: {
    const char* s[] = { "pPb", "Pbp", 0 };         
    UShort_t    e[] = { 5023, 0 };                 
    const char* t[] = { "CENTV0M", "CENTZNX", "CENTV0C", 0 };
    Drawer::DrawAll(s, e, t, 0x1);
  }
    break;
  }
}
