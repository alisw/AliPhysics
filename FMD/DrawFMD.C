/** @file    DrawFMD.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:21 2006
    @brief   Script to draw the FMD3 - obsolete
*/
void DrawFMD(const char* file="geometry.root", const char* option="ogl")
{
#if 1
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliGeomManager::LoadGeometry(file);
  TGeoVolume* top = gGeoManager->GetTopVolume();
  TGeoIterator next(top);
  TGeoNode* node = 0;
  while ((node = next())) {
    TGeoVolume* v = node->GetVolume();
    v->SetVisibility(kFALSE);
    v->InvisibleAll(kTRUE);
    v->VisibleDaughters(kFALSE);

    TString name(v->GetName());
    if (name[0] == 'F') { 
      if (name[1] == '1' || name[1] == '2' || name[1] == '3' ||
	  name[1] == 'I' || name[1] == 'O' || name[1] == 'M') {
	 v->SetVisibility(kTRUE);
	 v->InvisibleAll(kFALSE);
	 v->VisibleDaughters(kTRUE);
	 std::cout << "Making " << name << " visible" << std::endl;
	 continue;
      }
    }
    if (name[0] == 'C') { 
      if ((name[1] == 'P' || name[1] == 'p') && 
	  (name[2] == '1' || name[2] == '2' || name[3] == '3')) {
	v->SetVisibility(kTRUE);
	v->InvisibleAll(kFALSE);
	v->VisibleDaughters(kTRUE);
	std::cout << "Making " << name << " visible" << std::endl;
	continue;	
      }
  }
  // for (
  gGeoManager->SetVisLevel(5);
  top->InvisibleAll(kFALSE);
  top->Draw(option);

#else
  // gSystem->Load("/usr/lib/libshift");
  // gSystem->Load("/usr/lib/libgfortran");
  gSystem->Load("libgeant321");
  gMC = new TGeant3TGeo;
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  gROOT->LoadMacro("FMD/ViewFMD.C");
  gInterpreter->ProcessLine("ViewFMD()");
  gROOT->LoadMacro("ITS/ViewITS.C");
  gInterpreter->ProcessLine("ViewITS()");
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 12, .055, .055);
  gMC->Gdhead(1111, "Forward Multiplicity Detector");
  gMC->Gdman(16, 10, "MAN");
  gPad->Modified();
  gPad->cd();
#endif
}
