/** @file    DrawFMD.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:18:21 2006
    @brief   Script to draw the FMD3 - obsolete
*/
void DrawFMD(const char* file="geometry.root", const char* option="ogl")
{
#if 1
  AliGeomManager::LoadGeometry(file);
  TGeoVolume* top = gGeoManager->GetTopVolume();
  top->InvisibleAll(kTRUE);
  for (Int_t i = 1; i <= 3; i++) { 
    for (Int_t j = 0; j < 2; j++) { 
      TString name(Form("F%dM%c", i, (j == 0 ? 'T' : 'B')));
      TGeoVolume* v = gGeoManager->FindVolumeFast(name.Data());
      if (!v) { 
	std::cerr << "FMD" << i << " " 
		  << (j == 0 ? "top" : "bottom") 
		  << " Volume " << name << " not found" << std::endl;
	continue;
      }
      v->InvisibleAll(kFALSE);
      v->SetVisDaughters(kTRUE);
      v->SetVisLeaves(kTRUE);
    }
  }
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
