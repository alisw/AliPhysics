void SetViewVolumes(const char* opt){
    //
    if(!strstr(opt,"all")) return;
    //
}
//----------------------------------------------------------------------
void Displayv11(const char* filename=""){
    // Display AliITSv11 Geometry
    // Inputs:
    //    const char* filename output file with the display in it
    // Outputs:
    //    none.
    // Retrurn:
    //    none.

    gSystem->Load("libGeom");
    //
    TCanvas *c1 = new TCanvas("C1","ITS Simulation Geometry",400,400);
    //
    if(gGeoManager) delete gGeoManager;
    TGeoManager *mgr2 = gGeoManager = new TGeoManager("ITSGeometry",
				       " ITS Simulation Geometry Manager");
    //
    TGeoMaterial *vacmat = new TGeoMaterial("Vacume",0,0,0);
    TGeoMedium   *vacmed = new TGeoMedium("Vacume_med",1,vacmat);
    TGeoVolume *ALIC = mgr2->MakeBox("ALIC",vacmed,100.,100.,200.);
    mgr2->SetTopVolume(ALIC);
    //
    //AliITSv11 *its = new AliITSv11("ITS Simulation Volumes");
    AliITSv11 *its = new AliITSv11();
    its->SetDebug(1);
    its->CreateMaterials();
    its->CreateGeometry();
    //SetViewVolumes("all");
    //
    mgr2->CloseGeometry();
    //
    //mgr2->SetVisOption(0);
    //
    TView *view = gPad->GetView();
    if(view){
	view->RotateView(0,90);
	view->ShowAxis();
    } // end if
    ALIC->Draw();
    //
}
