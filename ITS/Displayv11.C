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

    Displayit();
}
//----------------------------------------------------------------------
void Displayit(){
    // Display AliITSv11 Geometry
    // Inputs:
    //    const char* filename output file with the display in it
    // Outputs:
    //    none.
    // Retrurn:
    //    none.

    gSystem->Load("libGeom");
    //
    TCanvas *c1 = new TCanvas("C1","ITS Simulation Geometry",900,900);
    c1->Divide(2,2);
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
    AliITSv11 *its = new AliITSv11();
    its->SetDebug(1);
    its->CreateMaterials();
    its->CreateGeometry();
    //SetViewVolumes("all");
    //
    mgr2->CloseGeometry();
    mgr2->SetNsegments(80);
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    //mgr2->SetPhiRange(0.0,180.0);
    //
    c1->cd(1);
    ALIC->Draw();
    TPad *p1 = c1->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetParralel();
        view1->Front();
        view1->ShowAxis();
    } // end if view2
    c1->cd(2);
    ALIC->Draw();
    TPad *p2 = c1->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetParralel();
        view2->RotateView(60.,30.);
        view2->ShowAxis();
    } // end if view2
    c1->cd(3);
    c1->SetPhi(90.0); c1->SetTheta(90.0);
    ALIC->Draw();
    TPad *p3 = c1->GetPad(3);
    TView *view3 = p3->GetView();
    if(view3){
        view3->SetParralel();
        view3->Top();
        view3->ShowAxis();
    } // end if view3
    c1->cd(4);
    ALIC->Draw();
    TPad *p4 = c1->GetPad(4);
    TView *view4 = p4->GetView();
    if(view4){
        view4->SetParralel();
        view4->Side();
        view4->ShowAxis();
    } // end if view4
    //
}
