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
    if(gGeoManager) delete gGeoManager;
    mgr2 = gGeoManager = new TGeoManager("ITSGeometry",
                                         " ITS Simulation Geometry Manager");
    //
    TGeoMaterial *vacmat = new TGeoMaterial("Vacume",0,0,0);
    TGeoMedium   *vacmed = new TGeoMedium("Vacume_med",1,vacmat);
    TGeoVolume *ALIC = mgr2->MakeBox("ALIC",vacmed,100.,100.,200.);
    mgr2->SetTopVolume(ALIC);
    //
    AliITSv11 *its = new AliITSv11();
    its->SetDebug(ISetits(0,-1));
    its->CreateMaterials();
    its->CreateGeometry();
    //
    mgr2->CloseGeometry();
    //
    TControlBar *bar=new TControlBar("vertical","ITS Geometry Display",10,10);
    bar->AddButton("Set ITS Debug level 1","ISetits(0,1)","Debug on");
    bar->AddButton("Set ITS Debug level 0","ISetits(0,0)","Debug off");
    bar->AddButton("Set Clipping on","ISetits(2,1)","Clipping on");
    bar->AddButton("Set Cllipping off","ISetits(2,0)","Clipping off");
    bar->AddButton("Set axis on","ISetits(3,1)","Show Axis on");
    bar->AddButton("Set axis off","ISetits(3,0)","Show Axis off");
    bar->AddButton("Set perspective on","ISetits(4,1)","Perspective on");
    bar->AddButton("Set perspective off","ISetits(4,0)","Perspective off");
    bar->AddButton("Set RayTrace on","ISetits(5,1)","Perspective on");
    bar->AddButton("Set RayTrace off","ISetits(5,0)","Perspective off");
    bar->AddButton("Set circle/80","ISetits(1,80)","circles ~ by 80 lines");
    bar->AddButton("Display Geometry","Displayit()","Run Displayit");
    bar->AddButton("Display SPD Thermal Sheald","EngineeringSPDThS()",
                   "Run EngineeringSPDThS");
    bar->AddButton("Display SDD Cone","EngineeringSDDCone()",
                   "Run EngineeringSDDCone");
    bar->AddButton("Display SDD Centeral Cylinder","EngineeringSDDCylinder()",
                   "Run EngineeringSDDCylinder");
    bar->AddButton("Display SSD Cone","EngineeringSSDCone()",
                   "Run EngineeringSSDCone");
    bar->AddButton("Display SSD Centeral Cylinder","EngineeringSSDCylinder()",
                   "Run EngineeringSSDCylinder");
    bar->AddButton("Display SUP RB24 side","EngineeringSupRB24()",
                   "Run EngineeringSDDCylinder");
    bar->AddButton("Display SUP RB26 side","EngineeringSupRB26()",
                   "Run EngineeringSupRB26");
    bar->AddButton("Quit/Exit",".q","Exit");
    bar->Show();
    gROOT->SaveContext();
         //Displayit();
}
//----------------------------------------------------------------------
Int_t ISetits(Int_t t,Int_t v){
    static Int_t itsdebug=1,nsegments=80,cut=0,axis=1,perspective=0,ray=0;

    switch (t) {
    case 0:
        if(v<0) return itsdebug;
        itsdebug = v;
        break;
    case 1:
        if(v<0) return nsegments;
        nsegments= v;
        break;
    case 2:
        if(v<0) return cut;
        cut = v;
        break;
    case 3:
        if(v<0) return axis;
        axis = v;
        break;
    case 4:
        if(v<0) return perspective;
        perspective = v;
        break;
    case 5:
        if(v<0) return ray;
        ray = v;
        break;
    }// end switch
    return 0;
}
//----------------------------------------------------------------------
Double_t DSetits(Int_t t,Double_t v){
    static Double_t phimincut=0.0,phimaxcut=180.0;
    static Double_t longitude=90.0,latitude=0.0,psi=0.0;

    switch (t) {
    case 0:
        if(v<0.) return phimincut;
        phimincut = v;
        break;
    case 1:
        if(v<0.) return phimaxcut;
        phimaxcut = v;
        break;
    case 2:
        if(v<0.) return longitude;
        longitude = v;
        break;
    case 3:
        if(v<0.) return latitude;
        latitude = v;
        break;
    case 4:
        if(v<0.) return latitude;
        latitude = v;
        break;
    }// end switch
    return 0;
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
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c1;
    if(!(c1 = (TCanvas*)gROOT->FindObject("C1")))
        c1 = new TCanvas("C1","ITS Simulation Geometry",900,900);
    c1->Divide(2,2);
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c1->cd(1);
    ALIC->Draw();
    TPad *p1 = c1->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) ALIC->Raytrace();
    c1->cd(2);
    ALIC->Draw();
    TPad *p2 = c1->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->RotateView(60.,30.);
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) ALIC->Raytrace();
    c1->cd(3);
    ALIC->Draw();
    c1->SetPhi(90.0); c1->SetTheta(90.0);
    TPad *p3 = c1->GetPad(3);
    TView *view3 = p3->GetView();
    if(view3){
        view3->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view3->SetParralel();
        else  view3->SetPerspective();
        view3->Top();
        if(ISetits(3,-1)!=0) view3->ShowAxis();
    } // end if view3
    if(ISetits(5,-1)==1) ALIC->Raytrace();
    c1->cd(4);
    ALIC->Draw();
    TPad *p4 = c1->GetPad(4);
    TView *view4 = p4->GetView();
    if(view4){
        view4->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view4->SetParralel();
        else  view4->SetPerspective();
        view4->Side();
        if(ISetits(3,-1)!=0) view4->ShowAxis();
    } // end if view4
    if(ISetits(5,-1)==1) ALIC->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSPDThS(){
    // Display SPD Thermal Sheald Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c4;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Cylinder Geometry",900,450);
    c4->Divide(2,1);
    TGeoVolume *ITS,*SPDThS=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSspdThermalSheald_1");
    SPDThS = node->GetVolume();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c4->cd(1);
    SPDThS->Draw();
    TPad *p1 = c4->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SPDThS->Raytrace();
    //
    c4->cd(2);
    SPDThS->Draw();
    TPad *p2 = c4->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SPDThS->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSDDCone(){
    // Display SDD Cone Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c2;
    if(!(c2 = (TCanvas*)gROOT->FindObject("C2")))
        c2 = new TCanvas("C2","ITS SDD Cone Geometry",900,450);
    c2->Divide(2,1);
    TGeoVolume *ITS,*SDD=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsddConeL_1");
    SDD = node->GetVolume();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c2->cd(1);
    SDD->Draw();
    TPad *p1 = c2->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDD->Raytrace();
    //
    c2->cd(2);
    SDD->Draw();
    TPad *p2 = c2->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDD->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSDDCylinder(){
    // Display SDD Cylinder Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c3;
    if(!(c3 = (TCanvas*)gROOT->FindObject("C3")))
        c3 = new TCanvas("C3","ITS SDD Cylinder Geometry",900,450);
    c3->Divide(2,1);
    TGeoVolume *ITS,*SDD=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsddCentCylCF_1");
    SDD = node->GetVolume();
    Double_t Rmin = ((TGeoTube*)(SDD->GetShape()))->GetRmin();
    Double_t Rmax = ((TGeoTube*)(SDD->GetShape()))->GetRmax();
    Double_t Dz   = ((TGeoTube*)(SDD->GetShape()))->GetDz();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c3->cd(1);
    SDD->Draw();
    TPad *p1 = c3->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDD->Raytrace();
    //arrow->DrawArrow(1.01*Rmax,-Dz,1.01*Rmax,+Dz);
    //
    c3->cd(2);
    SDD->Draw();
    TPad *p2 = c3->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SDD->Raytrace();
    //arrow->DrawArrow(Rmax,0.0,Rmax,0.0);
    //Double_t s = TMath::Sin(0.7),c = TMath::Cos(0.7);
    //arrow->DrawArrow(-Rmin*c,-Rmin*s,Rmin*c,Rmin*s);
    //
}
//----------------------------------------------------------------------
void EngineeringSSDCone(){
    // Display SSD Cone Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c2;
    if(!(c2 = (TCanvas*)gROOT->FindObject("C2")))
        c2 = new TCanvas("C2","ITS SSD Cone Geometry",900,450);
    c2->Divide(2,1);
    TGeoVolume *ITS,*SSD=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSssdConeA_1");
    SSD = node->GetVolume();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c2->cd(1);
    SSD->Draw();
    TPad *p1 = c2->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Top();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SSD->Raytrace();
    //
    c2->cd(2);
    SSD->Draw();
    TPad *p2 = c2->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Front();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SSD->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSSDCylinder(){
    // Display SSD Cylinder Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c3;
    if(!(c3 = (TCanvas*)gROOT->FindObject("C3")))
        c3 = new TCanvas("C3","ITS SSD Cylinder Geometry",900,450);
    c3->Divide(2,1);
    TGeoVolume *ITS,*SSD=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSssdCentCylCA_1");
    SSD = node->GetVolume();
    Double_t Rmin = ((TGeoTube*)(SSD->GetShape()))->GetRmin();
    Double_t Rmax = ((TGeoTube*)(SSD->GetShape()))->GetRmax();
    Double_t Dz   = ((TGeoTube*)(SSD->GetShape()))->GetDz();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c3->cd(1);
    SSD->Draw();
    TPad *p1 = c3->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SSD->Raytrace();
    //arrow->DrawArrow(1.01*Rmax,-Dz,1.01*Rmax,+Dz);
    //
    c3->cd(2);
    SSD->Draw();
    TPad *p2 = c3->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SSD->Raytrace();
    //arrow->DrawArrow(Rmax,0.0,Rmax,0.0);
    //Double_t s = TMath::Sin(0.7),c = TMath::Cos(0.7);
    //arrow->DrawArrow(-Rmin*c,-Rmin*s,Rmin*c,Rmin*s);
    //
}
//----------------------------------------------------------------------
void EngineeringSupRB24(){
    // Display  RB 24 side cable tray support structure Geometry
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c4;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Cylinder Geometry",900,450);
    c4->Divide(2,1);
    TGeoVolume *ITS,*SUPRB24=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsupFrameM24_1");
    SUPRB24 = node->GetVolume();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c4->cd(1);
    SUPRB24->Draw();
    TPad *p1 = c4->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SUPRB24->Raytrace();
    //
    c4->cd(2);
    SUPRB24->Draw();
    TPad *p2 = c4->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SUPRB24->Raytrace();
    //
}

//----------------------------------------------------------------------
void EngineeringSupRB26(){
    // Display RB 26 side cable tray support structure
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoVolume *ALIC = mgr2->GetTopVolume();
    TCanvas *c5;
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Cylinder Geometry",900,450);
    c5->Divide(2,1);
    TGeoVolume *ITS,*SUPRB26=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsupFrameM26_1");
    SUPRB26 = node->GetVolume();
    //
    mgr2->SetNsegments(ISetits(1,-1));
    //
    mgr2->SetVisLevel(6);
    mgr2->SetVisOption(0);
    //mgr2->CheckOverlaps(0.01);
    //mgr2->PrintOverlaps();
    if(ISetits(2,-1)==1){
        TGeoShape *clip = new TGeoTubeSeg(0, 1000, 2000, 45, 90);
        mgr2->SetClippingShape(clip);
    } // end if
    mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    if(ISetits(2,-1)!=0) mgr2->SetPhiRange(DSetits(1,-1.),DSetits(0,-1.));
    //
    c5->cd(1);
    SUPRB26->Draw();
    TPad *p1 = c5->GetPad(1);
    TView *view1 = p1->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParralel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SUPRB26->Raytrace();
    //
    c5->cd(2);
    SUPRB26->Draw();
    TPad *p2 = c5->GetPad(2);
    TView *view2 = p2->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParralel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SUPRB26->Raytrace();
    //
}
