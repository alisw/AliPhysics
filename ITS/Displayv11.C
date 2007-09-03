//----------------------------------------------------------------------

AliITSv11GeometrySPD *gspd;
AliITSv11GeometrySDD *gsdd;
AliITSv11GeometrySupport *gsupp;
AliITSv11GeometrySSD *gssd;
//
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
    TGeoVolume *ALIC = mgr2->MakeBox("ALIC",vacmed,1000.,1000.,2000.);
    mgr2->SetTopVolume(ALIC);
    TGeoVolume *ITS = mgr2->MakeBox("ITSV",vacmed,990.,990.,1990.);
    TGeoVolumeAssembly *ITSspd = new TGeoVolumeAssembly("ITSspd");
    ITS->AddNode(ITSspd,1);
    ALIC->AddNode(ITS,1);
    //
    /*
    AliITSv11 *its = new AliITSv11(0,3);
    its->SetDebug(ISetits(0,-1));
    its->GetSPDGeometry()->SetDebug(ISetits(0,-1));
    its->GetSupGeometry()->SetDebug(ISetits(0,-1));
    its->CreateMaterials();
    its->CreateGeometry();
    */
    gspd = new AliITSv11GeometrySPD(0);
    //gsdd = new AliITSv11GeometrySDD();
    gsupp = new AliITSv11GeometrySupport(0);
    //gssd = new AliITSv11GeometrySSD();
    //
    Int_t imat=1,imed=1,ireturn=0;
    ireturn = gspd->CreateSPDCentralMaterials(imed,imat);
    gspd->SPDSector(ITSspd,mgr2);
    gsupp->SPDCone(ITS,mgr2);
    gsupp->SetDebug(0);
    gsupp->SDDCone(ITS,mgr2);
    //gsdd->Layer3(ITS);
    //gsdd->Layer4(ITS);
    gsupp->SSDCone(ITS,mgr2);
    gsupp->ServicesCableSupport(ITS);
    //
    mgr2->CloseGeometry();
    //
    TControlBar *bar=new TControlBar("vertical","ITS Geometry Display",10,10);
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
    bar->AddButton("Display SPD Sector Volume","EngineeringSPDSector()",
                   "Run EngineeringSPDSector");
    bar->AddButton("Display SPD Ceneral Volume","EngineeringSPDCenter()",
                   "Run EngineeringSPDCenter");
    bar->AddButton("Display SPD Thermal Sheald","EngineeringSPDThS()",
                   "Run EngineeringSPDThS");
    bar->AddButton("Display SPD Sector Cross Sections","EngineeringSPDSCS()",
                   "Run EngineeringSPDSCS");
    bar->AddButton("Display SDD Layer 3","EngineeringSDDLayer3()",
                   "Run EngineeringSDDLayer3");
    bar->AddButton("Display SDD Layer 4","EngineeringSDDLayer4()",
                   "Run EngineeringSDDLayer4");
    bar->AddButton("Display SDD Cone","EngineeringSDDCone()",
                   "Run EngineeringSDDCone");
    bar->AddButton("Display SDD Centeral Cylinder","EngineeringSDDCylinder()",
                   "Run EngineeringSDDCylinder");
    bar->AddButton("Display SSD Cone","EngineeringSSDCone()",
                   "Run EngineeringSSDCone");
    bar->AddButton("Display SSD Centeral Cylinder","EngineeringSSDCylinder()",
                   "Run EngineeringSSDCylinder");
    bar->AddButton("Display SUP RB24 side","EngineeringSupRB24()",
                   "Run EngineeringSupRB24");
    bar->AddButton("Display Cable Trays RB24 side","EngineeringSupTrayRB24()",
                   "Run EngineeringSupTrayRB24");
    bar->AddButton("Display SUP RB26 side","EngineeringSupRB26()",
                   "Run EngineeringSupRB26");
    bar->AddButton("Save Geometry to File","ExportToFile()",
                   "Run ExportToFile");
    bar->AddButton("Quit/Exit",".q","Exit");
    bar->Show();
    gROOT->SaveContext();
         //Displayit();
}
//----------------------------------------------------------------------
void ExportToFile(){
    // Quirry file name and write geometry to a root file.
    // Inputs:
    //    const char* filename output file with the display in it
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Char_t filename[100];

    printf("Eneter File name:");
    scanf("%s",filename);

    gGeoManager->Export(filename);
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
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
        if(ISetits(4,-1)==0) view3->SetParallel();
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
        if(ISetits(4,-1)==0) view4->SetParallel();
        else  view4->SetPerspective();
        view4->Side();
        if(ISetits(3,-1)!=0) view4->ShowAxis();
    } // end if view4
    if(ISetits(5,-1)==1) ALIC->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSPDSCS(){
    // Display SPD Carbon Fiber Sector Cross sections A and B
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    TPolyLine *a0,*a1,*b0,*b1;
    TPolyMarker *p;
    TCanvas *cSPDSCS=0;
    Int_t i;
    Double_t max=0.0;

    a0 = new TPolyLine();
    a1 = new TPolyLine();
    b0 = new TPolyLine();
    b1 = new TPolyLine();
    p = new TPolyMarker();
    a0->SetLineColor(1);
    a1->SetLineColor(4);
    b0->SetLineColor(3);
    b0->SetLineStyle(2); // dashed
    b1->SetLineColor(6);
    b1->SetLineStyle(2); // dashed
    p->SetMarkerColor(2);
    p->SetMarkerStyle(5);
    if(gspd->Make2DcrossSections(*a0,*a1,*b0,*b1,*p)==kFALSE) return;
    for(i=0;i<a0->GetN();i++) {
      if(TMath::Abs(a0->GetX()[i])>max) max = TMath::Abs(a0->GetX()[i]);
      if(TMath::Abs(a0->GetY()[i])>max) max = TMath::Abs(a0->GetY()[i]);
    } // end for i
    for(i=0;i<a1->GetN();i++) {
      if(TMath::Abs(a1->GetX()[i])>max) max = TMath::Abs(a0->GetX()[i]);
      if(TMath::Abs(a1->GetY()[i])>max) max = TMath::Abs(a1->GetY()[i]);
    } // end for i
    for(i=0;i<b0->GetN();i++) {
      if(TMath::Abs(b0->GetX()[i])>max) max = TMath::Abs(b0->GetX()[i]);
      if(TMath::Abs(b0->GetY()[i])>max) max = TMath::Abs(b0->GetY()[i]);
    } // end for i
    for(i=0;i<b1->GetN();i++) {
      if(TMath::Abs(b1->GetX()[i])>max) max = TMath::Abs(b1->GetX()[i]);
      if(TMath::Abs(b1->GetY()[i])>max) max = TMath::Abs(b1->GetY()[i]);
    } // end for i
    max *= 1.05;
    cSPDSCS = gROOT->FindObject("cSPDSCS");
    if(cSPDSCS==0) delete cSPDSCS;
    cSPDSCS = new TCanvas("cSPDSCS","SPD Carbon Fiber Sector Cross sections",2);
    cSPDSCS->Range(-max,-max,max,max);
    //cSPDSCS->RangeAxis();
    cSPDSCS->SetFixedAspectRatio(kTRUE);
    //
    a0->Draw("");
    a1->Draw("same");
    p->Draw("same");
    b0->Draw("same");
    b1->Draw("same");
    cSPDSCS->Update();
    return;
}
//----------------------------------------------------------------------
void EngineeringSPDSector(){
    // Display SPD Sector Geometry
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
        c4 = new TCanvas("C4","ITS SPD Layer Geometry Side View",500,500);
    TGeoVolume *ITS,*ITSspd,*SPDLay=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSspd_1");
    ITSspd = node->GetVolume();
    node = ITSspd->FindNode("ITSSPDCarbonFiberSectorV_1");
    //node = ITSspd->FindNode("ITSSPDTempSPDMotherVolume_1");
    SPDLay = node->GetVolume();
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
    SPDLay->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SPDLay->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SPD Sector Geometry End View",500,500);
    SPDLay->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SPDLay->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSPDCenter(){
    // Display SPD Centeral Geometry
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
        c4 = new TCanvas("C4","ITS SPD Layer Geometry Side View",500,500);
    TGeoVolume *ITS,*ITSspd,*SPDLay=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSspd_1");
    ITSspd = node->GetVolume();
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
    ITSspd->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SPDLay->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SPD Centeral Geometry End View",500,500);
    ITSspd->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) ITSspd->Raytrace();
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
    node = ITS->FindNode("ITSspdThermalShealdMother_1");
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SPDThS->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSDDLayer3(){
    // Display SDD Layer 3 Geometry
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
        c4 = new TCanvas("C4","ITS SDD Layer 3 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay3=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsddLayer3_1");
    SDDLay3 = node->GetVolume();
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
    SDDLay3->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDDLay3->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Layer 3 Geometry End View",500,500);
    SDDLay3->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SDDLay3->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSDDLayer4(){
    // Display SDD Layer 4 Geometry
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
        c4 = new TCanvas("C4","ITS SDD Layer 4 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay4=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsddLayer4_1");
    SDDLay4 = node->GetVolume();
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
    SDDLay4->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDDLay4->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Layer 4 Geometry End View",500,500);
    SDDLay4->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SDDLay4->Raytrace();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
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
    //Double_t Rmin = ((TGeoPcon*)(SSD->GetShape()))->GetRmin();
    //Double_t Rmax = ((TGeoPcon*)(SSD->GetShape()))->GetRmax();
    //Double_t Dz   = ((TGeoPcon*)(SSD->GetShape()))->GetDz();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SUPRB24->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSupTrayRB24(){
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
    TCanvas *c4,*c5;
    //if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS RB24 Cable Trays and Patch Pannels",500,500);
    //c4->Divide(2,1);
    TGeoVolume *ITS,*SUPRB24=0;
    TGeoNode *node;
    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSsupCableTrayMotherMT24_1");
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
    //TPad *p1 = c4->GetPad(1);
    //TView *view1 = p1->GetView();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SUPRB24->Raytrace();
    //
    //c4->cd(2);
        c5 = new TCanvas("C5","ITS RB24 Cable Trays and Patch Pannels",500,500);
    c5->cd(1);
    SUPRB24->Draw();
    //TPad *p2 = c5->GetPad(1);
    //TView *view2 = p2->GetView();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
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
        if(ISetits(4,-1)==0) view1->SetParallel();
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
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SUPRB26->Raytrace();
    //
}