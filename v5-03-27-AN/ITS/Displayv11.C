//----------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TPolyLine.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TGeoManager.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TControlBar.h>
#include <TGeoTube.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TView.h>
#include "AliITSgeom.h"
#include "AliITSInitGeometry.h"
#include "AliITSv11GeometrySPD.h"
#include "AliITSv11GeometrySDD.h"
#include "AliITSv11GeometrySSD.h"
#include "AliITSv11GeometrySupport.h"
#endif


static AliITSv11GeometrySPD *gspd;
static AliITSv11GeometrySDD *gsdd;
static AliITSv11GeometrySSD *gssd;
static AliITSv11GeometrySupport *gsupp;
static AliITSgeom *geom;

void CreateMaterialsITS();

void AliMaterial(Int_t imat, const char* name, Float_t a, 
		 Float_t z, Float_t dens, Float_t radl,
		 Float_t absl);
void AliMedium(Int_t numed, const char *name, Int_t nmat,
	       Int_t isvol, Int_t ifield, Float_t fieldm,
	       Float_t tmaxfd, Float_t stemax, Float_t deemax,
	       Float_t epsil, Float_t stmin);
void AliMixture(Int_t kmat, const char* name, Float_t* a, Float_t* z,
		Double_t dens, Int_t nlmat=0, Float_t* wmat=0);
void Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z,
	     Double_t dens, Int_t nlmat, Double_t* wmat);
Double_t* CreateDoubleArray(Float_t* array, Int_t size);

Bool_t Make2DCrossSections(TPolyLine &a0,TPolyLine &a1,
			   TPolyLine &b0,TPolyLine &b1,TPolyMarker &p);

//
//----------------------------------------------------------------------
void Displayv11(){
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
    gGeoManager = new TGeoManager("ITSGeometry",
				  " ITS Simulation Geometry Manager");
    TGeoManager *mgr2 = gGeoManager;
    //
    const AliITSVersion_t kv11=(AliITSVersion_t)110;
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";
    const Int_t kLength=100;
    Char_t vstrng[kLength];
    AliITSInitGeometry initgeom(kv11,1);
    //
    TGeoMaterial *vacmat = new TGeoMaterial("Vacume",0,0,0);
    TGeoMedium   *vacmed = new TGeoMedium("Vacume_med",1,vacmat);
    TGeoVolume *ALIC = mgr2->MakeBox("ALIC",vacmed,1000.,1000.,2000.);
    mgr2->SetTopVolume(ALIC);
    TGeoVolume *ITS = mgr2->MakeBox("ITSV",vacmed,990.,990.,1990.);
    if(initgeom.WriteVersionString(vstrng,kLength,kv11,1,cvsDate,cvsRevision))
        ITS->SetTitle(vstrng);
    //TGeoVolumeAssembly *ITSSPD = new TGeoVolumeAssembly("ITSSPD");
    //ITS->AddNode(ITSSPD,1);
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
    gspd  = new AliITSv11GeometrySPD(0);
    gsdd  = new AliITSv11GeometrySDD();
    gsupp = new AliITSv11GeometrySupport(0);
    gssd  = new AliITSv11GeometrySSD();
    //
    CreateMaterialsITS();
    gspd->SPDSector(ITS,mgr2);
    gsupp->SPDCone(ITS,mgr2);
    gsupp->SetDebug(0);
    gsupp->SDDCone(ITS,mgr2);
    gsdd->Layer3(ITS);
    gsdd->Layer4(ITS);
    gsdd->ForwardLayer3(ITS);// in Hybrid its in IS02
    gsdd->ForwardLayer4(ITS);// in Hybrid its in IS02
    gssd->Layer5(ITS);
    gssd->Layer6(ITS);
    gssd->LadderSupportLayer5(ITS);
    gssd->LadderSupportLayer6(ITS);
    gssd->EndCapSupportSystemLayer6(ITS);
    gssd->EndCapSupportSystemLayer5(ITS);
    gsupp->SSDCone(ITS,mgr2);
    gsupp->ServicesCableSupport(ITS);
    //
    mgr2->CloseGeometry();
    //
    geom = new AliITSgeom();
    initgeom.InitAliITSgeom(geom);
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
    bar->AddButton("Print SPD Sector Volume data xfig","PrintSPDSectorData()",
                   "Run PrintSPDSectorData");
    bar->AddButton("Display SPD General Volume","EngineeringSPDCenter()",
                   "Run EngineeringSPDCenter");
    bar->AddButton("Display SPD Thermal Shield","EngineeringSPDThS()",
                   "Run EngineeringSPDThS");
    bar->AddButton("Display SPD Sector Cross Sections","EngineeringSPDSCS()",
                   "Run EngineeringSPDSCS");
    bar->AddButton("Display SDD Layer 3","EngineeringSDDLayer3()",
                   "Run EngineeringSDDLayer3");
    bar->AddButton("Display SDD Layer 4","EngineeringSDDLayer4()",
                   "Run EngineeringSDDLayer4");
    bar->AddButton("Display SDD Cone","EngineeringSDDCone()",
                   "Run EngineeringSDDCone");
    bar->AddButton("Display SDD Central Cylinder","EngineeringSDDCylinder()",
                   "Run EngineeringSDDCylinder");
    bar->AddButton("Display SSD Layer 5","EngineeringSSDLayer5()",
                   "Run EngineeringSSDLayer5");
    bar->AddButton("Display SSD Layer 6","EngineeringSSDLayer6()",
                   "Run EngineeringSSDLayer6");
    bar->AddButton("Display SSD Cone","EngineeringSSDCone()",
                   "Run EngineeringSSDCone");
    bar->AddButton("Display SSD Central Cylinder","EngineeringSSDCylinder()",
                   "Run EngineeringSSDCylinder");
//    bar->AddButton("Display SUP RB24 side","EngineeringSupRB24()",
//                   "Run EngineeringSupRB24");
//    bar->AddButton("Display Cable Trays RB24 side","EngineeringSupTrayRB24()",
//                   "Run EngineeringSupTrayRB24");
//    bar->AddButton("Display SUP RB26 side","EngineeringSupRB26()",
//                   "Run EngineeringSupRB26");
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

    printf("Enter File name:");
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
    static Double_t longitude=90.0,latitude=0.0;

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
    //c1->Divide(2,2);
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
    TPad *p1 = (TPad*)c1->GetPad(1);
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
    /*c1->cd(2);
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
    *///
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
    if(gspd==0||Make2DCrossSections(*a0,*a1,*b0,*b1,*p)==kFALSE)return;
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
    cSPDSCS = (TCanvas*)gROOT->FindObject("cSPDSCS");
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SPD Layer Geometry Side View",500,500);
    TGeoVolume *ITS=0,*ITSSPD=0,*SPDLay=0;
    TGeoNode *node=0;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSSPD_1");
    ITSSPD = node->GetVolume();
    node = ITSSPD->FindNode("ITSSPDCarbonFiberSectorV_1");
    if(node==0)Error("EngineeringSPDSector","could not find node %s",
                     "ITSSPDCarbonFiberSectorV_1");
    SPDLay = node->GetVolume();
    if(SPDLay==0)Error("EngineeringSPDSector","could not find volume SPDLay");
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
void PrintSPDSectorData(){
    // Print SPD Sector Data
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Int_t irr,i;
    //
    TGeoManager *mgr2 = gGeoManager;
    TGeoXtru * sA0;
    TGeoVolume *vA0=0;

    //mgr2->PushPath();
    //mgr2->cd("ITSSPDCarbonFiberSupportSectorA0_1");
    vA0 = mgr2->FindVolumeFast("ITSSPDCarbonFiberSupportSectorA0");
    sA0 = (TGeoXtru*) vA0->GetShape();
    irr = sA0->GetNvert();
    Double_t x,y;
    cout <<endl;
    cout <<"2 3 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 "<<irr;
    for(i=0;i<irr;i++){
        x = sA0->GetX(i)+2.5;
        y = sA0->GetY(i)+2.5;
        if(!(i%6)) { cout << endl; cout <<"        ";}
        cout<<" "<<TMath::Nint(x*450.)<<" "<<TMath::Nint(y*450);
        //cout<<" "<<x<<" "<<y;
    } // end for i
    x = sA0->GetX(0)+2.5;
    y = sA0->GetY(0)+2.5;
    if(!(i%6)) { cout << endl; cout <<"        ";}
    cout<<" "<<TMath::Nint(x*450.)<<" "<<TMath::Nint(y*450);
    //cout<<" "<<x<<" "<<y;
    cout << endl;
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SPD Layer Geometry Side View",500,500);
    TGeoVolume *ITS,*ITSSPD,*SPDLay=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSSPD_1");
    ITSSPD = node->GetVolume();
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
    ITSSPD->Draw();
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
    ITSSPD->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) ITSSPD->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSPDThS(){
    // Display SPD Thermal Shield Geometry
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
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSspdThermalShield_1");
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
    TPad *p1 = (TPad*)c4->GetPad(1);
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
    TPad *p2 = (TPad*)c4->GetPad(2);
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Layer 3 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay3=0;
    TGeoNode *node;
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Layer 4 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay4=0;
    TGeoNode *node;
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
    TCanvas *c2=0;
    if(!(c2 = (TCanvas*)gROOT->FindObject("C2")))
        c2 = new TCanvas("C2","ITS SDD Cone Geometry",900,450);
    c2->Divide(2,1);
    TGeoVolume *ITS,*SDD=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("SDDCarbonFiberCone_1");
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
    TPad *p1 = (TPad*)c2->GetPad(1);
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
    TPad *p2 = (TPad*)c2->GetPad(2);
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
    TCanvas *c3=0;
    if(!(c3 = (TCanvas*)gROOT->FindObject("C3")))
        c3 = new TCanvas("C3","ITS SDD Cylinder Geometry",900,450);
    c3->Divide(2,1);
    TGeoVolume *ITS,*SDD=0;
    TGeoNode *node;
//    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("SDDCarbonFiberCylinder_1");
    SDD = node->GetVolume();
//    Double_t Rmin = ((TGeoTube*)(SDD->GetShape()))->GetRmin();
//    Double_t Rmax = ((TGeoTube*)(SDD->GetShape()))->GetRmax();
//    Double_t Dz   = ((TGeoTube*)(SDD->GetShape()))->GetDz();
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
    TPad *p1 = (TPad*)c3->GetPad(1);
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
    TPad *p2 = (TPad*)c3->GetPad(2);
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
void EngineeringSSDLayer5(){
    // Display SSD Layer 5 Geometry
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Layer 5 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay5=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSssdLayer5_1");
    SDDLay5 = node->GetVolume();
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
    SDDLay5->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDDLay5->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Layer 5 Geometry End View",500,500);
    SDDLay5->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SDDLay5->Raytrace();
    //
}
//----------------------------------------------------------------------
void EngineeringSSDLayer6(){
    // Display SSD Layer 6 Geometry
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
    TCanvas *c4=0, *c5=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Layer 6 Geometry Side View",500,500);
    TGeoVolume *ITS,*SDDLay6=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSssdLayer6_1");
    SDDLay6 = node->GetVolume();
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
    SDDLay6->Draw();
    TView *view1 = c4->GetView();
    if(view1){
        view1->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view1->SetParallel();
        else  view1->SetPerspective();
        view1->Front();
        if(ISetits(3,-1)!=0) view1->ShowAxis();
    } // end if view1
    if(ISetits(5,-1)==1) SDDLay6->Raytrace();
    //
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Layer 6 Geometry End View",500,500);
    SDDLay6->Draw();
    TView *view2 = c5->GetView();
    if(view2){
        view2->SetView(DSetits(2,-1.),DSetits(3,-1.),DSetits(4,-1.),irr);
        if(irr) cout <<"error="<<irr<<endl;
        if(ISetits(4,-1)==0) view2->SetParallel();
        else  view2->SetPerspective();
        view2->Top();
        if(ISetits(3,-1)!=0) view2->ShowAxis();
    } // end if view2
    if(ISetits(5,-1)==1) SDDLay6->Raytrace();
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
    TCanvas *c2=0;
    if(!(c2 = (TCanvas*)gROOT->FindObject("C2")))
        c2 = new TCanvas("C2","ITS SSD Cone Geometry",900,450);
    c2->Divide(2,1);
    TGeoVolume *ITS,*SSD=0;
    TGeoNode *node;
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("ITSssdCone_1");
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
    TPad *p1 = (TPad*)c2->GetPad(1);
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
    TPad *p2 = (TPad*)c2->GetPad(2);
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
    TCanvas *c3=0;
    if(!(c3 = (TCanvas*)gROOT->FindObject("C3")))
        c3 = new TCanvas("C3","ITS SSD Cylinder Geometry",900,450);
    c3->Divide(2,1);
    TGeoVolume *ITS,*SSD=0;
    TGeoNode *node;
//    TArrow *arrow=new TArrow();
    //
    node = ALIC->FindNode("ITSV_1");
    ITS = node->GetVolume();
    node = ITS->FindNode("SSDexternalcylinder_1");
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
    TPad *p1 = (TPad*)c3->GetPad(1);
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
    TPad *p2 = (TPad*)c3->GetPad(2);
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
    TCanvas *c4=0;
    if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS SDD Cylinder Geometry",900,450);
    c4->Divide(2,1);
    TGeoVolume *ITS,*SUPRB24=0;
    TGeoNode *node;
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
    TPad *p1 = (TPad*)c4->GetPad(1);
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
    TPad *p2 = (TPad*)c4->GetPad(2);
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
    TCanvas *c4=0,*c5=0;
    //if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
        c4 = new TCanvas("C4","ITS RB24 Cable Trays and Patch Pannels",500,500);
    //c4->Divide(2,1);
    TGeoVolume *ITS,*SUPRB24=0;
    TGeoNode *node;
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
    TCanvas *c5=0;
    if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
        c5 = new TCanvas("C5","ITS SDD Cylinder Geometry",900,450);
    c5->Divide(2,1);
    TGeoVolume *ITS,*SUPRB26=0;
    TGeoNode *node;
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
    TPad *p1 = (TPad*)c5->GetPad(1);
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
    TPad *p2 = (TPad*)c5->GetPad(2);
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

// Local functions to replace "regular" ones not available
// in the environment where the macro runs

//______________________________________________________________________
void CreateMaterialsITS(){
    // Create ITS materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSv11Hybrid.
    // In general it is automatically replaced by
    // the CreateMaterials routine defined in AliITSv?. Should the function
    // CreateMaterials not exist for the geometry version you are using this
    // one is used. See the definition found in AliITSv5 or the other routine
    // for a complete definition.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    Bool_t fByThick=kTRUE;    // Flag to use services materials by thickness
                              // ture, or mass false.

//    Int_t   ifield = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
//    Float_t fieldm = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
    Int_t   ifield = 2;
    Float_t fieldm = 10.;

    Float_t tmaxfd = 0.1; // 1.0; // Degree
    Float_t stemax = 1.0; // cm
    Float_t deemax = 0.1; // 30.0; // Fraction of particle's energy 0<deemax<=1
    Float_t epsil  = 1.0E-4; // 1.0; // cm
    Float_t stmin  = 0.0; // cm "Default value used"

    Float_t tmaxfdSi = 0.1; // .10000E+01; // Degree
    Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
    Float_t deemaxSi = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;// .10000E+01;
    Float_t stminSi  = 0.0; // cm "Default value used"

    Float_t tmaxfdAir = 0.1; // .10000E+01; // Degree
    Float_t stemaxAir = .10000E+01; // cm
    Float_t deemaxAir = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;// .10000E+01;
    Float_t stminAir  = 0.0; // cm "Default value used"

    Float_t tmaxfdServ = 1.0; // 10.0; // Degree
    Float_t stemaxServ = 1.0; // 0.01; // cm
    Float_t deemaxServ = 0.5; // 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilServ  = 1.0E-3; // 0.003; // cm
    Float_t stminServ  = 0.0; //0.003; // cm "Default value used"

    // Freon PerFluorobuthane C4F10 see 
    // http://st-support-cooling-electronics.web.cern.ch/
    //        st-support-cooling-electronics/default.htm
    Float_t afre[2]  = { 12.011,18.9984032 };
    Float_t zfre[2]  = { 6., 9. };
    Float_t wfre[2]  = { 4.,10. };
    Float_t densfre  = 1.52;


    //CM55J

    Float_t aCM55J[4]={12.0107,14.0067,15.9994,1.00794};
    Float_t zCM55J[4]={6.,7.,8.,1.};
    Float_t wCM55J[4]={0.908508078,0.010387573,0.055957585,0.025146765};
    Float_t dCM55J = 1.8;

    //ALCM55J

    Float_t aALCM55J[5]={12.0107,14.0067,15.9994,1.00794,26.981538};
    Float_t zALCM55J[5]={6.,7.,8.,1.,13.};
    Float_t wALCM55J[5]={0.817657902,0.0093488157,0.0503618265,0.0226320885,0.1};
    Float_t dALCM55J = 1.9866;

    //Si Chips

    Float_t aSICHIP[6]={12.0107,14.0067,15.9994,1.00794,28.0855,107.8682};
    Float_t zSICHIP[6]={6.,7.,8.,1.,14., 47.};
    Float_t wSICHIP[6]={0.039730642,0.001396798,0.01169634,0.004367771,0.844665,0.09814344903};
    Float_t dSICHIP = 2.36436;

    //Inox
    
    Float_t aINOX[9]={12.0107,54.9380, 28.0855,30.9738,32.066,58.6928,51.9961,95.94,55.845};
    Float_t zINOX[9]={6.,25.,14.,15.,16., 28.,24.,42.,26.};
    Float_t wINOX[9]={0.0003,0.02,0.01,0.00045,0.0003,0.12,0.17,0.025,0.654};
    Float_t dINOX = 8.03;

    //SDD HV microcable

    Float_t aHVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHVm[5]={6.,1.,7.,8.,13.};
    Float_t wHVm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
    Float_t dHVm = 1.6087;

    //SDD LV+signal cable

    Float_t aLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zLVm[5]={6.,1.,7.,8.,13.};
    Float_t wLVm[5]={0.21722436468,0.0082859922,0.023028867,0.06574077612, 0.68572};
    Float_t dLVm = 2.1035;

    //SDD hybrid microcab

    Float_t aHLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHLVm[5]={6.,1.,7.,8.,13.};
    Float_t wHLVm[5]={0.24281879711,0.00926228815,0.02574224025,0.07348667449, 0.64869};
    Float_t dHLVm = 2.0502;

    //SDD anode microcab

    Float_t aALVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zALVm[5]={6.,1.,7.,8.,13.};
    Float_t wALVm[5]={0.392653705471,0.0128595919215,0.041626868025,0.118832707289, 0.431909};
    Float_t dALVm = 2.0502;

    //X7R capacitors

    Float_t aX7R[7]={137.327,47.867,15.9994,58.6928,63.5460,118.710,207.2};
    Float_t zX7R[7]={56.,22.,8.,28.,29.,50.,82.};
    Float_t wX7R[7]={0.251639432,0.084755042,0.085975822,0.038244751,0.009471271,0.321736471,0.2081768};
    Float_t dX7R = 7.14567;

    // AIR

    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;

    // Water

    Float_t aWater[2]={1.00794,15.9994};
    Float_t zWater[2]={1.,8.};
    Float_t wWater[2]={0.111894,0.888106};
    Float_t dWater   = 1.0;

    // CERAMICS
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3
    Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
    Float_t zcer[5]  = {       13.,     8.,    14.,     25.,    24. };
    Float_t wcer[5]  = {.4443408,.5213375,.0130872,.0178135,.003421};
    Float_t denscer  = 3.6;

    //G10FR4

    Float_t zG10FR4[14] = {14.00,	20.00,	13.00,	12.00,	5.00,	22.00,	11.00,	19.00,	26.00,	9.00,	8.00,	6.00,	7.00,	1.00};
    Float_t aG10FR4[14] = {28.0855000,40.0780000,26.9815380,24.3050000,10.8110000,47.8670000,22.9897700,39.0983000,55.8450000,18.9984000,15.9994000,12.0107000,14.0067000,1.0079400};
    Float_t wG10FR4[14] = {0.15144894,0.08147477,0.04128158,0.00904554,0.01397570,0.00287685,0.00445114,0.00498089,0.00209828,0.00420000,0.36043788,0.27529426,0.01415852,0.03427566};
    Float_t densG10FR4= 1.8;
    
     //--- EPOXY  --- C18 H19 O3
      Float_t aEpoxy[3] = {15.9994, 1.00794, 12.0107} ; 
      Float_t zEpoxy[3] = {     8.,      1.,      6.} ; 
      Float_t wEpoxy[3] = {     3.,     19.,     18.} ; 
      Float_t dEpoxy = 1.8 ;

      // rohacell: C9 H13 N1 O2
    Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
    Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
    Float_t wrohac[4] = { 14.,   10.,    2.,     6.};
    Float_t drohac    = 0.052;

    // If he/she means stainless steel (inox) + Aluminium and Zeff=15.3383 then
//
// %Al=81.6164 %inox=100-%Al

    Float_t aInAl[5] = {27., 55.847,51.9961,58.6934,28.0855 };
    Float_t zInAl[5] = {13., 26.,24.,28.,14. };
    Float_t wInAl[5] = {.816164, .131443,.0330906,.0183836,.000919182};
    Float_t dInAl    = 3.075;

    // Kapton

    Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
    Float_t zKapton[4]={1.,6.,7.,8.};
    Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
    Float_t dKapton   = 1.42;
    
    // Kapton + Cu (for Pixel Bus)

    Float_t aKaptonCu[5]={1.00794, 12.0107, 14.010, 15.9994, 63.5460};
    Float_t zKaptonCu[5]={1., 6., 7., 8., 29.};
    Float_t wKaptonCuBus[5];
    
    // Kapton + Cu (for Pixel MCM)

    Float_t wKaptonCuMCM[5];
    
    // Kapton + Cu (mix of two above)

    Float_t wKaptonCuMix[5];

    //SDD ruby sph.
    Float_t aAlOxide[2]  = { 26.981539,15.9994};
    Float_t zAlOxide[2]  = {       13.,     8.};
    Float_t wAlOxide[2]  = {0.4707, 0.5293};
    Float_t dAlOxide     = 3.97;

    // Silica for optical fibers: Si O2
    Float_t aoptfib[2] = { 28.0855, 15.9994};
    Float_t zoptfib[2] = { 14.,      8.    };
    Float_t woptfib[2] = {  1.,      2.    };
    Float_t doptfib    = 2.55;

    // Tetrafluorethylene-Perfluorpropylene (FEP) - 08 Mar 10
    Float_t aFEP[2] = { 12.0107, 18.9984};
    Float_t zFEP[2] = {  6.    ,  9.    };
    Float_t wFEP[2] = {  1.    ,  2.    };
    Float_t dFEP    = 2.15;

    // PVC (C2H3Cl)n - 08 Jul 10
    Float_t aPVC[3] = { 12.0107, 1.00794, 35.4527};
    Float_t zPVC[3] = {  6.    , 1.     , 35.   };
    Float_t wPVC[3] = {  2.    , 3.     ,  1.   };
    Float_t dPVC    = 1.3;

    //SSD NiSn capacitor ends
    Float_t aNiSn[2]  = { 56.6934,118.710};
    Float_t zNiSn[2]  = {     28.,     50.};
    Float_t wNiSn[2]  = {0.33, 0.67};
    Float_t dNiSn     = wNiSn[0]*8.908 + wNiSn[1]*7.310;

    AliMaterial(1,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(1,"SI$",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(2,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(2,"SPD SI CHIP$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(3,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(3,"SPD SI BUS$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(4,"C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(4,"C (M55J)$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(5,"AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(5,"AIR$",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(6,"GEN AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(6,"GEN AIR$",6,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(7,"SDD SI CHIP$",aSICHIP,zSICHIP,dSICHIP,6,wSICHIP);
    AliMedium(7,"SDD SI CHIP$",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(9,"SDD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(9,"SDD C (M55J)$",9,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(10,"SDD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(10,"SDD AIR$",10,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMaterial(11,"AL$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(11,"AL$",11,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(12, "Water$",aWater,zWater,dWater,2,wWater);
    AliMedium(12,"WATER$",12,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(13,"Freon$",afre,zfre,densfre,-2,wfre);
    AliMedium(13,"Freon$",13,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(14,"COPPER$",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
    AliMedium(14,"COPPER$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    AliMixture(15,"CERAMICS$",acer,zcer,denscer,5,wcer);
    AliMedium(15,"CERAMICS$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(20,"SSD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(20,"SSD C (M55J)$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(21,"SSD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(21,"SSD AIR$",21,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(25,"G10FR4$",aG10FR4,zG10FR4,densG10FR4,14,wG10FR4);
    AliMedium(25,"G10FR4$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(26,"GEN C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(26,"GEN C (M55J)$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(27,"GEN Air$",aAir,zAir,dAir,4,wAir);
    AliMedium(27,"GEN Air$",27,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(47,"PVC$",aPVC,zPVC,dPVC,-3,wPVC);
    AliMedium(47,"PVC$",47,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    Double_t cuFrac = 0.56;
    Double_t kFrac  = 1.0 - cuFrac;
    Double_t cuDens = 8.96;
    Float_t dKaptonCuBus   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuBus[j] = wKapton[j]*kFrac;
    wKaptonCuBus[4] = cuFrac;
    AliMixture(48, "SPD-BUS CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuBus, 5, wKaptonCuBus);
    AliMedium(48,"SPD-BUS CU KAPTON$",48,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
    cuFrac = 0.5;
    kFrac  = 1.0 - cuFrac;
    Float_t dKaptonCuMCM   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuMCM[j] = wKapton[j]*kFrac;
    wKaptonCuMCM[4] = cuFrac;
    AliMixture(49, "SPD-MCM CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuMCM, 5, wKaptonCuMCM);
    AliMedium(49,"SPD-MCM CU KAPTON$",49,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
    cuFrac = (0.56 + 0.5) / 2.0;
    kFrac  = 1.0 - cuFrac;
    Float_t dKaptonCuMix   = cuFrac * cuDens + kFrac * dKapton;
    for (Int_t j=0; j<4; j++)
      wKaptonCuMix[j] = wKapton[j]*kFrac;
    wKaptonCuMix[4] = cuFrac;
    AliMixture(50, "SPD-MIX CU KAPTON", aKaptonCu, zKaptonCu, dKaptonCuMix, 5, wKaptonCuMix);
    AliMedium(50,"SPD-MIX CU KAPTON$",50,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(51,"SPD SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(51,"SPD SI$",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(52,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(52,"SPD SI CHIP$",52,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(53,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(53,"SPD SI BUS$",53,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(54,"SPD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(54,"SPD C (M55J)$",54,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(55,"SPD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(55,"SPD AIR$",55,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(56, "SPD KAPTON(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(56,"SPD KAPTON(POLYCH2)$",56,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Gaseous Freon has same chemical composition but air density at 1.7 atm
    AliMixture(59,"GASEOUS FREON$",afre,zfre,1.7*dAir,-2,wfre);
    AliMedium(59,"GASEOUS FREON$",59,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(61,"EPOXY$",aEpoxy,zEpoxy,dEpoxy,-3,wEpoxy);
    AliMedium(61,"EPOXY$",61,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(62,"SILICON$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(62,"SILICON$",62,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(63, "KAPTONH(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(63,"KAPTONH(POLYCH2)$",63,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(64,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(64,"ALUMINUM$",64,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(65,"INOX$",aINOX,zINOX,dINOX,9,wINOX);
    AliMedium(65,"INOX$",65,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(66,"NiSn$",aNiSn,zNiSn,dNiSn,2,wNiSn);
    AliMedium(66,"NiSn$",66,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(67,"Sn$", 118.710, 50., 7.310, 1.206, 999.);
    AliMedium(67,"Sn$",67,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(68,"ROHACELL$",arohac,zrohac,drohac,-4,wrohac);
    AliMedium(68,"ROHACELL$",68,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(69,"SDD C AL (M55J)$",aALCM55J,zALCM55J,dALCM55J,5,wALCM55J);
    AliMedium(69,"SDD C AL (M55J)$",69,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
  
    AliMixture(70, "SDDKAPTON (POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(70,"SDDKAPTON (POLYCH2)$",70,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMaterial(71,"ITS SANDW A$",0.12011E+02,0.60000E+01,0.2115E+00,0.17479E+03,0.99900E+03);
    AliMedium(71,"ITS SANDW A$",71,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(72,"ITS SANDW B$",0.12011E+02,0.60000E+01,0.27000E+00,0.18956E+03,0.99900E+03);
    AliMedium(72,"ITS SANDW B$",72,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(73,"ITS SANDW C$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    AliMedium(73,"ITS SANDW C$",73,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(74,"HEAT COND GLUE$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
    AliMedium(74,"HEAT COND GLUE$",74,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(75,"ELASTO SIL$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(75,"ELASTO SIL$",75,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // SPD bus (data from Petra Riedler)
    Float_t aSPDbus[5] = {1.00794,12.0107,14.01,15.9994,26.982 };
    Float_t zSPDbus[5] = {1.,6.,7.,8.,13.};
    Float_t wSPDbus[5] = {0.023523,0.318053,0.009776,0.078057,0.570591};
    Float_t dSPDbus    = 2.128505;

    //   AliMaterial(76,"SPDBUS(AL+KPT+EPOX)$",0.19509E+02,0.96502E+01,0.19060E+01,0.15413E+02,0.99900E+03);
    AliMixture(76,"SPDBUS(AL+KPT+EPOX)$",aSPDbus,zSPDbus,dSPDbus,5,wSPDbus);
    AliMedium(76,"SPDBUS(AL+KPT+EPOX)$",76,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               
    AliMixture(77,"SDD X7R capacitors$",aX7R,zX7R,dX7R,7,wX7R);
    AliMedium(77,"SDD X7R capacitors$",77,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(78,"SDD ruby sph. Al2O3$",aAlOxide,zAlOxide,dAlOxide,2,wAlOxide);
    AliMedium(78,"SDD ruby sph. Al2O3$",78,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(79,"SDD SI insensitive$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(79,"SDD SI insensitive$",79,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(80,"SDD HV microcable$",aHVm,zHVm,dHVm,5,wHVm);
    AliMedium(80,"SDD HV microcable$",80,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(81,"SDD LV+signal cable$",aLVm,zLVm,dLVm,5,wLVm);
    AliMedium(81,"SDD LV+signal cable$",81,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(82,"SDD hybrid microcab$",aHLVm, zHLVm,dHLVm,5,wHLVm);
    AliMedium(82,"SDD hybrid microcab$",82,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(83,"SDD anode microcab$",aALVm,zALVm,dALVm,5,wALVm);
    AliMedium(83,"SDD anode microcab$",83,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    Float_t aDSring[4]={12.0107,      1.00794,     14.0067,      15.9994};
    Float_t zDSring[4]={ 6.,          1.,           7.,           8.};
    Float_t wDSring[4]={ 0.854323888, 0.026408778,  0.023050265,  0.096217069};
    Float_t dDSring = 0.2875;
    AliMixture(84,"SDD/SSD rings$",aDSring,zDSring,dDSring,4,wDSring);
    AliMedium(84,"SDD/SSD rings$",84,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(85,"inox/alum$",aInAl,zInAl,dInAl,5,wInAl);
    AliMedium(85,"inox/alum$",85,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // special media to take into account services in the SDD and SSD 
    // cones for the FMD
    //Begin_Html
    /*
      <A HREF="http://www.Physics.ohio-state.edu/~nilsen/ITS/ITS_MatBudget_4B.xls">
      </pre>
      <br clear=left>
      <font size=+2 color=blue>
      <p> The Exel spread sheet from which these density number come from.
      </font></A>
    */
    //End_Html

    //  AliMaterial(86,"AIRFMDSDD$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
    Float_t aA[13],zZ[13],wW[13],den;
    // From Pierluigi Barberis calculations of 2SPD+1SDD October 2 2002.
    zZ[0] = 1.0; aA[0] = 1.00794; // Hydrogen
    zZ[1] = 6.0; aA[1] = 12.011; // Carbon
    zZ[2] = 7.0; aA[2] = 14.00674; // Nitrogen
    zZ[3] = 8.0; aA[3] = 15.9994; // Oxigen
    zZ[4] = 14.0; aA[4] = 28.0855; // Silicon
    zZ[5] = 24.0; aA[5] = 51.9961; //Cromium
    zZ[6] = 25.0; aA[6] = 54.938049; // Manganese
    zZ[7] = 26.0; aA[7] = 55.845; // Iron
    zZ[8] = 28.0; aA[8] = 58.6934; // Nickle
    zZ[9] = 29.0; aA[9] = 63.546; // Copper
    zZ[10] = 13.0; aA[10] = 26.981539; // Alulminum
    zZ[11] = 47.0; aA[11] = 107.8682; // Silver
    zZ[12] = 27.0; aA[12] = 58.9332; // Cobolt
    wW[0] = 0.019965;
    wW[1] = 0.340961;
    wW[2] = 0.041225;
    wW[3] = 0.200352;
    wW[4] = 0.000386;
    wW[5] = 0.001467;
    wW[6] = 0.000155;
    wW[7] = 0.005113;
    wW[8] = 0.000993;
    wW[9] = 0.381262;
    wW[10] = 0.008121;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.5253276; // g/cm^3  Cell O370
    }else{
	den = 2.58423412; // g/cm^3 Cell L370
    } // end if fByThick
    //den = 6161.7/(3671.58978);//g/cm^3 Volume does not exclude holes
    AliMixture(86,"AIRFMDSDD$",aA,zZ,den,+11,wW);
    AliMedium(86,"AIRFMDSDD$",86,0,ifield,fieldm,tmaxfdAir,stemaxAir,
	      deemaxAir,epsilAir,stminAir);

    //AliMaterial(87,"AIRFMDSSD$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
    // From Pierluigi Barberis calculations of SSD October 2 2002.
    wW[0] = 0.019777;
    wW[1] = 0.325901;
    wW[2] = 0.031848;
    wW[3] = 0.147668;
    wW[4] = 0.030609;
    wW[5] = 0.013993;
    wW[6] = 0.001479;
    wW[7] = 0.048792;
    wW[8] = 0.009477;
    wW[9] = 0.350697;
    wW[10] = 0.014546;
    wW[11] = 0.005213;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.2464275; // g/cm^3   Cell O403
    }else{
	den = 1.28134409; // g/cm^3  Cell L403
    } // end if fByThick
    //den = 7666.3/(9753.553259); // volume does not exclude holes
    AliMixture(87,"AIRFMDSSD$",aA,zZ,den,+12,wW); 
    AliMedium(87,"AIRFMDSSD$",87,0,ifield,fieldm,tmaxfdAir,stemaxAir,
	      deemaxAir,epsilAir,stminAir);

    //AliMaterial(88,"ITS SANDW CFMDSDD$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of 1SDD+Carbon fiber October 2 2002
    wW[0] = 0.016302;
    wW[1] = 0.461870;
    wW[2] = 0.033662;
    wW[3] = 0.163595;
    wW[4] = 0.000315;
    wW[5] = 0.001197;
    wW[6] = 0.000127;
    wW[7] = 0.004175;
    wW[8] = 0.000811;
    wW[9] = 0.311315;
    wW[10] = 0.006631;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.9353276; // g/cm^3  Cell N370
    }else{
	den = 3.2788626; // g/cm^3 Cell F370
    } // end if fByThick
    //den = 7667.1/(3671.58978); // Volume does not excludeholes
    AliMixture(88,"ITS SANDW CFMDSDD$",aA,zZ,den,+11,wW); 
    AliMedium(88,"ITS SANDW CFMDSDD$",88,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //AliMaterial(89,"ITS SANDW CFMDSSD$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of SSD+Carbon fiber October 2 2002.
    wW[0] = 0.014065;
    wW[1] = 0.520598;
    wW[2] = 0.022650;
    wW[3] = 0.105018;
    wW[4] = 0.021768;
    wW[5] = 0.009952;
    wW[6] = 0.001051;
    wW[7] = 0.034700;
    wW[8] = 0.006740;
    wW[9] = 0.249406;
    wW[10] = 0.010345;
    wW[11] = 0.0003707;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.6564275; // g/cm^3  Cell N304
    }else{
	den = 1.7028296; // g/cm^3  Cell F304
    } // end if fByThick
    //den = 1166.5/(3671.58978); // Volume does not exclude holes
    AliMixture(89,"ITS SANDW CFMDSSD$",aA,zZ,den,+12,wW); 
    AliMedium(89,"ITS SANDW CFMDSSD$",89,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //AliMaterial(97,"SPD SERVICES$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of 1SPD October 2 2002.
    wW[0] = 0.005970;
    wW[1] = 0.304704;
    wW[2] = 0.042510;
    wW[3] = 0.121715;
    wW[4] = 0.001118;
    wW[5] = 0.030948;
    wW[6] = 0.003270;
    wW[7] = 0.107910;
    wW[8] = 0.020960;
    wW[9] = 0.360895;
    wW[10] = 0.000000;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 80.31136576; // g/cm^3 Cell H329
    }else{
	den = 87.13062; // g/cm^3  Cell G329
    } // end if fByThick
    //den = 1251.3/(0.05*2.0*TMath::Pi()*(7.75*7.75 - 3.7*3.7)); // g/cm^3
    AliMixture(97,"SPD SERVICES$",aA,zZ,den,+10,wW); 
    AliMedium(97,"SPD SERVICES$",97,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);


    // Special media

    AliMaterial(90,"SPD shield$", 12.011, 6., 1.93 , 22.36, 999);
    AliMedium(90,"SPD shield$",90,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    // SPD End Ladder (data from Petra Riedler)
    Float_t aSPDel[5] = {1.00794,12.0107,14.01,15.9994,63.54 };
    Float_t zSPDel[5] = {1.,6.,7.,8.,29.};
    Float_t wSPDel[5] = {0.004092,0.107274,0.011438,0.032476,0.844719};
    Float_t dSPDel    = 3.903403;

    //   AliMaterial(91, "SPD End ladder$", 47.0447, 21.7963, 3.6374, 4.4711, 999); 
    AliMixture(91,"SPD End ladder$",aSPDel,zSPDel,dSPDel,5,wSPDel);
    AliMedium(91,"SPD End ladder$",91,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    AliMaterial(92, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);    
    AliMedium(92,"SPD cone$",92,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /*  Material with fractional Z not actually used
    AliMaterial(93, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999);
    AliMedium(93,"SDD End ladder$",93,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(94, "SDD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(94,"SDD cone$",94,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /* Material with fractional Z not actually used
    AliMaterial(95, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
    AliMedium(95,"SSD End ladder$",95,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(96, "SSD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(96,"SSD cone$",96,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    AliMixture(98,"SDD OPTICFIB$",aoptfib,zoptfib,doptfib,-2,woptfib);
    AliMedium(98,"SDD OPTICFIB$",98,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(95,"SSD FEP$",aFEP,zFEP,dFEP,-2,wFEP);
    AliMedium(95,"SSD FEP$",95,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // Mean material for low-voltage cables on SPD trays Side A
    // (Copper + PolyEthylene (C2-H4)) (D.Elia for cable number and
    // cross-section area, M.Sitta for elemental computation) - 26 Feb 10
    wW[0] = 0.323024;//H
    wW[2] = 0.515464;//Cu
    wW[1] = 0.161512;//C
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 5.078866;
    AliMixture(60,"SPD_LOWCABLES$",aA,zZ,den,+3,wW);
    AliMedium(60,"SPD_LOWCABLES$",60,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    // Mean material for high-voltage cables on SPD trays Side A & C
    // (Copper + HD PolyEthylene (C2-H2)) (D.Elia for cable number and
    // cross-section area, M.Sitta for elemental computation) - 10 Jun 10
    wW[0] = 0.083766;//H
    wW[2] = 0.417136;//Cu
    wW[1] = 0.499098;//C
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.514930;
    AliMixture(58,"SPD_HICABLES$",aA,zZ,den,+3,wW);
    AliMedium(58,"SPD_HICABLES$",58,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    // PolyUrethane [C25-H42-N2-O6] - 07 Mar 10
    zZ[2] =  7.0; aA[2] =  14.0067; // Nitrogen - From Root TGeoElementTable

    wW[0] = 0.090724;//H
    wW[2] = 0.060035;//N
    wW[1] = 0.643513;//C
    wW[3] = 0.205728;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.158910;
    AliMixture(67,"POLYURETHANE$",aA,zZ,den,+4,wW);
    AliMedium(67,"POLYURETHANE$",67,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //  POM (Polyoxymethylene = (CH2O)n ) - 02 May 10
    zZ[2] =  8.0; aA[2] =  15.9994; // Oxigen

    wW[0] = 0.067137;//H
    wW[1] = 0.400016;//C
    wW[2] = 0.532847;//O
    wW[3] = 0.000000;//O
    wW[4] = 0.000000;//S
    wW[5] = 0.000000;//F
    wW[6] = 0.000000;//Sn
    wW[7] = 0.000000;//Pb
    wW[8] = 0.000000;//Cr
    wW[9] = 0.000000;//Si
    wW[10] = 0.000000;//Ni
    wW[11] = 0.000000;//Ca

    den = 1.4200;
    AliMixture(57,"POLYOXYMETHYLENE$",aA,zZ,den,+3,wW);
    AliMedium(57,"POLYOXYMETHYLENE$",57,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);


    // Anticorodal: Aliminum alloy for tray ring support on Side A
    den = 2.710301;
    AliMaterial(93,"ANTICORODAL$",0.26982E+02,0.13000E+02,den,0.89000E+01,0.99900E+03);
    AliMedium(93,"ANTICORODAL$",93,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    
}

//_____________________________________________________________________________
void AliMaterial(Int_t imat, const char* name, Float_t a, 
		 Float_t z, Float_t dens, Float_t radl,
		 Float_t absl) {

  TString uniquename = "ITS";
  uniquename.Append("_");
  uniquename.Append(name);

  gGeoManager->Material(uniquename.Data(), a, z, dens, imat, radl, absl);

}

//_____________________________________________________________________________
void AliMedium(Int_t numed, const char *name, Int_t nmat,
	       Int_t isvol, Int_t ifield, Float_t fieldm,
	       Float_t tmaxfd, Float_t stemax, Float_t deemax,
	       Float_t epsil, Float_t stmin){

  TString uniquename = "ITS";
  uniquename.Append("_");
  uniquename.Append(name);

  gGeoManager->Medium(uniquename.Data(),numed,nmat, isvol, ifield, fieldm, 
		      tmaxfd, stemax,deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliMixture(Int_t kmat, const char* name, Float_t* a, Float_t* z,
                    Double_t dens, Int_t nlmat, Float_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture.
  //
  // If nlmat < 0 then WMAT contains the number of atoms
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //

  TString uniquename = "ITS";
  uniquename.Append("_");
  uniquename.Append(name);


   Double_t* da = CreateDoubleArray(a, TMath::Abs(nlmat));
   Double_t* dz = CreateDoubleArray(z, TMath::Abs(nlmat));
   Double_t* dwmat = CreateDoubleArray(wmat, TMath::Abs(nlmat));

   Mixture(kmat, uniquename.Data(), da, dz, dens, nlmat, dwmat);
   for (Int_t i=0; i<nlmat; i++) {
      a[i] = da[i]; z[i] = dz[i]; wmat[i] = dwmat[i];
   }

   delete [] da;
   delete [] dz;
   delete [] dwmat;
}

//_____________________________________________________________________________
void Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z,
	     Double_t dens, Int_t nlmat, Double_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture.
  //
  // If nlmat < 0 then WMAT contains the number of atoms
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //

   if (nlmat < 0) {
      nlmat = - nlmat;
      Double_t amol = 0;
      Int_t i;
      for (i=0;i<nlmat;i++) {
         amol += a[i]*wmat[i];
      }
      for (i=0;i<nlmat;i++) {
         wmat[i] *= a[i]/amol;
      }
   }
   gGeoManager->Mixture(name, a, z, dens, nlmat, wmat, kmat);
}

//_____________________________________________________________________________
Double_t* CreateDoubleArray(Float_t* array, Int_t size)
{
// Converts Float_t* array to Double_t*,
// !! The new array has to be deleted by user.
// ---

   Double_t* doubleArray;
   if (size>0) {
      doubleArray = new Double_t[size];
      for (Int_t i=0; i<size; i++) doubleArray[i] = array[i];
   } else {
      //doubleArray = 0;
      doubleArray = new Double_t[1];
   }
   return doubleArray;
}

//______________________________________________________________________
Bool_t Make2DCrossSections(TPolyLine &a0,TPolyLine &a1,
			   TPolyLine &b0,TPolyLine &b1,TPolyMarker &p)
{
    //
    // Fill the objects with the points representing
    // a0 the outer carbon fiber SPD sector shape Cross Section A
    // a1 the inner carbon fiber SPD sector shape Cross Section A
    // b0 the outer carbon fiber SPD sector shape Cross Section B
    // b1 the inner carbon fiber SPD sector shape Cross Section B
    //
    // Inputs:
    //   TPolyLine &a0   The outer carbon fiber SPD sector shape
    //   TPolyLine &a1   The Inner carbon fiber SPD sector shape
    //   TPolyLine &b0   The outer carbon fiber SPD sector shape
    //   TPolyLine &b1   The Inner carbon fiber SPD sector shape
    //   TPolyMarker &p  The points where the ladders are to be placed
    // Outputs:
    //   TPolyLine &a0   The shape filled with the points
    //   TPolyLine &a1   The shape filled with the points
    //   TPolyLine &b0   The shape filled with the points
    //   TPolyLine &b1   The shape filled with the points
    //   TPolyMarker &p  The filled array of points
    // Return:
    //     An error flag.
    //
    Int_t n0,n1,i;
    Double_t x,y;
    TGeoVolume *a0V,*a1V,*b0V,*b1V;
    TGeoXtru *a0S,*a1S,*b0S,*b1S;
    TGeoManager *mgr = gGeoManager;

    a0V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorA0");
    a0S = dynamic_cast<TGeoXtru*>(a0V->GetShape());
    n0 = a0S->GetNvert();
    a0.SetPolyLine(n0+1);
    //for(i=0;i<fSPDsectorPoints0.GetSize();i++)
    //  printf("%d %d %d\n",i,fSPDsectorPoints0[i],fSPDsectorPoints1[i]);
    for(i=0;i<n0;i++){
        x = a0S->GetX(i);
          y = a0S->GetY(i);
          //printf("%d %g %g\n",i,x,y);
        a0.SetPoint(i,x,y);
          if(i==0) a0.SetPoint(n0,x,y);
    } // end for i
    a1V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorAirA1");
    a1S = dynamic_cast<TGeoXtru*>(a1V->GetShape());
    n1 = a1S->GetNvert();
    a1.SetPolyLine(n1+1);
    for(i=0;i<n1;i++){
        x = a1S->GetX(i);
          y = a1S->GetY(i);
        a1.SetPoint(i,x,y);
          if(i==0) a1.SetPoint(n1,x,y);
    } // end for i
    // Cross Section B
    b0V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndB0");
    b0S = dynamic_cast<TGeoXtru*>(b0V->GetShape());
    n0 = b0S->GetNvert();
    b0.SetPolyLine(n0+1);
    for(i=0;i<n0;i++){
        x = b0S->GetX(i);
          y = b0S->GetY(i);
        b0.SetPoint(i,x,y);
          if(i==0) b0.SetPoint(n0,x,y);
    } // end for i
    b1V = mgr->GetVolume("ITSSPDCarbonFiberSupportSectorEndAirB1");
    b1S = dynamic_cast<TGeoXtru*>(b1V->GetShape());
    n1 = b1S->GetNvert();
    b1.SetPolyLine(n1+1);
    for(i=0;i<n1;i++){
        x = b1S->GetX(i);
          y = b1S->GetY(i);
        b1.SetPoint(i,x,y);
          if(i==0) b1.SetPoint(n1,x,y);
    } // end for i
    //
    Double_t x0,y0,x1,y1;
    p.SetPolyMarker(2*gspd->GetSPDsectorX0Size());
    for(i=0;i<gspd->GetSPDsectorX0Size();i++){
          gspd->GetSectorMountingPoints(i,x0,y0,x1,y1);
          p.SetPoint(2*i,x0,y0);
          p.SetPoint(2*i+1,x1,y1);
    } // end for i
    return kTRUE;
}