void AliITSPlotGeom(const char *filename="galice.root"){
/////////////////////////////////////////////////////////////////////////
//   This macro displays part of the ITS sensitive volume as defined
// in the root file.
//   
//     Root > .L AliITSPlotGeom.C         //this loads the macro in memory
//     Root > AliITSPlotGeom();           //by default process first event   
//     Root > AliITSPlotGeom("galice2.root"); //process third event from 
//                                              galice2.root file.
//Begin_Html
/*
<img src="picts/AliITSplotgeom.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////
    if(gAlice){
	delete gAlice;
	gAlice=0;
    }else{
	// Dynamically link some shared libs
	if(gClassTable->GetID("AliRun") < 0) {
	    gROOT->LoadMacro("loadlibs.C");
	    loadlibs();
	} // end if
    } // end if gAlice
    // Connect the Root Galice file containing Geometry, Kine and Hits
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if(!file) file = new TFile(filename);

    // Get AliRun object from file or create it if not on file
    if(!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if(gAlice) printf("AliRun object found on file\n");
	if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    } // end if !gAlice
      

// Get pointers to the ITS Alice detectors and Hits containers
    AliITS    *ITS    = (AliITS*)    gAlice->GetDetector("ITS");
    if(ITS==0){
	cout << "ITS not found. Exiting." << endl;
	return;
    } // end if ITS==0
    AliITSgeom *itsgeom = ITS->GetITSgeom();
    if(itsgeom==0){
	cout << "ITSgeom not defined. Exiting." << endl;
	return;
    } // end if itsgeom==0
    Int_t is,ie,in;
    is = itsgeom->GetStartSPD();
    ie = itsgeom->GetLastSPD()+1;
    in = ie-is;

    // Transformations.
//    Float_t *x = new Float_t[in];
//    Float_t *y = new Float_t[in];
//    Float_t *z = new Float_t[in];
    Float_t x,y,z;
    TRotMatrix **r = new TRotMatrix*[in];
    Int_t i,j;
    Double_t m[9];
    char name[10],title[20],name2[10],title2[20];
//    cout << in << endl;
    for(i=0;i<in;i++){
	itsgeom->GetRotMatrix(i+is,m);
	sprintf(name,"ROT%d",i+is);
	sprintf(title,"ROT%d",i+is);
	r[i] = new TRotMatrix(name,title,m);
//	cout << name << title << endl;
//	itsgeom->GetTrans(i+is,x[i],y[i],z[i]);
//	cout << i << " " << x[i] << " " << y[i] << " " << z[i] <<endl;
    } // end for i

    // Sensitive volumes.
    if(itsgeom->IsShapeDefined(0)){
	TBRIK *spds = (TShape*) (itsgeom->GetShape(0));
    }else{
	TBRIK *spds = new TBRIK("SPD","SPD","void",0.64,0.015,3.48);
    } // end if
//    cout << spds << endl;
    if(itsgeom->IsShapeDefined(1)){
	TBRIK *sdds = (TShape*) (itsgeom->GetShape(240));
    }else{
	TBRIK *sdds = new TBRIK("SDD","SDD","void",3.50850,0.01499,3.76320);
    } // end if
//    cout << sdds << endl;
    if(itsgeom->IsShapeDefined(2)){
	TBRIK *ssds = (TShape*) (itsgeom->GetShape(1000));
    }else{
	TBRIK *ssds = new TBRIK("SSD","SSD","void",3.65,0.015,2.00);
    } // end if
//    cout << ssds << endl;

    // Set up display.

    TCanvas *c1 = new TCanvas("c1","ITS geometry",10,10,700,700);
//  TPad *p1 = new TPad("p1","p1",0.01,0.01,0.99,0.99,46,3,1);
//  p1->Draw();
//  p1->cd();
//  TView *view = new TView(1);
//  view->SetRange(-5,-5,-5,5,5,5);
      TShape  *mother = new TBRIK("Mother","Mother Volume","void",10,10,10);
//      TShape  *mother = new TBRIK("Mother","Mother Volume","void",30,30,30);
//      TShape  *mother = new TBRIK("Mother","Mother Volume","void",50,50,50);
      TNode *node = new TNode("node","Mother Node",mother);
      node->cd();

    // Set up nodes
    TNode **itsn = new TNode*[in];
    TPolyLine3D **pl = new TPolyLine3D*[in];
/*    Double_t p[5*3],axis[5*3]={1.,0.,0.,
                               0.,0.,0.,
			       0.,1.,0.,
			       0.,0.,0.,
			       0.,0.,1.}; */
    Double_t p[19*3],axis[19*3]={-0.25,0.,1.25,  // 1
				 +0.00,0.,1.25,  // 2
				 -0.25,0.,1.00,  // 3
				 +0.00,0.,1.00,  // 4
				 +0.00,0.,0.00,  // 5
				 +1.00,0.,0.00,  // 6
				 +1.25,0.,0.25,  // 7
				 +1.125,0.,0.125,  // 8
				 +1.00,0.,0.25,  // 9
				 +1.125,0.,0.125,  // 10
				 +1.25,0.,0.00,  // 11
				 +1.125,0.,0.125,  // 12
				 +1.00,0.,0.00,  // 13
				 +0.00,0.,0.00,  // 14
				 +0.00,1.,0.00,  // 15
				 +0.00,1.,0.125,  // 16
				 +0.25,1.,0.25,  // 17
				 +0.00,1.,0.125,  // 18
				 -0.25,1.,0.25,  // 19
    };
    for(i=0;i<in;i++){
	itsgeom->GetTrans(i+is,x,y,z);
/*	switch (itsgeom->GetGeomMatrix(i+is)->GetDetectorIndex())
        case 0:
	    sprintf(name,"SPD%d",i+is);
	    sprintf(title,"SPD%d",i+is);
	    sprintf(name2,"BSPD%d",i+is);
	    sprintf(title2,"BSPD%d",i+is);
	    itsn[i] = new TNode(name,title,new TBRIK(name2,title2,"void",
                      spds->GetDx(),spds->GetDy(),spds->GetDz()),x,y,z,r[i]);
	    break;
        case 1:
	    sprintf(name,"SDD%d",i+is);
	    sprintf(title,"SDD%d",i+is);
	    sprintf(name2,"BSDD%d",i+is);
	    sprintf(title2,"BSDD%d",i+is);
	    itsn[i] = new TNode(name,title,new TBRIK(name2,title2,"void",
                      sdds->GetDx(),sdds->GetDy(),sdds->GetDz()),x,y,z,r[i]);
	    break;
        case 2:
	    sprintf(name,"SSD%d",i+is);
	    sprintf(title,"SSD%d",i+is);
	    sprintf(name2,"BSSD%d",i+is);
	    sprintf(title2,"BSSD%d",i+is);
	    itsn[i] = new TNode(name,title,new TBRIK(name2,title2,"void",
                      ssds->GetDx(),ssds->GetDy(),ssds->GetDz()),x,y,z,r[i]);
	    break; */
	for(j=0;j<19;j++) itsgeom->LtoG(i+is,(Double_t*)&axis[3*j],
                                            (Double_t*)&p[3*j]);
	pl[i] = new TPolyLine3D(19,p);
    } // end for i

    // display it
    node->cd();
    node->Draw();
//    for(i=0;i<in;i++) itsn[i]->Draw();
//    node->Draw();
    for(i=0;i<in;i++) pl[i]->Draw();
    c1->Update();

    // clean up
//    delete[] x;
//    delete[] y;
//    delete[] z;
//    for(i=0;i<itsgeom->GetIndexMax();i++) delete[] r[i];
//    delete[] r;
}

