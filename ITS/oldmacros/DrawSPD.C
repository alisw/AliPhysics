void DrawSPD(TString view="Prospective" ){
    AliITS *ITS = (AliITS*) gAlice->GetModule("ITS");
    if(!ITS) {
	cout << "You must initilize AliRoot with an ITS before plotting it"<<endl;
	return;
    } // end if
    Int_t version1 = ITS->GetMajorVersion();
    Int_t version2 = ITS->GetMinorVersion();
    gMC->Gsatt("*", "seen", -1);
    gMC->Gsatt("*", "fill", 7);
    gMC->Gsatt("alic", "seen", 0);
    gROOT->LoadMacro("ViewITSSPD.C");
    ViewITSSPD(version1,version2);
    //gInterpreter->ProcessLine("ViewSPD()");
    gMC->Gdopt("hide", "on");
    gMC->Gdopt("shad", "on");
    gMC->SetClipBox(".");
    if(view.Contains("Prospective")){
	//SetClipBox(volume name,xmin,xmax,ymin,ymax,zmin,zmax)
	gMC->SetClipBox("*",0.0,100.0,-1000.0,1000.0,0.0,1000.0);
	gMC->DefaultRange();
	//Gdraw(volume name,theta,phi,psi,axis originX,axis originY,scaleX,scaleY) angles in degrees
	gMC->Gdraw("alic", 40.0, 30.0, 0.0, 11, 10, .4, .4);
	gMC->Gdhead(1111, "Inner Tracking System");
	//gMC->Gdman(16,6,"MAN");
    } // end if
    if(view.Contains("EndView")){
	//SetClipBox(name,xmin,xmax,ymin,ymax,zmin,zmax)
	gMC->SetClipBox("*",-100.0,100.0,-1000.0,1000.0,1.0,1000.0);
	gMC->DefaultRange();
	//Gdraw(volume name,theta,phi,psi,axis originX,axis originY,scaleX,scaleY) angles in degrees
	gMC->Gdraw("alic",0.0,0.0,0.0,10.0,10.0,1.0,1.0);
	gMC->Gdhead(1111, "Inner Tracking System");
	//gMC->Gdman(16, 6, "MAN");
    } // end if
    if(view.Contains("SideView")){
	//SetClipBox(name,xmin,xmax,ymin,ymax,zmin,zmax)
	gMC->SetClipBox("*",0.0,100.0,-1000.0,1000.0,0.0,1000.0);
	gMC->DefaultRange();
	//Gdraw(volume name,theta,phi,psi,axis originX,axis originY,scaleX,scaleY) angles in degrees
	gMC->Gdraw("alic",90.0,0.0,0.0,10.0,10.0,0.35,0.35);
	gMC->Gdhead(1111, "Inner Tracking System");
	//gMC->Gdman(16, 6, "MAN");
    } // end if
}
