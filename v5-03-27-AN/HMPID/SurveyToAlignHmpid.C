TVector3 fFM[28]; //array of global coordinates for 28 fiducial marks
Int_t sNch, oNch; // survey and offline chamber's number


TGeoHMatrix GetResSurvAlign(Int_t survNch, Int_t& offNch);

void SurveyToAlignHmpid(const char* filename="Survey_781282_HMPID.txt"){
	// Open file with AliSurveyPoints for the 7 HMPID chambers
	// Produce the corresponding alignment objects

	AliSurveyObj *so = new AliSurveyObj();

	Int_t size = so->GetEntries();
	printf("-> %d\n", size);

	so->FillFromLocalFile(filename);
	size = so->GetEntries();
	printf("--> %d\n", size);

	TObjArray *points = so->GetData();
	
	// We retrieve and open the ideal geometry
	AliCDBManager* cdbman = AliCDBManager::Instance();
	if(!cdbman->IsDefaultStorageSet()){
		cdbman->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	}else{
		cdbman->SetSpecificStorage("GRP/Geometry/*","local://$ALICE_ROOT/OCDB");
	}
	cdbman->SetRun(0);
	AliCDBEntry* cdbe = (AliCDBEntry*) cdbman->Get("GRP/Geometry/Data"); 


	for (Int_t i = 0; i < points->GetEntries(); ++i)
	{
		AliSurveyPoint *p=(AliSurveyPoint *) points->At(i);
		fFM[i].SetXYZ(p->GetX()*100.,p->GetY()*100.,p->GetZ()*100.);
	}

	TString chbasename("/HMPID/Chamber");
	for(Int_t sNch=0; sNch<7; sNch++){
		TGeoHMatrix mtx = GetResSurvAlign(sNch,oNch); //get global matrix from survey points
	
		TString chsymname = chbasename;
		chsymname += oNch;
		printf("getting global matrix for the alignable volume %s\n",chsymname.Data());
		TGeoHMatrix *gm = AliGeomManager::GetMatrix(chsymname.Data());

		if(!gm){
			printf("unable to get global matrix for the alignable volume %s\n",chsymname.Data());
			continue;
		}
		TGeoHMatrix gdelta = gm->Inverse();
		gdelta.MultiplyLeft(&mtx);

		//gdelta.Print();

		AliAlignObjMatrix* mobj = new
			AliAlignObjMatrix(AliGeomManager::SymName(AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,oNch)),
					AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,oNch),gdelta,kTRUE);
		/*
		   cout<<"\n************* obtained   AliAlignObjMatrix************\n";
		   mobj->Print();
		   cout<<""<<endl;

		   TGeoHMatrix pa=gdelta*g0;

		   pa.Print();
		 */
	}
}


TGeoHMatrix GetResSurvAlign(Int_t survNch, Int_t& offNch)
{
	// For a given chamber identified by survey chamber number 'survNch',
	// return the global matrix inferred from the survey points of its
	// 4 fiducial marks and set the offline chamber number 'offNch'
	//
	Int_t ChSrv2Off[7] = {4,3,5,1,6,2,0};
	//cout<<"  *********  Chamber Numbers  ******"<<endl;
	//cout<<"  ****  Survey   ****  Offline *****"<<endl;
	//for(Int_t ch=0; ch<7; ch++){
	//	cout<<"           "<<ch<<"              "<<ChSrv2Off[ch]<<endl;
	//}
	
	offNch=ChSrv2Off[survNch];

	Double_t ab[3], bc[3], n[3];
	Double_t plane[4], s;
	Double_t ngA[3]={fFM[0+4*survNch].X(),fFM[0+4*survNch].Y(),fFM[0+4*survNch].Z()};
	Double_t ngB[3]={fFM[1+4*survNch].X(),fFM[1+4*survNch].Y(),fFM[1+4*survNch].Z()};
	Double_t ngC[3]={fFM[2+4*survNch].X(),fFM[2+4*survNch].Y(),fFM[2+4*survNch].Z()};
	Double_t ngD[3]={fFM[3+4*survNch].X(),fFM[3+4*survNch].Y(),fFM[3+4*survNch].Z()};
	if(survNch>4)
	{
		// first vector on the plane of the fiducial marks
		for(Int_t i=0;i<3;i++){
			ab[i] = ngB[i] - ngA[i];
		}

		// second vector on the plane of the fiducial marks
		for(Int_t i=0;i<3;i++){
			bc[i] = ngC[i] - ngB[i];
		}
	}
	else{
		// first vector on the plane of the fiducial marks
		for(Int_t i=0;i<3;i++){
			ab[i] = ngB[i] - ngA[i];
		}

		// second vector on the plane of the fiducial marks
		for(Int_t i=0;i<3;i++){
			bc[i] = ngD[i] - ngB[i];
		}

	}
	// vector normal to the plane of the fiducial marks obtained
	// as cross product of the two vectors on the plane d0^d1
	n[0] = ab[1] * bc[2] - ab[2] * bc[1];
	n[1] = ab[2] * bc[0] - ab[0] * bc[2];
	n[2] = ab[0] * bc[1] - ab[1] * bc[0];

	Double_t sizen = TMath::Sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
	if(sizen>1.e-8){
		s = Double_t(1.)/sizen ; //normalization factor
	}else{
		return 0;
	}

	// plane expressed in the hessian normal form, see:
	// http://mathworld.wolfram.com/HessianNormalForm.html
	// the first three are the coordinates of the orthonormal vector
	// the fourth coordinate is equal to the distance from the origin

	for(i=0;i<3;i++){
		plane[i] = n[i] * s;
	}
	plane[3] = -( plane[0] * ngA[0] + plane[1] * ngA[1] + plane[2] * ngA[2] );
	cout<<"normal to plane and distance from IP: "<<plane[0]<<"  "<<plane[1]<<"  "<<plane[2]<<"  "<<plane[3]<<"  "<<endl;

	// The center of the square with fiducial marks as corners
	// as the middle point of one diagonal - md
	// Used below to get the center - orig - of the surveyed box
	Double_t orig[3], md[3];

	if(survNch>4){
		for(i=0;i<3;i++){
			md[i] = (ngA[i] + ngC[i]) * 0.5;//modified!!!!!!!!!
		}

	}

	else {
		for(i=0;i<3;i++){
			md[i] = (ngA[i] + ngD[i]) * 0.5;//modified!!!!!!!!!
		}
	}
	cout<<"The center of the box from Survey data: "<<md[0]<<"  "<<md[1]<<"  "<<md[2]<<endl;
	const Double_t zdepth=-0.9-4.85; //the fiducial marks are down the radiator (behind the honeycomb structure). They
	//lay on 4 cylinders whose height is 9 mm.

	// The center of the box
	for(i=0;i<1;i++){
		orig[i] = md[i] - (-plane[i])*(zdepth+plane[3]);
	}
	orig[1] = md[1] - (-plane[1])*(zdepth+plane[3]);
	orig[2] = md[2] - (-plane[2])*(zdepth+plane[3]);

	cout<<"The origin of the box: "<<orig[0]<<"  "<<orig[1]<<"  "<<orig[2]<<endl;

	// get x,y local directions needed to write the global rotation matrix
	// for the surveyed volume by normalising vectors ab and bc
	Double_t sx = TMath::Sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
	if(sx>1.e-8){
		for(i=0;i<3;i++){
			ab[i] /= sx;
		}
		cout<<"x "<<ab[0]<<"  "<<ab[1]<<"  "<<ab[2]<<endl;
	}
	Double_t sy = TMath::Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
	if(sy>1.e-8){
		for(i=0;i<3;i++){
			bc[i] /= sy;
		}
		cout<<"y "<<bc[0]<<"  "<<bc[1]<<"  "<<bc[2]<<endl;
	}


	// the global matrix for the surveyed volume - ng
        TVector3 v1;
        v1.SetXYZ(md[0],md[1],md[2]);

        TVector3 w=v1.Unit();
        Double_t chamberCenter[3];
        chamberCenter[0]=-w.X()*(zdepth-v1.Mag());
        chamberCenter[1]=-w.Y()*(zdepth-v1.Mag());
        chamberCenter[2]=-w.Z()*(zdepth-v1.Mag());

	Double_t rot[9] = {-ab[0],bc[0],-plane[0],-ab[1],bc[1],-plane[1],-ab[2],bc[2],-plane[2]};
	TGeoHMatrix ng;
	ng.SetTranslation(md);
	ng.SetRotation(rot);

	cout<<"\n********* global matrix inferred from surveyed fiducial marks for chamber"<<survNch<<"***********\n";
	ng.Print();


	return ng;

}




