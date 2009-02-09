void VZEROASurveyToAlignment(){

  // Macro to convert survey data into alignment data. 
  // The position of four fiducial marks, sticked on the 
  // entrance face of the V0A box is converted into the 
  // global position of the box. Positions given by surveyers 
  // are extracted from Survey Data Base. Thanks to Brigitte Cheynis
  // for providing this macro witch just had to be modified in order 
  // to obtain the desired results. 

  if(!gGeoManager) TGeoManager::Import("geometry.root");

//  TClonesArray *array = new TClonesArray("AliAlignObjMatrix",10);
  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &mobj = *array;
 
  Double_t l_vect[3]={0.,0.,0.}; // local vector (the origin)
  Double_t g_vect[3];            // vector corresp. to it in global RS
  Double_t m_vect[3];            // vector corresp. to it in mother RS
 
  // ************* get global matrix g3 *******************
  //  TGeoHMatrix *g3 = AliGeomManager::GetMatrix("VZERO/V0A");
  TGeoHMatrix *g3 = gGeoManager->GetCurrentMatrix();
  // this is used below as the IDEAL global matrix

  // ************* get local matrix l3 *******************
  TGeoNode* n3 = gGeoManager->GetCurrentNode();
  TGeoHMatrix *l3 = n3->GetMatrix(); 
  
 // point coordinates in the global RS
  g3->LocalToMaster(l_vect,g_vect);
  cout<<endl<<"Point coordinates in the global RS: "
      <<g_vect[0]<<" "<<g_vect[1]<<" "<<g_vect[2];

 // point coordinates in the mother volume RS
  l3->LocalToMaster(l_vect,m_vect);
  cout<<endl<<"Point coordinates in the mother's volume RS: \n"
      <<m_vect[0]<<" "<<m_vect[1]<<" "<<m_vect[2]<<" "<<endl;

 // Hereafter are the four ideal fiducial marks on the V0A box, 
 // expressed in local coordinates and in cms - hard coded. 

  const Double_t xside   = 22.627;
  const Double_t yside   = 22.627;
  const Double_t zsize   = 2.5;
  const Double_t zoffset = 0.3;
  
  const Double_t zdepth  = zsize+zoffset;
  Double_t A[3]={-xside,-yside,zdepth};
  Double_t B[3]={xside,-yside,zdepth};
  Double_t C[3]={xside,yside,zdepth};
  Double_t D[3]={-xside,yside,zdepth};

  TGeoTranslation* Atr = new TGeoTranslation("Atr",-xside,-yside,zdepth);
  TGeoTranslation* Btr = new TGeoTranslation("Btr",xside,-yside,zdepth);
  TGeoTranslation* Ctr = new TGeoTranslation("Ctr",xside,yside,zdepth);
  TGeoTranslation* Dtr = new TGeoTranslation("Dtr",-xside,yside,zdepth);

  //                    ^ local y
  //                    |
  //      D-------------|-------------C
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //  ------------------|------------------> local x
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      |             |             |
  //      A-------------|-------------B
  //
  // local z exiting the plane of the screen
   
  Double_t gA[3], gB[3], gC[3], gD[3];
  g3->LocalToMaster(A,gA);
  g3->LocalToMaster(B,gB);
  g3->LocalToMaster(C,gC);
  g3->LocalToMaster(D,gD);
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS: \n"
      <<"A "<<gA[0]<<" "<<gA[1]<<" "<<gA[2]<<" "<<endl
      <<"B "<<gB[0]<<" "<<gB[1]<<" "<<gB[2]<<" "<<endl
      <<"C "<<gC[0]<<" "<<gC[1]<<" "<<gC[2]<<" "<<endl
      <<"D "<<gD[0]<<" "<<gD[1]<<" "<<gD[2]<<" "<<endl;
  cout<<endl;  
    
// Retrieval of REAL survey data from ALICE Survey Data Depot : 

  AliSurveyObj *so = new AliSurveyObj();
 
  so->FillFromLocalFile("Survey_943928_V0.txt");
  Int_t size = so->GetEntries();

  Printf("Title: \"%s\"", so->GetReportTitle().Data());
  Printf("Date: \"%s\"", so->GetReportDate().Data());
  Printf("Detector: \"%s\"", so->GetDetector().Data());
  Printf("URL: \"%s\"", so->GetURL().Data());
  Printf("Number: \"%d\"", so->GetReportNumber());
  Printf("Version: \"%d\"", so->GetReportVersion());
  Printf("Observations: \"%s\"", so->GetObservations().Data());
  Printf("Coordinate System: \"%s\"", so->GetCoordSys().Data());
  Printf("Measurement Units: \"%s\"", so->GetUnits().Data());
  Printf("Nr Columns: \"%d\" \n", so->GetNrColumns());
  
  TObjArray *colNames = so->GetColumnNames();
  
  TObjArray   *points = so->GetData();
  const char  namePoint[4]  = "6001";
  Double_t    coordinates[4][3];
//  Printf("  ******* %c ******* \n\n ", namePoint[0]); 
  Printf("Relevant points to be used for alignment procedure (in mm):"); 
  for (Int_t i = 0; i < points->GetEntries(); ++i) {
    if(((AliSurveyPoint *) points->At(i))->GetPointName()[0] == namePoint[0]) {
           Printf("Point %d --> \"%s\" %f %f %f ", i, 
           ((AliSurveyPoint *) points->At(i))->GetPointName().Data(),
	   ((AliSurveyPoint *) points->At(i))->GetX(),
	   ((AliSurveyPoint *) points->At(i))->GetY(),
	   ((AliSurveyPoint *) points->At(i))->GetZ() ); 
	   if(i > 17){
	   coordinates[i-18][0]  = (AliSurveyPoint *) points->At(i))->GetX();
	   coordinates[i-18][1]  = (AliSurveyPoint *) points->At(i))->GetY();
	   coordinates[i-18][2]  = (AliSurveyPoint *) points->At(i))->GetZ(); } 
     }
  }   
    
   Double_t ngA[3], ngB[3], ngC[3], ngD[3];

   for(Int_t i=0; i<3; i++) 
     { ngA[i]  = coordinates[0][i] / 10.0 ; 
       ngD[i]  = coordinates[1][i] / 10.0 ;
       ngB[i]  = coordinates[2][i] / 10.0 ; 
       ngC[i]  = coordinates[3][i] / 10.0 ; }
          
   cout<<endl<<"Fiducial marks coordinates in the global RS given by surveyers: \n"
       <<"A "<<ngA[0]<<" "<<ngA[1]<<" "<<ngA[2]<<" "<<endl
       <<"B "<<ngB[0]<<" "<<ngB[1]<<" "<<ngB[2]<<" "<<endl
       <<"C "<<ngC[0]<<" "<<ngC[1]<<" "<<ngC[2]<<" "<<endl
       <<"D "<<ngD[0]<<" "<<ngD[1]<<" "<<ngD[2]<<" "<<endl;
    
  // From the new fiducial marks coordinates derive back the new global position
  // of the surveyed volume
  //*** What follows is the actual survey-to-alignment procedure which assumes,
  //*** as is the case of the present example, 4 fiducial marks
  //*** at the corners of a square lying on a plane parallel to a surface
  //*** of the surveyed box at a certain offset and with
  //*** x and y sides parallel to the box's x and y axes.
  //*** If the code below is placed in a separate class or method, it needs
  //*** as input the four points and the offset from the origin (zdepth)
  //*** The algorithm can be easily modified for different placement
  //*** and/or cardinality of the fiducial marks.
  
  Double_t ab[3], bc[3], n[3];
  Double_t plane[4], s;

  // first vector on the plane of the fiducial marks
  for(i=0;i<3;i++){
    ab[i] = ngB[i] - ngA[i];
  }

  // second vector on the plane of the fiducial marks
  for(i=0;i<3;i++){
    bc[i] = ngC[i] - ngB[i];
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
//  cout<<plane[0]<<"  "<<plane[1]<<"  "<<plane[2]<<"  "<<plane[3]<<"  "<<endl;

  // The center of the square with fiducial marks as corners
  // as the middle point of one diagonal - md
  // Used below to get the center - orig - of the surveyed box
  
  Double_t orig[3], md[3];
  for(i=0;i<3;i++){
    md[i] = (ngA[i] + ngC[i]) * 0.5;
  }

  //  center of the box
  for(i=0;i<3;i++){
    orig[i] = md[i] - plane[i]*zdepth;
  }

  cout<<endl<<"Center of the box: "<<orig[0]<<"  "<<orig[1]<<"  "<<orig[2]<<endl;

  // get x,y local directions needed to write the global rotation matrix
  // for the surveyed volume by normalising vectors ab and bc
  
  Double_t sx = TMath::Sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
  if(sx>1.e-8){
     for(i=0;i<3;i++){
          ab[i] /= sx;
     }
     cout<<"x direction "<<ab[0]<<"  "<<ab[1]<<"  "<<ab[2]<<endl;
  }
  Double_t sy = TMath::Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
  if(sy>1.e-8){
    for(i=0;i<3;i++){
          bc[i] /= sy;
    }
    cout<<"y direction "<<bc[0]<<"  "<<bc[1]<<"  "<<bc[2]<<endl;
  }

  // the global matrix for the surveyed volume - ng
  Double_t rot[9] = {ab[0],bc[0],plane[0],ab[1],bc[1],plane[1],ab[2],bc[2],plane[2]};
  TGeoHMatrix ng;
  ng.SetTranslation(orig);
  ng.SetRotation(rot);

//  cout<<"\n********* global matrix inferred from surveyed fiducial marks ***********\n";
//  ng.Print();

 // To produce the alignment object for the given volume you would
 // then do something like this:
 // Calculate the global delta transformation as ng * g3^-1
 
 TGeoHMatrix gdelta = g3->Inverse(); //now equal to the inverse of g3
 gdelta.MultiplyLeft(&ng);
 Int_t index = 0;
 
 // if the volume is in the look-up table use something like this instead:
 // AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,i); 
 
 //AliAlignObjMatrix* mobj[0] = new AliAlignObjMatrix("VZERO/V0A",index,gdelta,kTRUE);
 //  new(mobj[0]) AliAlignObjMatrix("VZERO/V0C",index,gdelta,kTRUE);

  new(mobj[0]) AliAlignObjParams("VZERO/V0A",index,gdelta,kTRUE);
  
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
     TFile f("V0ASurvey.root","RECREATE");
     if(!f) cerr<<"cannot open file for output\n";
     f.cd();
     f.WriteObject(array,"V0ASurveyObjs ","kSingleKey");
     f.Close();
  }else{
    // save in CDB storage
     AliCDBManager* cdb = AliCDBManager::Instance();
     AliCDBStorage* storage = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
     AliCDBMetaData* mda = new AliCDBMetaData();
     mda->SetResponsible("Lizardo Valencia");
     mda->SetComment("Alignment objects for V0A survey");
     mda->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
     AliCDBId id("VZERO/Align/Data",0,9999999);
     storage->Put(array,id,mda);
  }
  
  cout<<"\n********* Alignment constants contained in alignment object ***********\n";
  cout<<"*************** deduced from surveyed fiducial marks : ****************\n";
  array->Print();
  
  AliAlignObjParams* itsalobj = (AliAlignObjParams*) mobj.UncheckedAt(0);
  itsalobj->ApplyToGeometry();	
  
  array->Delete();

}

