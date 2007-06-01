void SurveyToAlignmentExample(){
  // Macro to show an example of conversion of survey data into alignment
  // data. The position of four fiducial marks, sticked above one surface
  // of a box is converted into the global position of the box.
  // 
  gSystem->Load("libGeom");
  TGeoManager *mgr = new TGeoManager("Geom","survey to alignment toy");
  TGeoMedium *medium = 0;
  TGeoVolume *top = mgr->MakeBox("TOP",medium,250,250,250);
  mgr->SetTopVolume(top);
  // make shape components
  // ******** red outermost box ***************
  TGeoBBox *sbox0  = new TGeoBBox(200,200,50);
  TGeoVolume* box0 = new TGeoVolume("B0",sbox0);
  box0->SetVisDaughters();
  box0->SetLineColor(2); //red
  top->AddNode(box0,1);
  // ******** green middle box ***************
  TGeoBBox *sbox1  = new TGeoBBox(180,180,40);
  TGeoVolume* box1 = new TGeoVolume("B1",sbox1);
  box1->SetLineColor(3);//green
  TGeoTranslation* tr = new TGeoTranslation("tr",10,0,0);
  box0->AddNode(box1,1,tr);
  // ******** bleu inner box ***************
  TGeoBBox *sbox2  = new TGeoBBox(160,160,30);
  TGeoVolume* box2 = new TGeoVolume("B2",sbox2);
  box2->SetLineColor(4);//bleu
  box1->AddNode(box2,1,tr);
  // ******** violet innermost box ***************
  Double_t zsize = 20.;
  TGeoBBox *sbox3  = new TGeoBBox(140,140,zsize);
  TGeoVolume* box3 = new TGeoVolume("B3",sbox3);
  box3->SetLineColor(6);//violet
  box2->AddNode(box3,1,tr);

  // Four fiducial marks on the box3, expressed in local coordinates
  // We imagine they are at 2mm above the upper surface of the volume
  // at the corners of a square of 200 cm side
  const Double_t xside = 100;
  const Double_t yside = 100;
  const Double_t zoffset = 0.2;
  const Double_t zdepth = zsize+zoffset;
  Double_t A[3]={-xside,-yside,zdepth};
  Double_t B[3]={xside,-yside,zdepth};
  Double_t C[3]={xside,yside,zdepth};
  Double_t D[3]={-xside,yside,zdepth};

  TGeoBBox *fmbox  = new TGeoBBox(1,1,1);
  TGeoVolume* fm = new TGeoVolume("FM",fmbox);
  fm->SetLineColor(7);//color
  TGeoTranslation* Atr = new TGeoTranslation("Atr",-xside,-yside,zdepth);
  TGeoTranslation* Btr = new TGeoTranslation("Btr",xside,-yside,zdepth);
  TGeoTranslation* Ctr = new TGeoTranslation("Ctr",xside,yside,zdepth);
  TGeoTranslation* Dtr = new TGeoTranslation("Dtr",-xside,yside,zdepth);
  
  box3->AddNode(fm,1,Atr);
  box3->AddNode(fm,2,Btr);
  box3->AddNode(fm,3,Ctr);
  box3->AddNode(fm,4,Dtr);

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
  
  mgr->CloseGeometry();
  mgr->GetTopVolume()->Draw();
  mgr->SetVisOption(0);
  mgr->SetVisLevel(6);

  Int_t i;
  // ************* get ideal global matrix *******************
  mgr->cd("TOP_1/B0_1/B1_1/B2_1/B3_1");
  TGeoHMatrix g3 = *mgr->GetCurrentMatrix(); // !!don't declare g3
  // as a pointer to mgr->GetCurrentMatrix(), mgr->cd("...")
  // would eventually change the content pointed by g3 behind your back  

  // ************* get ideal local matrix *******************
  TGeoNode* n3 = mgr->GetCurrentNode();
  TGeoMatrix* l3 = n3->GetMatrix(); 

  Double_t gA[3], gB[3], gC[3], gD[3]; // point coordinates in the global RS
  g3.LocalToMaster(A,gA);
  g3.LocalToMaster(B,gB);
  g3.LocalToMaster(C,gC);
  g3.LocalToMaster(D,gD);
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS:\n"<<
    "A "<<gA[0]<<" "<<gA[1]<<" "<<gA[2]<<" "<<endl<<
    "B "<<gB[0]<<" "<<gB[1]<<" "<<gB[2]<<" "<<endl<<
    "C "<<gC[0]<<" "<<gC[1]<<" "<<gC[2]<<" "<<endl<<
    "D "<<gD[0]<<" "<<gD[1]<<" "<<gD[2]<<" "<<endl;
  
  // We apply a delta transformation to the surveyed vol box3 to represent
  // its real position, given below by ng3 nl3, which differs from its
  // ideal position saved above in g3 and l3
  TGeoPhysicalNode* pn3 = mgr->MakePhysicalNode("TOP_1/B0_1/B1_1/B2_1/B3_1");
  Double_t dphi = 3; // tilt by 3 degrees around z
  Double_t dz = 5; // shift by 5 cm along z
  TGeoRotation* rrot = new TGeoRotation("rot",dphi,0.,0.);
  TGeoCombiTrans localdelta = *(new TGeoCombiTrans(0.,0.,dz, rrot));
  // new local matrix, representing real position
  TGeoHMatrix nlocal = *l3 * localdelta;
  TGeoHMatrix* nl3 = new TGeoHMatrix(nlocal);
  pn3->Align(nl3);

  //Let's get the global matrix for later comparison
  TGeoHMatrix* ng3 = pn3->GetMatrix(); //"real" global matrix, what survey sees 
  printf("\n\n************  real global matrix **************\n");
  ng3->Print();
  Double_t ngA[3], ngB[3], ngC[3], ngD[3];
  ng3->LocalToMaster(A,ngA);
  ng3->LocalToMaster(B,ngB);
  ng3->LocalToMaster(C,ngC);
  ng3->LocalToMaster(D,ngD);
    
  cout<<endl<<"Fiducial marks coordinates in the global RS given by survey:\n"<<
    "A "<<ngA[0]<<" "<<ngA[1]<<" "<<ngA[2]<<" "<<endl<<
    "B "<<ngB[0]<<" "<<ngB[1]<<" "<<ngB[2]<<" "<<endl<<
    "C "<<ngC[0]<<" "<<ngC[1]<<" "<<ngC[2]<<" "<<endl<<
    "D "<<ngD[0]<<" "<<ngD[1]<<" "<<ngD[2]<<" "<<endl;


  // From the new fiducial marks coordinates derive back the
  // new global position of the surveyed volume
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
  cout<<plane[0]<<"  "<<plane[1]<<"  "<<plane[2]<<"  "<<plane[3]<<"  "<<endl;

  // The center of the square with fiducial marks as corners
  // as the middle point of one diagonal - md
  // Used below to get the center - orig - of the surveyed box
  Double_t orig[3], md[3];
  for(i=0;i<3;i++){
    md[i] = (ngA[i] + ngC[i]) * 0.5;
  }

  // The center of the box
  for(i=0;i<3;i++){
    orig[i] = md[i] - plane[i]*zdepth;
  }
  orig[1] = md[1] - plane[1]*zdepth;
  orig[2] = md[2] - plane[2]*zdepth;
  cout<<endl<<"The origin of the box: "<<orig[0]<<"  "<<orig[1]<<"  "<<orig[2]<<endl;

  // get x,y local directions needed to write the global rotation matrix
  // for the surveyed volume by normalising vectors ab and bc
  Double_t sx = TMath::Sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
  if(sx>1.e-8){
     for(i=0;i<3;i++){
          ab[i] /= sx;
     }
     cout<<endl<<"x "<<ab[0]<<"  "<<ab[1]<<"  "<<ab[2]<<endl;
  }
  Double_t sy = TMath::Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
  if(sy>1.e-8){
    for(i=0;i<3;i++){
          bc[i] /= sy;
    }
    cout<<endl<<"y "<<bc[0]<<"  "<<bc[1]<<"  "<<bc[2]<<endl;
  }


  // the global matrix for the surveyed volume - ng
  Double_t rot[9] = {ab[0],bc[0],plane[0],ab[1],bc[1],plane[1],ab[2],bc[2],plane[2]};
  TGeoHMatrix ng;
  ng.SetTranslation(orig);
  ng.SetRotation(rot);

  cout<<"\n********* global matrix inferred from surveyed fiducial marks ***********\n";
  ng.Print();

//  // To produce the alignment object for the given volume you would
//  // then do something like this:
//  // Calculate the global delta transformation as ng * g3-1
//  TGeoHMatrix gdelta = g3->Inverse(); //now equal to the inverse of g3
//  gdelta.MultiplyLeft(&ng);
//  Int_t index = 0;
//  // if the volume is in the look-up table use something like this instead:
//  // AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,i); 
//  AliAlignObjMatrix* mobj = new AliAlignObjMatrix("symname",index,gdelta,kTRUE);
 
}
