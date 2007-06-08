
TVector3 v[28];
Int_t nCh;


TGeoHMatrix GetResSurvAlign(Int_t survNch);

void SurveyToAlignHmpid(){


 AliSurveyObj *so = new AliSurveyObj();


 Int_t size = so->GetEntries();
  printf("-> %d\n", size);
  
   so->FillFromLocalFile("Survey_781282_HMPID.txt");
  size = so->GetEntries();
  printf("--> %d\n", size);
  
  
   TObjArray *points = so->GetData();
//   TVector3 v[28];
   
  for (Int_t i = 0; i < points->GetEntries(); ++i)
   {
   AliSurveyPoint *p=(AliSurveyPoint *) points->At(i);
   v[i].SetXYZ(p->GetX()*100.,p->GetY()*100.,p->GetZ()*100.);
   }


//  // To produce the alignment object for the given volume you would
//  // then do something like this:
//  // Calculate the global delta transformation as ng * g3-1
//  TGeoHMatrix gdelta = g3->Inverse(); //now equal to the inverse of g3
//  gdelta.MultiplyLeft(&ng);
//  Int_t index = 0;
//  // if the volume is in the look-up table use something like this instead:
//  // AliGeomManager::LayerToVolUID(AliGeomManager::kTOF,i); 
//  AliAlignObjMatrix* mobj = new AliAlignObjMatrix("symname",index,gdelta,kTRUE);


TGeoHMatrix mtx = GetResSurvAlign(5);

TGeoManager::Import("/home/mserio/tstesdtrk/geometry.root");
gGeoManager->cd(Form("ALIC_1/Hmp_%1i",nCh));
TGeoHMatrix g0 = *gGeoManager->GetCurrentMatrix();
cout<<"\n\n*********Ideal Matrix (chamber "<<nCh<<")*********"<<endl;
g0.Print();
TGeoHMatrix gdelta = g0.Inverse();
gdelta.MultiplyLeft(&mtx);

//gdelta.Print();

AliAlignObjMatrix* mobj = new
AliAlignObjMatrix(AliGeomManager::SymName(AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,nCh)),
                                          AliGeomManager::LayerToVolUID(AliGeomManager::kHMPID,nCh),gdelta,kTRUE);
/*
cout<<"\n************* obtained   AliAlignObjMatrix************\n";
mobj->Print();
cout<<""<<endl;

TGeoHMatrix pa=gdelta*g0;

pa.Print();
*/
}


TGeoHMatrix GetResSurvAlign(Int_t survNch)
{
cout<<"    ************Survey numbering********Offline Numbering**********"<<endl;
cout<<"\nChamber No            0                      4                   "<<endl;
cout<<"Chamber No            1                      3                   "<<endl;
cout<<"Chamber No            2                      5                   "<<endl;
cout<<"Chamber No            3                      1                   "<<endl;
cout<<"Chamber No            4                      6                   "<<endl;
cout<<"Chamber No            5                      2                   "<<endl;
cout<<"Chamber No            6                      0                   "<<endl;


  // From the new fiducial marks coordinates derive back the
  // new global position of the surveyed volume
  //*** The 4 fiducial marks are assumed on a rectangle
  //*** parallel to a surface of the Hmp (main volume)
  //***  at a certain offset from the origin (zdepth) and with
  //*** x and y sides parallel to the box's x and y axes.

if(survNch==0) nCh=4;
if(survNch==1) nCh=3;
if(survNch==2) nCh=5;
if(survNch==3) nCh=1;
if(survNch==4) nCh=6;
if(survNch==5) nCh=2;
if(survNch==6) nCh=0;

  Double_t ab[3], bc[3], n[3];
  Double_t plane[4], s;
  Double_t ngA[3]={v[0+4*survNch].X(),v[0+4*survNch].Y(),v[0+4*survNch].Z()};
  Double_t ngB[3]={v[1+4*survNch].X(),v[1+4*survNch].Y(),v[1+4*survNch].Z()};
  Double_t ngC[3]={v[2+4*survNch].X(),v[2+4*survNch].Y(),v[2+4*survNch].Z()};
  Double_t ngD[3]={v[3+4*survNch].X(),v[3+4*survNch].Y(),v[3+4*survNch].Z()};
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
   cout<<endl<<"The center of the box from Survey data: "<<md[0]<<"  "<<md[1]<<"  "<<md[2]<<endl;
  const Double_t zdepth=-0.9-4.85; //the survey data are down the radiator (behind the honeycomb structure). They
                                   //lay on 4 cylinders whose height is 9 mm.

  // The center of the box
  for(i=0;i<1;i++){
    orig[i] = md[i] - (-plane[i])*(zdepth+plane[3]);
  }
  orig[1] = md[1] - (-plane[1])*(zdepth+plane[3]);
  orig[2] = md[2] - (-plane[2])*(zdepth+plane[3]);
 
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
  Double_t rot[9] = {-ab[0],bc[0],-plane[0],-ab[1],bc[1],-plane[1],-ab[2],bc[2],-plane[2]};
  TGeoHMatrix ng;
  ng.SetTranslation(md);
  ng.SetRotation(rot);

  cout<<"\n********* global matrix inferred from surveyed fiducial marks for chamber"<<survNch<<"***********\n";
  ng.Print();


return ng;

}




