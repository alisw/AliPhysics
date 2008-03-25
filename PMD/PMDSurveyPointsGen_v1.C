void PMDSurveyPointsGen_v1(){
  
  // Macro to generate survey points in the format of AliSurveyObj 
  // by giving some finite rotation and translation.
  
  if(!gGeoManager) AliGeomManager::LoadGeometry("geometry.root");
  
  TClonesArray *array = new TClonesArray("AliAlignObjMatrix",10);
  TClonesArray &mobj = *array;
  
  Double_t l_vect[3]={0.,0.,0.}; // a local vector (the origin)
  
  Double_t g_vect_1[3];    // vector corresp. to it in global RS for sector-1
  Double_t g_vect_2[3];    // vector corresp. to it in global RS for sector-2
  Double_t g_vect_3[3];    // vector corresp. to it in global RS for sector-3
  Double_t g_vect_4[3];    // vector corresp. to it in global RS for sector-4
  
  //**************** get global matrix *******************
  
  TGeoHMatrix *g3_1 = AliGeomManager::GetMatrix("PMD/Sector1");
  TGeoNode* n3_1 = gGeoManager->GetCurrentNode();
  TGeoHMatrix* l3_1 = n3_1->GetMatrix();// to get local matrix  
  g3_1->LocalToMaster(l_vect,g_vect_1); // point coordinates in the global RS
  
  cout<<endl<<"Origin of sector-1 in the global RS: "<<
    g_vect_1[0]<<" "<<g_vect_1[1]<<" "<<g_vect_1[2]<<" "<<endl;
  

  TGeoHMatrix *g3_2 = AliGeomManager::GetMatrix("PMD/Sector2");
  TGeoNode* n3_2 = gGeoManager->GetCurrentNode();
  TGeoHMatrix* l3_2 = n3_2->GetMatrix(); 
  g3_2->LocalToMaster(l_vect,g_vect_2);
  
  cout<<endl<<"Origin of sector-2 in the global RS: "<<
    g_vect_2[0]<<" "<<g_vect_2[1]<<" "<<g_vect_2[2]<<" "<<endl;
  
  
  TGeoHMatrix *g3_3 = AliGeomManager::GetMatrix("PMD/Sector3");
  TGeoNode* n3_3 = gGeoManager->GetCurrentNode();
  TGeoHMatrix* l3_3 = n3_3->GetMatrix(); 
  g3_3->LocalToMaster(l_vect,g_vect_3);
  
  cout<<endl<<"Origin of the sector-3 in the global RS: "<<
    g_vect_3[0]<<" "<<g_vect_3[1]<<" "<<g_vect_3[2]<<" "<<endl;
  
  
  TGeoHMatrix *g3_4 = AliGeomManager::GetMatrix("PMD/Sector4");
  TGeoNode* n3_4 = gGeoManager->GetCurrentNode();
  TGeoHMatrix* l3_4 = n3_4->GetMatrix();
  g3_4->LocalToMaster(l_vect,g_vect_4);
  
  cout<<endl<<"Origin of the sector-4 in the global RS: "<<
    g_vect_4[0]<<" "<<g_vect_4[1]<<" "<<g_vect_4[2]<<" "<<endl;
  
  
  
  // Hereafter the four ideal fiducial marks on the PMD 
  // are expressed in local coordinates. 
  //All the coordinates are expressed in cm here.
  const Double_t zsize   = 2.3;
  const Double_t zoffset = 0.0;
  const Double_t zdepth  = zsize + zoffset + g_vect_1[2];
  
  // Coordinates of four Ideal fudicial marks on the sector-1 of the PMD in 
  //global frame of reference.
  
  Double_t A1[3] = {12.025, 37.85, zdepth};// Fiducial mark # 19
  Double_t B1[3] = {74.275, 37.85, zdepth};// Fiducial mark # 20
  Double_t C1[3] = {74.275, 87.05, zdepth};// Fiducial mark # 18
  Double_t D1[3] = {12.025, 87.05, zdepth};// Fiducial mark # 17
     
  // Coordinates of four Ideal fudicial marks on the sector-2 of the PMD in 
  //global frame of reference.
    
  Double_t A2[3] = {-74.275, -87.05, zdepth};// Fiducial mark # 15
  Double_t B2[3] = {-12.025, -87.05, zdepth};// Fiducial mark # 16
  Double_t C2[3] = {-12.025, -37.85, zdepth};// Fiducial mark # 14
  Double_t D2[3] = {-74.275, -37.85, zdepth};// Fiducial mark # 13
  
  // Coordinates of four Ideal fudicial marks on the sector-3 of the PMD in 
  //global frame of reference.
  
  Double_t A3[3] = {-74.275, 11.15, zdepth};// Fiducial mark # 11
  Double_t B3[3] = {7.725, 11.15, zdepth};  // Fiducial mark # 12
  Double_t C3[3] = {7.725, 87.05, zdepth};  // Fiducial mark # 06
  Double_t D3[3] = {-74.25, 87.05, zdepth}; // Fiducial mark # 05
  
  // Coordinates of four Ideal fudicial marks on the sector-4 of the PMD in 
  //global frame of reference.
  
  Double_t A4[3] = {-7.725, -87.05, zdepth};// Fiducial mark # 27
  Double_t B4[3] = {74.275, -87.05, zdepth};// Fiducial mark # 28
  Double_t C4[3] = {74.275, -11.15, zdepth};// Fiducial mark # 22
  Double_t D4[3] = {-7.725, -11.15, zdepth};// Fiducial mark # 21
  
  // Let's  misalign the sector-1 for testing purpose only.
  for (Int_t i = 0; i < 3; i++) 
    {
      // To convert the coordinates of fiducial marks in local RS as those are
      // given  above in global RS.
      A1[i] = A1[i] - g_vect_1[i];
      B1[i] = B1[i] - g_vect_1[i];
      C1[i] = C1[i] - g_vect_1[i];
      D1[i] = D1[i] - g_vect_1[i];
      
    }
  
  
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
   
  Double_t gA1[3], gB1[3], gC1[3], gD1[3];
  
  g3_1->LocalToMaster(A1,gA1);
  g3_1->LocalToMaster(B1,gB1);
  g3_1->LocalToMaster(C1,gC1);
  g3_1->LocalToMaster(D1,gD1);

  /*
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS:1\n"<<
    "A1 "<<gA1[0]<<" "<<gA1[1]<<" "<<gA1[2]<<" "<<endl<<
    "B1 "<<gB1[0]<<" "<<gB1[1]<<" "<<gB1[2]<<" "<<endl<<
    "C1 "<<gC1[0]<<" "<<gC1[1]<<" "<<gC1[2]<<" "<<endl<<
    "D1 "<<gD1[0]<<" "<<gD1[1]<<" "<<gD1[2]<<" "<<endl;
  */
  
  // We apply a delta transformation to the surveyed vol. to represent
  // its real position, given below by ng3 and nl3, which differs from its
  // ideal position saved above in g3 and l3
  
  TGeoPhysicalNode* pn3_1 = gGeoManager->MakePhysicalNode("ALIC_1/EPM1_1");
  
  Double_t dphi_1   = 0.; //  tilt around Z
  Double_t dtheta_1 = 0.; //  tilt around X
  Double_t dpsi_1   = 0.; //  tilt around new Z
  
  Double_t dx_1 = 2.; // shift along X
  Double_t dy_1 = 3.; // shift along Y
  Double_t dz_1 = 4.; // shift along Z
  
  TGeoRotation* rrot_1 = new TGeoRotation("rot",dphi_1,dtheta_1,dpsi_1);
  TGeoCombiTrans localdelta_1 = *(new TGeoCombiTrans(dx_1,dy_1,dz_1, rrot_1));
    
  // new local matrix, representing real position
  TGeoHMatrix nlocal_1 = *l3_1 * localdelta_1;
  TGeoHMatrix* nl3_1 = new TGeoHMatrix(nlocal_1);
  pn3_1->Align(nl3_1);
  
  // Let's get the global matrix for later comparison
  TGeoHMatrix* ng3_1 = pn3_1->GetMatrix(); //"real" global matrix,what survey sees 
  printf("\n\n************  real global matrix for Sector-1 **************\n");
  ng3_1->Print();
  
  Double_t ngA1[3], ngB1[3], ngC1[3], ngD1[3];
  ng3_1->LocalToMaster(A1,ngA1);
  ng3_1->LocalToMaster(B1,ngB1);
  ng3_1->LocalToMaster(C1,ngC1);
  ng3_1->LocalToMaster(D1,ngD1);
  
  cout<<endl<<"Fiducial marks coordinates in the global RS given by survey:1\n"<<
    "A1 "<<ngA1[0]<<" "<<ngA1[1]<<" "<<ngA1[2]<<" "<<endl<<
    "B1 "<<ngB1[0]<<" "<<ngB1[1]<<" "<<ngB1[2]<<" "<<endl<<
    "C1 "<<ngC1[0]<<" "<<ngC1[1]<<" "<<ngC1[2]<<" "<<endl<<
    "D1 "<<ngD1[0]<<" "<<ngD1[1]<<" "<<ngD1[2]<<" "<<endl;
  
  
  // Let's  misalign the sector-2.
  for (Int_t i = 0; i < 3; i++) 
    {
      // To convert the coordinates of fiducial marks in local RS as those are
      // given  above in global RS.
      A2[i] = A2[i] - g_vect_2[i];
      B2[i] = B2[i] - g_vect_2[i];
      C2[i] = C2[i] - g_vect_2[i];
      D2[i] = D2[i] - g_vect_2[i];
      
    }
  
  Double_t gA2[3], gB2[3], gC2[3], gD2[3];
  
  g3_2->LocalToMaster(A2,gA2);
  g3_2->LocalToMaster(B2,gB2);
  g3_2->LocalToMaster(C2,gC2);
  g3_2->LocalToMaster(D2,gD2);
  
  /*
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS:2\n"<<
    "A2 "<<gA2[0]<<" "<<gA2[1]<<" "<<gA2[2]<<" "<<endl<<
    "B2 "<<gB2[0]<<" "<<gB2[1]<<" "<<gB2[2]<<" "<<endl<<
    "C2 "<<gC2[0]<<" "<<gC2[1]<<" "<<gC2[2]<<" "<<endl<<
    "D2 "<<gD2[0]<<" "<<gD2[1]<<" "<<gD2[2]<<" "<<endl;
  */
  
  // We apply a delta transformation to the surveyed vol. to represent
  // its real position, given below by ng3 and nl3, which differs from its
  // ideal position saved above in g3 and l3
  
  TGeoPhysicalNode* pn3_2 = gGeoManager->MakePhysicalNode("ALIC_1/EPM2_1");
  
  Double_t dphi_2   = 0.; // tilt around Z
  Double_t dtheta_2 = 0.; // tilt around X
  Double_t dpsi_2   = 0.; // tilt around new Z
  
  Double_t dx_2   = 2.; // shift along X
  Double_t dy_2   = 3.; // shift along Y
  Double_t dz_2   = 4.; // shift along Z
  
  TGeoRotation* rrot_2 = new TGeoRotation("rot",dphi_2,dtheta_2,dpsi_2);
  TGeoCombiTrans localdelta_2 = *(new TGeoCombiTrans(dx_2,dy_2,dz_2, rrot_2));
  
  // new local matrix, representing real position
  TGeoHMatrix nlocal_2 = *l3_2 * localdelta_2;
  TGeoHMatrix* nl3_2 = new TGeoHMatrix(nlocal_2);
  pn3_2->Align(nl3_2);
  
  // Let's get the global matrix for later comparison
  TGeoHMatrix* ng3_2 = pn3_2->GetMatrix(); //real global matrix 
  printf("\n\n************  real global matrix for Sector-2 **************\n");
  ng3_2->Print();
  
  Double_t ngA2[3], ngB2[3], ngC2[3], ngD2[3];
  ng3_2->LocalToMaster(A2,ngA2);
  ng3_2->LocalToMaster(B2,ngB2);
  ng3_2->LocalToMaster(C2,ngC2);
  ng3_2->LocalToMaster(D2,ngD2);
  
  cout<<endl<<"Fiducial marks coordinates in the global RS given by survey:2\n"<<
    "A2 "<<ngA2[0]<<" "<<ngA2[1]<<" "<<ngA2[2]<<" "<<endl<<
    "B2 "<<ngB2[0]<<" "<<ngB2[1]<<" "<<ngB2[2]<<" "<<endl<<
    "C2 "<<ngC2[0]<<" "<<ngC2[1]<<" "<<ngC2[2]<<" "<<endl<<
    "D2 "<<ngD2[0]<<" "<<ngD2[1]<<" "<<ngD2[2]<<" "<<endl;
  
// Let's  misalign the sector-3.
  for (Int_t i = 0; i < 3; i++) 
    {
      // To convert the coordinates of fiducial marks in local RS as those are
      // given  above in global RS.
      A3[i] = A3[i] - g_vect_3[i];
      B3[i] = B3[i] - g_vect_3[i];
      C3[i] = C3[i] - g_vect_3[i];
      D3[i] = D3[i] - g_vect_3[i];
    }
  
  Double_t gA3[3], gB3[3], gC3[3], gD3[3];
  
  g3_3->LocalToMaster(A3,gA3);
  g3_3->LocalToMaster(B3,gB3);
  g3_3->LocalToMaster(C3,gC3);
  g3_3->LocalToMaster(D3,gD3);
  
  /*
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS:3\n"<<
    "A3 "<<gA3[0]<<" "<<gA3[1]<<" "<<gA3[2]<<" "<<endl<<
    "B3 "<<gB3[0]<<" "<<gB3[1]<<" "<<gB3[2]<<" "<<endl<<
    "C3 "<<gC3[0]<<" "<<gC3[1]<<" "<<gC3[2]<<" "<<endl<<
    "D3 "<<gD3[0]<<" "<<gD3[1]<<" "<<gD3[2]<<" "<<endl;
  */

  // We apply a delta transformation to the surveyed vol. to represent
  // its real position, given below by ng3 and nl3, which differs from its
  // ideal position saved above in g3 and l3
  
  TGeoPhysicalNode* pn3_3 = gGeoManager->MakePhysicalNode("ALIC_1/EPM3_1");
  
  Double_t dphi_3   = 0.; //  tilt around Z
  Double_t dtheta_3 = 0.; //  tilt around X
  Double_t dpsi_3   = 0.; //  tilt around new Z
  
  Double_t dx_3   = 2.; // shift along X
  Double_t dy_3   = 3.; // shift along Y
  Double_t dz_3   = 4 ; // shift along Z
    
  TGeoRotation* rrot_3 = new TGeoRotation("rot",dphi_3,dtheta_3,dpsi_3);
  TGeoCombiTrans localdelta_3 = *(new TGeoCombiTrans(dx_3,dy_3,dz_3, rrot_3));
    
  // new local matrix, representing real position
  TGeoHMatrix nlocal_3 = *l3_3 * localdelta_3;
  TGeoHMatrix* nl3_3 = new TGeoHMatrix(nlocal_3);
  pn3_3->Align(nl3_3);
  
  // Let's get the global matrix for later comparison
  TGeoHMatrix* ng3_3 = pn3_3->GetMatrix(); //"real" global matrix,what survey sees 
  printf("\n\n************  real global matrix for Sector-3 **************\n");
  ng3_3->Print();
  
  Double_t ngA3[3], ngB3[3], ngC3[3], ngD3[3];
  ng3_3->LocalToMaster(A3,ngA3);
  ng3_3->LocalToMaster(B3,ngB3);
  ng3_3->LocalToMaster(C3,ngC3);
  ng3_3->LocalToMaster(D3,ngD3);
  
  cout<<endl<<"Fiducial marks coordinates in the global RS given by survey:3\n"<<
    "A3 "<<ngA3[0]<<" "<<ngA3[1]<<" "<<ngA3[2]<<" "<<endl<<
    "B3 "<<ngB3[0]<<" "<<ngB3[1]<<" "<<ngB3[2]<<" "<<endl<<
    "C3 "<<ngC3[0]<<" "<<ngC3[1]<<" "<<ngC3[2]<<" "<<endl<<
    "D3 "<<ngD3[0]<<" "<<ngD3[1]<<" "<<ngD3[2]<<" "<<endl;
 
  // Let's  misalign the sector-4.
  for (Int_t i = 0; i < 3; i++) 
    {
      // To convert the coordinates of fiducial marks in local RS as those are
      // given  above in global RS.
      A4[i] = A4[i] - g_vect_4[i];
      B4[i] = B4[i] - g_vect_4[i];
      C4[i] = C4[i] - g_vect_4[i];
      D4[i] = D4[i] - g_vect_4[i];
      
    }
  
  Double_t gA4[3], gB4[3], gC4[3], gD4[3];
  
  g3_4->LocalToMaster(A4,gA4);
  g3_4->LocalToMaster(B4,gB4);
  g3_4->LocalToMaster(C4,gC4);
  g3_4->LocalToMaster(D4,gD4);
  
  /*
  cout<<endl<<"Ideal fiducial marks coordinates in the global RS:4\n"<<
    "A4 "<<gA4[0]<<" "<<gA4[1]<<" "<<gA4[2]<<" "<<endl<<
    "B4 "<<gB4[0]<<" "<<gB4[1]<<" "<<gB4[2]<<" "<<endl<<
    "C4 "<<gC4[0]<<" "<<gC4[1]<<" "<<gC4[2]<<" "<<endl<<
    "D4 "<<gD4[0]<<" "<<gD4[1]<<" "<<gD4[2]<<" "<<endl;
  */

  // We apply a delta transformation to the surveyed vol. to represent
  // its real position, given below by ng3 and nl3, which differs from its
  // ideal position saved above in g3 and l3
  
  TGeoPhysicalNode* pn3_4 = gGeoManager->MakePhysicalNode("ALIC_1/EPM4_1");
  
  Double_t dphi_4   = 0.; //  tilt around Z
  Double_t dtheta_4 = 0.; //  tilt around X
  Double_t dpsi_4   = 0.; //  tilt around new Z
  
  Double_t dx_4   = 2.; // shift along X
  Double_t dy_4   = 3.; // shift along Y
  Double_t dz_4   = 4.; // shift along Z
  
    
  TGeoRotation* rrot_4 = new TGeoRotation("rot",dphi_4,dtheta_4,dpsi_4);
  TGeoCombiTrans localdelta_4 = *(new TGeoCombiTrans(dx_4,dy_4,dz_4, rrot_4));
    
  // new local matrix, representing real position
  TGeoHMatrix nlocal_4 = *l3_4 * localdelta_4;
  TGeoHMatrix* nl3_4 = new TGeoHMatrix(nlocal_4);
  pn3_4->Align(nl3_4);
  
  // Let's get the global matrix for later comparison
  TGeoHMatrix* ng3_4 = pn3_4->GetMatrix(); //"real" global matrix,what survey sees 
  printf("\n\n************  real global matrix for Sector-4 **************\n");
  ng3_4->Print();
  
  Double_t ngA4[3], ngB4[3], ngC4[3], ngD4[3];
  ng3_4->LocalToMaster(A4,ngA4);
  ng3_4->LocalToMaster(B4,ngB4);
  ng3_4->LocalToMaster(C4,ngC4);
  ng3_4->LocalToMaster(D4,ngD4);
  
  cout<<endl<<"Fiducial marks coordinates in the global RS given by survey:4\n"<<
    "A4 "<<ngA4[0]<<" "<<ngA4[1]<<" "<<ngA4[2]<<" "<<endl<<
    "B4 "<<ngB4[0]<<" "<<ngB4[1]<<" "<<ngB4[2]<<" "<<endl<<
    "C4 "<<ngC4[0]<<" "<<ngC4[1]<<" "<<ngC4[2]<<" "<<endl<<
    "D4 "<<ngD4[0]<<" "<<ngD4[1]<<" "<<ngD4[2]<<" "<<endl;
 

  Double_t ngP1[4][3], ngP2[4][3],ngP3[4][3], ngP4[4][3];
   
  for(Int_t i=0;i<4;i++)
    { 
      for (Int_t j=0; j<3; j++)
	{
	  ngP1[0][j] = ngA1[j];
	  ngP1[1][j] = ngB1[j];
	  ngP1[2][j] = ngC1[j];
	  ngP1[3][j] = ngD1[j];

	  ngP2[0][j] = ngA2[j];
	  ngP2[1][j] = ngB2[j];
	  ngP2[2][j] = ngC2[j];
	  ngP2[3][j] = ngD2[j];
	  
	  ngP3[0][j] = ngA3[j];
	  ngP3[1][j] = ngB3[j];
	  ngP3[2][j] = ngC3[j];
	  ngP3[3][j] = ngD3[j];
	  
	  ngP4[0][j] = ngA4[j];
	  ngP4[1][j] = ngB4[j];
	  ngP4[2][j] = ngC4[j];
	  ngP4[3][j] = ngD4[j];
  
	}
    }
  
  Float_t xx, yy, zz;
  Char_t seqno[6], PType[1], TUsed[1];
  Int_t precision;
  
  ofstream ftw;
  ftw.open("PMDGenSurveyPoints_v1.txt");

  ifstream ftr;
  ftr.open("Survey_Points_PMD01.txt",ios::in);
    
  Int_t nline = 0;
  const int Nrow=70;
  TString CommentLine[Nrow];
  TString tmp, tmp1;
  char line[80];
  Int_t surseq1[4] = {19, 20, 18, 17};
  Int_t surseq2[4] = {15, 16, 14, 13};
  Int_t surseq3[4] = {11, 12, 06, 05};
  Int_t surseq4[4] = {27, 28, 22, 21};
  
  while (nline<Nrow) 
    {
      ftr.getline(line,80, '\n');
      tmp = line;
      tmp1 = tmp(0,60);
      CommentLine[nline] = Strip(tmp1);
      
      ftw<<CommentLine[nline]<<endl;
      
      if (nline >= 40 && nline < 68)
	{ 
	  ftr>>seqno>>xx>>yy>>zz>>PType>>TUsed>>precision;
	  
	  // cout<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<"\t"<<TUsed<<
	  //"\t"<<precision<<endl;
	  
	  Int_t i= (nline-40)+1;
	  for (Int_t j=0;j<4;j++)
	    {	       
	      if (surseq1[j] == i)
		{
		  xx = ngP1[j][0]*10.0;
		  yy = ngP1[j][1]*10.0;
		  zz = ngP1[j][2]*10.0;
		  cout<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<
		    "\t"<<TUsed<< "\t"<<precision<<endl;   
		}
	     else if (surseq2[j] == i)
		{
		  xx = ngP2[j][0]*10.0;
		  yy = ngP2[j][1]*10.0;
		  zz = ngP2[j][2]*10.0;
		  cout<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<
		    "\t"<<TUsed<< "\t"<<precision<<endl;   
		}
	      else if (surseq3[j] == i)
		{
		  xx = ngP3[j][0]*10.0;
		  yy = ngP3[j][1]*10.0;
		  zz = ngP3[j][2]*10.0;
		  cout<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<
		    "\t"<<TUsed<< "\t"<<precision<<endl;   
		}
	     else if (surseq4[j] == i)
		{
		  xx = ngP4[j][0]*10.0;
		  yy = ngP4[j][1]*10.0;
		  zz = ngP4[j][2]*10.0;
		  cout<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<
		    "\t"<<TUsed<< "\t"<<precision<<endl;   
		}
	    }
	  ftw<<seqno<<"\t"<<xx<<"\t"<<yy<<"\t"<<zz<<"\t"<<PType<<"\t"<<TUsed<<
	    "\t"<<precision<<endl; 
	  
	}
      
      //cout << "CommentLine [" << nline<< "]: " << CommentLine[nline] <<endl;
      nline++;
    }
}
