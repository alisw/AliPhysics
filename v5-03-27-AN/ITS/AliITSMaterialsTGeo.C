void AliITSMaterialsTGeo(TString gfile="geometry.root"){
  // Macro to print out the ITS material definitions as found
  // in the TGeo geometry file.

  // retrives geometry 
  if(!gGeoManager) gGeoManager = new TGeoManager();
  TGeoManager::Import(gfile.Data());
  if (!gGeoManager) {
    cout<<"geometry not found\n";
    return;
  } // end if

  TList *medlist=gGeoManager->GetListOfMedia();
  TGeoMedium *med;
  TGeoMaterial *mat;
  Int_t imed,nmed,i;
  printf("imed  Id       Med_Name             Mat_Name       ");
  for(i=0;i<20;i++) printf("   par[%2d]   ",i);
  printf("\n");
  imed=0;
  do{
    med = (TGeoMedium*)(medlist->At(imed));
    if(!med) continue;
    /*if((((med->GetName())[0]=='I')&& // Only ITS.
        ((med->GetName())[1]=='T')&&
        ((med->GetName())[2]=='S')&&
	((med->GetName())[3]=='_')))*/{
    mat = med->GetMaterial();
    if(mat)
      printf("%4d %4d %30s %30s",imed,med->GetId(),med->GetName(),mat->GetName());
    else
      printf("%4d %4d %30s %30s",imed,med->GetId(),med->GetName(),"No Material");
    for(i=0;i<20;i++) printf(" %12g",med->GetParam(i));
    printf("\n");
    imed++;
    }
  }while(med!=medlist->Last());
}
