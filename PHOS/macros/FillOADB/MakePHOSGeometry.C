void MakePHOSGeometry(){
   //Fill OADB object with PHOS rotation-alignment matrixes
   //To fill this file you need any AliESD.root from official reconstruction
   //download it to the current directory
   
  AliOADBContainer * geomContainer = new AliOADBContainer("PHOSRotationMatrixes");
  TObjArray *matrixes = new TObjArray(5) ;
  matrixes->SetName("PHOSRotationMatrixes");

  //Read them from ESDs
  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    Error("MakePHOSGeometry", "Can not open ESD file, please download any from official reconstruction.");
    return ;
  }
  AliESDEvent * event = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("MakePHOSGeometry", "no ESD tree found");
    return ;
  }
  event->ReadFromTree(tree);
  tree->GetEvent(0) ;

  for(Int_t mod=0; mod<5; mod++) {
    TGeoHMatrix * m = event->GetPHOSMatrix(mod) ;
    if(!m)
      matrixes->AddAt(0x0,mod);
    else
      matrixes->AddAt(new TGeoHMatrix(*m), mod) ;
  }
  esdFile->Close() ;

  //Controll
  for(Int_t mod=0; mod<5; mod++) {
     printf("Module %d: \n", mod) ;
     if(matrixes->At(mod))
       ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
     else
       printf("Maxrix=null \n") ;
  }
  //Fill 
  geomContainer->AppendObject(matrixes,90000,AliCDBRunRange::Infinity()) ;
  geomContainer->AddDefaultObject(matrixes) ;
  geomContainer->WriteToFile("PHOSGeometry.root");
  
}
