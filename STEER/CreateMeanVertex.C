void CreateMeanVertex(Double_t xmed=0., Double_t ymed=0., Double_t sigx=0.005,Double_t sigy=0.005, Double_t sigz=5.3){
  Double_t resolx=35./10000.;
  Double_t resoly=35./10000.;
  Double_t sigma[3],position[3];
  position[0]=xmed;
  position[1]=ymed;
  position[2]=0.;
  sigma[0]=TMath::Sqrt(sigx*sigx+resolx*resolx);
  sigma[1]=TMath::Sqrt(sigy*sigy+resoly*resoly);
  sigma[2]=sigz;
  AliESDVertex *vave=new AliESDVertex(position,sigma,"vtxmean");
  vave->PrintStatus();

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliCDBId id("GRP/Calib/MeanVertex", 0, AliCDBRunRange::Infinity());
  AliCDBMetaData md;
  
  md.SetResponsible("prino@to.infn.it");
  md.SetComment("Default mean vertex position");
  md.SetAliRootVersion("Default mean vertex position");
  
  man->Put(vave, id, &md);
  
  man->Destroy();

//  TFile *outf=new TFile("AliESDMeanVertex.root","recreate");
// outf->cd();
//  vave->Write();
//  outf->Close();
}
