void AliTRDreadRaw2Clusters(const char *fname = "raw.root", const char *fnameGeom = "./geometry.root")
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  TFile *fGeometryFile = TFile::Open(fnameGeom);
  TGeoManager *fGeoManager = 0;
  if (fGeometryFile)
    {
      fGeoManager = (TGeoManager *)fGeometryFile->Get("Geometry");
    }

  if (fGeoManager == 0)
    {
      cout << "Geo Manager init failed." << endl;
    }

  AliTRDdigitsManager manR;
  manR.CreateArrays();

  AliRawReaderRoot reader(fname, 0);
  reader.Select("TRD");

  Int_t ievent = 0;
  while (reader.NextEvent())
    {
      TTree *treeR = new TTree(sdir, "TRD clusters");
      AliTRDclusterizer clusterizer("TRDclusterizer", "TRDclusterizer");
      clusterizer.OpenOutput(treeR);
      Int_t ir = clusterizer.Raw2ClustersChamber(&reader);
      
      cout << "Clusterizer returned " << ir << endl;

      // do something witht he clusters...

      ievent++;
    }

}
