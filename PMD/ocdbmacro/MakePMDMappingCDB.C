/**************************************************************************
 * To Create PMD Mapping CDB 
 * Input is taken from old ascii files
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 *                     
 **************************************************************************/

void MakePMDMappingCDB(){

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");	
  AliPMDMappingData *mapda = new AliPMDMappingData();


  Int_t moduleno, totPatchBus, bPatchBus, ePatchBus;
  Int_t ibus, totmcm, rows, rowe, cols, cole;
  Int_t moduleDDL[6] = {12,12,0,0,12,12};

  ifstream infile;



  for(Int_t iddl = 0; iddl < 6; iddl++)
    {
      if(iddl == 0) infile.open("PMD_Mapping_ddl0.dat"); // ascii file
      if(iddl == 1) infile.open("PMD_Mapping_ddl1.dat"); // ascii file
      if(iddl == 2) infile.open("PMD_Mapping_ddl2.dat"); // ascii file
      if(iddl == 3) infile.open("PMD_Mapping_ddl3.dat"); // ascii file
      if(iddl == 4) infile.open("PMD_Mapping_ddl4.dat"); // ascii file
      if(iddl == 5) infile.open("PMD_Mapping_ddl5.dat"); // ascii file
      
      if(!infile)
	Error("Could not read the mapping file for DDL No = 0");

      Int_t modulePerDDL = moduleDDL[iddl];

      for (Int_t im = 0; im < modulePerDDL; im++)
	{
	  infile >> moduleno;
	  infile >> totPatchBus >> bPatchBus >> ePatchBus;
	  mapda->SetPatchBus(iddl,moduleno,bPatchBus,ePatchBus);
	  if (totPatchBus == 0) continue;
	  
	  for(Int_t i=0; i<totPatchBus; i++)
	    {
	      infile >> ibus >> totmcm >> rows >> rowe >> cols >> cole;
	      
	      printf("%d %d %d %d %d %d %d %d \n",moduleno,totPatchBus,
		     ibus,totmcm,rows,rowe,cols,cole);

	      mapda->SetModuleNo(iddl,ibus,moduleno);
	      mapda->SetMcmperBus(iddl,ibus,totmcm);
	      mapda->SetRowBus(iddl,ibus,rows,rowe);
	      mapda->SetColBus(iddl,ibus,cols,cole);


	    }
	  
	}
      
      infile.close();

    }

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Basanta Nandi");
  md->SetComment("Patchbus Mapping info for PMD");

  AliCDBId id("PMD/Calib/Mapping",0,AliCDBRunRange::Infinity());

  man->GetDefaultStorage()->Put(mapda,id, md);

}
