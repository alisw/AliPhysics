/**************************************************************************
 * To Convert Ascii DDL file into DDL-CDB object
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 *                     
 **************************************************************************/
void MakePMDDDLinfoCDB()
{

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");	
  //  man->SetDefaultStorage("local:///Users/basanta/ALISOFT/PMD/VarInit/OCDB");	
  AliPMDddlinfoData *mapda = new AliPMDddlinfoData();

  Int_t ddlno;
  Int_t modno, totmod;
  Int_t serialmodno[12];
  //  Int_t moduleDDL[6] = {12,12,0,0,12,12};

  ifstream infile;
  infile.open("PMD_ddl_info.dat"); // ascii file
  if(!infile) Error("Could not open the DDL info file");


  for(Int_t iddl = 0; iddl < 6; iddl++)
    {
      infile >> ddlno >> totmod;

      mapda->SetNoOfModulePerDdl(iddl,totmod);

      //printf("%d %d\n", ddlno, totmod);
      if(totmod == 0) continue;

      for(Int_t imod = 0; imod < 12; imod++)
	{
	  infile >> modno;
	  serialmodno[imod] = modno;
	  //printf("%d \n", modno);
	}

      mapda->SetModuleNoPerDdl(iddl,serialmodno);
    }

  infile.close();
  
  infile.open("PMD_removed_chains.dat"); // ascii file
  if(!infile) Error("Could not open the DDL info file");

  Int_t det, smn;
  Int_t rows1, rowe1, cols1, cole1;
  Int_t rows2, rowe2, cols2, cole2;

  Int_t srowa[2][24];
  Int_t erowa[2][24];
  Int_t scola[2][24];
  Int_t ecola[2][24];
  Int_t srowb[2][24];
  Int_t erowb[2][24];
  Int_t scolb[2][24];
  Int_t ecolb[2][24];


  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  infile >> det >> smn >> rows1 >> rowe1 >> cols1 >> cole1
		     >> rows2 >> rowe2 >> cols2 >> cole2;
      
	  srowa[idet][ismn] = rows1;
	  erowa[idet][ismn] = rowe1;
	  scola[idet][ismn] = cols1;
	  ecola[idet][ismn] = cole1;
	  srowb[idet][ismn] = rows2;
	  erowb[idet][ismn] = rowe2;
	  scolb[idet][ismn] = cols2;
	  ecolb[idet][ismn] = cole2;

	}
    }

  mapda->SetStartRowA(srowa);
  mapda->SetStartRowB(srowb);
  mapda->SetStartColA(scola);
  mapda->SetStartColB(scolb);
  mapda->SetEndRowA(erowa);
  mapda->SetEndRowB(erowb);
  mapda->SetEndColA(ecola);
  mapda->SetEndColB(ecolb);


  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Basanta Nandi");
  md->SetComment("DDL info for PMD");
  
  AliCDBId id("PMD/Calib/Ddlinfo",0,AliCDBRunRange::Infinity());

  man->GetDefaultStorage()->Put(mapda,id, md);

}
