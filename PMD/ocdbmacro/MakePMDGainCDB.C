/**************************************************************************
 * To Create PMD Default Gain Object. 
 * sjena@cern.ch
 * Mon Nov 22 19:54:27 CET 2010
 * OCDB/PMD/Calib/Gain                    
 **************************************************************************/

void MakePMDGainCDB(TString type="IDEAL"){

	AliCDBManager* man = AliCDBManager::Instance();
	
	if(type == "IDEAL"){
		man->SetDefaultStorage("local://CDB_IDEAL");
	}else if (type == "DECALIB"){
		man->SetDefaultStorage("local://CDB_DECALIB");
	}else{
		cout << "Not a valid type!" << endl;
		break;
	
	}
	
	AliPMDCalibData *calibda = new AliPMDCalibData();
	
	TRandom random;
	AliCDBId id("PMD/Calib/Gain",0,0);
	
	const Int_t kDet = 2;
	const Int_t kMod = 24;
	const Int_t kRow = 48;
	const Int_t kCol = 96;

	if(type == "IDEAL"){
		// SET 1 (IDEAL)	
		for(int a=0;a<kDet;a++) 
			for(int b=0;b<kMod;b++) 
				for(int c=0;c<kRow;c++) 
					for(int d=0;d<kCol;d++) 
						calibda->SetGainFact(a, b, c, d, random.Gaus(15,2));
		id.SetRunRange(0,50);
		
		
		
	} else if (type == "DECALIB"){
		// SET 2 (DECALIB)	
		for(int a=0;a<kDet;a++) 
			for(int b=0;b<kMod;b++) 
				for(int c=0;c<kRow;c++) 
					for(int d=0;d<kCol;d++) 
						calibda->SetGainFact(a, b, c, d, TMath::Abs(random.Gaus(5,0.2)));
		id.SetRunRange(0,25);
	}


	
	
	AliCDBMetaData md;
	md.SetResponsible("Satyajit Jena");
	md.SetComment("Default Gain CDB");
	man->Put(calibda, id, &md);


}
