
void DBStorageSurvey(){

AliCDBManager *man = AliCDBManager::Instance();

AliCDBStorage *storLoc;
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

// Setting fake survey values :

AliVZEROSurveyData *surveyda = new AliVZEROSurveyData("Survey");

Float_t ngA[3] = { -22.327, -22.427, -84.7 };
surveyda->SetPointA(ngA);
Float_t ngB[3] = { 22.927, -22.427, -84.7 };
surveyda->SetPointB(ngB);
Float_t ngC[3] = { 22.927, 22.827, -84.7 };
surveyda->SetPointC(ngC);
Float_t ngD[3] = { -22.327, 22.827, -84.7 };
surveyda->SetPointD(ngD);

// Creation of the object VZERO Survey as a MetaData
	
TObjString str("VZERO Survey");      // object that will be stored

AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

AliCDBId id("VZERO/Survey/Data",0,9999999);

//md->SetObjectClassName("VZERO survey parameters"); automatically 
//set to AliVZEROSurveyData by the CDB classes during storage 

md->SetResponsible("Brigitte Cheynis");
md->SetBeamPeriod(0);
md->SetAliRootVersion("March2007");
md->SetComment("Prototype");
md->PrintMetaData();

storLoc = man->GetDefaultStorage();
storLoc->Put(surveyda, id, md);

storLoc->Delete();
delete md;

}
