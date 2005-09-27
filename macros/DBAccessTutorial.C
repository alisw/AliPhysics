AliCDBStorage *storLoc, *storGrid;

AliCDBEntry *entry=0;
TObjString *objstr=0;

void DBAccessTutorial(){

AliCDBManager *man = AliCDBManager::Instance();
// pointer man points to the single instance of AliCDBManager. 
// This will let us avoid typing AliCDBManager::Instance every time... 


printf("\n<< TUTORIAL >> Activating Grid storage...\n");
storGrid = man->GetStorage("alien://aliendb4.cern.ch:9000;colla;DBGrid;ALICE::CERN::se01"); // replace "colla" with your username! 
// The default storage is automatically set to this. 
// One can access default storage with:
// (AliCDBStorage*) AliCDBManager::Instance()->GetDefaultStorage()
// To check the activation of the default storage:
// (Bool_t) AliCDBManager::Instance()->IsDefaultStorageSet()

printf("\n<< TUTORIAL >> Activating Local storage...\n");
storLoc = man->GetStorage("local://DBLocal");
// To set the default storage to this one:
// AliCDBMAnager::Instance()->SetDefaultStorage("local://DBLocal") or
// AliCDBMAnager::Instance()->SetDefaultStorage(storLoc) 

/////////////////////////////////////////////
// Step 0: Write/read in Local storage     //
/////////////////////////////////////////////
printf("\n<< TUTORIAL >> ********************************************\n");
printf(  "<< TUTORIAL >> **** Step 0: write/read in local storage ***\n");
printf(  "<< TUTORIAL >> ********************************************\n");


//create the new object

TObjString str1("This is step zero"); // object that will be stored

AliCDBId id1("ZDC/Calib/Gain",0,10); // Id of the object: AliCDBId("name", firstRun, lastRun)

AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
md1->SetObjectClassName("TObjString");
md1->SetResponsible("Alberto Colla");
md1->SetBeamPeriod(1);
md1->SetAliRootVersion("05-04-00"); //root version
md1->SetComment("This is a test");
TObjString str("test");
md1->SetProperty("key1",&str);


// Store the object into local storage
printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(&str1,id1, md1); // filename: DBLocal/ZDC/Calib/Gain/Run0_10_v0_s0.root


// read, update, store again
printf("\n<< TUTORIAL >> Retrieve object from local storage...\n");
entry = storLoc->Get("ZDC/Calib/Gain", 5);

objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

objstr->SetString("This is step 0.1: slightly better!"); // update object

printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(entry); // store into local: filename = Run0_10_v0_s1.root

// read, update, store again
printf("\n<< TUTORIAL >> Retrieve object from local storage...\n");
entry = storLoc->Get("ZDC/Calib/Gain", 5);

objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

objstr -> SetString("This is step 0.2: much better!");

printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(entry); // store into local: filename = Run0_10_v0_s2.root

///////////////////////////////////////////////////////////////////////
// Step 1: read from Local, update, store locally and into Grid      //
///////////////////////////////////////////////////////////////////////
printf("\n<< TUTORIAL >> ********************************************\n");
printf(  "<< TUTORIAL >> **** Step 1: write/read in Grid storage  ***\n");
printf(  "<< TUTORIAL >> ********************************************\n");

// read from local
printf("\n<< TUTORIAL >> Retrieve object from local storage...\n");
entry = storLoc->Get("ZDC/Calib/Gain", 5);
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());


objstr -> SetString("This is step 1: stored into Local and into Grid!");

printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(entry); // store into local: filename =  Run0_10_v0_s3.root

printf("\n<< TUTORIAL >> Storing object into Grid storage...\n");
storGrid->Put(entry); // store into grid: filename =  DBGrid/ZDC/Calib/Gain/Run0_10_v1.root



// step 2: read from Grid, update, store again (into Grid)

printf("\n<< TUTORIAL >> Retrieve object from Grid storage...\n");
entry = storGrid->Get("ZDC/Calib/Gain", 5);

objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

objstr -> SetString("This is step 2: update and store into Grid!");

printf("\n<< TUTORIAL >> Storing object into Grid storage...\n");
storGrid->Put(entry); // store into grid: filename =   Run0_10_v2.root

// step 3: read, update, store again (into Grid)
printf("\n<< TUTORIAL >> Retrieve object from Grid storage...\n");
entry = storGrid->Get("ZDC/Calib/Gain", 5);
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

objstr = (TObjString*) entry->GetObject();
objstr -> SetString("This is step 3: update and store into Grid!");

printf("\n<< TUTORIAL >> Storing object into Grid storage...\n");
storGrid->Put(entry); // store into grid: filename =   Run0_10_v3.root

  ////////////////////////////////////////////////
 // Step 3.0: read from Grid, store locally!   //
////////////////////////////////////////////////
printf("\n<< TUTORIAL >> **********************************************\n");
printf(  "<< TUTORIAL >> **** Step 3: read from Grid, store locally ***\n");
printf(  "<< TUTORIAL >> **********************************************\n");

printf("\n<< TUTORIAL >> Retrieve object from Grid storage...\n");
entry = storGrid->Get("ZDC/Calib/Gain", 5);
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(entry); // local: Run0_10_v3_s0.root

// Step 3.1: read from Local, update, store again into Local
printf("\n<< TUTORIAL >> Retrieve object from local storage...\n");
entry = storLoc->Get("ZDC/Calib/Gain", 5);
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

objstr->SetString("This is step 3.1: updated locally!");

printf("\n<< TUTORIAL >> Storing object into local storage...\n");
storLoc->Put(entry); // local: Run0_10_v3_s1.root

/////////////////////////////////////////////////////////////
// Step 3.2: read again from Grid version 3, store locally //
//         -> ERROR, local update already present!!        //
/////////////////////////////////////////////////////////////
printf("\n<< TUTORIAL >> **********************************************\n");
printf(  "<< TUTORIAL >> **** Step 3.2: error test                  ***\n");
printf(  "<< TUTORIAL >> **********************************************\n");

printf("\n<< TUTORIAL >> Retrieve object from Grid storage...\n");
entry = (AliCDBEntry*) storGrid->Get("ZDC/Calib/Gain",5);
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

printf("\n<< TUTORIAL >> Trying to store object into local storage...\n");
storLoc->Put(entry); // ERROR message!

/////////////////////////////////////////////////////////////
// Step 4: read from local, DRAIN to a dump storage: DBDrain.root! //
/////////////////////////////////////////////////////////////
printf("\n<< TUTORIAL >> ************************************************************\n");
printf(  "<< TUTORIAL >> **** Step 4: Read from local and DRAIN to a dump storage ***\n");
printf(  "<< TUTORIAL >> ************************************************************\n");

printf("\n<< TUTORIAL >> Setting Drain storage ...\n");
AliCDBManager::Instance()->SetDrain("dump://DBDrain.root"); //setting Drain storage

// Testing default storage behavior: let's set default storage to Local storage
AliCDBManager::Instance()->SetDefaultStorage(storLoc);

// read from local (default) storage. The object is automatically drained into the drain storage!
printf("\n<< TUTORIAL >> Retrieve object from local storage...\n");
entry = man->GetDefaultStorage()->Get("ZDC/Calib/Gain",5); 
objstr = (TObjString*) entry->GetObject();
printf("\n<< TUTORIAL >> Object string: %s\n", objstr->GetName());

/////////////////////////////////////////////////////////////
// Step 5: READ AND DRAIN multiple objects (with GetAll)   //
/////////////////////////////////////////////////////////////
printf("\n<< TUTORIAL >> *******************************************************\n");
printf(  "<< TUTORIAL >> **** Step 5: Read and Drain multiple objects        ***\n");
printf(  "<< TUTORIAL >> *******************************************************\n");


// Step 5.1: Store an object into four different Grid databases 

TObjString str2("This is the TPC/Calib/Gain object valid for runs 0 to 10.");
AliCDBId id2("TPC/Calib/Gain",0,10);

TObjString str3("This is the TPC/Calib/Drift object valid for runs 0 to 20.");
AliCDBId id3("TPC/Calib/Drift",0,20);

TObjString str4("This is the TPC/Align/Angles object valid for runs 0 to 15.");
AliCDBId id4("TPC/Align/Angles",0,15);

TObjString str5("This is the TPC/Align/Position object valid for runs 0 to 8.");
AliCDBId id5("TPC/Align/Positions",0,8);

printf("\n<< TUTORIAL >> Storing more objects into Grid storage...\n");
storGrid->Put(&str2,id2,md1);
storGrid->Put(&str3,id3,md1);
storGrid->Put(&str4,id4,md1);
storGrid->Put(&str5,id5,md1);

// Step 5.2: Read all the TPC objects with GetAll and drain into DBDrain.root

printf("\n<< TUTORIAL >> Retrieve more objects from Grid storage and drain them into Dump ...\n");
TList *list = (TList*)storGrid->GetAll("TPC/*",5);

// That's all folks! Delete AliCDBManager instance and metaData object

printf(  "<< TUTORIAL >> **** That's all folks!!        ***\n");


AliCDBManager::Instance()->Destroy();
delete entry;
delete md1;

}
