
// This macro transfers OCDB data of one single ("detName") or all ALICE detectors from one location to another, 
// It is possible to set new run range, path etc... 

void CDBToGrid(const char *detName="", const char* fromUri="local://$ALICE_ROOT", 
	const char* toUri="local://newOCDB")
{

AliCDBManager *man = AliCDBManager::Instance();

man->SetDefaultStorage(fromUri);

AliCDBStorage *dest = man->GetStorage(toUri);

man->SetRun(0);

TList* calib=0;
TList* align=0;
TList* recpar=0;
TList* config=0;

if(detName == ""){
// drain everything
	//calib = man->GetAll("*/Calib/*");
	//align = man->GetAll("*/Align/Data");
	//recpar = man->GetAll("*/RecParam/*");
	//config = man->GetAll("*/Config/*");
	
// drain everything, really!
	calib = man->GetAll("*");

} else {
// drain calibration and alignment for detector "detName"

	TString calibPath = detName;
	TString alignPath = calibPath;
	calibPath+="/Calib/*";	
	alignPath+="/Align/Data";

	//calib = man->GetAll(calibPath);
	//align = man->GetAll(alignPath);

// drain everything, really!
	calib = man->GetAll(Form("%s/*",detName.Data()));
}

AliCDBEntry *entry;

Int_t ok=0; TString failed="";
if(calib){
	ok=0; failed="";
	for(int i=0;i<calib->GetEntries();i++){
 
      	entry = (AliCDBEntry*) calib->At(i);
      	entry->GetId().SetRunRange(0,999999999);
       
      	TString path=entry->GetId().GetPath();
	
	cout << path.Data() << endl;
      
      	if (path == "ITS/Resp/RespSDD") entry->GetId().SetPath("ITS/Calib/RespSDD"); // bug in ITS/Resp/RespSDD
	
	if(path == "TOF/Calib/Par" || path == "TOF/Calib/SimPar" || path.Contains("TOF/CDB")) continue;
 
      	if (dest->Put(entry)) {
      		ok++;
      	} else {
      		failed += path.Data(); failed += " "; 
      	}
	}
	printf("\n************ CALIB *********** \n");
	printf("************ Stored %d objects over %d *********** \n", ok, calib->GetEntries());
	if(failed != ""){
		printf("***** List of failed objects: %s \n", failed.Data());
	}
}


if(align){
ok=0; failed="";
for(int i=0;i<align->GetEntries();i++){

      entry = (AliCDBEntry*) align->At(i);
      entry->GetId().SetRunRange(0,999999999);
      TString path=entry->GetId().GetPath();

      	if (dest->Put(entry)) {
             	ok++;
      	} else {
             	failed += path.Data(); failed += " "; 
      	}
	}
	printf("\n************ ALIGN *********** \n");
	printf("************ Stored %d objects over %d *********** \n", ok, align->GetEntries());
	if(failed != ""){
		printf("***** List of failed objects: %s \n", failed.Data());
	}
}

if(recpar){
	ok=0; failed="";
	for(int i=0;i<recpar->GetEntries();i++){

     	entry = (AliCDBEntry*) recpar->At(i);
     	entry->GetId().SetRunRange(0,999999999);
     	TString path=entry->GetId().GetPath();

      	if (dest->Put(entry)) {
		     ok++;
      	} else {
	     failed += path.Data(); failed += " "; 
      	}
	}
	printf("\n************ RECPAR *********** \n");
	printf("************ Stored %d objects over %d *********** \n", ok, recpar->GetEntries());
	if(failed != ""){
		printf("***** List of failed objects: %s \n", failed.Data());
	}
}

if(config){
	ok=0; failed = "";
	for(int i=0;i<config->GetEntries();i++){

     	entry = (AliCDBEntry*) config->At(i);
     	entry->GetId().SetRunRange(0,999999999);
     	TString path=entry->GetId().GetPath();

      	if (dest->Put(entry)) {
		     ok++;
      	} else {
		     failed += path.Data(); failed += " "; 
      	}
	}
	printf("\n************ CONFIG *********** \n");
	printf("************ Stored %d objects over %d *********** \n", ok, config->GetEntries());
	if(failed != ""){
		printf("***** List of failed objects: %s \n", failed.Data());
	}
}

man->Destroy();

}
