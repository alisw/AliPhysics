// Simple example how to get CDB information stored in the ESD:
// - The list of parameters got from OCDB during reconstruction (cdbList)
// - The map of storages (default + specific) activated (cdbMap)
//
// author: alberto.colla@cern.ch

void GetCDBInfoFromESD(){

	TFile * f = new TFile("AliESDs.root");

	TTree* tree = f->Get("esdTree");
	TList* l = tree->GetUserInfo();

	TList* ids = l->FindObject("cdbList");
	ids->Print();

	TMap* storages = l->FindObject("cdbMap");
	storages->Print();
}
