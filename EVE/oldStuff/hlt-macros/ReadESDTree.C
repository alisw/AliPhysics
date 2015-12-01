int ReadESDTree(TString  fileName = "ESD_TPC.root")
{
	TFile *esdFile = TFile::Open(fileName.Data());
	AliESDEvent *esdEvent = new AliESDEvent();
	TTree *tr = (TTree*)esdFile->Get("esdTree");
	esdEvent->ReadFromTree(tr);
	tr->StartViewer();
 
	return 0;
}
