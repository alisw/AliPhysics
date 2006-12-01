#include "Riostream.h"
#include "TChain.h"

void CreateXML(const char *name, Int_t from, Int_t to, const char *path)
{
	const char *fixed1 = "<file name=\"AliESDs.root\" aclId=\"\" ctime=\"2006-08-14 17:17:43\" dir=\"2736575\" entryId=\"2736577\" expiretime=\"\" gowner=\"aliprod\" guid=\"0096DCE6-F62B-14DE-97BF-5519E793BEEF\"";
	const char *fixed2 = "md5=\"aaa\" owner=\"pulvir\" perm=\"755\" replicated=\"0\" seStringlist=\",72,\" size=\"19050411\"";
	
	char fileName[100];
	sprintf(fileName, "%s.xml", name);
	fstream file(fileName, ios::out);

	// header
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<alien>" << endl;
	file << "\t<collection name=\"" << name << "\">" << endl;
	
	// events
	for (int i = from; i <= to; i++) {
		file << "\t\t<event name=\"" << i << "\">" << endl;
		file << "\t\t\t" << fixed1 << " ";
		file << "lfn=\"" << path << "/" << i << "/AliESDs.root\" ";
		file << fixed2 << " turl=\"" << path << "/" << i << "/AliESDs.root\" type=\"f\" />" << endl;
		file << "\t\t</event>" << endl;
	}
	
	// footer
	file << "\t</collection>" << endl;
	file << "</alien>" << endl;
	file.close();
}

Bool_t AliRsnSelectRL
(
	Int_t        first_event = 1,
	Int_t        last_event = 10,
	const char * path = "/home/pulvir/events/head-2006-11-22/hijing",
    const char * chainName = "esdTree",
    Long64_t     nentries = TChain::kBigNumber
)
{
	gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_ROOT/PWG0/ -I$ROOTSYS/include -I$ALICE_ROOT/PWG2/RESONANCES");
	gSystem->Load("libXMLIO.so");
	gSystem->Load("libRAliEn.so");
	gSystem->Load("libPWG2.so");
	
	// Create automatically the XML collection
	char collection[255], collectionFile[255];
	sprintf(collection, "local_%d-%d", first_event, last_event);
	sprintf(collectionFile, "local_%d-%d.xml", first_event, last_event);
	CreateXML(collection, first_event, last_event, path);
	
	Info("Run", Form("Creating the collection from %s", collectionFile));
	TAlienCollection * myCollection = new TAlienCollection(collectionFile);
	if (!myCollection) {
		Error("Run", Form("Cannot create an AliEn collection from %s", collectionFile));
		return kFALSE;
	}
	
	Info("Run", Form("Creating the analysis chain %s", chainName));
	TChain* analysisChain = new TChain(chainName);
	
	Info("Run", "Preparing the file list");
	myCollection->Reset();
	while ( myCollection->Next() ) {
		char esdFile[255];
		sprintf(esdFile, "%s", myCollection->GetTURL(""));
		Info("Run", Form("Adding %s", esdFile));
		analysisChain->Add(esdFile);
	}
	
	Info("Run", "Setting parameters to selector...");
	AliRsnSelectorRL* selector = new AliRsnSelectorRL;
	selector->SetDebugFlag(kFALSE);
	selector->SetStoreKineInfo(1);
	selector->SetCheckITSRefit(1);
	selector->SetRejectFakes(0);
	selector->SetCopyMomentum(0);
	
	// Define prior probabilities (only for ESD PID)
	selector->SetPriorProbability(AliPID::kElectron, 0.042);
	selector->SetPriorProbability(AliPID::kMuon,     0.124);
	selector->SetPriorProbability(AliPID::kPion,     0.655);
	selector->SetPriorProbability(AliPID::kKaon,     0.041);
	selector->SetPriorProbability(AliPID::kProton,   0.130);
	
	// Define PID method
	selector->SetPIDMethod(AliRsnSelectorRL::kPerfectPID);  // kESDPID or kPerfectPID
	
	// Set output file name
	selector->SetOutputFile(Form("selection_%d-%d.root", first_event, last_event));
	
	// Process the chain
	Info("Run", "Processing the chain...");
	analysisChain->Process(selector, chainName, nentries, 0);
	
	return kTRUE;
}
