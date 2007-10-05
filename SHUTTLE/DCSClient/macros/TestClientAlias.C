TMap* GetValues(const char* host, Int_t port, const char* request,
	UInt_t startTime, UInt_t endTime, Int_t multiSplit) 
{
	AliDCSClient client(host, port, 1000, 20, multiSplit);
	// The 5th parameter switches from single alias to multi aliases!

	//Int_t result;

	TTimeStamp currentTime;
	//TMap values;

	TString rString(request);

	TObjArray* requests = rString.Tokenize(",");

	cout<<"Requests: "<<requests->GetEntries()<<endl;

	TMap* values=0;
	
	TStopwatch sw;
	sw.Start();
	
	if(requests->GetEntries() == 0) return NULL;

// 		TIter iter(requests);
// 		TObjString* aString;
// 		TObjArray* valueSet;
// 		while ((aString = (TObjString*) iter.Next())) {
// 		      cout<<"  Querying: "<<aString->GetName()<<endl;
// 		      valueSet = new TObjArray();
// 		      valueSet->SetOwner(1);
// 
// 			result = client.GetAliasValues(aString->GetName(), startTime,
// 				endTime, valueSet);
// 			values.Add(aString->Clone(), valueSet);
// 		}
 
	values = client.GetAliasValues(requests, startTime, endTime);

	if (!values) {
		cout<<"Query failed! Result error: "<<
			client.GetErrorString(client.GetResultErrorCode()) <<endl;
		if(client.GetResultErrorCode() == AliDCSClient::fgkServerError)	
			cout<<"Server error: "<<
				client.GetServerError().Data() <<endl;
		return NULL;
	}
	
	sw.Stop();
	cout<<"Elapsed time: "<<sw.RealTime()<<endl;

	cout<<"Time per alias: "<<sw.RealTime()/requests->GetEntries()<<endl;

	cout<<"Received values: "<<endl;

	Int_t nValues=0;

	TIter iter(values);
	TObjString* aRequest;
	while ((aRequest = (TObjString*) iter.Next())) {

		TObjArray* valueSet = (TObjArray*) values->GetValue(aRequest);
		
		cout<<" '"<<aRequest->String()<<"' values: " 
			<<valueSet->GetEntriesFast()<<endl;

		TIter valIter(valueSet);
		AliDCSValue* aValue;
		while ((aValue = (AliDCSValue*) valIter.Next())) {
			cout<<aValue->ToString()<<endl;
			nValues++;
		} 
	}
	
	cout<<"Number of received values: "<< nValues <<endl;
	

/*
	TFile file("dump.root", "UPDATE");
	file.cd();
	values.Write("DCSAliasMap", TObject::kSingleKey);
	file.Close(); 
*/

	//values.DeleteAll();
	//delete requests;

	cout<<"All values returned in runrange:  "<<endl;
	cout<<"StartTime: "<<TTimeStamp(startTime).AsString()<<endl;
	cout<<"EndTime: "<<TTimeStamp(endTime).AsString()<<endl;
	
	return values;
}

TMap* TestClientAlias(const char* host, Int_t port, const char* request,
	UInt_t startShift, UInt_t endShift, UInt_t multiSplit) {

	gSystem->Load("$ALICE_ROOT/SHUTTLE/DCSClient/AliDCSClient");

//	AliLog::EnableDebug(kFALSE);
//	AliLog::SetGlobalDebugLevel(3);

	TTimeStamp currentTime;

	
/*	
	TMap* values = GetValues(host, port, request,
		currentTime.GetSec() - startShift, 
		currentTime.GetSec() - endShift, multiSplit);
*/
		
// SHUTTLE query interval 
	
	TMap* values = GetValues(host, port, request,
		1181300060, 1181307260, multiSplit);

	if(values) values->Print();

	cout << endl;
	cout <<"Client done"<<endl;
	
	return values;
}
