void GetValues(const char* host, Int_t port, const char* request,
	UInt_t startTime, UInt_t endTime) 
{
	AliDCSClient client(host, port, 1000, 20);

	Int_t result;

	TTimeStamp currentTime;
	TMap values;

	TString rString(request);

	TObjArray* requests = rString.Tokenize(",");

	cout<<"Requests: "<<requests->GetEntries()<<endl;

	TStopwatch sw;
	sw.Start();

	if (requests->GetEntries() > 0) {
		TIter iter(requests);
		TObjString* aString;
		TObjArray* valueSet;
		while ((aString = (TObjString*) iter.Next())) {
			cout<<"  Querying: "<<aString->GetName()<<endl;
			valueSet = new TObjArray();
			valueSet->SetOwner(1);

			result = client.GetAliasValues(aString->GetName(), startTime,
				endTime, valueSet);
			values.Add(aString->Clone(), valueSet);
		}
	}

	if (result < 0) {
		cout<<"Communication failure: "<<
			AliDCSClient::GetErrorString(result)<<endl;

		if (result == AliDCSClient::fgkServerError) {
			cout<<"Server error code: "<<
				client.GetServerErrorCode()<<endl;
			cout<<client.GetServerError()<<endl;
		}
	}
	
	sw.Stop();
	cout<<"Elapsed time: "<<sw.RealTime()<<endl;
	if (result > 0) {
		cout<<"Time per value: "<<sw.RealTime()/result<<endl;
	}
	cout<<"Received values: "<<result<<endl;

	TIter iter(&values);
	TObjString* aRequest;
	while ((aRequest = (TObjString*) iter.Next())) {

		TObjArray* valueSet = (TObjArray*) values.GetValue(aRequest);
		
		cout<<" '"<<aRequest->String()<<"' values: " 
			<<valueSet->GetEntriesFast()<<endl;

		TIter valIter(valueSet);
		AliDCSValue* aValue;
		while ((aValue = (AliDCSValue*) valIter.Next())) {
			cout<<aValue->ToString()<<endl;
		} 
	}

/*
	TFile file("dump.root", "UPDATE");
	file.cd();
	values.Write("DCSAliasMap", TObject::kSingleKey);
	file.Close(); 
*/

	values.DeleteAll();
	delete requests;

	cout<<"All values returned in runrange:  "<<endl;
	cout<<"StartTime: "<<TTimeStamp(startTime).AsString()<<endl;
	cout<<"EndTime: "<<TTimeStamp(endTime).AsString()<<endl;
}

void TestClientAlias(const char* host, Int_t port, const char* request,
	UInt_t startShift, UInt_t endShift) {

	gSystem->Load("$ALICE_ROOT/SHUTTLE/DCSClient/AliDCSClient");

//	AliLog::EnableDebug(kFALSE);
//	AliLog::SetGlobalDebugLevel(3);

	TTimeStamp currentTime;

	GetValues(host, port, request,
		currentTime.GetSec() - startShift, 
		currentTime.GetSec() - endShift);

	cout<<"Client done"<<endl;
}
