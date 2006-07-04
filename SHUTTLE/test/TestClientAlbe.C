Int_t GetValues(const char* host, Int_t port, const char* request,
	Long_t startTime, Long_t endTime) 
{

	AliDCSClient client(host, port, 1000, 5);

	Int_t result;

	TTimeStamp currentTime;
	TMap values;

	TString rString(request);

	TObjArray* requests = rString.Tokenize(",");

	cout<<"Requests: "<<requests->GetEntries()<<endl;

	TStopwatch sw;
	sw.Start();

	if (requests->GetEntries() > 1) {
		
		TIter iter(requests);
		TObjString* aString;
		while ((aString = (TObjString*) iter.Next())) {
			values.Add(new TObjString(aString->String()), NULL);
		}	
		
		result = client.GetDPValues(startTime, endTime, values);

	} else {
		TObjArray* valueSet = new TObjArray();
		valueSet->SetOwner(1);

		values.Add(new TObjString(request), valueSet);

		result = client.GetAliasValues(request, startTime,
				endTime, *valueSet);
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


/* 
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

 */
/*	TFile file("dump.root", "UPDATE");
	file.WriteTObject(values, "values");
	file.Close(); */

	values.DeleteAll();
	delete requests;

	cout<<"All values returned in runrange:  "<<endl;
	cout<<"StartTime: "<<TTimeStamp(startTime).AsString()<<endl;
	cout<<"EndTime: "<<TTimeStamp(endTime).AsString()<<endl;
	
	return result;
}

void TestClientAlbe(const char* host, Int_t port, const char* request,
	UInt_t startShift, UInt_t endShift) {

	gSystem->Load("libSHUTTLE");

//	AliLog::EnableDebug(kFALSE);
//	AliLog::SetGlobalDebugLevel(3);

	TH1F *histo=new TH1F("h","h",36,-6,25);
	Int_t result=-1;
	for(int i=0;i<100;i++){

		TTimeStamp currentTime;

		result=GetValues(host, port, request,
			currentTime.GetSec() - startShift, 
			currentTime.GetSec() - endShift);
			if(result != 22) printf("\n\n\n NO 22!!!!!!!! %d\n\n",result);
		histo->Fill((Float_t)result);
		if(i != 99)gSystem->Sleep(10000);
	}

	cout<<"Client done"<<endl;
	histo->Draw();
}

