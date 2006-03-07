void Compare(TObjArray& a, TObjArray& b) {

	Assert(a.GetSize() == b.GetSize());

	TIter iterA(&a);
	TIter iterB(&b);
	
	AliDCSValue* valA;
	AliDCSValue* valB;
 
	while ((valA = (AliDCSValue*) iterA.Next())) {
		valB = (AliDCSValue*) iterB.Next();
		
		cout<<valA->ToString()<<endl;
		cout<<valB->ToString()<<endl;

		Assert(valA->GetSimpleValue() == valB->GetSimpleValue());
		Assert(valA->GetTimeStamp() == valB->GetTimeStamp());
	}
	
}

void CompareStr(TObjArray& a, TObjArray& b) {

	cout<<"A size: "<<a.GetSize()<<endl;
	cout<<"B size: "<<b.GetSize()<<endl;

        Assert(a.GetSize() == b.GetSize());

        TIter iterA(&a);
        TIter iterB(&b);

        TObjString* valA;
        TObjString* valB;

        while ((valA = (TObjString*) iterA.Next())) {
                valB = (TObjString*) iterB.Next();
         
                cout<<valA->String()<<endl;
                cout<<valB->String()<<endl;

		Assert(valA->String() == valB->String());
        }
}

void PrintBuffer(const char* buffer, UInt_t size) {
	
	cout<<"BUFFER START"<<endl;

	for (UInt_t k = 0; k < size; k ++) {
		cout<<" "<<(UInt_t) (UChar_t) buffer[k]<<" ";
		if (!(k + 1 % 50)) {
			cout<<endl;
		}
	}
	
	cout<<endl<<"BUFFER END"<<endl;
}

void TestMessage() {

	gSystem->Load("libSHUTTLE");

	TTimeStamp currentTime;
	TString requestString("value.int.something");

	// Request Message
	AliDCSMessage requestMsg;
	requestMsg.CreateRequestMessage(AliDCSMessage::kAlias, 
		currentTime.GetSec(), currentTime.GetSec() + 10, 
		requestString.Data());
	Assert(requestMsg.GetType() == AliDCSMessage::kRequest);
	Assert(requestMsg.GetRequestType() == AliDCSMessage::kAlias);
	Assert(requestMsg.GetStartTime() == currentTime.GetSec());
	Assert(requestMsg.GetEndTime() == currentTime.GetSec() + 10);
	Assert(requestMsg.GetRequestString() == requestString); 
	requestMsg.StoreToBuffer();

	AliDCSMessage newRequestMsg(requestMsg.GetMessage(),
		requestMsg.GetMessageSize());
	newRequestMsg.LoadFromBuffer();

	Assert(newRequestMsg.GetType() == AliDCSMessage::kRequest);
        Assert(newRequestMsg.GetRequestType() == AliDCSMessage::kAlias);
        Assert(newRequestMsg.GetStartTime() == currentTime.GetSec());
        Assert(newRequestMsg.GetEndTime() == currentTime.GetSec() + 10);
        Assert(newRequestMsg.GetRequestString() == requestString);


	// Count Message
	AliDCSMessage countMsg;
	countMsg.CreateCountMessage(200);
	
	Assert(countMsg.GetType() == AliDCSMessage::kCount);
	Assert(countMsg.GetCount() == 200);
	countMsg.StoreToBuffer();

	AliDCSMessage newCountMsg(countMsg.GetMessage(), 
		countMsg.GetMessageSize());
	newCountMsg.LoadFromBuffer();

	Assert(newCountMsg.GetType() == AliDCSMessage::kCount);
        Assert(newCountMsg.GetCount() == 200);

	// ResultSet Message

	Int_t maxCount = 3;

	// Bool
	AliDCSMessage rsMsg;
	rsMsg.CreateResultSetMessage(AliSimpleValue::kBool);

        Assert(rsMsg.GetType() == AliDCSMessage::kResultSet);

        TObjArray values;
	values.SetOwner(1);
        values.AddLast(new AliDCSValue(kFALSE, currentTime.GetSec()));
        values.AddLast(new AliDCSValue(kTRUE, currentTime.GetSec() + 1));
        values.AddLast(new AliDCSValue(kFALSE, currentTime.GetSec()+ 200));

        rsMsg.AddValue(*(AliDCSValue*) values.At(0));
        rsMsg.AddValue(*(AliDCSValue*) values.At(1));
        rsMsg.AddValue(*(AliDCSValue*) values.At(2));
        rsMsg.StoreToBuffer();

	PrintBuffer(rsMsg.GetMessage(), rsMsg.GetMessageSize());

        AliDCSMessage nrsMsg(rsMsg.GetMessage(), rsMsg.GetMessageSize());
        nrsMsg.LoadFromBuffer();

        Assert(nrsMsg.GetType() == AliDCSMessage::kResultSet);
	Assert(nrsMsg.GetSimpleValueType() == AliSimpleValue::kBool);

        TObjArray nvalues;
	nvalues.SetOwner(1);
        Assert(nrsMsg.GetValues(nvalues) == maxCount);

        Compare(values, nvalues);

	// Byte
        AliDCSMessage rsMsg;
	rsMsg.CreateResultSetMessage(AliSimpleValue::kByte);

        Assert(rsMsg.GetType() == AliDCSMessage::kResultSet);

        TObjArray values;
	values.SetOwner(1);
        values.AddLast(new AliDCSValue((Char_t) 25, currentTime.GetSec()));
        values.AddLast(new AliDCSValue((Char_t) 0, currentTime.GetSec() + 1));
        values.AddLast(new AliDCSValue((Char_t) 255, currentTime.GetSec()+ 200));

        rsMsg.AddValue(*(AliDCSValue*) values.At(0));
        rsMsg.AddValue(*(AliDCSValue*) values.At(1));
        rsMsg.AddValue(*(AliDCSValue*) values.At(2));
        rsMsg.StoreToBuffer();

        AliDCSMessage nrsMsg(rsMsg.GetMessage(), rsMsg.GetMessageSize());
        nrsMsg.LoadFromBuffer();

        Assert(nrsMsg.GetType() == AliDCSMessage::kResultSet);
        Assert(nrsMsg.GetSimpleValueType() == AliSimpleValue::kByte);

        TObjArray nvalues;
	nvalues.SetOwner(1);
        Assert(nrsMsg.GetValues(nvalues) == maxCount);

        Compare(values, nvalues);


	// Int
	AliDCSMessage rsMsg;
	rsMsg.CreateResultSetMessage(AliSimpleValue::kInt);

	Assert(rsMsg.GetType() == AliDCSMessage::kResultSet);

	TObjArray values;
	values.SetOwner(1);
	values.AddLast(new AliDCSValue(100, currentTime.GetSec()));
	values.AddLast(new AliDCSValue(10, currentTime.GetSec() + 1));
	values.AddLast(new AliDCSValue(666, currentTime.GetSec()+ 200));
	
	rsMsg.AddValue(*(AliDCSValue*) values.At(0));
	rsMsg.AddValue(*(AliDCSValue*) values.At(1));
	rsMsg.AddValue(*(AliDCSValue*) values.At(2));
	rsMsg.StoreToBuffer();

	AliDCSMessage nrsMsg(rsMsg.GetMessage(), rsMsg.GetMessageSize());
	nrsMsg.LoadFromBuffer();

	Assert(nrsMsg.GetType() == AliDCSMessage::kResultSet);
	Assert(nrsMsg.GetSimpleValueType() == AliSimpleValue::kInt);

        TObjArray nvalues;
	nvalues.SetOwner(1);
	Assert(nrsMsg.GetValues(nvalues) == maxCount);

	Compare(values, nvalues);
	
	// UInt
        AliDCSMessage rsMsg;
	rsMsg.CreateResultSetMessage(AliSimpleValue::kUInt);

        Assert(rsMsg.GetType() == AliDCSMessage::kResultSet);

        TObjArray values;
	values.SetOwner(1);
        values.AddLast(new AliDCSValue((UInt_t) 1000, currentTime.GetSec()));
        values.AddLast(new AliDCSValue((UInt_t) 104, currentTime.GetSec() + 1));
        values.AddLast(new AliDCSValue((UInt_t) 6665, currentTime.GetSec()+ 200));

        rsMsg.AddValue(*(AliDCSValue*) values.At(0));
        rsMsg.AddValue(*(AliDCSValue*) values.At(1));
        rsMsg.AddValue(*(AliDCSValue*) values.At(2));
        rsMsg.StoreToBuffer();

	AliDCSMessage nrsMsg(rsMsg.GetMessage(), rsMsg.GetMessageSize());
        nrsMsg.LoadFromBuffer();

        Assert(nrsMsg.GetType() == AliDCSMessage::kResultSet);
        Assert(nrsMsg.GetSimpleValueType() == AliSimpleValue::kUInt);

        TObjArray nvalues;
	nvalues.SetOwner(1);
        Assert(nrsMsg.GetValues(nvalues) == maxCount);

	Compare(values, nvalues); 

	// Float
        AliDCSMessage rsMsg;
	rsMsg.CreateResultSetMessage(AliSimpleValue::kFloat);

        Assert(rsMsg.GetType() == AliDCSMessage::kResultSet);

        TObjArray values;
	values.SetOwner(1);
        values.AddLast(new AliDCSValue((Float_t) 1000.55, currentTime.GetSec()));
        values.AddLast(new AliDCSValue((Float_t) 10.4, currentTime.GetSec() + 1));
        values.AddLast(new AliDCSValue((Float_t) 6665.1, currentTime.GetSec()+ 200));

        rsMsg.AddValue(*(AliDCSValue*) values.At(0));
        rsMsg.AddValue(*(AliDCSValue*) values.At(1));
        rsMsg.AddValue(*(AliDCSValue*) values.At(2));
        rsMsg.StoreToBuffer();

        AliDCSMessage nrsMsg(rsMsg.GetMessage(), rsMsg.GetMessageSize());
        nrsMsg.LoadFromBuffer();

        Assert(nrsMsg.GetType() == AliDCSMessage::kResultSet);
        Assert(nrsMsg.GetSimpleValueType() == AliSimpleValue::kFloat);

        TObjArray nvalues;
	nvalues.SetOwner(1);
        Assert(nrsMsg.GetValues(nvalues) == maxCount);

        Compare(values, nvalues);

	
	// Error Message
	TString errorString("This is an error string");

	AliDCSMessage errorMsg;
	errorMsg.CreateErrorMessage(AliDCSMessage::kUnknownAliasDPName, 
			errorString);		
	Assert(errorMsg.GetType() == AliDCSMessage::kError);
	Assert(errorMsg.GetErrorCode() == AliDCSMessage::kUnknownAliasDPName);
	Assert(errorMsg.GetErrorString() == errorString); 
	errorMsg.StoreToBuffer();

	AliDCSMessage newErrorMsg(errorMsg.GetMessage(),
		errorMsg.GetMessageSize());
	newErrorMsg.LoadFromBuffer();

        Assert(newErrorMsg.GetType() == AliDCSMessage::kError);
        Assert(newErrorMsg.GetErrorCode() 
		== AliDCSMessage::kUnknownAliasDPName);
        Assert(newErrorMsg.GetErrorString() == errorString); 

	// MultiRequest Message
	TObjArray requests;
	requests.SetOwner(1);
	requests.AddLast(new TObjString("alias1.test"));
	requests.AddLast(new TObjString("alias2.test"));
	requests.AddLast(new TObjString("alias3.test"));

	AliDCSMessage multiMsg;
	multiMsg.CreateMultiRequestMessage(AliDCSMessage::kAlias, 
			currentTime.GetSec(),
			currentTime.GetSec() + 10);
	multiMsg.AddRequestString("alias1.test");	
	multiMsg.AddRequestString("alias2.test");	
	multiMsg.AddRequestString("alias3.test");	
	
	Assert(multiMsg.GetType() == AliDCSMessage::kMultiRequest);
        Assert(multiMsg.GetRequestType() == AliDCSMessage::kAlias);
        Assert(multiMsg.GetStartTime() == currentTime.GetSec());
        Assert(multiMsg.GetEndTime() == currentTime.GetSec() + 10);
	multiMsg.StoreToBuffer();	

	PrintBuffer(multiMsg.GetMessage(), multiMsg.GetMessageSize());
	
	AliDCSMessage nmultiMsg(multiMsg.GetMessage(), 
		multiMsg.GetMessageSize());
	nmultiMsg.LoadFromBuffer();
	
	Assert(nmultiMsg.GetType() == AliDCSMessage::kMultiRequest);
        Assert(nmultiMsg.GetRequestType() == AliDCSMessage::kAlias);
        Assert(nmultiMsg.GetStartTime() == currentTime.GetSec());
        Assert(nmultiMsg.GetEndTime() == currentTime.GetSec() + 10); 

	PrintBuffer(nmultiMsg.GetMessage(), nmultiMsg.GetMessageSize());
	
	TObjArray nrequests;
	nrequests.SetOwner(1);
	nmultiMsg.GetRequestStrings(nrequests);
	
	CompareStr(requests, nrequests);
}

