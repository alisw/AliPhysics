void TestAliSimpleValue() {

	gSystem->Load("libSHUTTLE");

	// Bool
	AliSimpleValue a(AliSimpleValue::kBool);
	Assert(a.GetType() == AliSimpleValue::kBool);
	a.SetBool(kTRUE);
	Assert(a.GetBool() == kTRUE);	

	// equal
	AliSimpleValue b(kTRUE);
	Assert(a == b);

	// copy constructor
	AliSimpleValue c(a);
	Assert(b == c);

	// set
	a = AliSimpleValue(kFALSE);
	Assert(a.GetBool() == kFALSE);

	
	// Byte
        AliSimpleValue a(AliSimpleValue::kByte);
        Assert(a.GetType() == AliSimpleValue::kByte);
        a.SetByte(30);
        Assert(a.GetByte() == 30);

	// equal
        AliSimpleValue b((Char_t) 30);
	Assert(b.GetType() == AliSimpleValue::kByte);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(a);
        Assert(b == c);

        // set
        a = AliSimpleValue((Char_t) 60);
        Assert(a.GetByte() == 60);

	
	// Int
        AliSimpleValue a(AliSimpleValue::kInt);
        Assert(a.GetType() == AliSimpleValue::kInt);
        a.SetInt(30);
        Assert(a.GetInt() == 30);

        // equal
        AliSimpleValue b(30);
        Assert(b.GetType() == AliSimpleValue::kInt);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(a);
        Assert(b == c);

        // set
        a = AliSimpleValue(60);
        Assert(a.GetInt() == 60);


	// UInt
        AliSimpleValue a(AliSimpleValue::kUInt);
        Assert(a.GetType() == AliSimpleValue::kUInt);
        a.SetUInt(30000);
        Assert(a.GetUInt() == 30000);

        // equal
        AliSimpleValue b((UInt_t) 30000);
        Assert(b.GetType() == AliSimpleValue::kUInt);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(a);
        Assert(b == c);

        // set
        a = AliSimpleValue((UInt_t) 60);
        Assert(a.GetUInt() == 60);
	

	// Float
        AliSimpleValue a(AliSimpleValue::kFloat);
        Assert(a.GetType() == AliSimpleValue::kFloat);
        a.SetFloat(100.53);
	Float_t delta = 1e-38;
	Assert((a.GetFloat() <= 100.53 + delta) || (a.GetFloat() >= 100.53 - delta));

        // equal
        AliSimpleValue b((Float_t) 100.53);
      	Assert(b.GetType() == AliSimpleValue::kFloat);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(a);
        Assert(b == c);

        // set
        a = AliSimpleValue((Float_t) 666.666);
        Assert(a.GetFloat() <= 666.666 + delta || a.GetFloat() >= 666.666 - delta);


	// DynBool
	Bool_t buf[] = {kTRUE, kFALSE, kTRUE};	
	AliSimpleValue a(3, buf);
	Assert(a.GetType() == AliSimpleValue::kDynBool);
	Assert(a.GetDynamicSize() == 3);
	Assert(a.GetDynBool(0) == kTRUE);	
	Assert(a.GetDynBool(1) == kFALSE);	
	Assert(a.GetDynBool(2) == kTRUE);	

	a.SetDynBool(1, kTRUE);
	Assert(a.GetDynBool(1) == kTRUE);	
	a.SetDynBool(1, kFALSE);

	// equal
	AliSimpleValue b(3, buf);
	Assert(a == b);
	
	// copy constructor
	AliSimpleValue c(b);
	Assert(c == b);

	// set
	a.SetDynBool(1, kTRUE);
	a = c;
	Assert(a == b);


	// DynByte
        Char_t cbuf[] = {10, 50, 60};
        AliSimpleValue a(3, cbuf);
        Assert(a.GetType() == AliSimpleValue::kDynByte);
        Assert(a.GetDynamicSize() == 3);
        Assert(a.GetDynByte(0) == 10);
        Assert(a.GetDynByte(1) == 50);
        Assert(a.GetDynByte(2) == 60);

        a.SetDynByte(1, 100);
        Assert(a.GetDynByte(1) == 100);
        a.SetDynByte(1, 50);

        // equal
        AliSimpleValue b(3, cbuf);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(b);
        Assert(c == b);

        // set
        a.SetDynByte(1, 200);
        a = c;
        Assert(a == b);


	// DynInt
        Int_t ibuf[] = {100, 500, -60};
        AliSimpleValue a(3, ibuf);
        Assert(a.GetType() == AliSimpleValue::kDynInt);
        Assert(a.GetDynamicSize() == 3);
        Assert(a.GetDynInt(0) == 100);
        Assert(a.GetDynInt(1) == 500);
        Assert(a.GetDynInt(2) == -60);

        a.SetDynInt(1, 100);
        Assert(a.GetDynInt(1) == 100);
        a.SetDynInt(1, 500);

        // equal
        AliSimpleValue b(3, ibuf);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(b);
        Assert(c == b);

        // set
        a.SetDynInt(1, 200);
        a = c;
        Assert(a == b);

        // DynUInt
        UInt_t ubuf[] = {100, 504, 7060};
        AliSimpleValue a(3, ubuf);
        Assert(a.GetType() == AliSimpleValue::kDynUInt);
        Assert(a.GetDynamicSize() == 3);
        Assert(a.GetDynUInt(0) == 100);
        Assert(a.GetDynUInt(1) == 504);
        Assert(a.GetDynUInt(2) == 7060);

        a.SetDynUInt(1, 100);
        Assert(a.GetDynUInt(1) == 100);
        a.SetDynUInt(1, 504);

        // equal
        AliSimpleValue b(3, ubuf);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(b);
        Assert(c == b);

        // set
        a.SetDynUInt(1, 200);
        a = c;
        Assert(a == b);


	// DynFloat
        Float_t fbuf[] = {100.1, 504.5, 7060};
        AliSimpleValue a(3, fbuf);
        Assert(a.GetType() == AliSimpleValue::kDynFloat);
        Assert(a.GetDynamicSize() == 3);
        Assert(a.GetDynFloat(0) <= 100.1 + delta || aGetDynFloat(0) >= 100.1 - delta);
        Assert(a.GetDynFloat(1) <= 504.5 + delta || aGetDynFloat(1) >= 504.5 - delta);
        Assert(a.GetDynFloat(2) <= 7060 + delta || aGetDynFloat(2) >= 7060 - delta);

        a.SetDynFloat(1, -200.3);
        Assert(a.GetDynFloat(1) <= -200.3 + delta || aGetDynFloat(1) >= -200.3 - delta);
        a.SetDynFloat(1, 504.5);


        // equal
        AliSimpleValue b(3, fbuf);
        Assert(a == b);

        // copy constructor
        AliSimpleValue c(b);
        Assert(c == b);

        // set
        a.SetDynFloat(1, 200);
        a = c;
        Assert(a == b);

	
}
