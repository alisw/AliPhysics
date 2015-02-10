/// \file TestAliTPCDigitsArray.C

TFile f("pokus.root","recreate")
arr.MakeTree()
TBrowser b;
arr.AddSegment(1)
arr[1]
arr[1]->Dump()
arr.AddSegment(2)
arr[2]->Dump()
dig.StoreSegment(1)
arr.StoreSegment(1)
arr.StoreSegment(2)
arr.SetClass("AliDigits")
arr.ConnectTree()
TBrowser b;
arr.SetClass("AliDigits")
arr.MakeArray(100)
arr.AddSegment(1)
arr.AddSegment(2)
arr.StoreSegment(1)
arr.StoreSegment(2)
arr.GetTree()->Write("pokus1")
arr.GetTree()->Close()
arr.SetClass("AliSimDigits")
arr.MakeArray(1000)
arr.ConnectTree("pokus1")
arr.LoadSegment(1)
arr[1]
arr[1]->Dump()
arr[2]
arr.LoadSegment(2)
arr[2]->Dump()
.q
.! ddd root.exe &
.q
.! ddd root.exe &
AliTPCDigitsArray arr
AliTPCDigitsArray brr
.q
AliTPCDigitsArray arr
AliTPCParam par;
arr.Setup(&par)
arr.GetRow(1,1)
.q
.x mkhtml.C
.q
AliTPCDigitsArray arr;
arr.CreateRow(0,4)
.q
AliTPCDigitsArray arr;
AliTPCParam par;
arr.Setup(&par)
arr.CreateRow(0,4)
arr.Getrow(9,4)
arr.GetRow(0,4)
arr.GetRow(0,4)->Dump()
arr.GetRow(0,4)->IsA()->GetName()
AliTPCParam par;
par.Update()
par.GetMaxTBin()
.q
AliTPCParam par;
par.GetMaxTBin()
AliTPCDigitsArray arr;
arr.Setup(&par)
arr.CreateRow(0,4)
arr.GetRow(0,4)->Dump()
arr.GetIndex(0,4)
par.GetIndex(0,4)
arr.GetParam()->GetIndex(0,4)
(AliTPCParam* arr.GetParam())->GetIndex(0,4)
(AliTPCParam* arr.GetParam()))->GetIndex(0,4)
( (AliTPCParam*) arr.GetParam())->GetIndex(0,4)
arr.GetRow(0,4)->GetSegmentID()
arr.GetRow(0,4)->GetID()
