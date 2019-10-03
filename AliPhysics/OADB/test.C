void test()
{
  gSystem->Load("libOADB");
  con = new AliOADBContainer("OADB");
  //
  obj1 = new TNamed("obj1", "");
  obj2 = new TNamed("obj2", "");
  obj3 = new TNamed("obj3", "");
  obj4 = new TNamed("obj4", "");
  //
  con->AppendObject(obj1,  1,  10);
  con->AppendObject(obj2, 11,  20);
  con->AppendObject(obj3, 21,  30);
  con->UpdateObject(1, obj4, 100, 101);
  con->RemoveObject(0);
  // 
  con->WriteToFile("test.root");
  //
  AliOADBContainer cont0("");
  cont0.InitFromFile("test.root", "OADB");
  cont0.Dump();
  cont0.List();
}
