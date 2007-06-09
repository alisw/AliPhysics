void TestMultialiases(){
// Simple test for the multi request query

gSystem->Load("AliDCSClient.so");
 AliDCSClient client("aldcs053.cern.ch",4242,1000,5);


//TObjArray *arr = new TObjArray();
//arr->SetOwner(1);
//client.GetAliasValues("tpc_PT_322.Temperature", 1180686465, 1180686575, arr);



TList list;
list.Add(new TObjString("tpc_PT_322.Temperature"));
list.Add(new TObjString("tpc_PT_323.Temperature"));
list.Add(new TObjString("tpc_PT_324.Temperature"));
list.Add(new TObjString("tpc_PT_325.Temperature"));
list.Add(new TObjString("tpc_PT_326.Temperature"));
list.Add(new TObjString("tpc_PT_327.Temperature"));
list.Add(new TObjString("tpc_PT_328.Temperature"));
list.Add(new TObjString("tpc_PT_329.Temperature"));
list.Add(new TObjString("tpc_PT_330.Temperature"));
list.Add(new TObjString("tpc_PT_331.Temperature"));

TMap *map = client.GetAliasValues(&list, 1180586575, 1180686575, 2, 4);

TIter iter(map);
TObjString *objstr=0;

while(objstr = dynamic_cast<TObjString*>(iter.Next())){

   cout << objstr->GetName() << endl;
   TObjArray *arr = map->GetValue(objstr->GetName());   

   cout << "N of values: " << arr->GetEntries() << endl;
   //arr->Print(); 
}


}

