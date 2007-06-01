void TestMultiAliases(){
// Simple test for the multi request query

gSystem->Load("AliDCSClient.so");
 AliDCSClient client("192.168.39.54",4242,1000,5);


//TObjArray *arr = new TObjArray();
//arr->SetOwner(1);
//client.GetAliasValues("tpc_PT_322.Temperature", 1180686465, 1180686575, arr);



TMap map;
map.Add(new TObjString("tpc_PT_322.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_323.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_324.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_325.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_326.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_327.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_328.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_329.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_330.Temperature"), new TObjArray());
map.Add(new TObjString("tpc_PT_331.Temperature"), new TObjArray());

client.GetAliasValues(1180586575, 1180686575, map);

TIter iter(&map);
TObjString *objstr=0;

while(objstr = dynamic_cast<TObjString*>(iter.Next())){

   cout << objstr->GetName() << endl;
   TObjArray *arr = map.GetValue(objstr->GetName());   

   cout << "N of values: " << arr->GetEntries() << endl;
   //arr->Print(); 
}


}

