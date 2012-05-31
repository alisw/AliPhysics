void TestMultialiases(){
// Simple test for the multi request query

gSystem->Load("$ALICE_ROOT/SHUTTLE/DCSClient/AliDCSClient.so");
 AliDCSClient client("aldcs053.cern.ch",4242,1000,5,100);

//TObjArray *arr = new TObjArray();
//arr->SetOwner(1);
//client.GetAliasValues("tpc_PT_322.Temperature", 1180686465, 1180686575, arr);



TList list;
list.Add(new TObjString("tof_hv_vp_00"));
list.Add(new TObjString("tof_hv_vp_01"));
list.Add(new TObjString("tof_hv_vp_02"));
list.Add(new TObjString("tof_hv_vp_03"));
list.Add(new TObjString("tof_hv_vp_04"));

TMap *map = client.GetAliasValues(&list, 1180586575, 1180686575, 0, 5);

TIter iter(map);
TObjString *objstr=0;

while(objstr = dynamic_cast<TObjString*>(iter.Next())){

   cout << objstr->GetName() << endl;
   TObjArray *arr = map->GetValue(objstr->GetName());   

   cout << "N of values: " << arr->GetEntries() << endl;
   //arr->Print(); 
}


}

