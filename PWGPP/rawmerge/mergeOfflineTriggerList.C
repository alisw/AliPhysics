/*
   .L $ALICE_PHYSICS/../src/PWGPP/rawmerge/mergeOfflineTriggerList.C
   mergeOfflineTriggerList("xxx109.xml");
   or ...
   aliroot -b -q $ALICE_PHYSICS/../src/PWGPP/rawmerge/mergeOfflineTriggerList.C\(\"xxx109.xml\"\)
*/

void mergeOfflineTriggerList(const char * inputList){
  AliOfflineTrigger trigger("default", 30,100000000);
  trigger.LoadTriggerList("filteredHighPtV0s.list");
  trigger.LoadTriggerList("filteredMult.list");
  trigger.LoadTriggerList("filteredHighPt.list");
  trigger.ExtractSelected(inputList,"" , "rawSelected[].root",1000000000, 1  );
}


