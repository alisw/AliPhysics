/********************************************************************
 
 Bad Chunks Checking code, 15th April 2013
 
 --- Here you just create an array of TStrings with dataset names and
 this function will automatically compile the checking method and process
 all the datasets in the TString array. 
 
 ********************************************************************/

void ProcessDatasets(){
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;
  //Compile Macro
  Int_t workedornot = gSystem->CompileMacro("ProcessBadChunks02.C","-kfo");
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  if( workedornot == 0 ){
    cout<<"*************************************"<<endl;
    cout<<" ProcessBadChunks02.C compilation failed! "<<endl;
    cout<<"*************************************"<<endl;
    return;
  }
  
  //Load Class
  gSystem->Load("ProcessBadChunks02_C.so");
  
  
  //Process Datasets, dataset list
  TString lDatasets[] = {
    "LHC11a10b_plus",
    "LHC10a13",
    "LHC10b4",
    "LHC10b5",
    "LHC10d5",
    "LHC11a2h",
    "LHC11a2j",
    "LHC11a6a",
    "LHC11a6b",
    "LHC11b2",
    "LHC11b5",
    "LHC11b6",
    "LHC11b10a",
    "LHC11b10b",
    "LHC11b10c",
    "LHC12a13a",
    "LHC12a13b"
  };
  Long_t lNDatasets = sizeof(lDatasets)/sizeof(lDatasets[0]);
  cout<<"Will process this many datasets: "<<lNDatasets<<endl;
  cout<<"---> Output lists will be stored in the \"output\" directory!"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;

  for(Int_t lInd = 0; lInd < lNDatasets; lInd++){
    cout<<"Process dataset: "<<lDatasets[lInd]<<endl;
    ProcessBadChunks02( lDatasets[lInd] );
  }
}