class AliAnalysisAlien;

void AddRunsPbPb(AliAnalysisAlien* plugin,bool onAOD=false) {
  //  plugin->SetGridDataDir("/alice/data/2010/LHC10h");
  plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
  if( onAOD ){
    plugin->SetDataPattern("*ESDs/pass1/AOD033/*/AliAOD.root");
  }else{
    plugin->SetDataPattern("*ESDs/pass2/*/AliESDs.root");
  }

  int FieldMM[41]={
    167807, 167808, 167813, 167902, 167903, 
    167915, 167920, 167988, 168068, 168069, 
    168076, 168103, 168104, 168105, 168107, 
    168108, 168115, 168171, 168172, 168175, 
    168181, 168204, 168205, 168206, 168207, 
    168311, 168342, 168464, 168984, 169045, 
    169144, 169238, 169420, 169506, 169515, 
    169553, 169555, 169584, 169587, 169590, 
    169591 
  };
  int FieldPP[45]={
    169628, 169835, 169837, 169838, 169855, 
    169858, 169859, 169918, 169919, 169920, 
    169922, 169923, 169924, 169956, 169965, 
    169975, 169981, 170027, 170036, 170038, 
    170081, 170083, 170085, 170088, 170152, 
    170155, 170163, 170195, 170204, 170205, 
    170207, 170208, 170228, 170230, 170270, 
    170308, 170309, 170311, 170312, 170315, 
    170390, 170546, 170556, 170572, 170593
  };

  // for(int ir=0; ir!=41; ++ir){
  //     plugin->AddRunNumber( FieldMM[ir] );
  //   }
  for(int ir=0; ir!=45; ++ir){
    plugin->AddRunNumber( FieldPP[ir] );
  }
  
  plugin->SetNrunsPerMaster(10);
  
}

