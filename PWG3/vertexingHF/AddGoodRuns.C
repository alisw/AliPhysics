Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,TString mcprod="") {
  //
  // Adds good runs from the Monalisa Run Condition Table
  //
  if(mcprod=="") plugin->SetRunPrefix("000"); // DATA

  Int_t firstrun=0,lastrun=9999999;
  Int_t nruns=0,ngoodruns=0;

  if(mcprod=="LHC10d3") {firstrun=117054;lastrun=117222;}
  if(mcprod=="LHC10d5") {firstrun=117086;lastrun=117222;}


  if(lhcPeriod=="LHC10b") {
    nruns=31;
    Int_t runlist[31]={117222, 117220, 117116, 117112, 117109, 117099, 117092, 117086, 117077, 117065, 117063, 117060, 117059, 117054, 117053, 117052, 117050, 117048, 116645, 116643, 116574, 116571, 116562, 116403, 116288, 116102, 115401, 115393, 115193, 115186, 114931};
   
    for(Int_t k=0;k<nruns;k++){
      if(runlist[k]<firstrun || runlist[k]>lastrun) continue;
      plugin->AddRunNumber(runlist[k]);
      ngoodruns++;
    }
    plugin->SetNrunsPerMaster(ngoodruns);
  }

  if(lhcPeriod=="LHC10c") { 
    nruns=36;
    Int_t runlist[36]={120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244, 120079, 120076, 120073, 120072, 120069, 120067, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119163, 119161, 119159};
   
    for(Int_t k=0;k<nruns;k++){
      if(runlist[k]<firstrun || runlist[k]>lastrun) continue;
      plugin->AddRunNumber(runlist[k]);
      ngoodruns++;
    }
    plugin->SetNrunsPerMaster(ngoodruns);
  }

  if(lhcPeriod=="LHC10dhighmu") { // only runs with high mu
    nruns=17;
    Int_t runlist[17]={124750, 124746, 124702, 124608, 124607, 124606, 124605, 124604, 124381, 124380, 124378, 124367, 124362, 124358, 124355, 124191, 124187};
   
    for(Int_t k=0;k<nruns;k++){
      if(runlist[k]<firstrun || runlist[k]>lastrun) continue;
      plugin->AddRunNumber(runlist[k]);
      ngoodruns++;
    }
    plugin->SetNrunsPerMaster(ngoodruns);
  }

  if(lhcPeriod=="LHC10d") { // runs with high mu excluded
    nruns=55;
    Int_t runlist[55]={126437, 126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168, 126167, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843, 125842, 125633, 125632, 125630, 125296, 125134, 125101, 125100, 125097, 125085, 125023, 124751, 122375, 122374};
   
    for(Int_t k=0;k<nruns;k++){
      if(runlist[k]<firstrun || runlist[k]>lastrun) continue;
      plugin->AddRunNumber(runlist[k]);
      ngoodruns++;
    }
    plugin->SetNrunsPerMaster(ngoodruns);
  }

  if(lhcPeriod=="LHC10h") {
    nruns=109;
    Int_t runlist[109]={139510, 139507, 139505, 139504, 139503, 139465, 139440, 139439, 139438, 139437, 139360, 139329, 139328, 139314, 139311, 139310, 139309, 139308, 139173, 139172, 139107, 139105, 139104, 139042, 139038, 139037, 139036, 139029, 139028, 138980, 138979, 138978, 138977, 138872, 138870, 138837, 138830, 138740, 138732, 138731, 138730, 138666, 138662, 138653, 138652, 138638, 138637, 138624, 138621, 138583, 138582, 138578, 138534, 138533, 138469, 138442, 138439, 138438, 138396, 138359, 138275, 138225, 138201, 138200, 138197, 138192, 138190, 138154, 138150, 138126, 138125, 137848, 137843, 137752, 137751, 137748, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137230};  
   
    for(Int_t k=0;k<nruns;k++){
      if(runlist[k]<firstrun || runlist[k]>lastrun) continue;
      plugin->AddRunNumber(runlist[k]);
      ngoodruns++;
    }
    plugin->SetNrunsPerMaster(ngoodruns);
  }



  return ngoodruns;
}
