AliAnalysisTaskEMCALTriggerQA * AddTaskEMCALTriggerQA(TString outputFile = ""){
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALTriggerQA", "This task requires an input event handler");
    return NULL;
  }
    
  AliAnalysisTaskEMCALTriggerQA * qatrigger = new AliAnalysisTaskEMCALTriggerQA("QATrigger");
  
  AliEMCALRecoUtils * reco = qatrigger->GetRecoUtils();
  reco->SwitchOnRejectExoticCluster();
  
  // Pass the bad channels, temporary, to fix
  TString fileName="$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root";
  AliOADBContainer *contBC=new AliOADBContainer("");
  contBC->InitFromFile((char*)fileName.Data(),"AliEMCALBadChannels"); 
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(162500); // LHC11e run
  if(arrayBC){
    TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject("pass1");
    if(arrayBCpass){
      
      reco->SwitchOnBadChannelsRemoval();
      printf("trigger REMOVE bad cells \n");
      
      for (Int_t i=0; i<10; ++i) {
        TH2I *hbm = reco->GetEMCALChannelStatusMap(i);
        if (hbm)
          delete hbm;
        hbm=(TH2I*)arrayBCpass->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
        
        if (!hbm) {
          AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
          continue;
        }
        
        hbm->SetDirectory(0);
        reco->SetEMCALChannelStatusMap(i,hbm);
      }
    } else printf("trigger AliEMCALRecoUtils ---Do NOT remove bad channels 1\n");
  }  else  printf("trigger AliEMCALRecoUtils ---Do NOT remove bad channels 2\n");
  
  // Extra channels, temporary, need fix
  
  Int_t badAbsID[]={103, 1263, 1275, 1860, 2117, 2298, 2776, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3675, 3676, 3677, 3678, 3679, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725, 3726, 3727, 3764, 6111, 6800, 6801, 6802, 6803, 6804, 6805, 6806, 6807, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6815, 7430, 7491, 8352, 8353, 8354, 8356, 8357, 8362, 8808, 8810, 8812, 8814, 9217, 9361, 9704, 9769, 9802, 9837, 9839, 9888, 9946, 10117, 11462,
    74, 152, 759, 1059, 1175, 1204, 1288, 1376, 1382, 1386, 1519, 1967, 2026, 2047, 2112, 2114, 2115, 2116, 2118, 2119, 2120, 2123, 2124, 2125, 2350, 2506, 2540, 2793, 2891, 2985, 3135, 3503, 4377, 4817, 5600, 5601, 5602, 5603, 5612, 5613, 5614, 5615, 5648, 5649, 5650, 5651, 5660, 5661, 5662, 5663, 5836, 6104, 6481, 7371, 7375, 7425, 7572, 7874, 9269, 9302, 9389, 9696, 9697, 9698, 9700, 9701, 9702, 9703, 9705, 9706, 9707, 9708, 9709, 9710, 9711, 9750, 9758, 9792, 9793, 9794, 9795, 9798, 9800, 9801, 9803, 9804, 9815, 9824, 9825, 9828, 9829, 9830, 9831, 9832, 9833, 9834, 9835, 9836, 9838, 9872, 9874, 9875, 9878, 9882, 9883, 9889, 9890, 9891, 9892, 9893, 9894, 9896, 9897, 9898, 9899, 9900, 9901, 9902, 9903, 9927, 9937, 9938, 9939, 9940, 9941, 9942, 9943, 9947, 9948, 9949, 9950, 9951, 10112, 10113, 10114, 10115, 10116, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10718, 10723, 10771, 11042, 11091, 11363,
    191, 1968, 2210, 2339, 2391, 6331, 7089, 10126, 10127,
    46, 97, 145, 167, 173, 189, 238, 241, 337, 401, 433, 480, 527, 574, 594, 719, 816, 960, 1008, 1056, 1057, 1102, 1150, 1151, 1152, 1199, 1249, 1423, 1679, 1681, 1822, 1824, 1873, 1918, 2017, 2056, 2110, 2158, 2296, 2332, 2340, 2348, 2351, 2353, 2361, 2381, 2389, 2400, 2434, 2463, 2469, 2481, 2490, 2496, 2519, 2534, 2555, 2590, 2640, 2675, 2687, 2827, 2828, 2830, 2841, 2843, 2869, 2874, 2877, 2878, 2883, 2973, 3022, 3024, 3035, 3039, 3060, 3066, 3070, 3071, 3102, 3109, 3169, 3263, 3275, 3378, 3407, 3412, 4366, 4417, 4558, 4560, 4753, 4797, 8355, 9024, 9025, 9819, 11043};
  
  Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
  
  AliEMCALGeometry* geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  //printf("Total bad %d\n",sizeof(badAbsID)/sizeof(Int_t));
  for(Int_t i=0;i < sizeof(badAbsID)/sizeof(Int_t); i++){
    geom->GetCellIndex(badAbsID[i],iSM,iMod,iIphi,iIeta); 
    // Gives SuperModule and Tower numbers
    geom->GetCellPhiEtaIndexInSModule(iSM,iMod,
                                      iIphi, iIeta,iRow,iCol);
    //printf("bad ID %d, col %d, row %d, sm %d, previous status %d\n",badAbsID[i],iCol,iRow,iSM, reco->GetEMCALChannelStatus(iSM , iCol, iRow));
    reco->SetEMCALChannelStatus(iSM , iCol, iRow,1);
  }
  
  reco->GetEMCALChannelStatus(8 , 27, 101-4*24);
  reco->GetEMCALChannelStatus(9 , 68-48-2, 111-4*24);
  reco->GetEMCALChannelStatus(3 , 74-48, 31-24);
  
  reco->GetEMCALChannelStatus(8 , 23, 98-4*24);
  reco->GetEMCALChannelStatus(8 , 30, 100-4*24);
  reco->GetEMCALChannelStatus(8 , 31, 100-4*24);
  reco->GetEMCALChannelStatus(8 , 31, 96-4*24);
  
  reco->GetEMCALChannelStatus(8 , 37, 103-4*24);
  reco->GetEMCALChannelStatus(8 , 38, 102-4*24);
  reco->GetEMCALChannelStatus(8 , 39, 102-4*24);
  reco->GetEMCALChannelStatus(8 , 37, 101-4*24);
  reco->GetEMCALChannelStatus(8 , 38, 101-4*24);
  reco->GetEMCALChannelStatus(8 , 39, 101-4*24);
  reco->GetEMCALChannelStatus(8 , 37, 99-4*24);
  reco->GetEMCALChannelStatus(8 , 38, 99-4*24);    
  reco->GetEMCALChannelStatus(8 , 39, 99-4*24);
  
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer("EMCALQATrigger", TList::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:EMCALQATrigger",outputFile.Data()));
  mgr->AddTask(qatrigger);
  mgr->ConnectInput  (qatrigger, 0, cinput1);
  mgr->ConnectOutput (qatrigger, 1, coutput);
  
  return qatrigger;
  
}
