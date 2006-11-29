createTriggerDescriptor_PbPb()
{
   // create Trigger Descriptor for Pb-Pb interactions


   AliTriggerDescriptor descrip( "Pb-Pb", "Default Pb-Pb Descriptor" );

   // Define a Cluster Detector
   //descrip.AddDetectorCluster( "ALL" );
   descrip.AddDetectorCluster( "ITS TRD PHOS EMCAL MUON ZDC T0 VZERO" ); // no CRT yet

   // Define the trigger conditions form Table 4.2 TDR DAQ, Trigger pag 59

   //    TRD_pre_L0 input is not generated in simulation as it is a no busy signal, so, it is not 
   //    include in the conditions.

   // 1
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & ZDC_1_L1",// & TRD_pre_L0
                         "MB",                    "Minimum Bias",
                         (ULong64_t)0x1  );
   // 2
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & ZDC_2_L1",// & TRD_pre_L0
                         "SC",                    "Semi Central",
                         (ULong64_t)0x1 << 1 );
   // 3
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & ZDC_3_L1",// & TRD_pre_L0
                         "CE",                    "Central",
                         (ULong64_t)0x1 << 2 );
   // 4
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Unlike_HPt_L0 & ZDC_1_L1",  // & TRD_pre_L0
                         "DMUnlikeHPt_TPC_MB",    "Di Muon Unlike High Pt TPC Minimum Bias",
                         (ULong64_t)0x1 << 3 );
   // 5
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & MUON_Unlike_HPt_L0 & ZDC_2_L1",  // & TRD_pre_L0
                         "DMUnlikeHPt_TPC_SC",    "Di Muon Unlike High Pt TPC Semi Central",
                         (ULong64_t)0x1 << 4 );
   // 6 same as 4
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Unlike_HPt_L0 & ZDC_1_L1",
                         "DMUnlikeHPt_noTPC_MB",  "Di Muon Unlike High Pt no TPC Minimum Bias",
                         (ULong64_t)0x1 << 5 );
   // 7
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Unlike_LPt_L0 & ZDC_1_L1",
                         "DMUnlikeLPt_noTPC_MB",   "Di Muon Unlike Low Pt no TPC Minimum Bias",
                         (ULong64_t)0x1 << 6 );
   // 8
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & MUON_Unlike_LPt_L0 & ZDC_2_L1",
                         "DMUnlikeLPt_noTPC_SC",   "Di Muon Unlike Low Pt no TPC Semi Central",
                         (ULong64_t)0x1 << 7 );
   // 9
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Like_HPt_L0 & ZDC_1_L1", // & TRD_pre_L0
                         "DMLikeHPt_TPC_MB",       "Di Muon Like Low Pt TPC Minimum Bias", 
                         (ULong64_t)0x1 << 8 );
   // 10
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & MUON_Like_HPt_L0 & ZDC_2_L1", // & TRD_pre_L0
                         "DMLikeHPt_TPC_SC",       "Di Muon Like High Pt TPC Semi Central",
                         (ULong64_t)0x1 << 9 );
   // 11 same as 9
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Like_HPt_L0 & ZDC_1_L1",
                         "DMLikeHPt_noTPC_MB",     "Di Muon Like High Pt no TPC Minimum Bias", 
                         (ULong64_t)0x1 << 10 );
   // 12
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Like_LPt_L0 & ZDC_1_L1",
                         "DMLikeLPt_noTPC_MB",     "Di Muon Like Low Pt no TPC Minimum Bias", 
                         (ULong64_t)0x1 << 11 );
   // 13
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & MUON_Like_LPt_L0 & ZDC_2_L1",
                         "DMLikeLPt_noTPC_SC",     "Di Muon Like Low Pt no TPC Semi Central",
                         (ULong64_t)0x1 << 12 );

   // the current implementation of MUON trigger class has to change the inputs to match the following conditions
   // now they have 9 inputs instead of 1 for MUON single
   //          "MUON_SPlus_LPt_L0",  "Single Plus Low Pt",  
   //          "MUON_SPlus_HPt_L0",  "Single Plus High Pt", 
   //          "MUON_SPlus_All_L0",  "Single Plus All",     
   //          "MUON_SMinus_LPt_L0", "Single Minus Low Pt",  
   //          "MUON_SMinus_HPt_L0", "Single Minus High Pt", 
   //          "MUON_SMinus_All_L0", "Single Minus All",     
   //          "MUON_SUndef_LPt_L0", "Single Undefined Low Pt", 
   //          "MUON_SUndef_HPt_L0", "Single Undefined High Pt",
   //          "MUON_SUndef_All_L0", "Single Undefined All",    

   // 14
//   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & MUON_Single_L0 & TRD_Electron_L1 & ZDC_1_L1",  // & TRD_pre_L0
//                         "DMSingleTDRe_MB",       "Di Muon Single TDR electron Minimum Bias",
//                         (ULong64_t)0x1 << 13 );
   // 15
//   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & MUON_Single_L0 & TRD_Electron_L1 & ZDC_2_L1",  // & TRD_pre_L0
//                         "DMSingleTDRe_SC",       "Di Muon Single TDR electron Semi Central",
//                         (ULong64_t)0x1 << 14 );
   // 16
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_Electron_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRe_MB",               "TDR electron Minimum Bias",
                         (ULong64_t)0x1 << 15 );
   // 17
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_HadrLPt_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRLPt_MB",             "TDR Low Pt Minimum Bias",
                         (ULong64_t)0x1 << 16 );
   // 18
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_HadrHPt_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRHPt_MB",             "TDR High Pt Minimum Bias",
                         (ULong64_t)0x1 << 17 );
   // 19
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_Unlike_EPair_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRUnlikeHPt_MB",       "TDR Unlike Electron Pair High Pt Minimum Bias",
                         (ULong64_t)0x1 << 18 );
   // 20
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & TRD_Unlike_EPair_L1 & ZDC_2_L1",  // & TRD_pre_L0
                         "TDRUnlikeHPt_SC",       "TDR Unlike Electron Pair High Pt Semi Central",
                         (ULong64_t)0x1 << 19 );
   // 21
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_Like_EPair_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRLikeHPt_MB",         "TDR Like Electron Pair High Pt Minimum Bias",
                         (ULong64_t)0x1 << 20 );
   // 22
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & TRD_Like_EPair_L1 & ZDC_2_L1",  // & TRD_pre_L0
                         "TDRLikeHPt_SC",         "TDR Like Electron Pair High Pt Semi Central",
                         (ULong64_t)0x1 << 21 );

   // 23
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & TRD_Jet_HPt_L1 & ZDC_2_L1",  // & TRD_pre_L0
                         "TDRJetHPt_SC",          "TDR Jet High Pt Semi Central",
                         (ULong64_t)0x1 << 22 );
   // 24
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & TRD_Jet_LPt_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "TDRJetLPt_MB",          "TDR Jet Low Pt Minimum Bias",
                         (ULong64_t)0x1 << 23 );
   // 25
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & TRD_Jet_LPt_L1 & ZDC_2_L1",  // & TRD_pre_L0
                         "TDRJetLPt_SC",          "TDR Jet Low Pt Semi Central",
                         (ULong64_t)0x1 << 24 );

   // 26
   descrip.AddCondition( "(START_Vertex_L0 & VZERO_AND) & PHOS_JetHPt_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "PHOSHPt_MB",            "PHOS High Pt Minimum Bias",
                         (ULong64_t)0x1 << 25 );
   // 27
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & PHOS_JetLPt_L1 & ZDC_1_L1",  // & TRD_pre_L0
                         "PHOSLPt_MB",            "PHOS Low Pt Minimum Bias",
                         (ULong64_t)0x1 << 26 );
   // 28
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & PHOS_JetLPt_L1 & ZDC_2_L1",  // & TRD_pre_L0
                         "PHOSLPt_SC",            "PHOS Low Pt Semi Central",
                         (ULong64_t)0x1 << 27 );
   // 29
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & PHOS_L0 & ZDC_1_L1",
                         "PHOSStand-along",       "PHOS Stand-along Minimum Bias",
                         (ULong64_t)0x1 << 28 );
   // 30
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & EMCAL_JetHPt_L1 & ZDC_1_L1",
                         "EMCALJetHPt_MB",        "EMCAL Jet High Pt Minimum Bias",
                         (ULong64_t)0x1 << 29 );
   // 31
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & EMCAL_JetMPt_L1 & ZDC_1_L1",
                         "EMCALJetMPt_MB",        "EMCAL Jet Medium Pt Minimum Bias",
                         (ULong64_t)0x1 << 30 );
   // 32
   descrip.AddCondition( "START_Vertex_L0 & VZERO_AND & EMCAL_JetLPt_L1 & ZDC_1_L1",
                         "EMCALJetLPt_MB",        "EMCAL Jet Low Pt Minimum Bias",
                         (ULong64_t)0x1 << 31 );
   // 33
   descrip.AddCondition( "START_Vertex_L0 & VZERO_OR & EMCAL_JetLPt_L1 & ZDC_2_L1",
                         "EMCALJetLPt_SC",        "EMCAL Jet Low Pt Semi Central",
                         (ULong64_t)0x1 << 32 );
   // 34
   descrip.AddCondition( "ZDC_EMD_L1", // CRT_L0
                         "ZDC_diss",              "ZDC EMD Event",
                         (ULong64_t)0x1 << 33 ); 
   // 35
//   descrip.AddCondition( "CRT_cosmic_L0",
//                         "CRT_cosmic",            "CRT cosmic telescope",
//                         (ULong64_t)0x1 << 34 );

   // 36
   descrip.AddCondition( "!START_Vertex_L0",
                         "Beam-gas",              "Beam Gas Event",
                         (ULong64_t)0x1 << 35 );


   cout << endl << endl;

   if( !descrip.CheckInputsConditions("Config.C") ) {
      cerr << "\n ERROR: There are some problems on descriptor definition..." << endl;
      return;
   }

   cout << endl << endl;

   cout << "************************************************************" << endl;
   cout << "New Trigger descriptor" << endl;
   cout << "************************************************************" << endl;
   descrip.Print();

   // save the descriptor to file
   descrip.WriteDescriptor();

   cout << endl << endl << endl;
   cout << "************************************************************" << endl;
   cout << "Available Trigger descriptors" << endl;
   cout << "************************************************************" << endl;

   // Get and print all available descriptors
   TObjArray* desc = AliTriggerDescriptor::GetAvailableDescriptors();

   if( !desc || !desc->GetEntriesFast() ) {
      cerr << "Not descriptors availables" << endl;
      return;
   }

   Int_t ndesc = desc->GetEntriesFast();
   for( Int_t j=0; j<ndesc; j++ ) {  
      AliTriggerDescriptor* de = (AliTriggerDescriptor*)desc->At(j); 
      de->Print();
   }

}

