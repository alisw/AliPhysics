createTriggerDescriptor_MUON()
{
  // create Trigger Descriptor for MUON standalone
  // macro inspired from STEER/createTriggerDescriptor_*.C
  // May 2006, Ch. Finck

   AliTriggerDescriptor descrip( "MUON", "Standalone Trigger for MUON" );

   // Define a Cluster Detector
   descrip.AddDetectorCluster( "MUON" ); // only MUON

   // 1
   descrip.AddCondition( "MUON_SPlus_LPt_L0", 
                         "MUON_SPlus_LPt_L0",    "Muon Plus Low Pt",
                         (ULong64_t)0x1 << 1 );

   // 2
   descrip.AddCondition( "MUON_SPlus_HPt_L0", 
                         "MUON_SPlus_HPt_L0",    "Muon Plus High Pt",
                         (ULong64_t)0x1 << 2 );

   // 3
   descrip.AddCondition( "MUON_SPlus_All_L0", 
                         "MUON_SPlus_All_L0",    "Muon Plus All Pt",
                         (ULong64_t)0x1 << 3 );

   // 4
   descrip.AddCondition( "MUON_SMinus_LPt_L0", 
                         "MUON_SMinus_LPt_L0",    "Muon Minus Low Pt",
                         (ULong64_t)0x1 << 4 );

   // 5
   descrip.AddCondition( "MUON_SMinus_HPt_L0", 
                         "MUON_SMinus_HPt_L0",    "Muon Minus High Pt",
                         (ULong64_t)0x1 << 5 );

   // 6
   descrip.AddCondition( "MUON_SUndef_LPt_L0", 
                         "MUON_SUndef_LPt_L0",    "Muon Undefined sign Low Pt",
                         (ULong64_t)0x1 << 6 );

   // 7
   descrip.AddCondition( "MUON_SUndef_HPt_L0", 
                         "MUON_SUndef_HPt_L0",    "Muon Undefined sign High Pt",
                         (ULong64_t)0x1 << 7 );

   // 8
   descrip.AddCondition( "MUON_SUndef_All_L0", 
                         "MUON_SUndef_All_L0",    "Muon Undefined sign All Pt",
                         (ULong64_t)0x1 << 8 );

   // 9
   descrip.AddCondition( "MUON_Unlike_LPt_L0", 
                         "MUON_Unlike_LPt_L0",    "Di Muon Unlike sign Low Pt",
                         (ULong64_t)0x1 << 9 );

   // 10
   descrip.AddCondition( "MUON_Unlike_HPt_L0", 
                         "MUON_Unlike_HPt_L0",    "Di Muon Unlike sign High Pt",
                         (ULong64_t)0x1 << 10 );

   // 11
   descrip.AddCondition( "MUON_Unlike_All_L0", 
                         "MUON_Unlike_All_L0",    "Di Muon Unlike sign All Pt",
                         (ULong64_t)0x1 << 11 );

   // 12
   descrip.AddCondition( "MUON_Like_LPt_L0", 
                         "MUON_Like_LPt_L0",    "Di Muon Like sign Low Pt",
                         (ULong64_t)0x1 << 12 );

   // 13
   descrip.AddCondition( "MUON_Like_HPt_L0", 
                         "MUON_Like_HPt_L0",    "Di Muon Like sign High Pt",
                         (ULong64_t)0x1 << 13 );

   // 14
   descrip.AddCondition( "MUON_Like_All_L0", 
                         "MUON_Like_All_L0",    "Di Muon Like sign All Pt",
                         (ULong64_t)0x1 << 14 );
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

