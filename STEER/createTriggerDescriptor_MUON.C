createTriggerDescriptor_MUON()
{
  // create Trigger Descriptor for MUON standalone
  // macro inspired from STEER/createTriggerDescriptor_*.C
  // May 2006, Ch. Finck

   AliTriggerDescriptor descrip( "MUON", "Standalone Trigger for MUON" );

   // Define a Cluster Detector
   descrip.AddDetectorCluster( "MUON" ); // only MUON

   // 1
   descrip.AddCondition( "MUON_Single_LPt_L0", 
                         "MUON_Single_LPt_L0",    "Muon Single Low Pt",
                         (ULong64_t)0x1 );

   // 10
   descrip.AddCondition( "MUON_Single_HPt_L0", 
                         "MUON_Single_HPt_L0",    "Muon Single High Pt",
                         (ULong64_t)0x1 << 1 );

   // 100
   descrip.AddCondition( "MUON_Like_LPt_L0", 
                         "MUON_Like_LPt_L0",    "Di Muon Like sign Low Pt",
                         (ULong64_t)0x1 << 2 );

   // 1000
   descrip.AddCondition( "MUON_Like_HPt_L0", 
                         "MUON_Like_HPt_L0",    "Di Muon Like sign High Pt",
                         (ULong64_t)0x1 << 3 );

   // 10000
   descrip.AddCondition( "MUON_Unlike_LPt_L0", 
                         "MUON_Unlike_LPt_L0",    "Di Muon Unlike sign Low Pt",
                         (ULong64_t)0x1 << 4 );

   // 100000
   descrip.AddCondition( "MUON_Unlike_HPt_L0", 
                         "MUON_Unlike_HPt_L0",    "Di Muon Unlike sign High Pt",
                         (ULong64_t)0x1 << 5 );



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


