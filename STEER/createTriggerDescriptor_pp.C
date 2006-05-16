createTriggerDescriptor_pp()
{
   // create Trigger Descriptor for p-p interactions


   AliTriggerDescriptor descrip( "p-p", "Default p-p Descriptor" );

   // Define a Cluster Detector
   //descrip.AddDetectorCluster( "ALL" );
   descrip.AddDetectorCluster( "ITS START VZERO TOF" ); // no CRT yet

   descrip.AddCondition( "VZERO_LEFT",      "VZERO_LEFT",      "VZERO A (Left)",            (ULong64_t)0x1  );
   descrip.AddCondition( "VZERO_RIGHT",     "VZERO_RIGHT",     "VZERO C (Right)",           (ULong64_t)0x1 << 1 );
   descrip.AddCondition( "VZERO_BEAMGAS",   "VZERO_BEAMGAS",   "VZERO beam gas rejection",  (ULong64_t)0x1 << 2 );
   descrip.AddCondition( "START_A_L0",      "START_A_L0",      "START A (Left)",            (ULong64_t)0x1 << 3 );
   descrip.AddCondition( "START_C_L0",      "START_C_L0",      "START C (Right)",           (ULong64_t)0x1 << 4 );
   descrip.AddCondition( "ITS_SPD_GFO_L0",  "ITS_SPD_GFO_L0",  "SPD global fast-or",        (ULong64_t)0x1 << 5 );
   descrip.AddCondition( "ITS_SPD_HMULT_L0","ITS_SPD_HMULT_L0","SPD high mult. 100 ",       (ULong64_t)0x1 << 6 );

   descrip.AddCondition( "ITS_SPD_GFO_L0 & VZERO_AND",
                         "MB",
                         "Minimum Bias",
                         (ULong64_t)0x1 << 7 );
   descrip.AddCondition( "ITS_SPD_GFO_L0 & VZERO_AND | TOF_pp_MB_L0",
                         "MB-TOF",
                         "Minimum Bias with TOF",
                         (ULong64_t)0x1 << 8 );

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

