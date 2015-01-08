ReadTrigger( TString inFile = "galice.root" )
{
   // Dynamically link some shared libs
   if( gClassTable->GetID( "AliRun" ) < 0 ) {
      gROOT->ProcessLine( ".x $(ALICE_ROOT)/macros/loadlibs.C" );
   } else if( gAlice ) {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
   }

   AliRunLoader* rl = AliRunLoader::Open( inFile.Data() );
   if( rl == 0x0 ) {
      cerr << "ReadTrigger.C : Can not open session RunLoader=NULL"
           << endl;
       return 3;
   }

   // Read and Print Trigger

   rl->LoadTrigger();
   AliCentralTrigger *aCTP = rl->GetTrigger();
   aCTP->Print();

   // Loop over event and print trigger info
   Int_t nevent = rl->GetNumberOfEvents();
   for( Int_t i=0; i<nevent; i++ ) {
      rl->GetEvent( i );
      cout << endl << "Event " << i
           << " Global Trigger Class Mask: 0x" << hex << aCTP->GetClassMask() << endl;

      // Read trigger inputs from detector. Example. ITS
      AliLoader * loader = rl->GetDetectorLoader( "ITS" );
      if( loader ) {
         AliDataLoader * dataLoader = loader->GetDigitsDataLoader();
         if( !dataLoader->IsFileOpen() ) 
            dataLoader->OpenFile( "READ" );
         AliTriggerDetector* trgdet = (AliTriggerDetector*)dataLoader->GetDirectory()->Get( "Trigger" );
         if( trgdet ) {
            trgdet->Print();
         } else {
            cerr << "There is not trigger object for " << loader->GetName() << endl;
         }
      }
   }
}
