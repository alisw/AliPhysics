//_________________________________________________________________________
// Macro to create the tag files from:
// i)  AliESDs that are stored locally.
// ii) AliESDs that are stored in the CAF.
// iii)AliESDs that are stored in alien.
// The tag files are stored according to the value set in the SetStorage method
// Use Case : 
//          SetStorage(0) --> store the tags locally
//          SetStorage(1) --> store the tags in the grid
//          else the program will abort!!!                                    
//                                             
// As a final step the user can create a single merged tag file.         
//_________________________________________________________________________
void CreateTags() {
  TStopwatch timer;
  timer.Start();

  //___________________________________//
  //__Connect to AliEn's API services__//
  //___________________________________//
  //TGrid::Connect("alien://pcapiserv01.cern.ch:10000","pchrista"); 
  //TGrid::Connect("alien://");
    
  //___________________________________//
  //__Create an AliTagCreator object___//
  //___________________________________//
  AliTagCreator *t = new AliTagCreator(); 
    
  //___________________________________//
  //_____ Storage of the tag files:____//
  //________0-->local, 1-->alien_______//
  //___________________________________//
  t->SetStorage(0);
  //Define the grid's path where the tag files will be stored
  //t->SetGridPath("Tags/output");
  //Define the SE where the tag files will be stored
  //t->SetSE("ALICE::CERN::se01");

  //___________________________________//
  //________Locally stored ESDs________//
  //___________________________________//
  t->ReadLocalCollection("/home/pchrist/ALICE/Alien/Local/Tags");

  //___________________________________//
  //__________CAF stored ESDs__________//
  //___________________________________//
  //t->ReadCAFCollection("ESD1.txt");
  
  //___________________________________//
  //_________Alien stored ESDs_________//
  //___________________________________//
  //XML collection
  //TAlienCollection *collection = TAlienCollection::Open("pp.xml");
  //TGridResult* result = collection->GetGridResult("");
  //Read the TGridResult, create the tags and store them
  //t->ReadGridCollection(result);

  //___________________________________//
  //___________Merge the tags__________//
  //___________________________________//
  //Merge the tags and store the merged file
  t->MergeTags();

  timer.Stop();
  timer.Print();
}

