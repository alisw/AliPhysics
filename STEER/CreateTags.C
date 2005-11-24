//_________________________________________________________________________
// Macro to create the tag files from an ESD collection.
// The tag files are stored according to the value set in the SetStorage method
// Use Case : 
//          SetStorage(0) --> store the tgs locally
//          SetStorage(1) --> store the tgs in the grid
//          else the program will abort!!!                                    
//                                             
// As a final step the user can create a single merged tag file.         
//_________________________________________________________________________
void CreateTags()
{
  TStopwatch timer;
  timer.Start();

  //connect to AliEn's API services
  TGrid::Connect("alien://");

  //create an AliTagCreator object
  AliTagCreator *t = new AliTagCreator(); 
  
  //Query the file catalog and get a TGridResult
  //TGridResult* result = gGrid->Query("/alice/cern.ch/user/p/peters/analysis/undcent1","*.root","","-l 150 -Opublicaccess=1");  
  TGridResult* result = gGrid->Query("/alice/cern.ch/user/p/pchrista/PDC05/pp/*","AliESDs.root","",""); 
  
  //Store the tag files in AliEn's file catalog
  t->SetStorage(1);
  //Define the SE where the tag files will be stored
  t->SetSE("ALICE::CERN::se01");
  //Define the grid's path where the tag files will be stored
  t->SetGridPath("PDC05/pp/Tags");
  //Read the TGridResult, create the tags and store them
  t->ReadESDCollection(result);
  //Merge the tags and store the merged file
  t->MergeTags();
 
  timer.Stop();
  timer.Print();
}
