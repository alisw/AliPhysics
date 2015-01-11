//_________________________________________________________________________
// Macro to merge the tag files from an ESD collection.
// The tag files are stored according to the value set in the SetStorage method
// Use Case : 
//          SetStorage(0) --> store the tgs locally
//          SetStorage(1) --> store the tgs in the grid
//          else the program will abort!!!                                    
//                                             
// As a final step the user can create a single merged tag file.         
//_________________________________________________________________________
void MergeTags(const char* fCollectionName) {
  //connect to AliEn's API services
  TGrid::Connect("alien://pcapiserv01.cern.ch:10000","pchrist");  
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

  //___________________________________//
  //__________Merge tag files__________//
  //___________________________________//
  //XML collection
  TAlienCollection *collection = TAlienCollection::Open(fCollectionName);
  TGridResult* result = collection->GetGridResult("");
  //Read the TGridResult and merge the tags
  t->MergeTags(result);
}
