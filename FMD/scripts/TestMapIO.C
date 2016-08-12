//____________________________________________________________________
//
// $Id$
//
// Test I/O of ALiFMDMap
//
/** @defgroup FMD_MAPIO_TEST Map I/O test
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** @ingroup FMD_MAPIO_test
 */
void 
WriteTree()
{
  TFile* file = TFile::Open("map.root", "RECREATE");
  TTree* tree = new TTree("T", "T");
  AliFMDFloatMap* m = new AliFMDFloatMap(1, 1, 1, 3);
  tree->Branch("map", "AliFMDFloatMap", &m);
  for (int i = 0; i < 3; i++) m->operator()(1,'I',0,i) = i + 1;
  tree->Fill();
  file->Write();
  file->Close();
}

//____________________________________________________________________
/** @ingroup FMD_MAPIO_test
 */
void 
ReadTree()
{
  TFile* file = TFile::Open("map.root", "READ");
  TTree* tree = static_cast<TTree*>(file->Get("T"));
  AliFMDFloatMap* m = 0;
  tree->SetBranchAddress("map", &m);
  tree->GetEntry(0);
  for (int i = 0; i < 3; i++) {
    std::cout << "Map(1,'I',0," << i << "): " << m->operator()(1,'I',0,i)
	      << std::endl;
  }
  file->Close();
}

  
//____________________________________________________________________
/** @ingroup FMD_MAPIO_test
 */
void
WriteMap() 
{
  TFile* file = TFile::Open("map.root", "RECREATE");
  AliFMDFloatMap* m = new AliFMDFloatMap(1, 1, 1, 3);
  for (int i = 0; i < 3; i++) m->operator()(1,'I',0,i) = i + 1;
  m.Write("map");
  file->Close();
}

//____________________________________________________________________
/** @ingroup FMD_MAPIO_test
    @return  */
void
ReadMap() 
{
  TFile* file = TFile::Open("map.root", "READ");
  AliFMDFloatMap* m = static_cast<AliFMDFloatMap*>(file->Get("map"));
  std::cout << "Got map " << map << std::endl;
  for (int i = 0; i < 3; i++) {
    std::cout << "Map(1,'I',0," << i << "): " << m->operator()(1,'I',0,i)
	      << std::endl;
  }
  file->Close();
}


//____________________________________________________________________
/** @ingroup FMD_MAPIO_test
 */
void
TestMapIO()
{
  WriteMap();
  ReadMap();
  WriteTree();
  ReadTree();
}

//____________________________________________________________________
//
// EOF
//
