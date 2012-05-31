//-*- Mode: C++ -*-
// $Id: HOMERManager.C  $
/**
 * @file HOMERManager.C
 * @brief Sample macro for the use of the HOMERManager
 *
 * Usage:
 * <pre>
 *   aliroot -l HOMERManager.C 
 * </pre>
 *
 * This macro illustrates the usage of the AliHLTHOMERManager in order
 * to on-line read events outside the HLT.
 *
 * This macro can be run inside the
 * <ul>
 *    <li>DAQ/DCS(ACR) network</li>
 *    <li>HLT network</li>
 *    <li>CERN GPN network</li>
 * </ul>
 * 
 * See AliHLTHOMERManager for detailed description.
 *
 * @author Jochen Thaeder 
 * @ingroup alihlt_tutorial
 * @ingroup alihlt_homer
 */


Int_t HOMERManager() {
  
  Int_t iResult = 0;

  // -- Create new hM object
  AliHLTHOMERManager *hM = new AliHLTHOMERManager();

  printf( "== INITIALIZE ==\n" );

  iResult = hM->Initialize();
  if (iResult) return iResult;

  printf( "== CREATE SOURCE LIST ==\n" );

  iResult = hM->CreateSourcesList();
  if (iResult) return iResult; 

  printf( "== CONNECT HOMER ==\n" );
  
  iResult = hM->ConnectHOMER(); 
  if (iResult) return iResult;
  
  printf( "== NEXT EVENT ==\n" );
  
  iResult = hM->NextEvent();
  if (iResult) return iResult;
  
  printf( "== LOOP OVER BLOCKS ==\n" );

  if ( hM->GetBlockList() == NULL )
    return -1;

  if ( hM->GetBlockList()->IsEmpty() )
    return -2;

  TIter next(hM->GetBlockList());
  AliHLTHOMERBlockDesc* block = NULL;

  AliESDEvent* esd = NULL;
  
  while ((block = (AliHLTHOMERBlockDesc*)next())) {
   
    //printf ( "Detector : %s\n" ,block->GetDetector().Data() );
    //printf ( "Datatype : %s\n" ,block->GetDataType().Data() );
    
    if ( block->IsTObject() ) { 
      TObject*     object = block->GetTObject();
      // printf("ClassName %s\n", block->GetClassName().Data() );
    }

    // ++ HLT BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( ! block->GetDetector().CompareTo("HLT") ) {

      // -- ESDTREE
      if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
	esd = static_cast<AliESDEvent*> ( block->GetTObject() );
	esd->GetStdContent();
	
	printf( "Number of ESD Tracks : %d \n", esd->GetNumberOfTracks());
      } //if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
    } // if ( ! block->GetDetector().CompareTo("HLT") ) {

    else if ( ! block->GetDetector().CompareTo("TPC") ) {
      // -- Process TPC Clusters
      if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {

	Int_t   slice = block->GetSubDetector();
	Int_t   patch = block->GetSubSubDetector();
	//	Float_t phi   = ( slice + 0.5 ) * TMath::Pi() / 9.0;  
	//	Float_t cos   = TMath::Cos( phi );
	//	Float_t sin   = TMath::Sin( phi );

	Int_t ddl = slice*6 + patch;

	// cout << slice << " -- " << patch << " -- " << ddl << endl;
	
	//	AliHLTTPCClusterData *cd = reinterpret_cast<AliHLTTPCClusterData*> (block->GetData());
	//	UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);

      } // if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {
    } //     else if ( ! block->GetDetector().CompareTo("TPC") ) {
    
  } // while ((block = (AliHLTHOMERBlockDesc*)next())) {
  
  // -- Destroy hM object
  if (hM)
    delete hM;
  hM = NULL;
  
  return iResult;
}
