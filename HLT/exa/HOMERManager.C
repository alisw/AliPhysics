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
  
  // Empty argument should connect to all
  //iResult = hM->ConnectHOMER();
    iResult = hM->ConnectHOMER("HLT");
  if (iResult) return iResult;
  
  printf( "== NEXT EVENT ==\n" );
  
  iResult = hM->NextEvent();
  if (iResult) return iResult;

  printf( "== LOOP OVER BLOCKS ==\n" );

  TObject * object =  NULL;  

  TIter next(hM->GetBlockList());
  AliHLTHOMERBlockDesc* block = NULL;

  while ((block = (AliHLTHOMERBlockDesc*)next())) {
   
    printf ( "Detector : %s\n" ,block->GetDetector().Data() );
    printf ( "Datatype : %s\n" ,block->GetDataType().Data() );
    
    if ( block->IsTObject() ) {
      object = block->GetTObject();
      
        printf("ClassName %s\n", block->GetClassName().Data() );
    }
  }

  // -- Destroy hM object
  if (hM)
    delete hM;
  hM = NULL;

  return iResult;;
}
