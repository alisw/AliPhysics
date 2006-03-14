// $Id$
// $MpId: testDE.C,v 1.1 2006/03/14 09:06:54 ivana Exp $
//
// Test AliMpDEIterator & AliMpSegFactory classes

void testDE() 
{
  AliMpSegFactory factory; 

  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
    cout << "In detection element: " << it.CurrentDE() << endl;

    // Create/get segmentation via factory
    AliMpVSegmentation* segmentation 
      = factory.CreateMpSegmentation(it.CurrentDE(), 0);
      
    // Print number of pads
   cout << "   number of pads: " << segmentation->NofPads() << endl;   
      
    // Find pad by indices in this DE
    AliMpIntPair indices(segmentation->MaxPadIndexX()/2 , 
                         segmentation->MaxPadIndexY()/2);
    AliMpPad pad = segmentation->PadByIndices(indices, false);
    
    cout << "   found pad: " << pad << endl << endl;
  }

  // Delete all created segmentations
  factory.DeleteSegmentations();
}
