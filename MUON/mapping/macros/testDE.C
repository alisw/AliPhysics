// $Id$
// $MpId: testDE.C,v 1.1 2006/03/14 09:06:54 ivana Exp $
//
// Test AliMpDEIterator & AliMpSegFactory classes

void testDE() 
{
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
    cout << "In detection element: " << it.CurrentDEId() << endl;

    // Create/get segmentation via factory
    AliMpVSegmentation* segmentation 
      = AliMpSegmentation::Instance()
          ->GetMpSegmentation(it.CurrentDEId(), AliMp::kCath0);
      
    // Print number of pads
   cout << "   number of pads: " << segmentation->NofPads() << endl;   
      
    // Find pad by indices in this DE
    AliMpIntPair indices(segmentation->MaxPadIndexX()/2 , 
                         segmentation->MaxPadIndexY()/2);
    AliMpPad pad = segmentation->PadByIndices(indices, false);
    
    cout << "   found pad: " << pad << endl << endl;
  }
}
