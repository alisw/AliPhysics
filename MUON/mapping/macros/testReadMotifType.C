// $Id$
//
// Test macro for reading motif type data.

void testReadMotifType(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane)
{
  AliMpReader r(station, plane);
  //r.SetVerboseLevel(2);

  TString names="ABCDEFGHIJKLMN";
  for (Int_t i=0;i<names.Length();++i){
     r.BuildMotifType(names[i])->Print("G");
  }
}
