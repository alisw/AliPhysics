// $Id$
//
// Test macro for reading motif type data.

void testReadMotifType(AliMpStationType station = kStation1,
                       AliMpPlaneType plane = kBendingPlane)
{
  AliMpReader r(station, plane);
  //r.SetVerboseLevel(2);

  TString names;
  TString names2;
  Int_t nv =0;
  if ( station == kStation1 )
    if ( plane == kBendingPlane ) 
      names ="ABCDEFGHI";
    else
      names = "ABCDEFGHIJKLMN";
  else if ( station == kStation2 ) 
    if ( plane == kBendingPlane ) {
      names ="ABCDEFGHIJKLMNOPQRSTUVWXY";
      names2 ="abcdefghimnptuv";
      nv = 5;
    }  
    else {
      names = "ABCDEFGHIIJKLMN";
      names2 ="abcdefgijklmnopqrstuwv";
      nv = 6;
    }  
    
  for (Int_t i=0;i<names.Length();++i){
     r.BuildMotifType(names[i])->Print("G");
  }
  
  // motifs a1, b1, ..., u1, v1
  for (Int_t i2=0;i2<names2.Length();++i2){
    TString mtName = names2[i2];
    mtName += "1";
    r.BuildMotifType(mtName)->Print("G");
  }
  
  // motifs v2, ..., v5, v6
  TString names4="v";
  for (Int_t i3=2;i3<nv+1;++i3) { 
      TString mtName = "v";
      mtName += i3;
      r.BuildMotifType(mtName)->Print("G");
  }
}
