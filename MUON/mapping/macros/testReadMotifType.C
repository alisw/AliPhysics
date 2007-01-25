// $Id$
// $MpId: testReadMotifType.C,v 1.12 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for reading motif type data.

void testReadMotifType(AliMp::StationType station = AliMp::kStation1,
                       AliMp::PlaneType plane = AliMp::kBendingPlane,
	     	       Bool_t rootInput = false)
{
  AliMpMotifReader r(station, plane);
  //r.SetVerboseLevel(2);

  TString names;
  TString names2;
  Int_t nv =0;
  if ( station == AliMp::kStation1 )
    if ( plane == AliMp::kBendingPlane ) 
      names ="ABCDEFGHI";
    else
      names = "ABCDEFGHIJKLMN";
  else if ( station == AliMp::kStation2 ) 
    if ( plane == AliMp::kBendingPlane ) {
      names ="ABCDEFGHIJKLMNOPQRSTUVWXY";
      names2 ="abcdefghimnptuv";
      nv = 5;
    }  
    else {
      names = "ABCEFGHIJKLMN";
      names2 ="abcdefgijklmnopqrstuwv";
      nv = 5;
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
