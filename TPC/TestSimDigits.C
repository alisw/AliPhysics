/// \file TestSimDigits.C

void TestSimDigits() {
AliSimDigits dig;
dig.Allocate(10,10);
dig.AllocateTrack(3);

dig.SetTrackIDFast(550,5,5,0);
dig.SetTrackIDFast(551,5,5,1);
dig.SetTrackIDFast(552,5,5,2);


dig.SetTrackIDFast(550,6,5,0);
dig.SetTrackIDFast(551,6,5,1);
dig.SetTrackIDFast(552,6,5,2);

dig.SetTrackIDFast(550,7,5,0);
dig.SetTrackIDFast(551,7,5,1);
dig.SetTrackIDFast(552,7,5,2);


dig.SetTrackIDFast(440,4,4,0);
dig.SetTrackIDFast(441,4,4,1);
dig.SetTrackIDFast(442,4,4,2);




dig.CompresTrackBuffer(1);

dig.GetTrackID(5,5,0);
dig.ExpandTrackBuffer() ;
dig.GetTrackIDFast(5,5,0);
dig.GetTrackIDFast(5,5,0);
dig.GetTrackIDFast(5,5,1);
dig.GetTrackIDFast(5,5,2);
}
