{
  // change AliFast default parameters

   gAliFast->SetSmearing(1);
   gAliFast->SetTrackFinding(1);
   gAliFast->MCMaker()->Save(2);     //give 2 to save all particles
   gAliFast->TrackMaker()->Save(1);  // give 1 to save tracks
}

