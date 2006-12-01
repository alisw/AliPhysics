//-------------------------------------------------------------------------
// Setting the cuts for the V0 and cascade finding
// The values of the cuts below are "reasonable" for central PbPb events
//------------------------------------------------------------------------- 
{
  AliReconstruction rec;


  Double_t cuts[]={33,  // max. allowed chi2
                  0.16, // min. allowed negative daughter's impact parameter 
                  0.05, // min. allowed positive daughter's impact parameter 
                  0.08, // max. allowed DCA between the daughter tracks
                  0.99, // max. allowed cosine of V0's pointing angle
                  0.9,  // min. radius of the fiducial volume
                  2.9   // max. radius of the fiducial volume
  };
  AliV0vertexer::SetDefaultCuts(cuts);

  Double_t cts[]={33.,  // max. allowed chi2
                 0.05,  // min. allowed V0 impact parameter 
                 0.008, // window around the Lambda mass 
                 0.035, // min. allowed bachelor's impact parameter 
                 0.10,  // max. allowed DCA between a V0 and a track
                 0.9985,//max. allowed cosine of the cascade pointing angle
                 0.9,   // min. radius of the fiducial volume
                 2.9    // max. radius of the fiducial volume
  };
  AliCascadeVertexer::SetDefaultCuts(cts);


  rec.Run();
}
