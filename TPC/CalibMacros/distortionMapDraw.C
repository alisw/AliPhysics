/// \file distortionMapDraw.C
///
/// This is an example macro for drawing of the distortion maps which were created during the CPass0
/// respectivally,  CPass1
///
/// Macro to be extended - resulting QA plots to be published on the web page

void distortionMapDraw(){
  // 1.) connect the grid

  TGrid::Connect("alien");
  //
  // 2.) Open calibration of interest. e.g
  //
  TFile *f  = TFile::Open("alien:///alice/data/2011/LHC11h/000170593/cpass0_HLT/OCDB/meanITSVertex.root");
  //
  // 3.) Get the content - distortion map
  //   TPC=Vertex, TPC-ITS, TPC-TRD, TPC-TOF
  //   f.ls()
  //   TAlienFile**            alien:///alice/data/2011/LHC11h/000170593/cpass0_HLT/OCDB/meanITSVertex.root             //   TAlienFile*            alien:///alice/data/2011/LHC11h/000170593/cpass0_HLT/OCDB/meanITSVertex.root
  //   KEY: TTree    ITSdy;1 ITSdy
  //   KEY: TTree    ITSdz;1 ITSdz
  //   KEY: TTree    ITSdsnp;1       ITSdsnp
  //   KEY: TTree    Vertexdy;1      Vertexdy
  //   KEY: TTree    Vertexdsnp;1    Vertexdsnp
  //   KEY: TTree    Vertexdz;1      Vertexdz
  //   KEY: TTree    TOFdy;1 TOFdy
  //   KEY: TTree    TRDdy;1 TRDdy
  //
  // 4.) To see the variables of the distrotion tree
  //     ITSdy.Print()
  //
  // 5. Example draw - Note that  for the example run we show residuals obtained in  CPass0 
  //    before run dependent alignment and ExB twist 
  // 
  ITSdy->SetMarkerStyle(25);
  Vertexdy->SetMarkerStyle(25);
  TOFdy->SetMarkerStyle(25);
  // Draw mean residual 
  ITSdy->Draw("mean","entries>100&&theta<0&&abs(snp)<0.1","");
  //
  // Draw mean residual/distortion as function of sector position - theta in color code 
  //
  ITSdy->Draw("mean:sector:abs(theta)","entries>100&&theta<0","colz");
  Vertexdy->Draw("mean:sector:abs(theta)","entries>100&&theta<0","colz")
  TOFdy->Draw("mean:sector:abs(theta)","entries>100&&theta<0","colz");
}
