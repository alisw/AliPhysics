void Config(Int_t version)
{
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libRICH");

  AliRICH *RICH = 0;
  switch (version) {
    case 0: RICH  = new AliRICHv0("RICH","normal RICH"); break;
  }  

//
// Version 0
// Default Segmentation
    AliRICHSegmentationV0* SegmentationV0 = new AliRICHSegmentationV0;
//
//  Segmentation parameters
    SegmentationV0->SetPadSize(0.84,0.80);
    SegmentationV0->SetDAnod(0.84/2);
//
//  Geometry parameters
    AliRICHGeometry* GeometryV0 = new AliRICHGeometryV0;
    GeometryV0->SetGapThickness(7.6);
    GeometryV0->SetProximityGapThickness(.4);
    GeometryV0->SetQuartzLength(131);
    GeometryV0->SetQuartzWidth(126.2);
    GeometryV0->SetQuartzThickness(.5);
    GeometryV0->SetOuterFreonLength(131);
    GeometryV0->SetOuterFreonWidth(40.3);
    GeometryV0->SetInnerFreonLength(131);
    GeometryV0->SetInnerFreonWidth(40.3);
    GeometryV0->SetFreonThickness(1);
//
//  Response parameters
    AliRICHResponseV0*  Rresponse0   = new AliRICHResponseV0;
    Rresponse0->SetSigmaIntegration(5.);
    Rresponse0->SetChargeSlope(41.);
    Rresponse0->SetChargeSpread(0.18, 0.18);
    Rresponse0->SetMaxAdc(1024);
    Rresponse0->SetAlphaFeedback(0.05);
    Rresponse0->SetEIonisation(26.e-9);
    Rresponse0->SetSqrtKx3(0.77459667);
    Rresponse0->SetKx2(0.962);
    Rresponse0->SetKx4(0.379);
    Rresponse0->SetSqrtKy3(0.77459667);
    Rresponse0->SetKy2(0.962);
    Rresponse0->SetKy4(0.379);
    Rresponse0->SetPitch(0.25);
//
//      
  for (Int_t i=0; i<7; i++) {
    RICH->SetGeometryModel(i,GeometryV0);
    RICH->SetSegmentationModel(i, SegmentationV0);
    RICH->SetResponseModel(i, Rresponse0);
    RICH->SetNsec(i,1);
  }
}

