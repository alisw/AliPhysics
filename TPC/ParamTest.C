TFile f("prf2d.root");
AliTPCPRF2D prf1;
AliTPCPRF2D prf2; 
AliTPCRF1D  rf(kTRUE);
rf.SetGauss(0.3,0.5,1)
rf.Update();
prf1.Read("prf_205035_Gati_062074_d03");
prf2.Read("prf_205035_Gati_062074_d03");
AliTPCParamSR param;       
param.SetInnerPRF(&prf1)
param.SetOuterPRF(&prf2)
param.SetTimeRF(&rf);
Float_t x[3]={5.1,0,100}
Int_t index[3]
param.CalcResponse(x,index)
