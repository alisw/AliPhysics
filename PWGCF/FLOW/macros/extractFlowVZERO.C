extractFlowVZERO(Int_t icentr,Int_t spec,Int_t arm=2,Bool_t isMC=kFALSE){
  // NUA correction currently are missing
  char name[100];
  snprintf(name,100,"outVZEROv%i.root",arm);
  TFile *fo = new TFile(name);
  snprintf(name,100,"contVZEROv%i",arm);
  TList *cont = (TList *) fo->Get(name);

  Float_t xMin[5] = {icentr,-1,0,-10,0};
  Float_t xMax[5] = {icentr,1,1,10,1.5};

  cont->ls();

  TProfile *p1 = cont->At(2);
  TProfile *p2 = cont->At(3);
  TProfile *p3 = cont->At(4);

  Float_t res1=0,res2=0,res3=0; 
  Float_t eres1=0,eres2=0,eres3=0; 

  Int_t i = icentr;
  if(p1->GetBinError(i+1)){
    eres1 += 1./p1->GetBinError(i+1)/p1->GetBinError(i+1);
    res1 += p1->GetBinContent(i+1)/p1->GetBinError(i+1)/p1->GetBinError(i+1);      
  }
  if(p2->GetBinError(i+1)){
    eres2 += 1./p2->GetBinError(i+1)/p2->GetBinError(i+1);
    res2 += p2->GetBinContent(i+1)/p2->GetBinError(i+1)/p2->GetBinError(i+1);      
    }
  if(p3->GetBinError(i+1)){
    eres3 += 1./p3->GetBinError(i+1)/p3->GetBinError(i+1);
    res3 += p3->GetBinContent(i+1)/p3->GetBinError(i+1)/p3->GetBinError(i+1);      
  }
  
  res1 /= eres1;
  res2 /= eres2;
  res3 /= eres3;
  
  AliFlowVZEROResults *a = (AliFlowVZEROResults *) cont->At(0);
  AliFlowVZEROResults *b = (AliFlowVZEROResults *) cont->At(1);
  TProfile *pp = a->GetV2(spec,xMin,xMax);
  TProfile *pp2 = b->GetV2(spec,xMin,xMax);
    
  
  Float_t scaling = sqrt(res1*res3/res2);
  pp->Scale(1./scaling);
  
  printf("resolution V0A = %f\n",scaling);
  Float_t err1_2 = eres1*eres1/res1/res1/4 +
    eres2*eres2/res2/res2/4 +
    eres3*eres3/res3/res3/4;
  Float_t err2_2 = err1_2;
  err1_2 /= scaling*scaling;
  scaling = sqrt(res2*res3/res1);
  err2_2 /= scaling*scaling;
  pp2->Scale(1./scaling);
  printf("resolution V0C =%f\n",scaling);

  pp->SetName("V0A");
  pp2->SetName("V0C");

  pp->Draw();
  pp2->Draw("SAME");

  if(arm == 2 && isMC){
    snprintf(name,100,"outVZEROmc.root");
    fo = new TFile(name);
    snprintf(name,100,"contVZEROmc");
    cont = (TList *) fo->Get(name);
    AliFlowVZEROResults *c = (AliFlowVZEROResults *) cont->At(0);
    c->GetV2(spec,xMin,xMax)->Draw("SAME");
  }
}
