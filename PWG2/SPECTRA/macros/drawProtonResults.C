void drawProtonResults(const char* esdFileName) {
  //Macro to visualize the proton ratio results
  //It also visualizes the QA plots
  drawQAPlots(esdFileName);
}

//___________________________________________________//
void drawQAPlots(const char* esdFileName) {
  //Draws the QA plots from the output of the analysis
  //=========================================================//
  //List of cuts
  TFile *fCutFile = TFile::Open("ListOfCuts.root");
  TCanvas *cListOfCuts = dynamic_cast<TCanvas *>(fCutFile->Get("cListOfCuts"));
  cListOfCuts->Draw();

  //=========================================================//
  //QA plots
  TFile *f = TFile::Open(esdFileName);
  TList *listQA = dynamic_cast<TList *>(f->Get("outputQAList"));
  TList *gListGlobalQA = dynamic_cast<TList *>(listQA->At(0));

  //================QA plots================//
  TList *fQA2DList = dynamic_cast<TList *>(gListGlobalQA->At(0));

  //Accepted protons
  TList *fQAProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQA->At(1));
  TH1F *gProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(0));
  TH1F *gProtonsChi2PerClusterITSPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(1));
  TH1F *gProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(2));
  TH1F *gProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(3));
  TH1F *gProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(4));
  TH1F *gProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(5));
  TH1F *gProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(6));
  TH1F *gProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(7));
  TH1F *gProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(8));
  TH1F *gProtonsSigmaToVertexPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(9));
  TH1F *gProtonsSigmaToVertexTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(10));
  TH1F *gProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(11));
  TH1F *gProtonsDCAXYTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(12));
  TH1F *gProtonsDCAZPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(13));
  TH1F *gProtonsDCAZTPCPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(14));
  TH1F *gProtonsConstrainChi2Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(15));
  TH1F *gProtonsITSRefitPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(16));
  TH1F *gProtonsTPCRefitPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(17));
  TH1F *gProtonsESDpidPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(18));
  TH1F *gProtonsTPCpidPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(19));
  TH1F *gProtonsPointOnITSLayer1Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(20));
  TH1F *gProtonsPointOnITSLayer2Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(21));
  TH1F *gProtonsPointOnITSLayer3Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(22));
  TH1F *gProtonsPointOnITSLayer4Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(23));
  TH1F *gProtonsPointOnITSLayer5Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(24));
  TH1F *gProtonsPointOnITSLayer6Pass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(25));
  TH1F *gProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAProtonsAcceptedList->At(26));

  //Rejected protons
  TList *fQAProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQA->At(2));
  TH1F *gProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(0));
  TH1F *gProtonsChi2PerClusterITSReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(1));
  TH1F *gProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(2));
  TH1F *gProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(3));
  TH1F *gProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(4));
  TH1F *gProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(5));
  TH1F *gProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(6));
  TH1F *gProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(7));
  TH1F *gProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(8));
  TH1F *gProtonsSigmaToVertexReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(9));
  TH1F *gProtonsSigmaToVertexTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(10));
  TH1F *gProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(11));
  TH1F *gProtonsDCAXYTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(12));
  TH1F *gProtonsDCAZReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(13));
  TH1F *gProtonsDCAZTPCReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(14));
  TH1F *gProtonsConstrainChi2Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(15));
  TH1F *gProtonsITSRefitReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(16));
  TH1F *gProtonsTPCRefitReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(17));
  TH1F *gProtonsESDpidReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(18));
  TH1F *gProtonsTPCpidReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(19));
  TH1F *gProtonsPointOnITSLayer1Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(20));
  TH1F *gProtonsPointOnITSLayer2Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(21));
  TH1F *gProtonsPointOnITSLayer3Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(22));
  TH1F *gProtonsPointOnITSLayer4Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(23));
  TH1F *gProtonsPointOnITSLayer5Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(24));
  TH1F *gProtonsPointOnITSLayer6Reject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(25));
  TH1F *gProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAProtonsRejectedList->At(26));

  //Accepted antiprotons
  TList *fQAAntiProtonsAcceptedList = dynamic_cast<TList *>(gListGlobalQA->At(3));
  TH1F *gAntiProtonsITSClustersPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(0));
  TH1F *gAntiProtonsChi2PerClusterITSPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(1));
  TH1F *gAntiProtonsTPCClustersPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(2));
  TH1F *gAntiProtonsChi2PerClusterTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(3));
  TH1F *gAntiProtonsExtCov11Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(4));
  TH1F *gAntiProtonsExtCov22Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(5));
  TH1F *gAntiProtonsExtCov33Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(6));
  TH1F *gAntiProtonsExtCov44Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(7));
  TH1F *gAntiProtonsExtCov55Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(8));
  TH1F *gAntiProtonsSigmaToVertexPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(9));
  TH1F *gAntiProtonsSigmaToVertexTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(10));
  TH1F *gAntiProtonsDCAXYPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(11));
  TH1F *gAntiProtonsDCAXYTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(12));
  TH1F *gAntiProtonsDCAZPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(13));
  TH1F *gAntiProtonsDCAZTPCPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(14));
  TH1F *gAntiProtonsConstrainChi2Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(15));
  TH1F *gAntiProtonsITSRefitPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(16));
  TH1F *gAntiProtonsTPCRefitPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(17));
  TH1F *gAntiProtonsESDpidPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(18));
  TH1F *gAntiProtonsTPCpidPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(19));
  TH1F *gAntiProtonsPointOnITSLayer1Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(20));
  TH1F *gAntiProtonsPointOnITSLayer2Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(21));
  TH1F *gAntiProtonsPointOnITSLayer3Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(22));
  TH1F *gAntiProtonsPointOnITSLayer4Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(23));
  TH1F *gAntiProtonsPointOnITSLayer5Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(24));
  TH1F *gAntiProtonsPointOnITSLayer6Pass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(25));
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsPass = dynamic_cast<TH1F *>(fQAAntiProtonsAcceptedList->At(26));

  //Rejected antiprotons
  TList *fQAAntiProtonsRejectedList = dynamic_cast<TList *>(gListGlobalQA->At(4));
  TH1F *gAntiProtonsITSClustersReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(0));
  TH1F *gAntiProtonsChi2PerClusterITSReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(1));
  TH1F *gAntiProtonsTPCClustersReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(2));
  TH1F *gAntiProtonsChi2PerClusterTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(3));
  TH1F *gAntiProtonsExtCov11Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(4));
  TH1F *gAntiProtonsExtCov22Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(5));
  TH1F *gAntiProtonsExtCov33Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(6));
  TH1F *gAntiProtonsExtCov44Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(7));
  TH1F *gAntiProtonsExtCov55Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(8));
  TH1F *gAntiProtonsSigmaToVertexReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(9));
  TH1F *gAntiProtonsSigmaToVertexTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(10));
  TH1F *gAntiProtonsDCAXYReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(11));
  TH1F *gAntiProtonsDCAXYTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(12));
  TH1F *gAntiProtonsDCAZReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(13));
  TH1F *gAntiProtonsDCAZTPCReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(14));
  TH1F *gAntiProtonsConstrainChi2Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(15));
  TH1F *gAntiProtonsITSRefitReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(16));
  TH1F *gAntiProtonsTPCRefitReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(17));
  TH1F *gAntiProtonsESDpidReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(18));
  TH1F *gAntiProtonsTPCpidReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(19));
  TH1F *gAntiProtonsPointOnITSLayer1Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(20));
  TH1F *gAntiProtonsPointOnITSLayer2Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(21));
  TH1F *gAntiProtonsPointOnITSLayer3Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(22));
  TH1F *gAntiProtonsPointOnITSLayer4Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(23));
  TH1F *gAntiProtonsPointOnITSLayer5Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(24));
  TH1F *gAntiProtonsPointOnITSLayer6Reject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(25));
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsReject = dynamic_cast<TH1F *>(fQAAntiProtonsRejectedList->At(26));

  //__________________________________________________//
  TCanvas *c1 = new TCanvas("c1","",0,0,600,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->Divide(2,1);
  c1->cd(1); gProtonsITSClustersPass->Draw(); 
  gProtonsITSClustersReject->Draw("same");
  c1->cd(2); gAntiProtonsITSClustersPass->Draw(); 
  gAntiProtonsITSClustersReject->Draw("same");

  TCanvas *c2 = new TCanvas("c2","",0,100,600,400);
  c2->SetFillColor(10); c2->SetHighLightColor(10);
  c2->Divide(2,1);
  c2->cd(1); gProtonsChi2PerClusterITSPass->Draw(); 
  gProtonsChi2PerClusterITSReject->Draw("same");
  c2->cd(2); gAntiProtonsChi2PerClusterITSPass->Draw(); 
  gAntiProtonsChi2PerClusterITSReject->Draw("same");

  TCanvas *c3 = new TCanvas("c3","",0,200,600,400);
  c3->SetFillColor(10); c3->SetHighLightColor(10);
  c3->Divide(2,1);
  c3->cd(1); gProtonsTPCClustersPass->Draw();
  gProtonsTPCClustersReject->Draw("same");
  c3->cd(2); gAntiProtonsTPCClustersPass->Draw();
  gAntiProtonsTPCClustersReject->Draw("same");

  TCanvas *c4 = new TCanvas("c4","",0,300,600,400);
  c4->SetFillColor(10); c4->SetHighLightColor(10);
  c4->Divide(2,1);
  c4->cd(1); gProtonsChi2PerClusterTPCPass->Draw(); 
  gProtonsChi2PerClusterTPCReject->Draw("same");
  c4->cd(2); gAntiProtonsChi2PerClusterTPCPass->Draw(); 
  gAntiProtonsChi2PerClusterTPCReject->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",0,400,600,400);
  c5->SetFillColor(10); c5->SetHighLightColor(10);
  c5->Divide(2,1);
  c5->cd(1); gProtonsExtCov11Pass->Draw(); 
  gProtonsExtCov11Reject->Draw("same");
  c5->cd(2); gAntiProtonsExtCov11Pass->Draw(); 
  gAntiProtonsExtCov11Reject->Draw("same");

  TCanvas *c6 = new TCanvas("c6","",0,500,600,400);
  c6->SetFillColor(10); c6->SetHighLightColor(10);
  c6->Divide(2,1);
  c6->cd(1); gProtonsExtCov22Pass->Draw(); 
  gProtonsExtCov22Reject->Draw("same");
  c6->cd(2); gAntiProtonsExtCov22Pass->Draw(); 
  gAntiProtonsExtCov22Reject->Draw("same");

  TCanvas *c7 = new TCanvas("c7","",600,0,600,400);
  c7->SetFillColor(10); c7->SetHighLightColor(10);
  c7->Divide(2,1);
  c7->cd(1); gProtonsExtCov33Pass->Draw(); 
  gProtonsExtCov33Reject->Draw("same");
  c7->cd(2); gAntiProtonsExtCov33Pass->Draw(); 
  gAntiProtonsExtCov33Reject->Draw("same");

  TCanvas *c8 = new TCanvas("c8","",600,100,600,400);
  c8->SetFillColor(10); c8->SetHighLightColor(10);
  c8->Divide(2,1);
  c8->cd(1); gProtonsExtCov44Pass->Draw(); 
  gProtonsExtCov44Reject->Draw("same");
  c8->cd(2); gAntiProtonsExtCov44Pass->Draw(); 
  gAntiProtonsExtCov44Reject->Draw("same");

  TCanvas *c9 = new TCanvas("c9","",600,200,600,400);
  c9->SetFillColor(10); c9->SetHighLightColor(10);
  c9->Divide(2,1);
  c9->cd(1); gProtonsExtCov55Pass->Draw(); 
  gProtonsExtCov55Reject->Draw("same");
  c9->cd(2); gAntiProtonsExtCov55Pass->Draw(); 
  gAntiProtonsExtCov55Reject->Draw("same");

  TCanvas *c10 = new TCanvas("c10","",600,300,600,400);
  c10->SetFillColor(10); c10->SetHighLightColor(10);
  c10->Divide(2,1);
  c10->cd(1); gProtonsSigmaToVertexPass->Draw(); 
  gProtonsSigmaToVertexReject->Draw("same");
  c10->cd(2); gAntiProtonsSigmaToVertexPass->Draw(); 
  gAntiProtonsSigmaToVertexReject->Draw("same");

  TCanvas *c11 = new TCanvas("c11","",600,400,600,400);
  c11->SetFillColor(10); c11->SetHighLightColor(10);
  c11->Divide(2,1);
  c11->cd(1); gProtonsSigmaToVertexTPCPass->Draw(); 
  gProtonsSigmaToVertexTPCReject->Draw("same");
  c11->cd(2); gAntiProtonsSigmaToVertexTPCPass->Draw(); 
  gAntiProtonsSigmaToVertexTPCReject->Draw("same");

  TCanvas *c12 = new TCanvas("c12","",600,500,600,400);
  c12->SetFillColor(10); c12->SetHighLightColor(10);
  c12->Divide(2,1);
  c12->cd(1); gProtonsDCAXYPass->Draw(); 
  gProtonsDCAXYReject->Draw("same");
  c12->cd(2); gAntiProtonsDCAXYPass->Draw(); 
  gAntiProtonsDCAXYReject->Draw("same");

  TCanvas *c13 = new TCanvas("c13","",1200,0,600,400);
  c13->SetFillColor(10); c13->SetHighLightColor(10);
  c13->Divide(2,1);
  c13->cd(1); gProtonsDCAXYTPCPass->Draw(); 
  gProtonsDCAXYTPCReject->Draw("same");
  c13->cd(2); gAntiProtonsDCAXYTPCPass->Draw(); 
  gAntiProtonsDCAXYTPCReject->Draw("same");

  TCanvas *c14 = new TCanvas("c14","",1200,100,600,400);
  c14->SetFillColor(10); c14->SetHighLightColor(10);
  c14->Divide(2,1);
  c14->cd(1); gProtonsDCAZPass->Draw(); 
  gProtonsDCAZReject->Draw("same");
  c14->cd(2); gAntiProtonsDCAZPass->Draw(); 
  gAntiProtonsDCAZReject->Draw("same");

  TCanvas *c15 = new TCanvas("c15","",1200,200,600,400);
  c15->SetFillColor(10); c15->SetHighLightColor(10);
  c15->Divide(2,1);
  c15->cd(1); gProtonsDCAZTPCPass->Draw(); 
  gProtonsDCAZTPCReject->Draw("same");
  c15->cd(2); gAntiProtonsDCAZTPCPass->Draw(); 
  gAntiProtonsDCAZTPCReject->Draw("same");

  TCanvas *c16 = new TCanvas("c16","",1200,300,600,400);
  c16->SetFillColor(10); c16->SetHighLightColor(10);
  c16->Divide(2,1);
  c16->cd(1); gProtonsNumberOfTPCdEdxPointsPass->Draw(); 
  gProtonsNumberOfTPCdEdxPointsReject->Draw("same");
  c16->cd(2); gAntiProtonsNumberOfTPCdEdxPointsPass->Draw(); 
  gAntiProtonsNumberOfTPCdEdxPointsReject->Draw("same");

  //================Vertex QA================//
  TList *gListVertexQA = dynamic_cast<TList *>(listQA->At(1));
  TH1F *gHistVx = dynamic_cast<TH1F *>(gListVertexQA->At(0));
  TH1F *gHistVxAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(1));
  gHistVxAccepted->SetFillColor(10);
  TH1F *gHistVy = dynamic_cast<TH1F *>(gListVertexQA->At(2));
  TH1F *gHistVyAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(3));
  gHistVyAccepted->SetFillColor(10);
  TH1F *gHistVz = dynamic_cast<TH1F *>(gListVertexQA->At(4));
  TH1F *gHistVzAccepted = dynamic_cast<TH1F *>(gListVertexQA->At(5));
  gHistVzAccepted->SetFillColor(10);

  TCanvas *cVertex = new TCanvas("cVertex","Veretx QA",0,0,900,400);
  cVertex->SetFillColor(10); cVertex->SetHighLightColor(10);
  cVertex->Divide(3,1);
  cVertex->cd(1); gHistVx->Draw(); gHistVxAccepted->Draw("same");
  cVertex->cd(2); gHistVy->Draw(); gHistVyAccepted->Draw("same");
  cVertex->cd(3); gHistVz->Draw(); gHistVzAccepted->Draw("same");
}
  
