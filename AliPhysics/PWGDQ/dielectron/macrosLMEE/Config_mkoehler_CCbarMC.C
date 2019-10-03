void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD);
void AddMCsignal(AliDielectron *die);




AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition);

TString namesDieleData=("CF;kAnyLeg;kMixLegs");



TObjArray *arrNamesDieleData=namesDieleData.Tokenize(";");

const Int_t nDie=arrNamesDieleData->GetEntries();

AliDielectron* Config_mkoehler_CCbarMC(Int_t cutDefinition, Bool_t isAOD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNamesDieleData->GetEntriesFast()){
    name=arrNamesDieleData->At(cutDefinition)->GetName();
  }
  AliDielectron *diele = new AliDielectron(Form("%s",name.Data()),
                             Form("Track cuts: %s",name.Data()));

  // cut setup
  SetupTrackCutsDieleData(diele, cutDefinition, isAOD);
  SetupPairCutsDieleData(diele, cutDefinition, isAOD);

  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielelectron framework histograms will be filled
  //
  InitHistogramsDieleData(diele, cutDefinition, isAOD);
  InitCFDieleData(diele, cutDefinition, isAOD);

  
  return diele;
}


//______________________________________________________________________________________
void SetupTrackCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setup the track cuts
  //

  AddMCsignal(diele);
  
  //ESD quality cuts DielectronTrackCuts

    diele->GetTrackFilter().AddCuts(SetupESDtrackCutsDieleData(cutDefinition));
  
}//SetupTrackCutsDieleData




//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCutsDieleData(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;

  // basic track quality cuts  (basicQ)
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0);

  esdTrackCuts->SetEtaRange( -0.9 , 0.9 );

  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);

  esdTrackCuts->SetPtRange(.8,1e30);

  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);

  // default SPD any
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

  return esdTrackCuts;
}

 

//______________________________________________________________________________________
void SetupPairCutsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  // Setup the pair cuts

  if(cutDefinition >=1){
  //rapidity selection
  AliDielectronVarCuts *pairCut=new AliDielectronVarCuts("pairCut","pairCut");
  pairCut->AddCut(AliDielectronVarManager::kY,-1.,1.);
  diele->GetPairFilter().AddCuts(pairCut);


  AliDielectronVarCuts *mycut = new AliDielectronVarCuts("CutEMCAL","cut for EMCal");
  mycut->AddCut(AliDielectronVarManager::kEMCALEoverP,0.75,1.25);
  mycut->AddCut(AliDielectronVarManager::kEMCALE,3.5,100.);

  

  AliDielectronPairLegCuts *varpair=new AliDielectronPairLegCuts();
  varpair->GetLeg1Filter().AddCuts(mycut);
  varpair->GetLeg2Filter().AddCuts(mycut);
  if(cutDefinition == 1)varpair->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
  if(cutDefinition == 2)varpair->SetCutType(AliDielectronPairLegCuts::kMixLegs);
  diele->GetPairFilter().AddCuts(varpair);
  }


}//SetupPairCutsDieleData

//______________________________________________________________________________________
void InitHistogramsDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(diele->GetName(),diele->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10

   for (Int_t i=0; i<3; ++i){
  	histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));

  }
  
  //legs from pair
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
  }
  //track rotation
  //histos->AddClass(Form("Pair_%s",PairClassName(AliDielectron::kEv1PMRot)));
  //histos->AddClass(Form("Track_Legs_%s",PairClassName(AliDielectron::kEv1PMRot)));
  

    //add histograms to event class

    histos->AddClass("Event");
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",500,-40.,40.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","VtxYxVtxZ","Vertexyz;Z[cm];Y[cm]",500,-40.,40.,400,-0.5,0.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserHistogram("Event","VtxXxVtxZ","Vertexxz;Z[cm];X[cm]",500,-40.,40.,400,-0.5,0.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kXvPrim);
    histos->UserHistogram("Event","VtxYxVtxX","Vertexxz;Z[cm];X[cm]",400,-0.5,0.5,400,-0.5,0.5,AliDielectronVarManager::kXvPrim,AliDielectronVarManager::kYvPrim);
    histos->UserHistogram("Event","MultV0A","MultV0A;multiplicity",1000,0.,1000.,AliDielectronVarManager::kMultV0A);
    histos->UserHistogram("Event","MultV0C","MultV0C;multiplicity",1000,0.,1000.,AliDielectronVarManager::kMultV0C);
    histos->UserHistogram("Event","MultV0","MultV0;multiplicity",1000,0.,1000.,AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","RefMult","RefMult;multiplicity",1000,0.,1000.,AliDielectronVarManager::kRefMult);   
    histos->UserHistogram("Event","RefMultTPConly","RefMultTPConly;multiplicity",1000,0.,1000.,AliDielectronVarManager::kRefMultTPConly);    
    histos->UserHistogram("Event","VZEROchMult","VZEROchMult;multiplicity",1000,0.,1000.,AliDielectronVarManager::kVZEROchMult);     
    histos->UserHistogram("Event","MixingBin","kMixingBin;",100,0.,100.,AliDielectronVarManager::kMixingBin);     
    
  
  
  
  //add histograms to Track classes
  histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",250,0,50.,AliDielectronVarManager::kPt,kTRUE);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC,kTRUE);
  histos->UserHistogram("Track","TPCchi2Cl","Chi-2/Clusters TPC;Chi2/ncls number clusters;#tracks",100,0,10,AliDielectronVarManager::kTPCchi2Cl,kTRUE);
  histos->UserHistogram("Track","TPCnFCls","Number of findable Clusters TPC;Number of findable Clusters TPC;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPC,kTRUE);
  histos->UserHistogram("Track","TPCnFClsfCross","fraction crossed rows/findable;fraction crossed rows/findable;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCfCross,kTRUE);
  histos->UserHistogram("Track","TPCnFClsr","Number of findable Clusters(crossed rows) TPC;Number of findable crossed rows TPC;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr,kTRUE);
  histos->UserHistogram("Track","TPCnFClsrFrac","Number of found/findable Clusters TPC;Number of found/findable Clusters TPC;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCrFrac,kTRUE);
  histos->UserHistogram("Track","TPCnFClsTPCfCross","Fraction of findable Clusters/Cr.rows TPC;Fraction of findable Clusters/Cr.rows TPC;#tracks",200,0.,2.,AliDielectronVarManager::kNFclsTPCfCross,kTRUE);  
  histos->UserHistogram("Track","TPCsignalN","Number of points for TPC Signal;TPC Npoints dEdx;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN,kTRUE);    
  
  histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",200,-1.5,1.5,AliDielectronVarManager::kImpactParXY,kTRUE);
  histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",200,-4.,4.,AliDielectronVarManager::kImpactParZ,kTRUE);
  
  histos->UserHistogram("Track","Eta_Phi","Eta vs Phi; Eta; Phi;#tracks",
                        100,-1.2,1.2,144,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi,kTRUE);

  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,800,20.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_Pt","dEdx;Pt [GeV];TPC signal (arb units);#tracks",
                        200,0.2,20.,800,20.,200.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_P","TPCnSigmaEle;P [GeV];TPCnSigmaEle;#tracks",
                        200,0.2,20.,800,-12.,12.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","TPCnSigmaEle_Pt","TPCnSigmaEle;Pt [GeV];TPCnSigmaEle;#tracks",
                        200,0.2,20.,800,-12.,12.,AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);

  histos->UserHistogram("Track","TPCnSigmaEle_Phi","TPCnSigmaEle;#phi [rad];TPCnSigmaEle;#tracks",
                        200,0.,2*TMath::Pi(),800,-12.,12.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
  
  histos->UserHistogram("Track","TPCnSigmaEle_Eta","TPCnSigmaEle;#eta;TPCnSigmaEle;#tracks",
                        200,-1.,1.,800,-12.,12.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    

  histos->UserHistogram("Track","dEdx_Phi","dEdx vs phi;#phi [rad];TPC signal (arb units);#tracks",
                        200,0.,2*TMath::Pi(),800,20.,200.,AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_Eta","dEdx vs eta;#eta;TPC signal (arb units);#tracks",
                        200,-1.,1.,800,20.,200.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,kTRUE);



  histos->UserHistogram("Track","dEdx_nSigmaEMCal","dEdx vs nSigmaEMCal;NsigmaEmcal;TPC signal (arb units);NSigmaEMCAL",
                        200,-5.,5.,800,20.,200.,AliDielectronVarManager::kEMCALnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_TPCnSigmaEle","dEdx vs TPCnSigmaEle;TPC signal electrons(arbunits);TPC number of sigmas Electrons;TPC signal (a.u.);#tracks",
                        100,-10.,10.,800,20.,200.,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal,kTRUE);

  histos->UserHistogram("Track","dEdx_EoverP","dEdx;EoverP;TPC signal (arbunits);E/P",100,0.,5.,800,20.,200.,AliDielectronVarManager::kEMCALEoverP,AliDielectronVarManager::kTPCsignal,kTRUE);
  
  histos->UserHistogram("Track","nSigmaEMCal_EoverP","NsigmaEmcal;EoverP;NSigmaEMCAL;E/P",100,0.,5.,200,-5.,5.,AliDielectronVarManager::kEMCALEoverP,AliDielectronVarManager::kEMCALnSigmaEle,kTRUE);
 
  histos->UserHistogram("Track","EMCal_E","EmcalE;Cluster Energy [GeV];#Clusters",200,0.,40.,AliDielectronVarManager::kEMCALE,kTRUE);

  histos->UserHistogram("Track","ITS_FirstCls","ITS First Layer;ITS First Layer;#Entries",6,0.,6.,AliDielectronVarManager::kITSLayerFirstCls,kTRUE);
 
  //add histograms to Pair classes
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                        251,-.01,5.01,AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","InvMass2D","Inv.Mass;Pt [GeV]; Inv. Mass [GeV]",
                        20,0.,20.,251,-.01,5.01,AliDielectronVarManager::kPt,AliDielectronVarManager::kM);
     
     
    histos->UserHistogram("Pair","InvMasslong","Inv.Mass;Inv. Mass [GeV];#pairs",
                        301,-.02,15.02,AliDielectronVarManager::kM);
    
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        50,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        50,0.,3.15,AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair","PseudoProperTime","Pseudoproper decay length; pseudoproper-decay-length[#mum];Entries/40#mum",
                          150,-0.3.,0.3,AliDielectronVarManager::kPseudoProperTime);
  
 histos->UserHistogram("Pair","Chi2/NDF","#Chi^{2}/NDF;#Chi^{2}/NDF",
                        100, 0., 20., AliDielectronVarManager::kChi2NDF);
  
  
  diele->SetHistogramManager(histos);
}


void InitCFDieleData(AliDielectron *diele, Int_t cutDefinition, Bool_t isAOD)
{
  //
  // Setupd the CF Manager if needed
  //
  
  AliDielectronCF *cf=new AliDielectronCF(diele->GetName(),diele->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kPt,"1.0,2.0,3.0,4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,15.0,16.0,17.0,18.0,19.0,20.0");
  cf->AddVariable(AliDielectronVarManager::kM,500,0.,10.);//also try variable bi sizes later...
  
    cf->AddVariable(AliDielectronVarManager::kPairType,12,0,12);
  cf->AddVariable(AliDielectronVarManager::kOpeningAngle,31,-0.15,3.15);
  cf->AddVariable(AliDielectronVarManager::kEta,40,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kY,40,-1.,1.);
  cf->AddVariable(AliDielectronVarManager::kPhi,20,0.,20*0.32);
//  cf->AddVariable(AliDielectronVarManager::kPseudoProperTime,300,-0.3,0.3);
//  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeErr,200,0.,0.1);
//  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimeResolution,400,-0.1,0.1);
//  cf->AddVariable(AliDielectronVarManager::kPseudoProperTimePull,400,-0.1,0.1); 
//  cf->AddVariable(AliDielectronVarManager::kChi2NDF,40, 0., 20.);

  //global leg variables
  cf->AddVariable(AliDielectronVarManager::kP,50,0.,5.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPt,"0.,0.5,0.75,0.9,1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0, 8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,30.0,50.0,100.0",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"65, 70, 75, 80, 85, 90, 95, 100, 120, 160",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCchi2Cl,100, 0., 10.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignalN,160,-0.5,159.5,kTRUE);   
  cf->AddVariable(AliDielectronVarManager::kEta,44,-1.2,1.2,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kPhi,64,0.,64*0.1,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kImpactParXY,200,-1.,1.,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kImpactParZ,600,-3.,3.,kTRUE);
  //TPC
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,"40.,50.,55.,60.,65.,68.,70.,72.,75.,80.,90.,100.,110.,200.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,20,-3.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,8,1.,4.5,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,8,0.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsSTPC,20,0.,1.,kTRUE);//shared cluster
    //ITS
  cf->AddVariable(AliDielectronVarManager::kNclsITS,6,0.,6.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,6,0.,6.,kTRUE);
  
  //TOF
//  cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,20,-3.5,4.,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kTOFnSigmaPio,8,0.,4.,kTRUE);
  //EMCal variables
//  cf->AddVariable(AliDielectronVarManager::kEMCALE,20,0.,20.,kTRUE); 
//  cf->AddVariable(AliDielectronVarManager::kEMCALnSigmaEle,50,-5.,5.,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kEMCALNCells,50,0,50,kTRUE);
//  cf->AddVariable(AliDielectronVarManager::kEMCALEoverP,"0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.8,2.0,4.0",kTRUE);
  
//    cf->AddVariable(AliDielectronVarManager::kMixingBin,100,0.,100.);
    cf->AddVariable(AliDielectronVarManager::kZvPrim,20,-20.,20.);

   
    /*
    if (cutDefinition == 0){
      
	    cf->SetStepForMCtruth();
	    
  }
    */
    
    
    
        cf->SetStepForMCtruth();
	
	
  //leg variables
  diele->SetCFManagerPair(cf);
  
}

void AddMCsignal(AliDielectron *die){
  
  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("diEleOpenCharm","di-electrons from open charm");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenCharm);

//    501 - open beauty mesons     
  AliDielectronSignalMC* diEleOpenBeauty = new AliDielectronSignalMC("diEleOpenBeauty","di-electrons from open beauty");  // dielectrons originating from open beauty hadrons
  diEleOpenBeauty->SetLegPDGs(11,-11);
  diEleOpenBeauty->SetMotherPDGs(501,501);
  diEleOpenBeauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenBeauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBeauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBeauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenBeauty);

}//void AddMCsignal(AliDielectron *die)




