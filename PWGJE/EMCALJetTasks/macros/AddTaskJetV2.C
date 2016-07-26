/* AddTask macro for Jet V2 task 
 * Redmer Alexander Bertens, rbertens@cern.ch
 * Utrecht University, Utrecht, Netherlands */

AliAnalysisTaskJetV2* AddTaskJetV2(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t   jetradius          = 0.2,
  Double_t   jetptcut           = 1,
  Double_t   jetareacut         = 0.557,
  const char* type              = "TPC",
  Int_t      leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskJetV2",
  UInt_t     runMode            = AliAnalysisTaskJetV2::kGrid,
  Bool_t     fillQA             = kTRUE,
  TString    fitOpts            = "WLQI",
  UInt_t     fitType            = AliAnalysisTaskJetV2::kCombined,
  TArrayD    *centralities      = 0x0,
  TRandom3   *randomizer        = 0x0,
  Double_t   trackptcut         = .15,
  Bool_t     LHC10h             = kFALSE,
  Bool_t     addEPweights       = kFALSE,
  Bool_t     baseClassHistos    = kTRUE,
  Float_t    minEta             = -.7,
  Float_t    maxEta             = .7)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (!strcmp(type, "TPC"))
    name += "_TPC";
  else if (!strcmp(type, "EMCAL"))
    name += "_EMCAL";
  else if (!strcmp(type, "USER")) 
    name += "_USER";

  // create instance of the object
  AliAnalysisTaskJetV2* jetTask = new AliAnalysisTaskJetV2(name, runMode, baseClassHistos);
  
  // create and connect data containers
  AliParticleContainer* partCont = jetTask->AddParticleContainer(ntracks);
  if(partCont) {
      partCont->SetName("Tracks");
      partCont->SetParticlePtCut(trackptcut);
      partCont->SetClassName("AliAODTrack");
  }
  TString tmp(nclusters);
  AliClusterContainer* clusterCont = 0x0;
  if(!tmp.IsNull()) {
      clusterCont = jetTask->AddClusterContainer(nclusters);
      jetTask->SetAnalysisType(AliAnalysisTaskJetV2::kFull);
  }
  AliJetContainer* jetCont = jetTask->AddJetContainer(njets, type, jetradius);
  if(jetCont) {
      jetCont->SetName("Jets");
      jetCont->SetPercAreaCut(jetareacut);
      jetCont->SetRhoName(nrho);
      if(minEta > -.7 || maxEta < 0.7) {
          jetCont->SetJetAcceptanceType(AliJetContainer::kUser);
          jetCont->SetJetEtaLimits(minEta, maxEta);
      }
      if(partCont)      jetCont->ConnectParticleContainer(partCont);
      if(clusterCont)   jetCont->ConnectClusterContainer(clusterCont);
      jetCont->PrintCuts();
  }

  // task specific setters
  jetTask->SetFillQAHistograms(fillQA);
  jetTask->SetModulationFitType(fitType);
  jetTask->SetModulationFitOptions(fitOpts);
  jetTask->SetModulationFitMinMaxP(.01, 1);
  // if centralities haven't been specified use defaults
  if(!centralities) {
      if(LHC10h) {
          Double_t c[] = {0., 5., 10., 30., 50., 70., 90.};
          jetTask->SetCentralityClasses(new TArrayD(sizeof(c)/sizeof(c[0]), c));
      } else {
          Double_t c[] = {0., 5., 10., 30., 50., 70., 90.};
          jetTask->SetCentralityClasses(new TArrayD(sizeof(c)/sizeof(c[0]), c));
      }
  }  else jetTask->SetCentralityClasses(centralities);
  // if a randomized hasn't specified use a safe default 
  if(!randomizer) jetTask->SetRandomSeed(new TRandom3(0));

  // pass the expected run lists to the task. the total list is used for QA plots which are stored per run-number, the semi-good list is used to change the phi acceptance of jets and pico trakcs, and - if an alternatie is provided - switch to a 'small rho' task, which also runs on limited acceptance
  Int_t totalRuns[] = {167813, 167988, 168066, 168068, 168069, 168076, 168104, 168212, 168311, 168322, 168325, 168341, 168361, 168362, 168458, 168460, 168461, 168992, 169091, 169094, 169138, 169143, 169167, 169417, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 170027, 170036, 170081, /* up till here original good TPC list */169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309, /* original semi-good tpc list */169415, 169411, 169035, 168988, 168984, 168826, 168777, 168512, 168511, 168467, 168464, 168342, 168310, 168115, 168108, 168107, 167987, 167915, 167903, /*new runs, good according to RCT */ 169238, 169160, 169156, 169148, 169145, 169144 /* run swith missing OROC 8 but seem ok in QA */};

  Int_t totalRuns10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137544, 137541, 137539, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};

  // set the runnumbers for either 10h or 11h
  (LHC10h) ? jetTask->SetExpectedRuns(new TArrayI(sizeof(totalRuns10h)/sizeof(totalRuns10h[0]), totalRuns10h)) : jetTask->SetExpectedRuns(new TArrayI(sizeof(totalRuns)/sizeof(totalRuns[0]), totalRuns));

  Int_t semiGoodRuns[] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};

  // set the semi-good runnumbers for 10h or 11h (10h has no semi-good numbers, so pass a NULL pointer)
  (LHC10h) ? jetTask->SetExpectedSemiGoodRuns(0x0) : jetTask->SetExpectedSemiGoodRuns(new TArrayI(sizeof(semiGoodRuns)/sizeof(semiGoodRuns[0]), semiGoodRuns));

  // and if 10h, pass this info to the task so acceptance isn't changed
  if(LHC10h) {
      jetTask->SetCollisionType(AliAnalysisTaskJetV2::kPbPb10h);
      jetTask->SetNeedEmcalGeom(kFALSE);
  }

  if(addEPweights) {
     // add a set of default event plane weights, FIXME should be more elegant
      TH1F *cen_0_2 = new TH1F("cen_0_2","cen_0_2",40,-1.5708,1.5708);
      cen_0_2->SetBinContent(1,94735);        cen_0_2->SetBinError(1,307.791);
      cen_0_2->SetBinContent(2,95146);        cen_0_2->SetBinError(2,308.457);
      cen_0_2->SetBinContent(3,95351);        cen_0_2->SetBinError(3,308.79);
      cen_0_2->SetBinContent(4,94992);        cen_0_2->SetBinError(4,308.208);
      cen_0_2->SetBinContent(5,95601);        cen_0_2->SetBinError(5,309.194);
      cen_0_2->SetBinContent(6,96160);        cen_0_2->SetBinError(6,310.097);
      cen_0_2->SetBinContent(7,96695);        cen_0_2->SetBinError(7,310.958);
      cen_0_2->SetBinContent(8,97219);        cen_0_2->SetBinError(8,311.8);
      cen_0_2->SetBinContent(9,97190);        cen_0_2->SetBinError(9,311.753);
      cen_0_2->SetBinContent(10,97791);       cen_0_2->SetBinError(10,312.716);
      cen_0_2->SetBinContent(11,97401);       cen_0_2->SetBinError(11,312.091);
      cen_0_2->SetBinContent(12,97236);       cen_0_2->SetBinError(12,311.827);
      cen_0_2->SetBinContent(13,96961);       cen_0_2->SetBinError(13,311.386);
      cen_0_2->SetBinContent(14,96768);       cen_0_2->SetBinError(14,311.076);
      cen_0_2->SetBinContent(15,95826);       cen_0_2->SetBinError(15,309.558);
      cen_0_2->SetBinContent(16,95398);       cen_0_2->SetBinError(16,308.866);
      cen_0_2->SetBinContent(17,95515);       cen_0_2->SetBinError(17,309.055);
      cen_0_2->SetBinContent(18,94947);       cen_0_2->SetBinError(18,308.135);
      cen_0_2->SetBinContent(19,95186);       cen_0_2->SetBinError(19,308.522);
      cen_0_2->SetBinContent(20,94636);       cen_0_2->SetBinError(20,307.63);
      cen_0_2->SetBinContent(21,94392);       cen_0_2->SetBinError(21,307.233);
      cen_0_2->SetBinContent(22,93569);       cen_0_2->SetBinError(22,305.891);
      cen_0_2->SetBinContent(23,93569);       cen_0_2->SetBinError(23,305.891);
      cen_0_2->SetBinContent(24,93071);       cen_0_2->SetBinError(24,305.075);
      cen_0_2->SetBinContent(25,93268);       cen_0_2->SetBinError(25,305.398);
      cen_0_2->SetBinContent(26,93554);       cen_0_2->SetBinError(26,305.866);
      cen_0_2->SetBinContent(27,93144);       cen_0_2->SetBinError(27,305.195);
      cen_0_2->SetBinContent(28,92574);       cen_0_2->SetBinError(28,304.26);
      cen_0_2->SetBinContent(29,93524);       cen_0_2->SetBinError(29,305.817);
      cen_0_2->SetBinContent(30,93150);       cen_0_2->SetBinError(30,305.205);
      cen_0_2->SetBinContent(31,93167);       cen_0_2->SetBinError(31,305.233);
      cen_0_2->SetBinContent(32,92736);       cen_0_2->SetBinError(32,304.526);
      cen_0_2->SetBinContent(33,92283);       cen_0_2->SetBinError(33,303.781);
      cen_0_2->SetBinContent(34,92465);       cen_0_2->SetBinError(34,304.081);
      cen_0_2->SetBinContent(35,91792);       cen_0_2->SetBinError(35,302.972);
      cen_0_2->SetBinContent(36,92733);       cen_0_2->SetBinError(36,304.521);
      cen_0_2->SetBinContent(37,92764);       cen_0_2->SetBinError(37,304.572);
      cen_0_2->SetBinContent(38,92983);       cen_0_2->SetBinError(38,304.931);
      cen_0_2->SetBinContent(39,93020);       cen_0_2->SetBinError(39,304.992);
      cen_0_2->SetBinContent(40,94158);       cen_0_2->SetBinError(40,306.852);
                                              cen_0_2->SetEntries(3.78267e+06);

      TH1F *cen_2_4 = new TH1F("cen_2_4","cen_2_4",40,-1.5708,1.5708);
      cen_2_4->SetBinContent(1,96710);        cen_2_4->SetBinError(1,310.982);
      cen_2_4->SetBinContent(2,96521);        cen_2_4->SetBinError(2,310.678);
      cen_2_4->SetBinContent(3,97002);        cen_2_4->SetBinError(3,311.451);
      cen_2_4->SetBinContent(4,96756);        cen_2_4->SetBinError(4,311.056);
      cen_2_4->SetBinContent(5,97302);        cen_2_4->SetBinError(5,311.933);
      cen_2_4->SetBinContent(6,97922);        cen_2_4->SetBinError(6,312.925);
      cen_2_4->SetBinContent(7,98138);        cen_2_4->SetBinError(7,313.27);
      cen_2_4->SetBinContent(8,97963);        cen_2_4->SetBinError(8,312.99);
      cen_2_4->SetBinContent(9,98238);        cen_2_4->SetBinError(9,313.429);
      cen_2_4->SetBinContent(10,99241);       cen_2_4->SetBinError(10,315.025);
      cen_2_4->SetBinContent(11,99001);       cen_2_4->SetBinError(11,314.644);
      cen_2_4->SetBinContent(12,98935);       cen_2_4->SetBinError(12,314.539);
      cen_2_4->SetBinContent(13,98333);       cen_2_4->SetBinError(13,313.581);
      cen_2_4->SetBinContent(14,98090);       cen_2_4->SetBinError(14,313.193);
      cen_2_4->SetBinContent(15,97929);       cen_2_4->SetBinError(15,312.936);
      cen_2_4->SetBinContent(16,97053);       cen_2_4->SetBinError(16,311.533);
      cen_2_4->SetBinContent(17,97195);       cen_2_4->SetBinError(17,311.761);
      cen_2_4->SetBinContent(18,97052);       cen_2_4->SetBinError(18,311.532);
      cen_2_4->SetBinContent(19,97234);       cen_2_4->SetBinError(19,311.824);
      cen_2_4->SetBinContent(20,97431);       cen_2_4->SetBinError(20,312.139);
      cen_2_4->SetBinContent(21,96637);       cen_2_4->SetBinError(21,310.865);
      cen_2_4->SetBinContent(22,96737);       cen_2_4->SetBinError(22,311.026);
      cen_2_4->SetBinContent(23,96031);       cen_2_4->SetBinError(23,309.889);
      cen_2_4->SetBinContent(24,95833);       cen_2_4->SetBinError(24,309.569);
      cen_2_4->SetBinContent(25,95999);       cen_2_4->SetBinError(25,309.837);
      cen_2_4->SetBinContent(26,95868);       cen_2_4->SetBinError(26,309.626);
      cen_2_4->SetBinContent(27,95085);       cen_2_4->SetBinError(27,308.359);
      cen_2_4->SetBinContent(28,95340);       cen_2_4->SetBinError(28,308.772);
      cen_2_4->SetBinContent(29,95021);       cen_2_4->SetBinError(29,308.255);
      cen_2_4->SetBinContent(30,95272);       cen_2_4->SetBinError(30,308.662);
      cen_2_4->SetBinContent(31,95196);       cen_2_4->SetBinError(31,308.538);
      cen_2_4->SetBinContent(32,94857);       cen_2_4->SetBinError(32,307.989);
      cen_2_4->SetBinContent(33,94366);       cen_2_4->SetBinError(33,307.19);
      cen_2_4->SetBinContent(34,93971);       cen_2_4->SetBinError(34,306.547);
      cen_2_4->SetBinContent(35,94704);       cen_2_4->SetBinError(35,307.74);
      cen_2_4->SetBinContent(36,93867);       cen_2_4->SetBinError(36,306.377);
      cen_2_4->SetBinContent(37,95382);       cen_2_4->SetBinError(37,308.84);
      cen_2_4->SetBinContent(38,94614);       cen_2_4->SetBinError(38,307.594);
      cen_2_4->SetBinContent(39,94825);       cen_2_4->SetBinError(39,307.937);
      cen_2_4->SetBinContent(40,95784);       cen_2_4->SetBinError(40,309.49);
                                              cen_2_4->SetEntries(3.85944e+06);

      TH1F *cen_4_6 = new TH1F("cen_4_6","cen_4_6",40,-1.5708,1.5708);
      cen_4_6->SetBinContent(1,93671);        cen_4_6->SetBinError(1,306.057);
      cen_4_6->SetBinContent(2,93879);        cen_4_6->SetBinError(2,306.397);
      cen_4_6->SetBinContent(3,94534);        cen_4_6->SetBinError(3,307.464);
      cen_4_6->SetBinContent(4,95317);        cen_4_6->SetBinError(4,308.735);
      cen_4_6->SetBinContent(5,95943);        cen_4_6->SetBinError(5,309.747);
      cen_4_6->SetBinContent(6,96615);        cen_4_6->SetBinError(6,310.83);
      cen_4_6->SetBinContent(7,96664);        cen_4_6->SetBinError(7,310.908);
      cen_4_6->SetBinContent(8,97577);        cen_4_6->SetBinError(8,312.373);
      cen_4_6->SetBinContent(9,98161);        cen_4_6->SetBinError(9,313.307);
      cen_4_6->SetBinContent(10,99558);       cen_4_6->SetBinError(10,315.528);
      cen_4_6->SetBinContent(11,99298);       cen_4_6->SetBinError(11,315.116);
      cen_4_6->SetBinContent(12,99623);       cen_4_6->SetBinError(12,315.631);
      cen_4_6->SetBinContent(13,98866);       cen_4_6->SetBinError(13,314.43);
      cen_4_6->SetBinContent(14,99483);       cen_4_6->SetBinError(14,315.409);
      cen_4_6->SetBinContent(15,99474);       cen_4_6->SetBinError(15,315.395);
      cen_4_6->SetBinContent(16,98786);       cen_4_6->SetBinError(16,314.302);
      cen_4_6->SetBinContent(17,98481);       cen_4_6->SetBinError(17,313.817);
      cen_4_6->SetBinContent(18,98584);       cen_4_6->SetBinError(18,313.981);
      cen_4_6->SetBinContent(19,97531);       cen_4_6->SetBinError(19,312.3);
      cen_4_6->SetBinContent(20,97096);       cen_4_6->SetBinError(20,311.602);
      cen_4_6->SetBinContent(21,97336);       cen_4_6->SetBinError(21,311.987);
      cen_4_6->SetBinContent(22,96302);       cen_4_6->SetBinError(22,310.326);
      cen_4_6->SetBinContent(23,94775);       cen_4_6->SetBinError(23,307.855);
      cen_4_6->SetBinContent(24,94466);       cen_4_6->SetBinError(24,307.353);
      cen_4_6->SetBinContent(25,93538);       cen_4_6->SetBinError(25,305.84);
      cen_4_6->SetBinContent(26,93687);       cen_4_6->SetBinError(26,306.083);
      cen_4_6->SetBinContent(27,92062);       cen_4_6->SetBinError(27,303.417);
      cen_4_6->SetBinContent(28,91618);       cen_4_6->SetBinError(28,302.685);
      cen_4_6->SetBinContent(29,91283);       cen_4_6->SetBinError(29,302.131);
      cen_4_6->SetBinContent(30,91370);       cen_4_6->SetBinError(30,302.275);
      cen_4_6->SetBinContent(31,91440);       cen_4_6->SetBinError(31,302.39);
      cen_4_6->SetBinContent(32,90818);       cen_4_6->SetBinError(32,301.36);
      cen_4_6->SetBinContent(33,90675);       cen_4_6->SetBinError(33,301.123);
      cen_4_6->SetBinContent(34,90517);       cen_4_6->SetBinError(34,300.86);
      cen_4_6->SetBinContent(35,91048);       cen_4_6->SetBinError(35,301.742);
      cen_4_6->SetBinContent(36,90588);       cen_4_6->SetBinError(36,300.978);
      cen_4_6->SetBinContent(37,91704);       cen_4_6->SetBinError(37,302.827);
      cen_4_6->SetBinContent(38,91461);       cen_4_6->SetBinError(38,302.425);
      cen_4_6->SetBinContent(39,92626);       cen_4_6->SetBinError(39,304.345);
      cen_4_6->SetBinContent(40,93072);       cen_4_6->SetBinError(40,305.077);
                                              cen_4_6->SetEntries(3.79953e+06);

      TH1F *cen_6_8 = new TH1F("cen_6_8","cen_6_8",40,-1.5708,1.5708);
      cen_6_8->SetBinContent(1,89711);        cen_6_8->SetBinError(1,299.518);
      cen_6_8->SetBinContent(2,90565);        cen_6_8->SetBinError(2,300.94);
      cen_6_8->SetBinContent(3,91592);        cen_6_8->SetBinError(3,302.642);
      cen_6_8->SetBinContent(4,92066);        cen_6_8->SetBinError(4,303.424);
      cen_6_8->SetBinContent(5,93309);        cen_6_8->SetBinError(5,305.465);
      cen_6_8->SetBinContent(6,94230);        cen_6_8->SetBinError(6,306.969);
      cen_6_8->SetBinContent(7,94570);        cen_6_8->SetBinError(7,307.522);
      cen_6_8->SetBinContent(8,95336);        cen_6_8->SetBinError(8,308.765);
      cen_6_8->SetBinContent(9,96106);        cen_6_8->SetBinError(9,310.01);
      cen_6_8->SetBinContent(10,97463);       cen_6_8->SetBinError(10,312.191);
      cen_6_8->SetBinContent(11,97315);       cen_6_8->SetBinError(11,311.954);
      cen_6_8->SetBinContent(12,97228);       cen_6_8->SetBinError(12,311.814);
      cen_6_8->SetBinContent(13,97833);       cen_6_8->SetBinError(13,312.783);
      cen_6_8->SetBinContent(14,97544);       cen_6_8->SetBinError(14,312.32);
      cen_6_8->SetBinContent(15,98113);       cen_6_8->SetBinError(15,313.23);
      cen_6_8->SetBinContent(16,98115);       cen_6_8->SetBinError(16,313.233);
      cen_6_8->SetBinContent(17,97034);       cen_6_8->SetBinError(17,311.503);
      cen_6_8->SetBinContent(18,96951);       cen_6_8->SetBinError(18,311.37);
      cen_6_8->SetBinContent(19,96210);       cen_6_8->SetBinError(19,310.177);
      cen_6_8->SetBinContent(20,96277);       cen_6_8->SetBinError(20,310.285);
      cen_6_8->SetBinContent(21,95433);       cen_6_8->SetBinError(21,308.922);
      cen_6_8->SetBinContent(22,93814);       cen_6_8->SetBinError(22,306.291);
      cen_6_8->SetBinContent(23,92606);       cen_6_8->SetBinError(23,304.312);
      cen_6_8->SetBinContent(24,92304);       cen_6_8->SetBinError(24,303.816);
      cen_6_8->SetBinContent(25,90901);       cen_6_8->SetBinError(25,301.498);
      cen_6_8->SetBinContent(26,90235);       cen_6_8->SetBinError(26,300.391);
      cen_6_8->SetBinContent(27,89045);       cen_6_8->SetBinError(27,298.404);
      cen_6_8->SetBinContent(28,88375);       cen_6_8->SetBinError(28,297.279);
      cen_6_8->SetBinContent(29,88009);       cen_6_8->SetBinError(29,296.663);
      cen_6_8->SetBinContent(30,87783);       cen_6_8->SetBinError(30,296.282);
      cen_6_8->SetBinContent(31,87114);       cen_6_8->SetBinError(31,295.151);
      cen_6_8->SetBinContent(32,87103);       cen_6_8->SetBinError(32,295.132);
      cen_6_8->SetBinContent(33,86535);       cen_6_8->SetBinError(33,294.168);
      cen_6_8->SetBinContent(34,86843);       cen_6_8->SetBinError(34,294.691);
      cen_6_8->SetBinContent(35,87077);       cen_6_8->SetBinError(35,295.088);
      cen_6_8->SetBinContent(36,87961);       cen_6_8->SetBinError(36,296.582);
      cen_6_8->SetBinContent(37,87711);       cen_6_8->SetBinError(37,296.16);
      cen_6_8->SetBinContent(38,88030);       cen_6_8->SetBinError(38,296.699);
      cen_6_8->SetBinContent(39,88614);       cen_6_8->SetBinError(39,297.681);
      cen_6_8->SetBinContent(40,89383);       cen_6_8->SetBinError(40,298.97);
                                              cen_6_8->SetEntries(3.69244e+06);

      TH1F *cen_8_10 = new TH1F("cen_8_10","cen_8_10",40,-1.5708,1.5708);
      cen_8_10->SetBinContent(1,53717);       cen_8_10->SetBinError(1,231.769);
      cen_8_10->SetBinContent(2,54275);       cen_8_10->SetBinError(2,232.97);
      cen_8_10->SetBinContent(3,54896);       cen_8_10->SetBinError(3,234.299);
      cen_8_10->SetBinContent(4,55150);       cen_8_10->SetBinError(4,234.84);
      cen_8_10->SetBinContent(5,55468);       cen_8_10->SetBinError(5,235.516);
      cen_8_10->SetBinContent(6,56135);       cen_8_10->SetBinError(6,236.928);
      cen_8_10->SetBinContent(7,56655);       cen_8_10->SetBinError(7,238.023);
      cen_8_10->SetBinContent(8,56675);       cen_8_10->SetBinError(8,238.065);
      cen_8_10->SetBinContent(9,57552);       cen_8_10->SetBinError(9,239.9);
      cen_8_10->SetBinContent(10,57453);      cen_8_10->SetBinError(10,239.694);
      cen_8_10->SetBinContent(11,57887);      cen_8_10->SetBinError(11,240.597);
      cen_8_10->SetBinContent(12,58079);      cen_8_10->SetBinError(12,240.996);
      cen_8_10->SetBinContent(13,58405);      cen_8_10->SetBinError(13,241.671);
      cen_8_10->SetBinContent(14,58433);      cen_8_10->SetBinError(14,241.729);
      cen_8_10->SetBinContent(15,58389);      cen_8_10->SetBinError(15,241.638);
      cen_8_10->SetBinContent(16,58376);      cen_8_10->SetBinError(16,241.611);
      cen_8_10->SetBinContent(17,58086);      cen_8_10->SetBinError(17,241.01);
      cen_8_10->SetBinContent(18,58042);      cen_8_10->SetBinError(18,240.919);
      cen_8_10->SetBinContent(19,57578);      cen_8_10->SetBinError(19,239.954);
      cen_8_10->SetBinContent(20,57797);      cen_8_10->SetBinError(20,240.41);
      cen_8_10->SetBinContent(21,57334);      cen_8_10->SetBinError(21,239.445);
      cen_8_10->SetBinContent(22,56602);      cen_8_10->SetBinError(22,237.912);
      cen_8_10->SetBinContent(23,56125);      cen_8_10->SetBinError(23,236.907);
      cen_8_10->SetBinContent(24,55528);      cen_8_10->SetBinError(24,235.644);
      cen_8_10->SetBinContent(25,54730);      cen_8_10->SetBinError(25,233.944);
      cen_8_10->SetBinContent(26,54620);      cen_8_10->SetBinError(26,233.709);
      cen_8_10->SetBinContent(27,53905);      cen_8_10->SetBinError(27,232.175);
      cen_8_10->SetBinContent(28,53114);      cen_8_10->SetBinError(28,230.465);
      cen_8_10->SetBinContent(29,52992);      cen_8_10->SetBinError(29,230.2);
      cen_8_10->SetBinContent(30,53162);      cen_8_10->SetBinError(30,230.569);
      cen_8_10->SetBinContent(31,52587);      cen_8_10->SetBinError(31,229.319);
      cen_8_10->SetBinContent(32,52432);      cen_8_10->SetBinError(32,228.98);
      cen_8_10->SetBinContent(33,52357);      cen_8_10->SetBinError(33,228.817);
      cen_8_10->SetBinContent(34,52614);      cen_8_10->SetBinError(34,229.377);
      cen_8_10->SetBinContent(35,52605);      cen_8_10->SetBinError(35,229.358);
      cen_8_10->SetBinContent(36,52974);      cen_8_10->SetBinError(36,230.161);
      cen_8_10->SetBinContent(37,52967);      cen_8_10->SetBinError(37,230.146);
      cen_8_10->SetBinContent(38,53153);      cen_8_10->SetBinError(38,230.549);
      cen_8_10->SetBinContent(39,53320);      cen_8_10->SetBinError(39,230.911);
      cen_8_10->SetBinContent(40,53744);      cen_8_10->SetBinError(40,231.828);
                                              cen_8_10->SetEntries(2.21591e+06);

      // pass the calibration data to the task
      jetTask->SetEventPlaneWeights(cen_0_2, 0);
      jetTask->SetEventPlaneWeights(cen_2_4, 1);
      jetTask->SetEventPlaneWeights(cen_4_6, 2);
      jetTask->SetEventPlaneWeights(cen_6_8, 3);
      jetTask->SetEventPlaneWeights(cen_8_10, 4);

      // a standard centrality binning is necessary (which corresponds to the binning of the ep weights)
      // so we overwrite possible existing binnings
      Double_t cw[] = {0., 2., 4., 6., 8., 10., 30., 50., 90.};
      jetTask->SetCentralityClasses(new TArrayD(sizeof(cw)/sizeof(cw[0]), cw));
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  jetTask->SetVzRange(-10., 10.);
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname+="_PWGJE";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
          TList::Class(),AliAnalysisManager::kOutputContainer,
	  Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  if(baseClassHistos) {
      contname+="BaseClassHistograms";
      AliAnalysisDataContainer *coutputBase = mgr->CreateContainer(contname.Data(), 
	  TList::Class(),AliAnalysisManager::kOutputContainer,
	  Form("%s", AliAnalysisManager::GetCommonFileName()));
      mgr->ConnectOutput(jetTask, 1, coutputBase);
  }
  mgr->ConnectOutput(jetTask, (baseClassHistos) ? 2 : 1, coutput1);

  switch (runMode) {
      case AliAnalysisTaskJetV2::kLocal : {
          gStyle->SetOptFit(1);
          AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("good_fits_%s", name.Data()), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
          AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("bad_fits_%s", name.Data()),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							     Form("%s", AliAnalysisManager::GetCommonFileName()));
          mgr->ConnectOutput (jetTask, (baseClassHistos) ? 3 : 2, coutput2);
          mgr->ConnectOutput (jetTask, (baseClassHistos) ? 4 : 3, coutput3);
      } break;
      default: break;
  }
  return jetTask;
}
