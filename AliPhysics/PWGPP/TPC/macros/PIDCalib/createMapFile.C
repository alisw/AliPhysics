Int_t createMapFile(Bool_t specialVersion = kFALSE, Bool_t commitVersion = kFALSE,
                    TString normalisation = "NoNormalisation"/*, TString pathMapPackage = "etaMaps",
                    TString fileNameMapPackage = "TPCetaMaps.root"*/)
{
  gROOT->LoadMacro("addMapToFile.C+");
  
  if (specialVersion) {
    addMapToFile("finalCuts/PbPb/2.76ATeV/LHC11h.pass2/uncorrected/LegoTrain_LuciasRunList/outputCheckShapeEtaTree_2014_10_06__21_13__NewSplines___2014_10_06__10_21.root", normalisation, "LHC11h", 2, kFALSE, kFALSE, "etaMaps/special", "TPCetaMaps_special.root");
    addMapToFile("finalCuts/PbPb/2.76ATeV/LHC11h.pass2/uncorrected/LegoTrain_LuciasRunList/outputCheckShapeEtaTree_2014_10_06__22_21__CorrectedWithMap_outputCheckShapeEtaTree_2014_10_06__21_13__NewSplines___2014_10_06__10_21___2014_10_06__21_28.root", normalisation, "LHC11h", 2, kFALSE, kTRUE, "etaMaps/special", "TPCetaMaps_special.root");
    
    printf("\nCreated maps for special version!\n");
    
    return 0;
  }
  
  // Data
  
  // 2012 pp, 8 TeV, pass2
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12a.pass2/outputCheckShapeEtaTree_2015_06_05__18_52__NewSplines___2015_05_29__12_14.root", normalisation, "LHC12a", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12a.pass2/outputCheckShapeEtaTree_2015_06_05__22_40__CorrectedWithMap_outputCheckShapeEtaTree_2015_06_05__18_52__NewSplines___2015_05_29__12_14___2015_06_05__22_39.root", normalisation, "LHC12a", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12b.pass2/outputCheckShapeEtaTree_2015_06_16__20_16__NewSplines___2015_05_29__12_55.root", normalisation, "LHC12b", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12b.pass2/outputCheckShapeEtaTree_2015_06_16__20_28__CorrectedWithMap_outputCheckShapeEtaTree_2015_06_16__20_16__NewSplines___2015_05_29__12_55___2015_06_16__20_27.root", normalisation, "LHC12b", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12c.pass2/outputCheckShapeEtaTree_2015_05_29__15_17__NewSplines___2015_05_29__15_15.root", normalisation, "LHC12c", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12c.pass2/outputCheckShapeEtaTree_2015_05_29__16_40__CorrectedWithMap_outputCheckShapeEtaTree_2015_05_29__15_17__NewSplines___2015_05_29__15_15___2015_05_29__15_46.root", normalisation, "LHC12c", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12d.pass2/outputCheckShapeEtaTree_2015_05_29__16_50__NewSplines___2015_05_29__16_49.root", normalisation, "LHC12d", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12d.pass2/outputCheckShapeEtaTree_2015_05_29__17_22__CorrectedWithMap_outputCheckShapeEtaTree_2015_05_29__16_50__NewSplines___2015_05_29__16_49___2015_05_29__17_20.root", normalisation, "LHC12d", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12e.pass2/outputCheckShapeEtaTree_2015_06_16__20_22__NewSplines___2015_05_29__16_47.root", normalisation, "LHC12e", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12e.pass2/outputCheckShapeEtaTree_2015_06_16__20_27__CorrectedWithMap_outputCheckShapeEtaTree_2015_06_16__20_22__NewSplines___2015_05_29__16_47___2015_06_16__20_25.root", normalisation, "LHC12e", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12f.pass2/outputCheckShapeEtaTree_2015_05_29__17_22__NewSplines___2015_05_29__17_21.root", normalisation, "LHC12f", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12f.pass2/outputCheckShapeEtaTree_2015_05_29__17_31__CorrectedWithMap_outputCheckShapeEtaTree_2015_05_29__17_22__NewSplines___2015_05_29__17_21___2015_05_29__17_29.root", normalisation, "LHC12f", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12g_pp.pass2/outputCheckShapeEtaTree_2015_05_29__17_31__NewSplines___2015_05_29__17_30.root", normalisation, "LHC12g_pp", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12g_pp.pass2/outputCheckShapeEtaTree_2015_05_29__17_37__CorrectedWithMap_outputCheckShapeEtaTree_2015_05_29__17_31__NewSplines___2015_05_29__17_30___2015_05_29__17_36.root", normalisation, "LHC12g_pp", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12h.pass2/outputCheckShapeEtaTree_2015_06_16__20_58__NewSplines___2015_05_29__17_37.root", normalisation, "LHC12h", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12h.pass2/outputCheckShapeEtaTree_2015_06_16__21_08__CorrectedWithMap_outputCheckShapeEtaTree_2015_06_16__20_58__NewSplines___2015_05_29__17_37___2015_06_16__21_07.root", normalisation, "LHC12h", 2, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12i.pass2/outputCheckShapeEtaTree_2015_06_16__20_34__NewSplines___2015_05_29__17_43.root", normalisation, "LHC12i", 2, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/8TeV/LHC12i.pass2/outputCheckShapeEtaTree_2015_06_16__20_40__CorrectedWithMap_outputCheckShapeEtaTree_2015_06_16__20_34__NewSplines___2015_05_29__17_43___2015_06_16__20_37.root", normalisation, "LHC12i", 2, kFALSE, kTRUE);
  
  
  
  
  
  // Re-reconstruction 2010
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10b.pass4/outputCheckShapeEtaTree_2014_11_17__15_17__NewSplines___2014_11_17__15_13.root", normalisation, "LHC10b", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10b.pass4/outputCheckShapeEtaTree_2014_11_17__15_31__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_17__15_17__NewSplines___2014_11_17__15_13___2014_11_17__15_29.root", normalisation, "LHC10b", 4, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10c.pass4/outputCheckShapeEtaTree_2014_11_17__15_42__NewSplines___2014_11_17__15_39.root", normalisation, "LHC10c", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10c.pass4/outputCheckShapeEtaTree_2014_11_20__10_49__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_17__15_42__NewSplines___2014_11_17__15_39___2014_11_17__15_49.root", normalisation, "LHC10c", 4, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10d.pass4/outputCheckShapeEtaTree_2014_11_17__19_17__NewSplines___2014_11_17__18_40.root", normalisation, "LHC10d", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10d.pass4/outputCheckShapeEtaTree_2014_11_20__11_03__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_17__19_17__NewSplines___2014_11_17__18_40___2014_11_20__10_44.root", normalisation, "LHC10d", 4, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10e.pass4/outputCheckShapeEtaTree_2014_11_20__11_27__NewSplines___2014_11_18__11_09.root", normalisation, "LHC10e", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10e.pass4/outputCheckShapeEtaTree_2014_11_20__11_08__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_18__11_14__NewSplines___2014_11_18__11_09___2014_11_19__11_10.root", normalisation, "LHC10e", 4, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10f.pass4/outputCheckShapeEtaTree_2014_11_19__11_19__NewSplines___2014_11_19__11_14.root", normalisation, "LHC10f", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10f.pass4/outputCheckShapeEtaTree_2014_11_20__11_44__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_19__11_19__NewSplines___2014_11_19__11_14___2014_11_19__14_12.root", normalisation, "LHC10f", 4, kFALSE, kTRUE);
  
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10g.pass4/outputCheckShapeEtaTree_2014_11_20__07_58__NewSplines___2014_11_20__07_52.root", normalisation, "LHC10g", 4, kFALSE, kFALSE);
  addMapToFile("etaMaps/gridOutput/pp/7TeV/LHC10g.pass4/outputCheckShapeEtaTree_2014_11_21__14_27__CorrectedWithMap_outputCheckShapeEtaTree_2014_11_20__07_58__NewSplines___2014_11_20__07_52___2014_11_21__14_26.root", normalisation, "LHC10g", 4, kFALSE, kTRUE);
  
  
  
  
  addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_16__NewSplines___2012_11_30__14_14.root", normalisation, "LHC11a_2.76TeV", 2, kFALSE, kFALSE);
  
  if (commitVersion) {
    addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_18__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_16__NewSplines___2012_11_30__14_14___2012_11_30__14_17.root", normalisation, "LHC11a_2.76TeV", 2, kFALSE, kTRUE);// original
  
    //TODO other maps with worse resolution in MIP region - seems to give worse results addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_03_20__12_41__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_16__NewSplines___2012_11_30__14_14___2012_11_30__14_17.root", normalisation, "LHC11a_2.76TeV", 2, kFALSE, kTRUE);
  }
  else {
    addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_07_07__20_02__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_16__NewSplines___2012_11_30__14_14___2012_11_30__14_17.root", normalisation, "LHC11a_2.76TeV", 2, kFALSE, kTRUE);
  }
  
  
  
  addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass4/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_27__NewSplines___2012_11_30__14_25.root", normalisation, "LHC11a_2.76TeV", 4, kFALSE, kFALSE);
  
  if (commitVersion) {
    addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass4/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_36__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_27__NewSplines___2012_11_30__14_25___2012_11_30__14_28.root", normalisation, "LHC11a_2.76TeV", 4, kFALSE, kTRUE);// original 
  }
  else {
    addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass4/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_07_07__20_23__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_27__NewSplines___2012_11_30__14_25___2012_11_30__14_28.root", normalisation, "LHC11a_2.76TeV", 4, kFALSE, kTRUE);
  
    //TODO other maps with worse resolution in MIP region - seems to give worse results addMapToFile("finalCuts/pp/2.76TeV/LHC11a_without_SDD.pass4/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_03_20__10_58__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_27__NewSplines___2012_11_30__14_25___2012_11_30__14_28.root", normalisation, "LHC11a_2.76TeV", 4, kFALSE, kTRUE);
  }
  
  
  addMapToFile("finalCuts/pp/7TeV/LHC10b.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_43__NewSplines___2012_11_30__13_41.root", normalisation, "LHC10b", 2, kFALSE, kFALSE);
  addMapToFile("finalCuts/pp/7TeV/LHC10b.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_46__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_43__NewSplines___2012_11_30__13_41___2012_11_30__13_44.root", normalisation, "LHC10b", 2, kFALSE, kTRUE);
  
  addMapToFile("finalCuts/pp/7TeV/LHC10c.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_52__NewSplines___2012_11_30__13_50.root", normalisation, "LHC10c", 2, kFALSE, kFALSE);
  addMapToFile("finalCuts/pp/7TeV/LHC10c.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_55__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_52__NewSplines___2012_11_30__13_50___2012_11_30__13_54.root", normalisation, "LHC10c", 2, kFALSE, kTRUE);
  
  addMapToFile("finalCuts/pp/7TeV/LHC10d.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_06__NewSplines___2012_11_30__12_55.root", normalisation, "LHC10d", 2, kFALSE, kFALSE);
  
  if (commitVersion) {
    addMapToFile("finalCuts/pp/7TeV/LHC10d.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_11__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_06__NewSplines___2012_11_30__12_55___2012_11_30__13_09.root", normalisation, "LHC10d", 2, kFALSE, kTRUE);// original
  }
  else {
    addMapToFile("finalCuts/pp/7TeV/LHC10d.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_07_09__14_59__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_06__NewSplines___2012_11_30__12_55___2012_11_30__13_09.root", normalisation, "LHC10d", 2, kFALSE, kTRUE);
  }
  
  
  addMapToFile("finalCuts/pp/7TeV/LHC10e.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_28__NewSplines___2012_11_30__13_19.root", normalisation, "LHC10e", 2, kFALSE, kFALSE);
  
  if (commitVersion) {
    addMapToFile("finalCuts/pp/7TeV/LHC10e.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__13_36__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_28__NewSplines___2012_11_30__13_19___2012_11_30__13_30.root", normalisation, "LHC10e", 2, -kFALSE, kTRUE);// original 
  }
  else {
    addMapToFile("finalCuts/pp/7TeV/LHC10e.pass2/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2013_07_07__17_06__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__13_28__NewSplines___2012_11_30__13_19___2012_11_30__13_30.root", normalisation, "LHC10e", 2, kFALSE, kTRUE);
    
    /*replaced by dzdr maps
     addMapToFile("finalCuts/pp/7TeV/LHC10e.pass2/uncorrected/dzdr/outputCheckShapeEtaTree_2013_04_29__08_50.root", normalisation, "LHC10e", 2, kFALSE, kFALSE);
     addMapToFile("finalCuts/pp/7TeV/LHC10e.pass2/uncorrected/dzdr/outputCheckShapeEtaTree_2013_04_29__08_57__CorrectedWithMap_outputCheckShapeEtaTree_2013_04_29__08_50___2013_04_29__08_53.root", normalisation, "LHC10e", 2, kFALSE, kTRUE);
     */
  }  
  
  addMapToFile("finalCuts/pp/7TeV/LHC11a.pass1/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_03__NewSplines___2012_11_30__14_02.root", normalisation, "LHC11a_7TeV", 1, kFALSE, kFALSE);
  addMapToFile("finalCuts/pp/7TeV/LHC11a.pass1/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_06__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_03__NewSplines___2012_11_30__14_02___2012_11_30__14_05.root", normalisation, "LHC11a_7TeV", 1, kFALSE, kTRUE);
  
  //TODO: Improve PbPb
  if (commitVersion) {
    // Without multiplicity correction:
    addMapToFile("finalCuts/old/uncorrected/0.6GeVcut/PbPb/2.76ATeV/LHC10h.pass2/newV0tree/outputCheckShapeEtaTree_2012_11_06__14_19.root", normalisation, "LHC10h", 2, kFALSE, kFALSE);
    addMapToFile("finalCuts/old/uncorrected/0.6GeVcut/PbPb/2.76ATeV/LHC10h.pass2/newV0tree/outputCheckShapeEtaTree_2012_11_07__07_58__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_06__14_19___2012_11_06__15_03.root", normalisation, "LHC10h", 2, kFALSE, kTRUE);
  }
  else {
    addMapToFile("finalCuts/PbPb/2.76ATeV/LHC10h.pass2/uncorrected/multiplicityCorrUndoneTotalTracks/outputCheckShapeEtaTree_2013_01_29__15_04__NewSplines___2013_01_28__09_09_mult0_1999.root", normalisation, "LHC10h", 2, kFALSE, kFALSE);
    addMapToFile("finalCuts/PbPb/2.76ATeV/LHC10h.pass2/uncorrected/multiplicityCorrUndoneTotalTracks/outputCheckShapeEtaTree_2013_01_31__08_37__multiplicityCorrected_CorrectedWithMap_outputCheckShapeEtaTree_2013_01_29__15_04__NewSplines___2013_01_28__09_09_mult0_1999___2013_01_30__13_56_mult0_1999.root", normalisation, "LHC10h", 2, kFALSE, kTRUE);
  }
  
  
  
  // pPb
  /*
  if (commitVersion) {
    // Without multiplicity correction:
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54.root", normalisation, "LHC13b", 2, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_14__CorrectedWithMap_outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54___2013_02_08__13_11.root", normalisation, "LHC13b", 2, kFALSE, kTRUE);
    // Use same maps for pass1 and pass3
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54.root", normalisation, "LHC13b", 1, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_14__CorrectedWithMap_outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54___2013_02_08__13_11.root", normalisation, "LHC13b", 1, kFALSE, kTRUE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54.root", normalisation, "LHC13b", 3, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass2/uncorrected/highStatistics/outputCheckShapeEtaTree_2013_02_08__13_14__CorrectedWithMap_outputCheckShapeEtaTree_2013_02_08__13_07__NewSplines___2013_02_08__12_54___2013_02_08__13_11.root", normalisation, "LHC13b", 3, kFALSE, kTRUE);
  }
  else {*/
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499.root", normalisation, "LHC13b", 2, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_07_12__14_51__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 2, kFALSE, kTRUE);
    //original addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_33__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 2, kFALSE, kTRUE);
    // Use same maps for pass1 and pass3
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499.root", normalisation, "LHC13b", 1, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_07_12__14_51__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 1, kFALSE, kTRUE);
    //original addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_33__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 1, kFALSE, kTRUE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499.root", normalisation, "LHC13b", 3, kFALSE, kFALSE);
    addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_07_12__14_51__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 3, kFALSE, kTRUE);
    //original addMapToFile("finalCuts/pPb/5.023ATeV/LHC13b.pass3/uncorrected/outputCheckShapeEtaTree_2013_05_17__10_33__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__10_10__mult0_499___2013_05_17__10_13_mult0_499.root", normalisation, "LHC13b", 3, kFALSE, kTRUE);    
  //}
  
  // MC_pp
  addMapToFile("finalCuts/MC_pp/2.76TeV/LHC12f1_merged/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_53__NewSplines___2012_11_30__14_50.root", normalisation, "LHC12f1", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pp/2.76TeV/LHC12f1_merged/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__14_55__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__14_53__NewSplines___2012_11_30__14_50___2012_11_30__14_54.root", normalisation, "LHC12f1", 1, kTRUE, kTRUE);

  addMapToFile("finalCuts/MC_pp/LHC11b_merged/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__15_05__NewSplines___2012_11_30__15_04.root", normalisation, "LHC11b2", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pp/LHC11b_merged/uncorrected/optimisedSplines/outputCheckShapeEtaTree_2012_11_30__15_07__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__15_05__NewSplines___2012_11_30__15_04___2012_11_30__15_06.root", normalisation, "LHC11b2", 1, kTRUE, kTRUE);
  
  addMapToFile("finalCuts/MC_pp/7TeV/LHC10d1/uncorrected/outputCheckShapeEtaTree_2012_11_30__15_14__NewSplines___2012_11_30__15_12.root", normalisation, "LHC10d1", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pp/7TeV/LHC10d1/uncorrected/outputCheckShapeEtaTree_2012_11_30__15_16__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__15_14__NewSplines___2012_11_30__15_12___2012_11_30__15_15.root", normalisation, "LHC10d1", 1, kTRUE, kTRUE);
  
  addMapToFile("finalCuts/MC_pp/7TeV/LHC10f6a/uncorrected/outputCheckShapeEtaTree_2012_11_30__15_25__NewSplines___2012_11_30__15_23.root", normalisation, "LHC10f6a", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pp/7TeV/LHC10f6a/uncorrected/outputCheckShapeEtaTree_2012_11_30__15_29__CorrectedWithMap_outputCheckShapeEtaTree_2012_11_30__15_25__NewSplines___2012_11_30__15_23___2012_11_30__15_27.root", normalisation, "LHC10f6a", 1, kTRUE, kTRUE);
  
  
  // MC_pPb
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC12g1/uncorrected/newSplines/outputCheckShapeEtaTree_2012_12_10__08_40__NewSplines___2012_12_10__08_26.root", normalisation, "LHC12g", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC12g1/uncorrected/newSplines/outputCheckShapeEtaTree_2012_12_10__08_44__CorrectedWithMap_outputCheckShapeEtaTree_2012_12_10__08_40__NewSplines___2012_12_10__08_26___2012_12_10__08_43.root", normalisation, "LHC12g", 1, kTRUE, kTRUE);
  
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_fix/uncorrected/outputCheckShapeEtaTree_2013_04_17__12_44__NewSplines___2013_04_17__12_29.root", normalisation, "LHC13b2_fix", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_fix/uncorrected/outputCheckShapeEtaTree_2013_04_17__13_37__CorrectedWithMap_outputCheckShapeEtaTree_2013_04_17__12_44__NewSplines___2013_04_17__12_29___2013_04_17__13_35.root", normalisation, "LHC13b2_fix", 1, kTRUE, kTRUE);
  
  // TODO Improve with higher statistics
  //addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_efix_p1/uncorrected/outputCheckShapeEtaTree_2013_06_10__11_11__NewSplines___2013_06_10__11_08.root", normalisation, "LHC13b2_fixn1", 1, kTRUE, kFALSE);
  //addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_efix_p1/uncorrected/outputCheckShapeEtaTree_2013_06_10__11_16__CorrectedWithMap_outputCheckShapeEtaTree_2013_06_10__11_11__NewSplines___2013_06_10__11_08___2013_06_10__11_12.root", normalisation, "LHC13b2_fixn1", 1, kTRUE, kTRUE);
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_efix_p1/uncorrected/higherStatistics/outputCheckShapeEtaTree_2013_06_17__08_02__NewSplines___2013_06_17__07_56.root", normalisation, "LHC13b2_fixn1", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_pPb/5.023ATeV/LHC13b2_efix_p1/uncorrected/higherStatistics/outputCheckShapeEtaTree_2013_06_17__08_29__CorrectedWithMap_outputCheckShapeEtaTree_2013_06_17__08_02__NewSplines___2013_06_17__07_56___2013_06_17__08_04.root", normalisation, "LHC13b2_fixn1", 1, kTRUE, kTRUE);
  
  
  //TODO: Improve MC_PbPb -> Take 11a10a for map, since it has better statistics at low p compared to 11a10b, which is the region of the map
  addMapToFile("finalCuts/MC_PbPb/2.76ATeV/LHC11a10a_bis/uncorrected/outputCheckShapeEtaTree_2013_05_17__13_37__NewSplines___2013_05_17__13_31.root", normalisation, "LHC11A10", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_PbPb/2.76ATeV/LHC11a10a_bis/uncorrected/outputCheckShapeEtaTree_2013_05_17__13_44__CorrectedWithMap_outputCheckShapeEtaTree_2013_05_17__13_37__NewSplines___2013_05_17__13_31___2013_05_17__13_41.root", normalisation, "LHC11A10", 1, kTRUE, kTRUE);  
  /*
  addMapToFile("finalCuts/MC_PbPb/2.76ATeV/LHC11a10a_bis/uncorrected/multiplicityCorrUndoneTotalTracks/outputCheckShapeEtaTree_2012_12_20__15_57__NewSplines___2012_12_20__15_54_mult0_1999.root", normalisation, "LHC11A10", 1, kTRUE, kFALSE);
  addMapToFile("finalCuts/MC_PbPb/2.76ATeV/LHC11a10a_bis/uncorrected/multiplicityCorrUndoneTotalTracks/outputCheckShapeEtaTree_2013_01_30__12_04__multiplicityCorrected_CorrectedWithMap_outputCheckShapeEtaTree_2012_12_20__15_57__NewSplines___2012_12_20__15_54_mult0_1999___2013_01_30__11_47_mult0_1999.root", normalisation, "LHC11A10", 1, kTRUE, kTRUE);  
  */
  
  printf("TODO Using new sigma map for 10d (log binning for map creation + new sigma fit w/o restriction)!\n");
  printf("TODO Using new sigma map for 10e (log binning for map creation + new sigma fit w/o restriction)!\n");
  printf("TODO Using new sigma map for 11a.p2 (log binning for map creation + new sigma fit w/o restriction)!\n");
  printf("TODO Using new sigma map for 11a.p4 (log binning for map creation + new sigma fit w/o restriction)!\n");
  printf("TODO Using new sigma map for 13b.p[1-3] (log binning for map creation + new sigma fit w/o restriction)!\n");
  
  if (commitVersion)
    printf("\nCreated maps for commit version!\n");
  else
    printf("\nCreated maps for private trunk version!\n");
    
  return 0;
}