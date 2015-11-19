
Int_t qaConfig(TTree* tree, TString* returnStrings)
{
  Float_t entryFrac=0.8, nsigmaOutlier=6., nsigmaWarning=3., epsilon=1.0e-6;

  //
  // specify all variables for which the status aliases shall be defined.
  // only for these variables the lines in the trending plots can be computed and plotted.
  //
  TString sTrendVars="meanTPCncl;meanTPCnclF;meanMIP;resolutionMIP;MIPattachSlopeA;MIPattachSlopeC;";
  sTrendVars+=";meanMIPele;resolutionMIPele;electroMIPSeparation;meanVertX;meanVertY;meanVertZ;meanMultPos;meanMultNeg;";
  sTrendVars+=";tpcItsMatchA;tpcItsMatchC;tpcItsMatchHighPtA;tpcItsMatchHighPtC;lambdaPull;ptPull;yPull;zPull;";
  sTrendVars+=";tpcConstrainPhiA;tpcConstrainPhiC;deltaPt;";
  sTrendVars+=";offsetdRA;offsetdZA;offsetdRC;offsetdZC;dcarAP0;dcarAP1;dcarCP0;dcarCP1;";
  sTrendVars+=";dcar_posA_0;dcar_posA_1;dcar_posA_2;dcaz_posA_0;dcaz_posA_1;dcaz_posA_2;";
  sTrendVars+=";dcar_posC_0;dcar_posC_1;dcar_posC_2;dcaz_posC_0;dcaz_posC_1;dcaz_posC_2;";
  sTrendVars+=";dcar_negA_0;dcar_negA_1;dcar_negA_2;dcaz_negA_0;dcaz_negA_1;dcaz_negA_2;";
  sTrendVars+=";dcar_negC_0;dcar_negC_1;dcar_negC_2;dcaz_negC_0;dcaz_negC_1;dcaz_negC_2;";

  //
  // combined variables
  // name them '..._combN' with N being the number of combined variables!
  //
  tree->SetAlias("meanMult_comb2"         , "((meanMultPos+meanMultNeg)/2.)");
  tree->SetAlias("tpcItsMatch_comb4"      , "((tpcItsMatchA+tpcItsMatchC+tpcItsMatchHighPtA+tpcItsMatchHighPtC)/4)"); // mean of all 4.
  tree->SetAlias("itsTpcPulls_comb4"      , "(TMath::Sqrt(lambdaPull**2+ptPull**2+yPull**2+zPull**2))");
  tree->SetAlias("tpcConstrainPhi_comb2"  , "(TMath::Sqrt(tpcConstrainPhiA**2+tpcConstrainPhiC**2))"); // sqrt of quadr. sum ok because it's a bias.
  tree->SetAlias("offsetd_comb4"          , "(TMath::Sqrt(offsetdRA**2+offsetdZA**2+offsetdRC**2+offsetdZC**2))");
  tree->SetAlias("dcarFitpar_comb4"       , "((dcarAP0+dcarAP1+dcarCP0+dcarCP1)/4)"); // mean of 4.
  tree->SetAlias("dcar0_comb4"            , "(TMath::Sqrt(dcar_posA_0**2+dcar_posC_0**2+dcar_negA_0**2+dcar_negC_0**2))");
  tree->SetAlias("dcar1_comb4"            , "(TMath::Sqrt(dcar_posA_1**2+dcar_posC_1**2+dcar_negA_1**2+dcar_negC_1**2))");
  tree->SetAlias("dcar2_comb4"            , "(TMath::Sqrt(dcar_posA_2**2+dcar_posC_2**2+dcar_negA_2**2+dcar_negC_2**2))");
  tree->SetAlias("dcaz0_comb4"            , "(TMath::Sqrt(dcaz_posA_0**2+dcaz_posC_0**2+dcaz_negA_0**2+dcaz_negC_0**2))");
  tree->SetAlias("dcaz1_comb4"            , "(TMath::Sqrt(dcaz_posA_1**2+dcaz_posC_1**2+dcaz_negA_1**2+dcaz_negC_1**2))");
  tree->SetAlias("dcaz2_comb4"            , "(TMath::Sqrt(dcaz_posA_2**2+dcaz_posC_2**2+dcaz_negA_2**2+dcaz_negC_2**2))");
  tree->SetAlias("MIPattachSlope_comb2"   , "((MIPattachSlopeA+MIPattachSlopeC*(-1))/2)");
  tree->SetAlias("PIDSepPow_comb2"        , "(meanMIPele-meanMIP)/(0.5*(resolutionMIP*meanMIP+resolutionMIPele*meanMIPele))");

  //
  // add all combined variables to sTrendVars!
  // only then the statistics aliases will be computed for them as well.
  //
  sTrendVars+=";meanMult_comb2;tpcItsMatch_comb4;itsTpcPulls_comb4;tpcConstrainPhi_comb2;";
  sTrendVars+=";offsetd_comb4;dcarFitpar_comb4;dcar0_comb4;dcar1_comb4;dcar2_comb4;dcaz0_comb4;dcaz1_comb4;dcaz2_comb4;";
  sTrendVars+=";MIPattachSlope_comb2;PIDSepPow_comb2;";

  //
  // specify criterion to mark runs with enough statistics.
  // these runs are the basis for the computation of the outlier criteria.
  // -> robust mean and rms are computed from given entry fraction (EF) of these runs!
  //
  tree->SetAlias("statisticOK", "(meanTPCncl>0)");

  //
  // creation of aliases for Outliers, Warnings, PhysicsAcceptable ...
  // upper and lower limits needed separately to retrieve the numerical values for the line positions.
  // for that reason the aliases are expected to be computable expressions.
  //
  TObjArray* oaTrendVars = sTrendVars.Tokenize(",;");
  for (Int_t vari=0; vari<oaTrendVars->GetEntriesFast(); vari++)
  {
    TString sVar( oaTrendVars->At(vari)->GetName() );

    // outliers, warnings and robust mean are set for all variables identically.
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_OutlierMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaOutlier, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_OutlierMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaOutlier, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_WarningMin:(MeanEF-%f*RMSEF-%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_WarningMax:(MeanEF+%f*RMSEF+%f):%f", nsigmaWarning, epsilon, entryFrac));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_RobustMean:(MeanEF+0):%f", entryFrac));

    // physics acceptable should be set for each type of variable individually.
    // some sets of variables are already set here, the rest should be set appropriately after the loop!
    Float_t combfac=1.;
    if (sVar.Contains("_comb")) {
      TString last = sVar(sVar.Last('b')+1, sVar.Length());
      combfac = TMath::Sqrt( atoi(last.Data()) );
    }
    if (sVar.Contains("dca")) {
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f+0)", 0.2*combfac)); // 2 mm around mean
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f+0)", 0.2*combfac)); // 2 mm around mean
    }
    else if (sVar.Contains("Pull")) {
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f+0)", 1.0*combfac)); // check them!
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f+0)", 1.0*combfac));
    }
    else { // other variables set to +- 5% of the mean as default to avoid crashes:
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.05*combfac, entryFrac));
      TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.05*combfac, entryFrac));
    }

  }

  //
  // all aliases can just be overwritten here...
  // PhysAcc should be set to relative or absolute values appropriately!
  //
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.05, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanTPCncl",    "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.05, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanMIP",       "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.01, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanMIP",       "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.01, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "resolutionMIP", "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.10, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "resolutionMIP", "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.10, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "MIPattachSlopeA",  "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.005, entryFrac)); // p0 + p1 * tan(theta) // showing p1
  TStatToolkit::SetStatusAlias(tree, "MIPattachSlopeA",  "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.005, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "MIPattachSlopeC",  "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.005, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "MIPattachSlopeC",  "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.005, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanMIPele",       "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.01, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanMIPele",       "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.01, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "resolutionMIPele", "statisticOK", Form("varname_PhysAccMin:(MeanEF-%f*MeanEF):%f", 0.10, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "resolutionMIPele", "statisticOK", Form("varname_PhysAccMax:(MeanEF+%f*MeanEF):%f", 0.10, entryFrac));
  TStatToolkit::SetStatusAlias(tree, "meanVertZ",     "statisticOK", Form("varname_PhysAccMin:(%f+0)",-1.0 ));
  TStatToolkit::SetStatusAlias(tree, "meanVertZ",     "statisticOK", Form("varname_PhysAccMax:(%f+0)", 1.0 ));
//  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "statisticOK", Form("varname_PhysAccMin:(%f+0)", 0.7 )); // TODO: set to +-3% after multiplicity correction (see wiki)
//  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchA",  "statisticOK", Form("varname_PhysAccMax:(%f+0)", 1.0 ));
//  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchC",  "statisticOK", Form("varname_PhysAccMin:(%f+0)", 0.7 ));
//  TStatToolkit::SetStatusAlias(tree, "tpcItsMatchC",  "statisticOK", Form("varname_PhysAccMax:(%f+0)", 1.0 ));

  //
  // now the actual criteria for each variable are set automatically.
  //
  for (Int_t vari=0; vari<oaTrendVars->GetEntriesFast(); vari++)
  {
    TString sVar( oaTrendVars->At(vari)->GetName() );
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "", Form("varname_Warning:(varname>varname_WarningMax||varname<varname_WarningMin)"));
    TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "", Form("varname_PhysAcc:(varname>varname_PhysAccMin&&varname<varname_PhysAccMax)"));
  }

  // TStatToolkit::SetStatusAlias(tree, sVar.Data(),    "", Form("varname_Outlier:(varname>varname_OutlierMax||varname<varname_OutlierMin)"));


  //
  // COMBINED CRITERIA
  //
  // Adding new subcriteria to the following combined criteria may change the descision based on the individual ones.
  // e.g.: Assume a meanMIP warning that was vetoed by a meanMIP physacc before. It will not be vetoed anymore, if the physacc of MIPattach is very small.
  //

  // combined MIP quality
  // tree->SetAlias("MIPquality_Outlier", "(meanMIP_Outlier||resolutionMIP_Outlier||MIPattachSlopeA_Outlier||MIPattachSlopeC_Outlier||PIDSepPow_comb2_Outlier)");
  // tree->SetAlias("MIPquality_Warning", "(meanMIP_Warning||resolutionMIP_Warning||MIPattachSlopeA_Warning||MIPattachSlopeC_Warning||PIDSepPow_comb2_Warning)");
  // tree->SetAlias("MIPquality_PhysAcc", "(meanMIP_PhysAcc&&resolutionMIP_PhysAcc&&MIPattachSlopeA_PhysAcc&&MIPattachSlopeC_PhysAcc&&PIDSepPow_comb2_PhysAcc)");
  Float_t MeanMip_limitmax = 53.;
  Float_t MeanMip_limitmin = 47.;
  tree->SetAlias("MIPquality_Outlier", Form("(meanMIP_Outlier||resolutionMIP_Outlier||PIDSepPow_comb2_Outlier||meanMIP<%f||meanMIP>%f)", MeanMip_limitmin, MeanMip_limitmax));
  tree->SetAlias("MIPquality_Warning", "(meanMIP_Warning||resolutionMIP_Warning||PIDSepPow_comb2_Warning)");
  tree->SetAlias("MIPquality_PhysAcc", Form("(meanMIP_PhysAcc&&resolutionMIP_PhysAcc&&PIDSepPow_comb2_PhysAcc&&meanMIP>%f&&meanMIP<%f)", MeanMip_limitmin, MeanMip_limitmax));

  // combined matching efficiency
  Float_t TPC_ITS_limit = 0.5;
  tree->SetAlias("tpcItsMatch_Outlier", Form("(tpcItsMatchA_Outlier||tpcItsMatchC_Outlier||tpcItsMatchHighPtA_Outlier||tpcItsMatchHighPtC_Outlier||tpcItsMatchA<%f||tpcItsMatchC<%f)", TPC_ITS_limit, TPC_ITS_limit));
  tree->SetAlias("tpcItsMatch_Warning", "(tpcItsMatchA_Warning||tpcItsMatchC_Warning||tpcItsMatchHighPtA_Warning||tpcItsMatchHighPtC_Warning)");
  tree->SetAlias("tpcItsMatch_PhysAcc", Form("(tpcItsMatchA_PhysAcc&&tpcItsMatchC_PhysAcc&&tpcItsMatchHighPtA_PhysAcc&&tpcItsMatchHighPtC_PhysAcc&&tpcItsMatchA>%f&&tpcItsMatchC>%f)", TPC_ITS_limit, TPC_ITS_limit));

  // combined matching quality (lambdaPull ptPull yPull zPull)
  tree->SetAlias("itsTpcPulls_Outlier", "(lambdaPull_Outlier||ptPull_Outlier||yPull_Outlier||zPull_Outlier)");
  tree->SetAlias("itsTpcPulls_Warning", "(lambdaPull_Warning||ptPull_Warning||yPull_Warning||zPull_Warning)");
  tree->SetAlias("itsTpcPulls_PhysAcc", "(lambdaPull_PhysAcc&&ptPull_PhysAcc&&yPull_PhysAcc&&zPull_PhysAcc)");

  // combined DCA R and Z
  tree->SetAlias("dcar0_Outlier", "(dcar_posA_0_Outlier||dcar_posC_0_Outlier||dcar_negA_0_Outlier||dcar_negC_0_Outlier)");
  tree->SetAlias("dcar1_Outlier", "(dcar_posA_1_Outlier||dcar_posC_1_Outlier||dcar_negA_1_Outlier||dcar_negC_1_Outlier)");
  tree->SetAlias("dcar2_Outlier", "(dcar_posA_2_Outlier||dcar_posC_2_Outlier||dcar_negA_2_Outlier||dcar_negC_2_Outlier)");
  tree->SetAlias("dcar_Outlier" , "(dcar0_Outlier||dcar1_Outlier||dcar2_Outlier)");
  tree->SetAlias("dcaz0_Outlier", "(dcaz_posA_0_Outlier||dcaz_posC_0_Outlier||dcaz_negA_0_Outlier||dcaz_negC_0_Outlier)");
  tree->SetAlias("dcaz1_Outlier", "(dcaz_posA_1_Outlier||dcaz_posC_1_Outlier||dcaz_negA_1_Outlier||dcaz_negC_1_Outlier)");
  tree->SetAlias("dcaz2_Outlier", "(dcaz_posA_2_Outlier||dcaz_posC_2_Outlier||dcaz_negA_2_Outlier||dcaz_negC_2_Outlier)");
  tree->SetAlias("dcaz_Outlier" , "(dcaz0_Outlier||dcaz1_Outlier||dcaz2_Outlier)");

  tree->SetAlias("dcar0_Warning", "(dcar_posA_0_Warning||dcar_posC_0_Warning||dcar_negA_0_Warning||dcar_negC_0_Warning)");
  tree->SetAlias("dcar1_Warning", "(dcar_posA_1_Warning||dcar_posC_1_Warning||dcar_negA_1_Warning||dcar_negC_1_Warning)");
  tree->SetAlias("dcar2_Warning", "(dcar_posA_2_Warning||dcar_posC_2_Warning||dcar_negA_2_Warning||dcar_negC_2_Warning)");
  tree->SetAlias("dcar_Warning" , "(dcar0_Warning||dcar1_Warning||dcar2_Warning)");
  tree->SetAlias("dcaz0_Warning", "(dcaz_posA_0_Warning||dcaz_posC_0_Warning||dcaz_negA_0_Warning||dcaz_negC_0_Warning)");
  tree->SetAlias("dcaz1_Warning", "(dcaz_posA_1_Warning||dcaz_posC_1_Warning||dcaz_negA_1_Warning||dcaz_negC_1_Warning)");
  tree->SetAlias("dcaz2_Warning", "(dcaz_posA_2_Warning||dcaz_posC_2_Warning||dcaz_negA_2_Warning||dcaz_negC_2_Warning)");
  tree->SetAlias("dcaz_Warning" , "(dcaz0_Warning||dcaz1_Warning||dcaz2_Warning)");

  tree->SetAlias("dcar0_PhysAcc", "(dcar_posA_0_PhysAcc&&dcar_posC_0_PhysAcc&&dcar_negA_0_PhysAcc&&dcar_negC_0_PhysAcc)");
  tree->SetAlias("dcar1_PhysAcc", "(dcar_posA_1_PhysAcc&&dcar_posC_1_PhysAcc&&dcar_negA_1_PhysAcc&&dcar_negC_1_PhysAcc)");
  tree->SetAlias("dcar2_PhysAcc", "(dcar_posA_2_PhysAcc&&dcar_posC_2_PhysAcc&&dcar_negA_2_PhysAcc&&dcar_negC_2_PhysAcc)");
  tree->SetAlias("dcar_PhysAcc" , "(dcar0_PhysAcc&&dcar1_PhysAcc&&dcar2_PhysAcc)");
  tree->SetAlias("dcaz0_PhysAcc", "(dcaz_posA_0_PhysAcc&&dcaz_posC_0_PhysAcc&&dcaz_negA_0_PhysAcc&&dcaz_negC_0_PhysAcc)");
  tree->SetAlias("dcaz1_PhysAcc", "(dcaz_posA_1_PhysAcc&&dcaz_posC_1_PhysAcc&&dcaz_negA_1_PhysAcc&&dcaz_negC_1_PhysAcc)");
  tree->SetAlias("dcaz2_PhysAcc", "(dcaz_posA_2_PhysAcc&&dcaz_posC_2_PhysAcc&&dcaz_negA_2_PhysAcc&&dcaz_negC_2_PhysAcc)");
  tree->SetAlias("dcaz_PhysAcc" , "(dcaz0_PhysAcc&&dcaz1_PhysAcc&&dcaz2_PhysAcc)");


  //
  // specify the variables which shall be included in the status bar.
  // give them meaningful names for the y axis of the status bar!
  //
  TString sStatusbarVars ("MIPquality;dcaz;dcar;tpcItsMatch;itsTpcPulls;meanTPCncl;global");
  TString sStatusbarNames("dE/dx (MIP);v_drift (DCAZ);space p.(DCAR);TPC-ITS m.eff.;ITS-TPC m.qual.;mean TPC Ncl;global result");

  //
  // global descision is made automatically from the other entries in the status bar.
  //
  TString sVarsNoGlobal = sStatusbarVars;
  sVarsNoGlobal.Remove(sVarsNoGlobal.Last(';'), sVarsNoGlobal.Length()); // remove last item (';global')
  // global PhysAcc
  TString sGlobalCriterion = sVarsNoGlobal;
  sGlobalCriterion.ReplaceAll(";","_PhysAcc&&"); sGlobalCriterion.Append("_PhysAcc");
  cout << "global PhysAcc:  " << sGlobalCriterion.Data() << endl;
  tree->SetAlias("global_PhysAcc" , sGlobalCriterion.Data());
  // global Outlier
  sGlobalCriterion="";
  TObjArray* oaVarsNoGlobal = sVarsNoGlobal.Tokenize(";");
  for (Int_t vari=0; vari<oaVarsNoGlobal->GetEntriesFast(); vari++)
  {
    TString sVar( oaVarsNoGlobal->At(vari)->GetName() );
    sGlobalCriterion+=Form("(%s_Outlier&&!%s_PhysAcc)", sVar.Data(), sVar.Data());
    if (vari<oaVarsNoGlobal->GetEntriesFast()-1) sGlobalCriterion+="||";
  }
  cout << "global Outlier:  " << sGlobalCriterion.Data() << endl;
  tree->SetAlias("global_Outlier" , sGlobalCriterion.Data());
  // global Warning
  sGlobalCriterion.ReplaceAll("_Outlier","_Warning");
  cout << "global Warning:  " << sGlobalCriterion.Data() << endl;
  tree->SetAlias("global_Warning" , sGlobalCriterion.Data());

  //
  // set Boolean criteria to be checked for the markers in the status bar:
  // separate them by colon ":"   (1) = true --> first marker always plotted
  //
  TString sCriteria("(1):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show vetos: (varname_PhysAcc&&varname_Warning)

  returnStrings[0] = sStatusbarVars;
  returnStrings[1] = sStatusbarNames;
  returnStrings[2] = sCriteria;

  return 1;
}
