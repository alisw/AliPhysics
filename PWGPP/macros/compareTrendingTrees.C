void compareTrendingTrees(const char* file1, const char* file2)
{
  TString minEntries = "1000";
  TString indexBranch = "run";
  TString friendPrefix = "b";

  TChain tree2("trending");
  TChain tree1("trending");

  //tree2.AddFile("/hera/alice/mkrzewic/benchmark_flatdev/output/full_NoAD_NoPHOS_6/QAplots/TPC/2015/LHC15f/cpass1/trending.root");
  //tree1.AddFile("/hera/alice/miranov/alice-tpc-notes/reconstruction/dataProductionPreparation/ATO-240/data/benchmarkFull/testATO240_07_31_15_Gitv5-06-35-228-gbf86468/testATO240_07_31_15_Gitv5-06-35-228-gbf86468/QAplots/TPC/data/2015/LHC15f/cpass1/trending.root");
  tree1.AddFile(file1);
  tree2.AddFile(file2);

  tree1.BuildIndex(indexBranch);
  tree2.BuildIndex(indexBranch);

  tree1.AddFriend(&tree2,friendPrefix.Data());
  friendPrefix.Append(".");

  TObjArray* listOfBranches = tree1.GetListOfBranches();
  TObjArray* listOfBranchesFriend = tree2.GetListOfBranches();
  for (int i=0; i<listOfBranches->GetEntries(); i++)
  {
    TString branchName = listOfBranches->At(i)->GetName();
    TObject* branchInFriendTree = listOfBranchesFriend->FindObject(branchName);
    if (!branchInFriendTree)
    {
      Printf("branch %s not preset in friend tree",branchInFriendTree->GetName());
    }
    if (branchName.BeginsWith("info"))
    {
      TString nentries = branchName+"fElements[0]";
      TString mean = branchName+"fElements[1]";
      TString meanError = branchName+"fElements[2]";
      TString rms = branchName+"fElements[3]";
      TString rmsError = branchName+"fElements[4]";
      
      TCanvas* canvas = new TCanvas();
      canvas->Divide(1,2);

      canvas->cd(1);
      TString mean1 = mean;
      TString mean2 = friendPrefix + mean;
      //TString error1 = rms + "/" + "sqrt(" + nentries + ")";
      //TString error2 = friendPrefix + rms + "/" + "sqrt(" + friendPrefix + nentries + ")";
      TString error1 = meanError;
      TString error2 = friendPrefix+meanError;
      //TString error1 = rms;
      //TString error2 = friendPrefix + rms;
      TString expression = "("+mean1+"-"+mean2+")/sqrt("+error1+"*"+error1+"+"+error2+"*"+error2+")";
      TString filterExpression = indexBranch + "==" + friendPrefix + indexBranch + "&&" + nentries+">"+minEntries+"&&"+friendPrefix+nentries+">"+minEntries;
      Printf("branch: %s", branchName.Data());
      Printf("expression: %s", expression.Data());
      Printf("  filter: %s", filterExpression.Data());
      Printf("");
      tree1.Draw(expression+":"+indexBranch,filterExpression,"*");
      
      canvas->cd(2);
      mean1 = mean;
      mean2 = friendPrefix + mean;
      //error1 = rms + "/" + "sqrt(" + nentries + ")";
      //error2 = friendPrefix + rms + "/" + "sqrt(" + friendPrefix + nentries + ")";
      //error1 = meanError;
      //error2 = friendPrefix+meanError;
      error1 = rms;
      error2 = friendPrefix + rms;
      expression = "("+mean1+"-"+mean2+")/sqrt("+error1+"*"+error1+"+"+error2+"*"+error2+")";
      filterExpression = indexBranch + "==" + friendPrefix + indexBranch + "&&" + nentries+">"+minEntries+"&&"+friendPrefix+nentries+">"+minEntries;
      Printf("branch: %s", branchName.Data());
      Printf("expression: %s", expression.Data());
      Printf("  filter: %s", filterExpression.Data());
      Printf("");
      tree1.Draw(expression+":"+indexBranch,filterExpression,"*");
      Printf("__________________________________________________________________");

      canvas->SaveAs(branchName+"pdf");
    }
  }
  

}

