void ExtractQAFill(const TString fillFileList="fillFileList.txt",
	       const TString fillList="fillList.txt")
{
  gROOT->LoadMacro("ExtractQA.C");
  ExtractQA(fillFileList, fillList, "outputQAFill.root");
}
