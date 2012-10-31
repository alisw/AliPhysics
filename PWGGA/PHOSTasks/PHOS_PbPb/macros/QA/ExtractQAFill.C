void ExtractQAFill(const TString fillFile="fillFile.txt")
{
  gROOT->LoadMacro("ExtractQA.C");
  ExtractQA(fillFile, "outputQAFill.root");
}
