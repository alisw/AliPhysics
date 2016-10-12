/*
  Append logbook metadata decribing tree structure,  and annotating varaibles
  
  Usage:
     1.) Metadata can be setup infoking macro:
         AliExternalInfo info;
         TTree * tree = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;QA.ITS;QA.EVS;Logbook.detector");             
         metadataMacro="$ALICE_ROOT/../src/STAT/Macros/qatpcAddMetadata.C";
         gROOT->ProcessLine(TString::Format(".x %s((TTree*)%p,0);",metadataMacro,tree).Data());
      
     2.) Macro can be executed automatically if proper configuation file leaded AliExternalInfo.cfg
     
     3.) Example query metadata:
           AliTreePlayer::selectMetadata(tree, "[class==\"DCAR&&!ERR&&!CHI2\"]",0)->Print();
	   AliTreePlayer::selectMetadata(tree, "[class==\"DCAR&&ERR\"]",0)->Print();
           AliTreePlayer::selectMetadata(tree, "[class==\"DCAZ\"]",0)->Print();

*/

void qatpcAddMetadata(TTree*tree, Int_t verbose=0){
  //
  // Set metadata infomation 
  //
  if (tree==NULL) {
    ::Error("qatpcAddMetadata","Start processing. Emtpy tree");
    return;
  }
  ::Info("qatpcAddMetadata","Start processing Tree %s",tree->GetName());
  TObjArray * branches=tree->GetListOfBranches();
  //
  // regular expression to defined automaticaly some variables following naming conventions - used to fefine classes and default description
  TPRegexp regError("err$");
  TPRegexp regChi2("chi2$");
  TPRegexp regAside("(A$|A.?side)");   // is A side   
  TPRegexp regCside("(A$|C.?side)");   // is C side 
  TPRegexp regPositive("_pos");        // QA for positive charge particle
  TPRegexp regNegative("_neg");        // positive charge particle
  TPRegexp regHighPt("ighPt");         // hightPt QA

  for (Int_t ibr=0; ibr<branches->GetEntriesFast(); ibr++){
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Description",branches->At(ibr)->GetName()).Data(),
			      TString::Format("TPC standard QA variables ").Data());
    TString brNameCase(branches->At(ibr)->GetName());
    brNameCase.ToLower();
    Bool_t isEff      =  brNameCase.Contains("tpcItsMatch",TString::kIgnoreCase);
    Bool_t isdEdx     =  brNameCase.Contains("MIP",TString::kIgnoreCase);
    Bool_t isdcar     =  brNameCase.Contains("dcar",TString::kIgnoreCase);
    Bool_t isdcaz     =  brNameCase.Contains("dcaz",TString::kIgnoreCase);
    Bool_t isErrorBar =  brNameCase.Contains(regError);
    Bool_t isChi2     =  brNameCase.Contains(regChi2);
    TString sside      =  "";    
    if (brNameCase.Contains(regAside)) sside="Side A";
    if (brNameCase.Contains(regCside)) sside="Side C";
    TString scharge      =  "";    
    if (brNameCase.Contains(regPositive)) scharge="Q>0";
    if (brNameCase.Contains(regNegative)) scharge="Q<0";
    TString shighPt="";
    if (brNameCase.Contains(regHighPt)) shighPt="Q>0";

    TString brClass="TPC";    
    TString braxisTitle="";
    TString brTitle="";
    TString brLegend="";

    if (isdEdx) {
      brClass+=" dEdx";
      TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				TString::Format("dEdx (MIP/50)").Data());            
    }

    if (isdcar)  {
      brClass+=" DCAR";
      if (!isErrorBar  && !isChi2) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("DCAr#phi (cm)").Data());      
      }
      if (isChi2) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("#chi2 DCAr#phi (cm)").Data());      
	brClass+=" CHI2";
      }
      if (isErrorBar) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("#sigma DCAr#phi (cm) fir error").Data()); 
	brClass+=" ERR";
      }
    }

    if (isdcaz)  {
      brClass+=" DCAZ";
      if (!isErrorBar  && !isChi2) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("DCAr#phi (cm)").Data());      
      }
      if (isChi2) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("#chi2 DCAr#phi (cm)").Data());      
	brClass+=" CHI2";
      }
      if (isErrorBar) {
	TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),
				  TString::Format("#sigma DCAr#phi (cm) fir error").Data());      
	brClass+=" ERR";

      }
    }    
    TStatToolkit::AddMetadata(tree,TString::Format("%s.class",branches->At(ibr)->GetName()).Data(),brClass.Data());    
  }

  // Fill Based and  custom metadata
  //
  // Index
  TStatToolkit::AddMetadata(tree,"run.class","Base Logbook Index");
  TStatToolkit::AddMetadata(tree,"meanMIP.class","Base ");
  //
  TList * mlist = (TList*)(tree->GetUserInfo()->FindObject("metaTable"));
  mlist->Sort();
  if ((verbose&1)>0){
    mlist->Print();
  }
  if ((verbose&2)>0){
    AliTreePlayer::selectMetadata(tree, "[class==\"DCAR\"]",0)->Print();
  }
  ::Info("qatpcAddMetadata","End");

}
