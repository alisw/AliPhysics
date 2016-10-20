/*
  Append QA TPC  metadata decribing tree structure,  and annotating branche variables. 
  Partialy inspired by CSS (https://de.wikipedia.org/wiki/Cascading_Style_Sheets) but not full functionality implemented
  
  Instead of the manual entering all metadatas, arrays of regular expressions defining  classes of variables
  --  kineVariablesClass
  --  qaVariableClass
  --  statClassClass
  --  categoryClass

  varName.class=<statClass>+<kineVariableClass>+<qaVariableClass>+<categoryClass>+<Class:>

  Examples: 
     dcar_negC_1_Err.class:          TPC Err FitGX DCAr   CSide Neg
     dcar_negC_1_Err.Title:          #sigma x_{G} DCAr   C side Q<0
     meanMIPvsSector..class:         TPC Mean dEdx   Sector Class:TVectorT<double>
     meanMIPvsSector.Title:         mean dEdx   Sector

  Work in progress: include other types of metadata
        -- automatic html metadata
        -- markers and colors for predefined variable (charge, side?)
        -- automatic parser of aliases ?

  Usage and debugging of metadata setting:

     1.) Metadata can be setup invoking macro:
         AliExternalInfo info;
         TTree * tree = info.GetTree("QA.TPC","LHC15o","cpass1_pass1","QA.TPC;QA.TRD;QA.TOF;QA.ITS;QA.EVS;Logbook.detector");             
	 .x $ALICE_ROOT/../src/STAT/Macros/qatpcAddMetadata.C+(tree,4)

            
     2.) Macro can be executed automatically if proper configuation file leaded AliExternalInfo.cfg - see line:
         QA.TPC.metadataMacro $ALICE_ROOT/STAT/Macros/qatpcAddMetadata.C+

     3.) Printing all metadata:
         a.) 
	 AliTreePlayer::selectMetadata(tree, "[class==\"\"]",0)->Print();

     3.) Example query particular info:
           AliTreePlayer::selectMetadata(tree, "[class==\"DCAR&&!ERR&&!CHI2\"]",0)->Print();
	   AliTreePlayer::selectMetadata(tree, "[class==\"DCAR&&ERR\"]",0)->Print();
           AliTreePlayer::selectMetadata(tree, "[class==\"DCAZ\"]",0)->Print();
	   
	   AliTreePlayer::selectMetadata(tree,"[class==DCAr&&!Pos&&!Neg",0).Print();
	   // alternative
	   ((TObjArray*)(tree->GetUserInfo()->FindObject("metaTable"))).Print("",TPRegexp("dca.*class"))


*/

#include "TTree.h"
#include "TPRegexp.h"
#include "TList.h"
#include "AliTreePlayer.h"
#include "TStatToolkit.h"

void qatpcAddMetadata(TTree*tree, Int_t verbose){
  //
  // Set metadata infomation 
  //
  if (tree==NULL) {
    ::Error("qatpcAddMetadata","Start processing. Emtpy tree");
    return;
  }
  ::Info("qatpcAddMetadata","Start processing Tree %s",tree->GetName());
  TObjArray * branches=tree->GetListOfBranches();
  // Clasigication of variables
  //   regular expression to defined automaticaly some variables following naming conventions - used to define classes/Axis/legends 
  //           default description
  //   regular expression used can be tested on site https://regex101.com/
  //           hovewer root (perl)  regular expression looks to be in some cases different - in some case double escape had to be used
  //           e.g to math c.  expression c\\. has to be used
  //   variables
  const Int_t nKine=11;
  const TString kineVariableClass[nKine]={"X", "Y","Z", "Phi", "Theta", "Pt", "QOverPt", "FitPhi","FitGX", "FitGY",""};
  const TString kineVariableAxisTitle[nKine]={"x(cm)", "y(cm)","z(cm)", "#phi", "#Theta", "p_{T}(Gev/c)", "q/p_{T}(c/GeV)", "#phi","x_{G}", "y_{G}", " "};
  const TString kineVariableLegend[nKine]={"x", "y","z", "#phi", "#Theta", "p_{T}", "q/p_{T}", "#phi","x_{G}", "y_{G}", " "};
  const TString kineVariableTitle[nKine]={"x", "y","z", "#phi", "#Theta", "p_{T}", "q/p_{T}", "#phi","x_{G}", "y_{G}", " "};
  TPRegexp regKineVariables[nKine];
  regKineVariables[0]=TPRegexp("^x|x$");               // X - varaible begining on X
  regKineVariables[1]=TPRegexp("(^y|^infoy|y$)");      // Y - 
  regKineVariables[2]=TPRegexp("(^z|^infoz|z$)");      // Z
  regKineVariables[3]=TPRegexp("(phi|infophi)");              // phi
  regKineVariables[4]=TPRegexp("(theta|lambda|infolambda)");   // theta
  regKineVariables[5]=TPRegexp("(^pt|^infopt|meanpt|deltapt)");               // pt
  regKineVariables[6]=TPRegexp("qoverpt");          // qoverPt
  regKineVariables[7]=TPRegexp("(_0^|_-_)");        // fit Phi
  regKineVariables[8]=TPRegexp("(_1^|_1_)");        // fit GlobalX      
  regKineVariables[9]=TPRegexp("(_2^|_2_)");        // fit GlobalY
  regKineVariables[10]=TPRegexp("(dummy)");        // dummy
  //
  //   QA variables  
  const Int_t nQA=8;
  const TString qaVariableClass[nQA]={"Ncl", "Eff","dEdx", "DCAr", "DCAz","Occ","Attach","Electron"};
  const TString qaVariableAxisTitle[nQA]={"Ncl(#)", "Eff(unit)","dEdx(MIP/50)", "DCAr(cm)", "DCAz(cm)", "Occupancy","Attach","Electron"};
  const TString qaVariableLegend[nQA]={"Ncl", "Eff","dEdx", "DCAr", "DCAz","Occ.","Attach","e-"};
  const TString qaVariableTitle[nQA]={"Ncl", "Eff","dEdx", "DCAr", "DCAz", "Occupancy","Attach","e-"};
  TPRegexp regQAVariable[nQA];
  regQAVariable[0]=TPRegexp("ncl");               // Ncl
  regQAVariable[1]=TPRegexp("tpcitsmatch");       // TPC->ITS eff
  regQAVariable[2]=TPRegexp("mip");               // dEdx
  regQAVariable[3]=TPRegexp("dcar");              // dcar
  regQAVariable[4]=TPRegexp("dcaz");              // dcaz
  regQAVariable[5]=TPRegexp("occ");               // occupancy
  regQAVariable[6]=TPRegexp("attach");            // attachement
  regQAVariable[7]=TPRegexp("(electron|ele$)");    // Electron char.

  //
  // statistic
  //
  const Int_t nStat=10;
  const TString statClass[nStat]={"Constrain", "Mean","Delta","Median", "RMS","Pull", "Err","Chi2", "StatInfo[]","FitInfo[]"};
  const TString statAxisTitle[nStat]={"Constrain", "mean","#Delta","med.", "rms","pull", "#sigma"," #chi2","stat[]","fit[]"};
  const TString statTitle[nStat]={"Constrain", "mean","#Delta","med.", "rms","pull", "#sigma"," #chi2","stat[]","fit[]" };
  TPRegexp regStat[nStat];
  regStat[0]=TPRegexp("constrain");        // constrain
  regStat[1]=TPRegexp("^mean");            // mean
  regStat[2]=TPRegexp("^delta");           // delta
  regStat[3]=TPRegexp("^median");          // median variable
  regStat[4]=TPRegexp("(^rms|resolution)");// rms resolution        
  regStat[5]=TPRegexp("pull");             // pull
  regStat[6]=TPRegexp("err");              // error
  regStat[7]=TPRegexp("chi2");             // chi2
  regStat[8]=TPRegexp("^info");            // stat info array 
  regStat[9]=TPRegexp("^fit");             // fit info  array
  //
  // 
  const Int_t nCategory=9;
  const TString categoryClass[nCategory]={"ASide","CSide","Pos", "Neg", "ROC", "IROC", "OROC", "Sector",  "HighPt"};
  const TString categoryLegend[nCategory]={"A side","C side","Q>0", "Q<0","ROC", "IROC", "OROC", "Sector", "high p_{T}"};
  const TString categoryTitle[nCategory]={"A side","C side","Q>0", "Q<0","ROC", "IROC", "OROC", "Sector", "high p_{T}"};
  // TString categoryAxisTitle[nCategory] - empty
  TPRegexp regCategory[nCategory];
  regCategory[0]=TPRegexp("(a$|a.?side|a\\.$|dcarap|a_(0|1|2))");   // is A side - options( a at the beigininng, a ... side, a. at the end)  
  regCategory[1]=TPRegexp("(c$|c.?side|c\\.$|dcarcp|c_(0|1|2))");   // is C side
  regCategory[2]=TPRegexp("(_pos|pos$)");           // QA for positive charge particle
  regCategory[3]=TPRegexp("(_neg|neg$)");           // positive charge particle
  regCategory[4]=TPRegexp("^roc");                  // ROC
  regCategory[5]=TPRegexp("iroc");                  // IROC
  regCategory[6]=TPRegexp("oroc");                  // OROC
  regCategory[7]=TPRegexp("sector");                // sector
  regCategory[8]=TPRegexp("highpt");                // hightPt QA
  


  for (Int_t ibr=0; ibr<branches->GetEntriesFast(); ibr++){
    TBranch * branch = (TBranch*)branches->At(ibr);
    TString brClass="";
    TString brAxisTitle="";
    TString brTitle="";
    TString brLegend="";
    //
    TString brNameCase(branches->At(ibr)->GetName());
    brNameCase.ToLower();
    //
    // define met
    brClass="TPC";
    brAxisTitle="";
    // stat
    for (Int_t ivar=0; ivar<nStat; ivar++) if  (brNameCase.Contains( regStat[ivar])) { 
      brClass+=" "+statClass[ivar];
      brTitle+=statTitle[ivar];
    }
    // kine variables
    for (Int_t ivar=0; ivar<nKine; ivar++) if  (brNameCase.Contains( regKineVariables[ivar])) {
      brClass+=" "+kineVariableClass[ivar];      
      brAxisTitle+=" "+kineVariableAxisTitle[ivar];
      brTitle+=" "+kineVariableTitle[ivar];
      brLegend+=" "+kineVariableLegend[ivar];
    }
    // QA variables
    for (Int_t ivar=0; ivar<nQA; ivar++) if  (brNameCase.Contains( regQAVariable[ivar])) {
      brClass+=" "+qaVariableClass[ivar];
      brAxisTitle+=" "+qaVariableAxisTitle[ivar];
      brTitle+=" "+qaVariableTitle[ivar];
      brLegend+=" "+qaVariableLegend[ivar];
    }
    // category
    for (Int_t ivar=0; ivar<nCategory; ivar++) if  (brNameCase.Contains(regCategory[ivar])) {
      brClass+=" "+categoryClass[ivar];
      brLegend+=" "+categoryLegend[ivar];
      brTitle+=" "+categoryTitle[ivar];
    }
    if (branch!=NULL && branch->GetClassName()!=NULL && strlen(branch->GetClassName())>0){
      brClass+=" Class:";
      brClass+=branch->GetClassName();
    }
    TPRegexp emptyB("^\ *");
    emptyB.Substitute(brClass,"");
    emptyB.Substitute(brTitle,"");
    emptyB.Substitute(brAxisTitle,"");
    emptyB.Substitute(brLegend,"");

    //
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Description",branches->At(ibr)->GetName()).Data(),
			      TString::Format("TPC standard QA variables.  Class %s", brClass.Data()).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("%s.class",branches->At(ibr)->GetName()).Data(),brClass.Data());    
    TStatToolkit::AddMetadata(tree,TString::Format("%s.AxisTitle",branches->At(ibr)->GetName()).Data(),brAxisTitle.Data());    
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Title",branches->At(ibr)->GetName()).Data(),brTitle.Data());    
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Legend",branches->At(ibr)->GetName()).Data(),brLegend.Data());    
    if (verbose&4) printf("Class %s: \t%s\n", branches->At(ibr)->GetName(),brClass.Data());
    if (verbose&8) printf("Axis title %s: \t%s\n", branches->At(ibr)->GetName(),brAxisTitle.Data());
    if (verbose&16) printf("Title %s: \t%s\n", branches->At(ibr)->GetName(),brTitle.Data());
    if (verbose&32) printf("Legend %s: \t%s\n", branches->At(ibr)->GetName(),brLegend.Data());

    
  }

  // Fill Based and  custom metadata
  //
  // Index
  TStatToolkit::AddMetadata(tree,"run.class","Base Index");
  TStatToolkit::AddMetadata(tree,"run.Title","run");
  TStatToolkit::AddMetadata(tree,"run.AxisTitle","run");
  //
  TList * mlist = (TList*)(tree->GetUserInfo()->FindObject("metaTable"));
  mlist->Sort();
  if ((verbose&1)>0){
    mlist->Print();
  }
  if ((verbose&2)>0){
    AliTreePlayer::selectMetadata(tree, "[class==\"\"]",0)->Print();
  }
  ::Info("qatpcAddMetadata","End");

}
