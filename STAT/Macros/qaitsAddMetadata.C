/*
  

  Append QA ITS  metadata decribing tree structure,  and annotating branche variables. 
  Partialy inspired by CSS (https://de.wikipedia.org/wiki/Cascading_Style_Sheets) but not full functionality implemented

  authors: 
  first version  - marian.ivanov@cern.ch


  
  Instead of the manual entering all metadatas, arrays of regular expressions defining  classes of variables
  --  kineVariablesClass
  --  qaVariableClass
  --  statClassClass
  --  categoryClass

  varName.class=<statClass>+<kineVariableClass>+<qaVariableClass>+<categoryClass>+<Class:>

  Examples: 
   
  Work in progress: include other types of metadata
        -- automatic html metadata
        -- markers and colors for predefined variable (charge, side?)
        -- automatic parser of aliases ? 
        -- automatic paerser of the error variables for varaible annotation

  Usage and debugging of metadata setting:

     1.) Metadata can be setup invoking macro:
         AliExternalInfo info;
         TTree * tree = info.GetTree("QA.ITS","LHC15o","pass1","QA.TPC;QA.TRD;QA.TOF;QA.EVS;QA.ITS;Logbook.detector");             
	 .x $ALICE_ROOT/../src/STAT/Macros/qaitsAddMetadata.C+(tree,4)

            
     2.) Macro can be executed automatically if proper configuation file leaded AliExternalInfo.cfg - see line:
         QA.ITS.metadataMacro $ALICE_ROOT/STAT/Macros/qaitsAddMetadata.C+

     3.) Printing all metadata:
         a.) 
	 AliTreePlayer::selectMetadata(tree, "[class==\"\"]",0)->Print();

     3.) Example query particular info:
	   
	   AliTreePlayer::selectMetadata(tree, "[class==\"spd&&eff&&!Err\"]",0)->Print();  // efficiency spd nor error
	   OBJ: TNamed    EffSPDPt02.class        ITS    eff        spd pt02      
	   OBJ: TNamed    EffSPDPt1.class ITS    eff        spd pt1      
	   OBJ: TNamed    EffSPDPt10.class        ITS    eff        spd pt10      
	   OBJ: TNamed    EffoneSPDPt02.class     ITS    effone        spd pt02      
	   OBJ: TNamed    EffoneSPDPt1.class      ITS    effone        spd pt1      
	   OBJ: TNamed    EffoneSPDPt10.class     ITS    effone        spd pt10      
		
	   AliTreePlayer::selectMetadata(tree, "[class==\"eff&&!Err&&pt10\"]",0)->Print();
	   Collection name='TObjArray', class='TObjArray', size=581
	   OBJ: TNamed    Eff2Pt10.class  ITS    eff2        pt10      
	   OBJ: TNamed    Eff3Pt10.class  ITS    eff3        pt10      
	   OBJ: TNamed    Eff4Pt10.class  ITS    eff4        pt10      
	   OBJ: TNamed    Eff5Pt10.class  ITS    eff5        pt10      
	   OBJ: TNamed    Eff6Pt10.class  ITS    eff6        pt10      
	   OBJ: TNamed    EffSPDPt10.class        ITS    eff        spd pt10      
	   OBJ: TNamed    EffTOTPt10.class        ITS    eff        pt10      
	   OBJ: TNamed    EffoneSPDPt10.class     ITS    effone        spd pt10      
	   ..
	   

*/
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TTree.h"
#include "TPRegexp.h"
#include "TList.h"
#include "AliTreePlayer.h"
#include "TStatToolkit.h"
#endif

void qaitsAddMetadata(TTree*tree, Int_t verbose){
  //
  // Set metadata infomation 
  //
  if (tree==NULL) {
    ::Error("qaitsAddMetadata","Start processing. Emtpy tree");
    return;
  }
  ::Info("qaitsAddMetadata","Start processing Tree %s",tree->GetName());
  TObjArray * branches=tree->GetListOfBranches();
  // Clasigication of variables
  //   regular expression to defined automaticaly some variables following naming conventions - used to define classes/Axis/legends 
  //           default description
  //   regular expression used can be tested on site https://regex101.com/
  //           hovewer root (perl)  regular expression looks to be in some cases different - in some case double escape had to be used
  //           e.g to math c.  expression c\\. has to be used
  //   variables
  const TString kineVariableClass[11]={"X", "Y","Z", "Phi", "Theta", "Pt", "QOverPt"};
  const TString kineVariableAxisTitle[11]={"x(cm)", "y(cm)","z(cm)", "#phi", "#Theta", "p_{T}", "q/p_{T}(c/GeV)"};
  const TString kineVariableLegend[11]={"x", "y","z", "#phi", "#Theta", "p_{T}", "q/p_{T}"};
  const TString kineVariableTitle[11]={"x", "y","z", "#phi", "#Theta", "p_{T}", "q/p_{T}"};
  TPRegexp regKineVariables[11];
  regKineVariables[0]=TPRegexp("^x|x$");               // X - varaible begining on X
  regKineVariables[1]=TPRegexp("(^y|^infoy|y$)");      // Y - 
  regKineVariables[2]=TPRegexp("(^z|^infoz|z$)");      // Z
  regKineVariables[3]=TPRegexp("(phi|infophi)");              // phi
  regKineVariables[4]=TPRegexp("(theta|lambda|infolambda)");   // theta
  regKineVariables[5]=TPRegexp("(^pt)");           // pt
  regKineVariables[6]=TPRegexp("qoverpt");          // qoverPt
  //
  //   QA variables  
  const TString qaVariableClass[11]={"Frac", "$","dEdx","ChargeRatio"};
  const TString qaVariableAxisTitle[11]={"Frac","$","dEdx","ChargeRatio"};
  const TString qaVariableLegend[11]={"Frac","$","dEdx","ChargeRatio"};
  const TString qaVariableTitle[11]={"Frac","$","dEdx","ChargeRatio"};
  TPRegexp regQAVariable[11];
  regQAVariable[0]=TPRegexp("frac");              // Frac
  regQAVariable[1]=TPRegexp("(eff[0-9]|effone|^eff)");          // eff layer?
  regQAVariable[2]=TPRegexp("(mpv|dedx)");         // dEdx 
  regQAVariable[3]=TPRegexp("chargeratio");        // ChargeRatio

 
  //
  // statistic
  //
  const TString statClass[10]={"Constrain", "Mean","Delta","Median", "RMS","Pull", "Err","Chi2", "StatInfo[]","FitInfo[]"};
  const TString statAxisTitle[10]={"Constrain", "mean","#Delta","med.", "rms","pull", "#sigma"," #chi2","stat[]","fit[]"};
  const TString statTitle[10]={"Constrain", "mean","#Delta","med.", "rms","pull", "#sigma"," #chi2","stat[]","fit[]" };
  TPRegexp regStat[10];
  regStat[0]=TPRegexp("constrain");        // constrain
  regStat[1]=TPRegexp("^mean");            // mean
  regStat[2]=TPRegexp("^delta");           // delta
  regStat[3]=TPRegexp("^median");          // median variable
  regStat[4]=TPRegexp("(^rms|resolution|sigma)");// rms resolution        
  regStat[5]=TPRegexp("pull");             // pull
  regStat[6]=TPRegexp("err");              // error
  regStat[7]=TPRegexp("chi2");             // chi2
  regStat[8]=TPRegexp("^info");            // stat info array 
  regStat[9]=TPRegexp("^fit");             // fit info  array
  //
  // 
  const TString categoryClass[10]={"Vertex","$1","$1"};
  const TString categoryLegend[10]={"Vertex","$1","$1"};
  const TString categoryTitle[10]={"Vertex","$1","$1"};
  TPRegexp regCategory[10];   // proper parsing of layer numbers to be added
  regCategory[0]=TPRegexp("(vertex|vtx)");             //
  regCategory[1]=TPRegexp("s(p|d|d)d[0-2]?");    // reg.exp match silical layers+number 
  regCategory[2]=TPRegexp("pt[0-9]+");         // pt bin
  


  for (Int_t ibr=0; ibr<branches->GetEntriesFast(); ibr++){
    TBranch * branch = (TBranch*)branches->At(ibr);
    TString matchClass="";   // class match
    TString brClass="";
    TString brAxisTitle="";
    TString brTitle="";
    TString brLegend="";
    //
    TString brNameCase(branches->At(ibr)->GetName());
    brNameCase.ToLower();
    //
    // define met
    brClass="ITS";
    brAxisTitle="";
    // stat
    for (Int_t ivar=0; ivar<11; ivar++) if  (brNameCase.Contains( regStat[ivar])) { 
      brClass+=" "+statClass[ivar];
      brTitle+=statTitle[ivar];
    }
    // kine variables
    for (Int_t ivar=0; ivar<7; ivar++) if  (brNameCase.Contains( regKineVariables[ivar])) {
      brClass+=" "+kineVariableClass[ivar];      
      brAxisTitle+=" "+kineVariableAxisTitle[ivar];
      brTitle+=" "+kineVariableTitle[ivar];
      brLegend+=" "+kineVariableLegend[ivar];
    }
    // QA variables
    for (Int_t ivar=0; ivar<5; ivar++) if  (brNameCase.Contains( regQAVariable[ivar])) {
      if ( qaVariableClass[ivar].Contains("$")==kFALSE){
	brClass+=" "+ qaVariableClass[ivar];
	brLegend+=" "+ qaVariableLegend[ivar];
	brTitle+=" "+ qaVariableTitle[ivar];
      }else{
	TObjArray *amatch=regQAVariable[ivar].MatchS(brNameCase);
	if (amatch){
	  TString match=amatch->At(0)->GetName();
	  brClass+=" "+match;
	  brLegend+=" "+match;
	  brTitle+=" "+match;
	}
      }
    }
    // category
    for (Int_t ivar=0; ivar<3; ivar++) if  (brNameCase.Contains(regCategory[ivar])) {
      if (categoryClass[ivar].Contains("$")==kFALSE){
	brClass+=" "+categoryClass[ivar];
	brLegend+=" "+categoryLegend[ivar];
	brTitle+=" "+categoryTitle[ivar];
      }else{
	TObjArray *amatch=regCategory[ivar].MatchS(brNameCase);
	if (amatch){
	  TString match=amatch->At(0)->GetName();
	  brClass+=" "+match;
	  brLegend+=" "+match;
	  brTitle+=" "+match;
	}
      }
    }
    if (branch!=NULL && branch->GetClassName()!=NULL && strlen(branch->GetClassName())>0){
      brClass+=" Class:";
      brClass+=branch->GetClassName();
    }
    
    //
    TStatToolkit::AddMetadata(tree,TString::Format("%s.Description",branches->At(ibr)->GetName()).Data(),
			      TString::Format("ITS standard QA variables.  Class %s", brClass.Data()).Data());
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
  ::Info("qaitsAddMetadata","End");

}
