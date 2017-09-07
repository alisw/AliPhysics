/*
  .x $NOTES/aux/rootlogon.C
  .L $ALICE_PHYSICS/../QA/ATO-177/code/dumpTable.C+
  //fname="/hera/alice/miranov/alice-tpc-notes/QA/ATO-102/data/production0604/output/TPC/data/2012/LHC12c/pass2/trending.root";
  fname="trending.root";
  treeName="tpcQA";
  tableName="runTable";
  dumpTable(fname, treeName)

  initTree(fname,treeName);
  initQuery();
  dumpWebRunTable(tree, tableName);

*/


#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TObjArray.h"
#include "AliXRDPROOFtoolkit.h"
#include "TError.h"
#include "AliTPCPerformanceSummary.h"
#include "AliExternalInfo.h"
#include "TPRegexp.h" 

TTree * tree=0;
TString query="";
TString queryDescriptor="";
TString queryFormat="";
TString queryLink="";

TString year="";
TString pass="";
TString anchprodname="";

void initTree(const char* fname, const char* treeName);
void initQuery();
void dumpWebRunTable(TTree* tree, const char* tableName );
void dumpWebPeriodTable(TTree* tree, const char* tableName, const char* MCper=NULL, const char* anchorper=NULL, const char* anchorpass=NULL);
void dumpWebPeriodTableMCRD(TTree* tree, const char* tableName,const char* MCper=NULL, const char* anchorper=NULL, const char* anchorpass=NULL);
void dumpTable(const char *fname, const char * tableName, const char* MCper=NULL, const char* anchorper=NULL, const char* anchorpass=NULL){
  //
  //
  //
  const char *treeName="tpcQA";
  initTree(fname,"tpcQA");
  if (!tree || tree->GetEntries()==0)  initTree(fname,"trending");
  initQuery();
  dumpWebRunTable(tree, tableName);
  dumpWebPeriodTable(tree, tableName, MCper, anchorper, anchorpass);
  if(anchorper!=NULL){
      dumpWebPeriodTableMCRD(tree, tableName, MCper, anchorper, anchorpass);
  }
}

void initTree(const char* fname, const char* treeName){
  //
  /*
   */
  if (TString(fname).Contains(".list")){
    tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(fname,treeName,0,-1);
  }else{
    TFile * fin = TFile::Open(fname);
    tree=(TTree*)fin->Get(treeName);
  }
}

void initQuery(){
  //
  //
  AliTPCPerformanceSummary::MakeMissingChambersAliases(tree);
  query="";
  queryDescriptor="";
  queryLink="";
  //
  query+="run:";
  queryDescriptor+="Run number:";
  queryFormat+="0:";
  queryLink+="0:";
  //
  query+="runType.GetName():";
  queryDescriptor+="Run type:";
  queryFormat+="-1:";
  queryLink+="0:";
  //
  query+="year:";
  queryDescriptor+="Year:";
  queryFormat+="0:";
  queryLink+="0:";
  //
  query+="period.GetName():";
  queryDescriptor+="Period:";
  queryFormat+="-1:";
  queryLink+="0:";
  //
  query+="pass.GetName():";
  queryDescriptor+="Pass:";
  queryFormat+="-1:";
  queryLink+="0:";
  //
  query+="entriesVertX:";
  queryDescriptor+="Events used:";
  queryFormat+="0:";
  queryLink+="TPC_event_info.png:";
  //
  query+="rawLowCounter75:";
  queryDescriptor+="#ROC RawQA <small>(Q<sub>max</sub><0.75&timesMed.) || (Occ<0.75&timesMed.)</small>  :";
  queryFormat+="0:";
  queryLink+="rawQAInformation.png:";
  //
  query+="ocdbStatusCounter:";
  queryDescriptor+="#ROC OCDB <small>(Non Active Map)</small>:";
  queryFormat+="0:";
  queryLink+="canvasROCStatusOCDB.png:";
  //
  query+="ocdbHVStatusCounter:";
  queryDescriptor+="#ROC Low Voltage <small>(HV disabled)</small>:";
  queryFormat+="0:";
  queryLink+="canvasROCStatusOCDB.png:";
  //
  query+="sectorNclMissing70:";
  queryDescriptor+="#ROC QA <small>(N<sub>cl</sub><0.70&timesMed.)</small>:";
  queryFormat+="0:";
  queryLink+="cluster_in_detail.png:";
  //
  query+="sectorNtrMissing70:";
  queryDescriptor+="#ROC QA <small>(N<sub>tr</sub><0.70&timesMed.)</small>:";
  queryFormat+="0:";
  queryLink+="eta_phi_pt.png:";
  //
  query+="qaClOccupancyCounter60:";
  queryDescriptor+="#ROC Cluster <small>(Occ<0.60&timesMed.)</small>:";
  queryFormat+="0:";
  queryLink+="/cluster_occupancy.png:";
  //
  query+="(MIPquality_Warning+10*MIPquality_Outlier+100*(wrongdEdxSectorCounter5>0)):";
  queryDescriptor+="MIP <small>status</small>:";
  queryFormat+="0:";
  queryLink+="TPC_dEdx_track_info.png:";
  //
  query+="(dcar_Warning+10*dcar_Outlier):";
  queryDescriptor+="DCAr <small>status</small>:";
  queryFormat+="0:";
  queryLink+="dca_and_phi.png:";
  //
  query+="(dcaz_Warning+10*dcaz_Outlier)";
  queryDescriptor+="DCAz <small>status</small>:";
  queryFormat+="0";
  queryLink+="dca_and_phi.png";
  //

}


void dumpWebRunTable(TTree* tree, const char * tableName  ){
  //
  // Dump tree content to the - HTML (DOM) sourced data
  //
  tree->SetEstimate(tree->GetEntries());
  Int_t entries=tree->Draw(query.Data(),"1","goffpara");
  TObjArray *arrayDescriptor = queryDescriptor.Tokenize(":");
  TObjArray *arrayFormat     = queryFormat.Tokenize(":");
  TObjArray *arrayLink       = queryLink.Tokenize(":");
  Int_t ncols=arrayDescriptor->GetEntries();
  if (arrayDescriptor->GetEntries()!=arrayFormat->GetEntries()){
    ::Error("dumpWebRunTable",TString::Format("Mismatch Descriptor(%d):Format(%d)",arrayDescriptor->GetEntries(), arrayFormat->GetEntries() ).Data());
  }
  if (arrayDescriptor->GetEntries()!=arrayLink->GetEntries()){
    ::Error("dumpWebRunTable",TString::Format("Mismatch Descriptor(%d):Link(%d)",arrayDescriptor->GetEntries(), arrayLink->GetEntries() ).Data());
  }
  if (arrayFormat->GetEntries()!=arrayLink->GetEntries()){
    ::Error("dumpWebRunTable",TString::Format("Mismatch Format(%d):Link(%d)",arrayFormat->GetEntries(), arrayLink->GetEntries() ).Data());
  }

  FILE * fp=0;
  fp = fopen ("treeRunTable.inc", "w");
  {
    fprintf(fp,"<table id=\"%s\" class=\"display\" cellspacing=\"0\" width=\"100%\">\n", tableName);
    //
    // 1.) make table header
    //
    fprintf(fp,"\t<thead class=\"header\">\n\t\t<tr>\n");
    for (Int_t iheader=0; iheader<ncols; iheader++){
      fprintf(fp,"\t\t\t<th class=\"tooltip\" data-tooltip=\"\">%s</th>\n",arrayDescriptor->At(iheader)->GetName());
    }
    fprintf(fp,"\t\t</tr>\n\t</thead>\n");
    fprintf(fp,"\t<tfoot>\n\t\t<tr>\n");
    for (Int_t iheader=0; iheader<arrayDescriptor->GetEntries(); iheader++){
      fprintf(fp,"\t\t\t<th>%s</th>\n",arrayDescriptor->At(iheader)->GetName());
    }
    fprintf(fp,"\t\t</tr>\n\t</tfoot>\n");
    //
    // 2.) write table content
    //
    fprintf(fp,"\t<tbody>\n");
    for (Int_t ientry=0; ientry<entries; ientry++){
      fprintf(fp,"\t\t<tr>\n");
      //tree->GetEntry(ientry);
      for (Int_t icol=0; icol<ncols; icol++){
	//char* buffer[20];
	//sprintf(buffer, "\t\t\t<td>\%%s</td>\n","d");
	Int_t width=atoi(arrayFormat->At(icol)->GetName());
	if (width>=0) {
	  if (icol==0) {
	    fprintf(fp,"\t\t\t<td><a href=\"000%.0f/index.html\">%.0f</a> </td>\n", tree->GetVal(icol)[ientry],  tree->GetVal(icol)[ientry]); // run numberspecial - another try
	  }else{
	    TString hlink= arrayLink->At(icol)->GetName();
	    if (hlink.Length()<2){
	      fprintf(fp,"\t\t\t<td>  %.*f  </td>\n",width, tree->GetVal(icol)[ientry]);
	    }else{
	      fprintf(fp,"\t\t\t<td>  <a href=\"000%.0f/%s\">%.*f</a>  </td>\n", tree->GetVal(0)[ientry], hlink.Data(), width, tree->GetVal(icol)[ientry]);
	      //<div class="box"><iframe src="http://www.jquery.com/" width = "500px" height = "500px"></iframe></div>
	      //fprintf(fp,"\t\t\t<td>  <a href=\"000%.0f/%s\"> <div class=\"box\"><iframe src=\"000%.0f/%s\" width = \"500px\" height = \"700px\"></iframe></div> %.*f</a>  </td>\n", tree->GetVal(0)[ientry], hlink.Data(),  tree->GetVal(0)[ientry], hlink.Data(), width, tree->GetVal(icol)[ientry]);
	    }
	  }
	}else{
	  TString strName((char *)tree->GetVar(icol)->EvalObject(ientry));
	  fprintf(fp,"\t\t\t<td>  %s</td>\n",  strName.Data());
	}
      }
      fprintf(fp,"\t\t</tr>\n");
    }
    fprintf(fp,"\t</tbody>\n");
  }
  fprintf(fp,"</table>\n");
  fclose(fp);
}

void dumpWebPeriodTable(TTree* tree, const char * tableName, const char* MCper, const char* anchorper, const char* anchorpass) {
  FILE * fp=0;
  fp = fopen ("treePeriodTable.inc", "w");
  fprintf(fp, "\
  <h2 style=\"margin-bottom:3px;margin-top:3px\">TPC QA trending</h2>\n \
  <p style=\"margin-left:10px;margin-top:0px\">\n \
  <b>Tracking</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"meanTPCncl_vs_run.png\">Mean Number of TPC Clusters</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"> 85%%\" href=\"meanTPCnclF_vs_run.png\">Found Clusters / Findable Clusters</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"pp 7 TeV ~ 6; pPb 5 TeV ~ 18;\" href=\"meanMult_vs_run.png\">Multiplicities of Primary Tracks</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"-0.002 < q/pt < +0.002\" href=\"1overPt_vs_run.png\">Delta q/pt</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"<= 0.1\" href=\"DCAOffset_vs_run.png\">DCA Residuals</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Has to be stable over the period\" href=\"pullPhiConstrain_vs_run.png\">Tracking parameter phi</a><br>\n \
  <b>PID</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"50 +/- 1\" href=\"meanMIP_vs_run.png\">Mean of MIPs</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"85 +/- 2\" href=\"meandEdxele_vs_run.png\">Mean of electron energy loss</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Run1: ~ 0.08, Run2: < Run1\" href=\"resolutionMIP_vs_run.png\">Resolution of MIPs</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Run1: ~ 0.065, Run2: < Run1\" href=\"resolutionMeandEdxEle_vs_run.png\">Resolution of mean electron energy loss</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"35 +/- 3\" href=\"ElectroMIPSeparation_vs_run.png\">Separation between electron and MIPs energy loss</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"SeparationPower_vs_run.png\">Separation Power between electron and MIPs energy loss</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Check with Marian\" href=\"MIPattachSlopes_vs_run.png\">Attachment parameter p1</a><br>\n \
  <b>Vertex</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Has to be stable over the period\" href=\"meanVertX_vs_run.png\">Mean of Vertex X</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"Has to be stable over the period\" href=\"meanVertY_vs_run.png\">Mean of Vertex Y</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"-1.5 < VertexZ < 1.5\" href=\"meanVertZ_vs_run.png\">Mean of Vertex Z</a><br>\n \
  <b>TPC ITS Matching</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"> 85%%\" href=\"TPC-ITS-matching-efficiency_vs_run.png\">TPC-ITS matching efficiency</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"~ 0.4; Has to be stable over the period\" href=\"ITS-TPC-matching-quality_vs_run.png\">TPC-ITS matching quality</a><br>\n \
  <b>DCA Fitting</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcar_fitting_run.png\">DCAr fitting parameters</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcar_0_vs_run.png\">DCAr fitting parameters (0)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcar_1_vs_run.png\">DCAr fitting parameters (1)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcar_2_vs_run.png\">DCAr fitting parameters (2)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcaz_0_vs_run.png\">DCAz fitting parameters (0)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcaz_1_vs_run.png\">DCAz fitting parameters (1)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"dcaz_2_vs_run.png\">DCAz fitting parameters (2)</a><br>\n \
  <b>MISC</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"==18\" href=\"occ_AC_Side_IROC_OROC_vs_run.png\">Chambers with lower gain (occupancy)</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"\" href=\"prodinfo\">Production information</a>\n \
  ");
  if(anchorper==NULL){

      fprintf(fp, "\
                  </p>\n \
                  ");
  }
  else{                 //remove hard coded mcrd_com directory/filename name!?  
  //Get directory of anchor production
  AliExternalInfo i;   

  TString AnchProdNamenPass = i.GetMCPassGuess(TString::Format("%s",MCper));
  TObjArray *subStrL;
  subStrL = TPRegexp("(?:[0-9])[0-9]{1}").MatchS(AnchProdNamenPass);
  year = ((TObjString*)subStrL->At(0))->GetString();
  
  subStrL = TPRegexp("pass(.*)").MatchS(AnchProdNamenPass);
  pass = ((TObjString*)subStrL->At(0))->GetString();
  
  subStrL = TPRegexp("^([^ ]+)").MatchS(AnchProdNamenPass);
  anchprodname = ((TObjString*)subStrL->At(0))->GetString();
  
  fprintf(fp, "\
  <br>\n <b>MCRD</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"==18\" href=\"mcrd_com/index_mcrd.html\">MC-RD comparison</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"==18\" href=\"../../../../data/20%s/%s/%s/index.html\">RD QA</a><br>\n \
                  </p>\n \
                  ",year.Data(),anchprodname.Data(),pass.Data());      
  
  }
   
}


void dumpWebPeriodTableMCRD(TTree* tree, const char * tableName,const char* MCper, const char* anchorper, const char* anchorpass) {
  FILE * fp=0;
  fp = fopen ("treePeriodTableMCRD.inc", "w");
  fprintf(fp, "\
  <h2 style=\"margin-bottom:3px;margin-top:3px\">TPC QA trending: MC - RD comparison</h2>\n \
  <p style=\"margin-left:10px;margin-top:0px\">\n \
  <b>MC - RD comparison</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"matchingTPC-ITSEffe.png\">TPC-ITS Matching Eff.</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"matchingTPC-ITSEffe_1.png\">TPC-ITS Matching Eff1</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"matchingTPC-ITSEffe_2.png\">TPC-ITS Matching Eff2</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"matchingTPC-ITSEffe3.png\">TPC-ITS Matching Eff3</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"meanTPCNclRatiMCtoAnchor.png\">meanTPCNclRatiMCtoAnchor</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"meanTPCNclMCtoAnchor.png\">meanTPCNclMCtoAnchor</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"rmsDCAAt1pt0MCtoAnchor.png\">rmsDCAAt1pt0MCtoAnchor</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\">> 120\" href=\"rmsDCAMultSpartMCtoAnchor.png\">rmsDCAMultSpartMCtoAnchor</a><br>\n \
  <br>\n <b>MC</b><br>\n \
  <a class=\"tooltip\" data-tooltip=\"==18\" href=\"../index.html\">MC QA</a><br>\n \
  <a class=\"tooltip\" data-tooltip=\"==18\" href=\"../../../../../data/20%s/%s/%s/index.html\">RD QA</a><br>\n \
                  </p>\n \
                  ",year.Data(),anchprodname.Data(),pass.Data());      
  
  }
   

