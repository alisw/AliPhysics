//
// This macro is a part of "Alice PPR fast flow analysis package"
// 
// The macro draw event plane resolution for a given 
// multiplicity and V2. As source data it is using
// "flowData" NTuple in "flowPicoEvent.root" file.
// Data source structure have to be created 
// by AliFlowCreate.C and filed by AliFlowReconstruction.C
// 
// INPUT PARAMETERS:
// type: type of plot:
// 0 - Event Plane Resolution
// 1 - Sub-Events Correlation
// 2 - V2 distribution 
//
// info: if !0 display histogram info
//
// if (plotName != 0) two files are created in "plots" directory  
// plotName.eps and plotMame.gif
//
// Sylwester Radomski, GSI
// mail: S.Radomski@gsi
// 23. Oct. 2002
//


AliFlowDist(int multiplicity, float trueV2, int type = 0, int info = 0,
              const char *plotName = 0) {
		  
   gROOT->SetStyle("Plain");
   gStyle->SetOptTitle(0);

   if (info) gStyle->SetOptStat(1110);
   else gStyle->SetOptStat(0);


   if (type > 2 || type < 0) {
     ::Error("AliFlowDrawSummary","Wrong Type [0-1] : %d", type);
     return;
   }
   
   const char *what[3] = {"Psi - truePsi>> htemp(40,-40,40)",
			    "PsiA - PsiB >> htemp(40,-40,40)",
			  "V2 >> htemp(40,0.0, 0.14)"  };
   
   const char *xTitle[3] = {"Event Plane Resolution [deg]",
			      "Sub-Events Correlation [deg]",
			    "V_{2}"};

   const char *where = "Mult == %d && (trueV2 > %f) && (trueV2 < %f)";  
   
   TFile *inFile = new TFile("flowPicoEvent.root");
   TNtuple *data = (TNtuple*) inFile->Get("flowData");		  
    
   char buff[80];
   sprintf(buff, where, multiplicity, trueV2-0.001, trueV2+0.001);
   
   data->Draw(what[type], buff, "");
    
   TH1 *h = (TH1*)gPad->FindObject("htemp");
   h->GetXaxis()->SetTitle(xTitle[type]);
   
   gPad->SetTicks(1,1);      
   //gPad->Update();	   	  

   if (plotName) {
     
     char buffer[60];
     
     sprintf(buffer, "plots/%s.eps", plotName);
     gPad->Print(buffer, "eps");

     sprintf(buffer, "plots/%s.gif", plotName);
     gPad->Print(buffer, "gif");
   }
}
