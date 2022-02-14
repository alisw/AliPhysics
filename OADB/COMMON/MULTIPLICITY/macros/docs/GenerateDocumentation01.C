#include "AliMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultSelection.h"

void GenerateDocumentation01() {
    //This macro auto-generates documentation (twiki-style)
    //based on all OADBs currently deployed in the AliPhysics
    //directory specified here:

    TString lPathToOADBs = "$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/";
    cout<<"Expanding path name..."<<endl;
    gSystem->ExpandPathName( lPathToOADBs );
    cout<< "Expanded to: "<<lPathToOADBs.Data()<<endl;

    //Start: determine list of available OADBs
    TSystemDirectory dir(lPathToOADBs.Data(), lPathToOADBs.Data());
    TList *files = dir.GetListOfFiles();

    //sorting...
    files->Sort();

    //Number of automatically found OADB files
    Long_t lNOADBs = 0;

    //Output text file
    TString lOutputDoc = "doc.txt";
    ofstream lDocs;
    lDocs.open ( lOutputDoc.Data() );

    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            TString lProdName = fname;
            TString lCalibRuns = "";
            TString lTime = "";
            
            TString lCommitMessage = "";
            TString lEvSels = "";
            TString lEstims = ""; 
	    TString lEstimPerRun = ""; 
            lProdName.ReplaceAll("OADB-","");
            lProdName.ReplaceAll(".root","");
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                cout<<"======================================================"<<endl;
                cout<<" OADB entry for period ["<<lProdName.Data() <<"] found!"<<endl;
                cout<<" -> Opening OADB and checking for information ... "<<endl;

                //Mod time
                //Long_t id,size,flags,mt;
                //gSystem->GetPathInfo(Form("%s%s",lPathToOADBs.Data(),fname.Data()),&id,&size,&flags,&mt);
                //TTimeStamp st(mt);
                //cout<<"Date : "<<st.GetDate()<<endl;
                //lTime.Append(Form("%i",st.GetDate()));
                
                //Modification time, as per git command:
                // git log -1 --format="%ai" -- OADB-blablabla.root
                lTime = gSystem->GetFromPipe(Form("cd %s; git log -1 --format=\"%%ai\" -- OADB-%s.root", lPathToOADBs.Data(), lProdName.Data()));
                lCommitMessage = gSystem->GetFromPipe(Form("cd %s; git log -1 --oneline -- OADB-%s.root", lPathToOADBs.Data(), lProdName.Data()));
                cout<<"======================================================"<<endl;
                cout<<lTime.Data()<<endl;
                cout<<"Commit Message Text: "<<lCommitMessage.Data()<<endl;
                //Will now append lCommitMessage to lTime
                lCommitMessage.Prepend(" <noautolink><element title=\"Last Commit Message: ");
                lCommitMessage.Append("\"> %ICONURL{help}% </element></noautolink>");
                
                lTime.Append(lCommitMessage);
                
                cout<<" lTime string, final: "<<lTime.Data()<<endl;
                
                cout<<"======================================================"<<endl;
                
                TFile * f = new TFile (Form("%s%s",lPathToOADBs.Data(),fname.Data()));
                AliOADBContainer * oadbContMS = (AliOADBContainer*) f->Get("MultSel");
                cout<<" ---> contains this many runs: "<<oadbContMS->GetNumberOfEntries()<<endl;
                AliOADBMultSelection * oadbMultSelection = 0x0;
                for(Long_t ir=0; ir<oadbContMS->GetNumberOfEntries(); ir++) {

                    //Acquire the MultSelection Object for this run and get the specific
                    //calibration information (Anchor point and anchor percentile) for each estimator
                    oadbMultSelection = (AliOADBMultSelection* )oadbContMS->GetObjectByIndex(ir);

                    //Get estimator info
                    lCalibRuns.Append("<element title=\"");
                    for(Int_t iEst = 0; iEst<oadbMultSelection->GetMultSelection()->GetNEstimators(); iEst++) {
                        TString def     = oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetName();
                        Float_t lAP     = oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetAnchorPoint();
                        Float_t lAPPerc = oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetAnchorPercentile();

                        lCalibRuns.Append(Form("%s",def.Data())); //Name
			if(oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetUseAnchor()) { 
                        lCalibRuns.Append(Form(": Raw Val %.3f, Percentile: %.3f",lAP,lAPPerc)); //Characteristics
			} else { 
			  lCalibRuns.Append(": Unanchored"); //Characteristics
			}
                        if ( iEst != oadbMultSelection->GetMultSelection()->GetNEstimators()-1 )
                            lCalibRuns.Append("&#10;"); //Name
                    }
                    lCalibRuns.Append(Form("\"> %i</element>",oadbContMS->LowerLimit(ir)));
                    if( ir != oadbContMS->GetNumberOfEntries()-1 ) lCalibRuns.Append(", ");
                }
                cout<<" ---> Detected calibrated runs: "<<endl;
                cout<<lCalibRuns.Data()<<endl;
                cout<<" ---> Acquiring example object by index... "<<endl;
                oadbMultSelection = (AliOADBMultSelection* )oadbContMS->GetObjectByIndex(0);
                cout<<" -> Adding: Event selection information ..."<<endl;
                if( oadbMultSelection ) {
                    lEvSels.Append(Form("Vertex position &#124;z&#124; < %.1f cm",oadbMultSelection->GetEventCuts()->GetVzCut() ) );
                    if( oadbMultSelection->GetEventCuts()->GetTriggerCut() ) lEvSels.Append(" + Physics selected (Minimum Bias: kMB/kINT7)");
                    if( oadbMultSelection->GetEventCuts()->GetINELgtZEROCut() ) lEvSels.Append(" + INEL&gt;0 with SPD tracklets");
                    if( oadbMultSelection->GetEventCuts()->GetTrackletsVsClustersCut() ) lEvSels.Append(" + tracklets versus clusters cut");
                    if( oadbMultSelection->GetEventCuts()->GetRejectPileupInMultBinsCut() ) lEvSels.Append(" + Pileup in Mult bins rejection");
                    if( oadbMultSelection->GetEventCuts()->GetVertexConsistencyCut() ) lEvSels.Append(" + SPD and Track Vertex Consistency OK");
                    if( oadbMultSelection->GetEventCuts()->GetNonZeroNContribs() ) lEvSels.Append(" + At least 1 contributor to PV");
                    oadbMultSelection->GetMultSelection()->PrintInfo();
                    for(Int_t iEst = 0; iEst<oadbMultSelection->GetMultSelection()->GetNEstimators(); iEst++) {
                        //syntax = <element title="test test test mom mom mom"> %ICONURL{help}% </div>
                        TString def = oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetDefinition();
                        cout<<"Def = "<<def.Data()<<endl;
			//Replace &gt; and &lt;s for html code
			def.ReplaceAll("<","&lt;");
			def.ReplaceAll(">","&gt;");
                        lEstims.Append("<element title=\"");
                        lEstims.Append(Form("%s",def.Data()));
                        lEstims.Append(Form("\"> !%s</element>",oadbMultSelection->GetMultSelection()->GetEstimator(iEst)->GetName()));
                        if(iEst!=oadbMultSelection->GetMultSelection()->GetNEstimators()-1) lEstims.Append(", ");
                    }
                }
                cout<<" ---> Event Selection: "<<lEvSels.Data()<<endl;
                //Appending event selection as a tooltip next to the dataset name
                lProdName.Prepend("*");
                lProdName.Append("*");
                lProdName.Append(Form(" <element title=\"Event Selection: %s\"",lEvSels.Data() ));
                lProdName.Append("> %ICONURL{help}% </element>");

                cout<<" ---> Estimator string: "<<lEstims.Data()<<endl;
                cout<<" -> Adding: estimator definition information ..."<<endl;
                cout<<"======================================================"<<endl;
                //Print this in this format:
                //| LHCxxx   |   RUNLIST   |   xx/xx/xxx   |   ESTIMATORS   |
                lDocs<<"| "<<lProdName.Data()<<" | "<<lCalibRuns.Data()<<" |   "<<lTime.Data()<<"   | "<<lEstims.Data()<<" |"<<endl;

            }
        }
    }

}
