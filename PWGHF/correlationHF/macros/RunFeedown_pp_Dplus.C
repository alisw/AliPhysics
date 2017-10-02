//*******************************************************************************/
//  MACRO TO PERFORM THE FD SUBTRACTION FOR D+-hadron correlations in pp
//   Jitendra.Kumar (jitendra.kumar@cern.ch)  
//
//*******************************************************************************/
TString inputfc = "./Dplus/HFPtSpectrum_pp.root"; // input fprompt
TString templatedir = "./Templates_pp/"; // template path
TString inputcorrelationDir = "./Input_Plots_pp/";// directory where input files are stored
TString inputfileroot="1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated";
TString fdsubtrmacrodir="";
TString strSystemFDtempl="none";
void SetFDmacroDirectory(TString macrodir){
  fdsubtrmacrodir=macrodir;
}
void SetFpromptInputFile(TString inputfcfile){
  inputfc=inputfcfile;
}
void SetTemplateDir(TString templdite){
  templatedir=templdite;
}
void SetDirectoryInputFiles(TString inputdir){
  inputcorrelationDir=inputdir;
}
void SetInputFileNameRoot(TString fileinputroot){
  inputfileroot=fileinputroot;
}

void RunFeedown_pp_Dplus(){
    GetEnvelopeForEachV2();
}
void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}

//_____________________________________________________________
void GetEnvelopeForEachV2(){
    

    Int_t collsyst = 0; // 0 is pp, 1 is p-Pb
    gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
    SetSystemStringForTemplateFDnames(strSystemFDtempl.Data());
      
    TString inputcorrelation = ""; // input data file (not needed here)
    TString outputfilename = ""; //  (not needed here)
    
    //    Double_t Dpt[] = {3,5,8,16}; // set D meson pt bins
    //    Double_t hadpt[]    = {0.3, 0.5, 0.3, 1};//03-05, 05-1, 03-1, 1-99
    //    Double_t hadptMax[] = {0.5, 1.0, 1.0, 99.0};//03-05, 05-1, 03-1, 1-99
    

    Double_t Dpt[] = {3,5,8,16}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1}; // set associated tracks pt bins (lower) //03-05, 05-1, 03-1, 1-99
    Double_t hadptMaxInput[] = {99.0,1.0,99.0}; // set associated tracks pt bins (upper) //03-05, 05-1, 03-1, 1-99
    Double_t hadptMax[] = {99.0,1.0,99.0}; // set associated tracks pt bins (upper) //03-05, 05-1, 03-1, 1-99
    Double_t purities[3];
    if(collsyst ==0){
      purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973; //03-05, 05-1, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; //03-05, 05-1, 03-1, 1-99
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");

    }
    else if(collsyst ==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");

    }
    else cout << " ! Purity Value is not system wise " << endl;

    
    Int_t systmode = 3; // mode for evaluation of the envelope - standard is set to 3
    for(Int_t iDpt = 0; iDpt <3; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <3; ihadpt++){

            cout << " " << endl;
            cout << "====================================" << endl;
            cout << "BFeed down subtraction              " << endl;
            cout << "Dmeson pt " << Dpt[iDpt] << "-" << Dpt[iDpt+1] << endl;
            cout << "had pt = " << hadpt[ihadpt] << " to " << hadptMax[ihadpt] << " - applied purity " << purities[ihadpt] << endl;
            cout << "====================================" << endl;
            cout << " " << endl;
            
            //For D+ input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
            if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 5;}
            if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 8;}
            if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 16;}
           
            // set correct paths
            inputcorrelation = Form("%s/%s%.0f_%.0f_%.1f_%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),Dptbin[0],Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have

            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            
            // set correct paths
            outputfilename = "./Final_Plots_pp/Singlev2Envelope/pp_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            //FDsubDstarPt3to5assoc0.3to1.0_v2D0.00_v2had0.00.root

            // first one - no v2
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode);
            if(collsyst) cout << "Check: This is not pp" << endl; // if pp, stop here
            
        }
    }
} // end function
