//*******************************************************************************/
//  MACRO TO PERFORM THE FD SUBTRACTION FOR D+-hadron correlations in pPb
//   Jitendra.Kumar (jitendra.kumar@cern.ch)  
//
//*******************************************************************************/
TString inputfc = "./Dplus/HFPtSpectrum_pp.root"; // input fprompt
TString templatedir = "./Templates_pp/"; // template path
TString inputcorrelationDir = "./Input_Plots_pp/";// directory where input files are stored
TString inputfileroot="1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated";
TString fdsubtrmacrodir="";
TString strSystemFDtempl="none";
Bool_t isLoadedFDsubtract=kFALSE;
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


void RunFeedown_pPb_Dplus(){
    GetEnvelopeForEachV2();
    GetTotalEnvelopeFromV2();
}

void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}

//_____________________________________________________________
void GetEnvelopeForEachV2(){
    
    //**********************************
    // This function loops on all the templates, creating 5 envelopes for different v2 values.
    // The output file contains also the modulated templates
    //***********************************

    Int_t collsyst = 1; // 0 is pp, 1 is p-Pb (note that if you run on pp, it will perform only the v2=0 feeddown
    if(!isLoadedFDsubtract){
      gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
      isLoadedFDsubtract=kTRUE;
    }
    SetSystemStringForTemplateFDnames(strSystemFDtempl.Data());

    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    TString inputcorrelation = ""; // input data file (not needed here)
  
    TString outputfilename = ""; //  (not needed here)

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
    for(Int_t iDpt = 1; iDpt <3; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <3; ihadpt++){

            cout << " " << endl;
            cout << "====================================" << endl;
            cout << "BFeed down subtraction              " << endl;
            cout << "Dmeson pt " << Dpt[iDpt] << "-" << Dpt[iDpt+1] << endl;
            cout << "had pt = " << hadpt[ihadpt] << " to " << hadptMax[ihadpt] << " - applied purity " << purities[ihadpt] << endl;
            cout << "====================================" << endl;
            cout << " " << endl;
            
            if(hadpt[ihadpt] < 0.6){v2hadmin = 0.05; v2hadmax = 0.07;}
            if(hadpt[ihadpt] > 0.6){v2hadmin = 0.08; v2hadmax = 0.13;}
 
            //For D+ input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
            if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 5;}
            if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 8;}
            if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 16;}
           
            // set correct paths
            inputcorrelation = Form("%s/%s%.0f_%.0f_%.1f_%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),Dptbin[0],Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have

            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            
            // set correct paths
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            // be careful... for this check, do not end the outputfilename with .root... the functions that are called will add a prefix based on the value of the v2 to be modulated
            
      
            // first one - no v2
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode);
            if(!collsyst) continue; // if pp, stop here
            
        
            //________________________________________________________________
            // v2 combination 1
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmax,systmode);
           
           
            //________________________________________________________________
            // v2 combination 2
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmax,systmode);
            
            
            //________________________________________________________________
            // v2 combination 3
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmin,systmode);
            
            //________________________________________________________________
            // v2 combination 4
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmin,systmode);
        }
        
    }
    
} // end function
//___________________________________________________________
void GetTotalEnvelopeFromV2(){
    // note - this functions has to be runned only on p-Pb... setting it by mistake to pp should abort the process
    Int_t collsyst = 1; // 0 is pp, 1 is p-Pb (note that if you run on pp, it will perform only the v2=0 feeddown
    

    gSystem->Exec("mkdir -p ./Final_Plots_pPb/TotalEnvelope/");
    gSystem->Sleep(100);
    if(!isLoadedFDsubtract){
      gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
      isLoadedFDsubtract=kTRUE;
    }
    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    TString inputcorrelation = ""; // input data file
    
    TString outputfilename = "";

    Double_t Dpt[] = {3,5,8,16}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1}; // set associated tracks pt bins (lower) //03-05, 05-1, 03-1, 1-99
    Double_t hadptMaxInput[] = {99.0,1.0,99.0}; // set associated tracks pt bins (upper) //03-05, 05-1, 03-1, 1-99
    Double_t hadptMax[] = {99.0,1.0,99.0}; // set associated tracks pt bins (upper) //03-05, 05-1, 03-1, 1-99
     
    Double_t purities[3];// = {0.964,0.967,0.968};
    purities[0] = 0.965;
    purities[1] = 0.965;
    purities[2] = 0.965;
    
    Int_t systmode = 3;
    for(Int_t iDpt = 1; iDpt <3; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <3; ihadpt++){
            
            cout << " " << endl;
            cout << "====================================" << endl;
            cout << "BFeed down subtraction              " << endl;
            cout << "Dmeson pt " << Dpt[iDpt] << "-" << Dpt[iDpt+1] << endl;
            cout << "had pt = " << hadpt[ihadpt] << " to " << hadptMax[ihadpt] << " - applied purity " << purities[ihadpt] << endl;
            cout << "====================================" << endl;
            cout << " " << endl;
            
            if(hadpt[ihadpt] < 0.6){v2hadmin = 0.05; v2hadmax = 0.07;}
            if(hadpt[ihadpt] > 0.6){v2hadmin = 0.08; v2hadmax = 0.13;}
            
            //For D0 input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
            if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 5;}
            if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 8;}
            if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 16;}
            
            // set correct paths
            inputcorrelation = Form("%s/%s%.0f_%.0f_%.1f_%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),Dptbin[0],Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            // set correct paths
            outputfilename = "./Final_Plots_pPb/TotalEnvelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f.root",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            // be careful... for this check, END the outputfilename with .root... the functions that are called will NOT change the name of your output anymore
            // GetTemplateFromFit(hFDtemplFile[0],hFDtempl[0],"cFitFD",0,v2D,v2Had);

            // running everything
            SubtractFDexploitingClassDplusv2Modulations(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),1,systmode); // collsyst is always set to 1!
            
        }
        
    }
    
    
}


// null
