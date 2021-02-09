//*******************************************************************************/
//  MACRO TO PERFORM THE FD SUBTRACTION FOR D0-hadron correlations in pp
//   Jitendra.Kumar (jitendra.kumar@cern.ch)  
//
//*******************************************************************************/
TString inputfc = "./Dzero/HFPtSpectrum_pp.root"; // input fprompt
TString templatedir = "./Templates_pp/"; // template path
TString inputcorrelationDir = "./Input_Plots_pp/";// directory where input files are stored
TString inputfileroot="1D_Signal_WithEMCorr_Normal_Charg_OriginSuper_Integrated";
TString strSystemFDtempl="none";
TString fdsubtrmacrodir="";
TString purityTempldir="";
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

void RunFeedown_pp_Dzero(Int_t collsyst, Bool_t subtrMCclos, Bool_t oldnames, Int_t centbin, Int_t purityOpt){
    GetEnvelopeForEachV2(collsyst,subtrMCclos,oldnames,centbin,purityOpt);
}
void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}
void SetPurityTemplateDir(TString purtempldite){
  purityTempldir=purtempldite;
}

//_____________________________________________________________
void GetEnvelopeForEachV2(Int_t collsyst, Bool_t subtrMCclos, Bool_t oldnames, Int_t centbin, Int_t purityOpt){
    
    //**********************************
    // This function loops on all the templates, creating 5 envelopes for different v2 values.
    // The output file contains also the modulated templates
    //***********************************

    gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
    SetSystemStringForTemplateFDnames(strSystemFDtempl.Data());

    TString inputcorrelation;
  
    TString outputfilename = ""; //  (not needed here)
    
    Double_t Dpt[] = {2,3,5,8,16,24}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1.0,2.0,3.0,1.0,2.0}; // set associated tracks pt bins (lower) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMaxInput[] = {99.0,0.99,99.0,99.0,99.0,1.99,2.99}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMax[] = {99.0,1.0,99.0,99.0,99.0,2.0,3.0}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t purities[7];
    if(collsyst ==0){
      purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973; //03-05, 05-1, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; //03-05, 05-1, 03-1, 1-99
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");

    }
    else if(collsyst ==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");

    }
    else if(collsyst ==3){ //pp@5TeV
        purities[0] = 0.956; purities[1] = 0.949; purities[2] = 0.967; purities[3] = 0.974; purities[4] = 0.980; purities[5] = 0.964; purities[6] = 0.973; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");

    }    
    else cout << " ! Purity Value is not system wise " << endl;

    Int_t dPtmax = 4;
    Int_t dhadPtmax = 3;
    Int_t dPtMin = 1;
    if(collsyst==3) {dPtmax = 5; dhadPtmax = 7; dPtMin = 0;printf("STO\n");}
    
    Int_t systmode = 3; // mode for evaluation of the envelope - standard is set to 3
    for(Int_t iDpt = dPtMin; iDpt <dPtmax; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <dhadPtmax; ihadpt++){

            cout << " " << endl;
            cout << "====================================" << endl;
            cout << "BFeed down subtraction              " << endl;
            cout << "Dmeson pt " << Dpt[iDpt] << "-" << Dpt[iDpt+1] << endl;
            cout << "had pt = " << hadpt[ihadpt] << " to " << hadptMaxInput[ihadpt] << " - applied purity " << purities[ihadpt] << endl;
            cout << "====================================" << endl;
            cout << " " << endl;
            
 
            //For D0 input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
            if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 6;}
            if(iDpt == 2) {Dptbin[0] = 7; Dptbin[1] = 9;}
            if(iDpt == 3) {Dptbin[0] = 10; Dptbin[1] = 11;}
            if(collsyst==3) {
              if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 3;}
              if(iDpt == 1) {Dptbin[0] = 4; Dptbin[1] = 5;}
              if(iDpt == 2) {Dptbin[0] = 6; Dptbin[1] = 8;}
              if(iDpt == 3) {Dptbin[0] = 9; Dptbin[1] = 11;}
              if(iDpt == 4) {Dptbin[0] = 12; Dptbin[1] = 13;}
            }

            // set correct paths
            inputcorrelation = Form("%s/%s%dto%d_Limits_2_4_TreshPt_%.1f_to_%.2f_Data.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            if(!oldnames) inputcorrelation = Form("%s/%s%dto%d_PoolInt_thr%.1fto%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            
            // set correct paths
            outputfilename = "./Final_Plots_pp/Singlev2Envelope/pp_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
                        
            TString binDIn = -1, binDOut = -1;
            TString binhIn = -1, binhOut = -1;
            if(iDpt == 0) {binDIn = "2"; binDOut = "3";}
            if(iDpt == 1) {binDIn = "3"; binDOut = "5";}
            if(iDpt == 2) {binDIn = "5"; binDOut = "8";}
            if(iDpt == 3) {binDIn = "8"; binDOut = "16";}
            if(iDpt == 4) {binDIn = "16"; binDOut = "24";} 
            if(ihadpt == 0) {binhIn = "03"; binhOut = "99";}
            if(ihadpt == 1) {binhIn = "03"; binhOut = "1";}
            if(ihadpt == 2) {binhIn = "1"; binhOut = "99";}
            if(ihadpt == 3) {binhIn = "2"; binhOut = "3";} //here I put 2-3 since this is what I have for templates (we don't use this bin) - it shall be 2-99
            if(ihadpt == 4) {binhIn = "3"; binhOut = "99";} 
            if(ihadpt == 5) {binhIn = "1"; binhOut = "2";}
            if(ihadpt == 6) {binhIn = "2"; binhOut = "3";}
            TString filepurity = Form("%s/DeltaPhi_%sto%s_%sto%s_RatioPrimOverAll.root",purityTempldir.Data(),binDIn.Data(),binDOut.Data(),binhIn.Data(),binhOut.Data());
            cout << "Purity option: " << purityOpt << " - FILE INPUT FOR PURITY TO BE LOADED: " << filepurity << endl;

            // first one - no v2
            if(!purityOpt) SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode,oldnames,subtrMCclos,centbin,filepurity);
            
           
        }
    }
    
} // end function
//=========================================================

