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
TString purityTempldir="";
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

void RunFeedown_pPb_Dplus(Int_t collsyst, Bool_t subtrMCclos, Bool_t oldnames, Int_t centbin, Int_t purityOpt){
  GetEnvelopeForEachV2(collsyst,subtrMCclos,oldnames,centbin,purityOpt);
  GetTotalEnvelopeFromV2(collsyst,subtrMCclos,oldnames,centbin,purityOpt);
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

    //Int_t collsyst = 2; // 0 is pp, 1 is p-Pb 2013, 2 is pPb 2016 (note that if you run on pp, it will perform only the v2=0 feeddown     
    if(!isLoadedFDsubtract){
      gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
      isLoadedFDsubtract=kTRUE;
    }
    SetSystemStringForTemplateFDnames(strSystemFDtempl.Data());

    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    TString inputcorrelation = ""; // input data file (not needed here)
  
    TString outputfilename = ""; //  (not needed here)

    Double_t Dpt[] = {3,5,8,16,24}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1.0,2.0,3.0,1.0,2.0}; // set associated tracks pt bins (lower) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMaxInput[] = {99.0,0.99,99.0,99.0,99.0,1.99,2.99}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMax[] = {99.0,1.0,99.0,99.0,99.0,2.0,3.0}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t purities[7];
    if(collsyst==0){
     purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973;  //03-99, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; 
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");
    }
    else if(collsyst==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; //03-99, 03-1, 1-99
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst==2 && centbin==0){ //0-100
        purities[0] = 0.952; purities[1] = 0.949; purities[2] = 0.962; purities[3] = 0.970; purities[4] = 0.977; purities[5] = 0.958; purities[6] = 0.967; //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
        //CORRECTION FOR RATIO OF EFFICIENCIES WITH ONLY pi,K,p,e,mu (CORRECT) OVER ALL TRACKS (WRONG AND USED FOR RESULTS)
 /* DEACTIVATED */ //        purities[0]/=1.015;  purities[1]/= 1.010; purities[2]/= 1.042; purities[3]/= 1.063; purities[4]/= 1.069; purities[5]/= 1.038; purities[6]/= 1.063; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst==2 && centbin==1){ //0-20
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.961; purities[3] = 0.970; purities[4] = 0.977; purities[5] = 0.958; purities[6] = 0.967; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else if(collsyst==2 && centbin==2){ //20-60
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.962; purities[3] = 0.971; purities[4] = 0.979; purities[5] = 0.959; purities[6] = 0.968; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else if(collsyst==2 && centbin==3){ //60-100
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.964; purities[3] = 0.973; purities[4] = 0.979; purities[5] = 0.961; purities[6] = 0.971; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else cout << " ! Purity Value is not system wise " << endl;

    
    Int_t systmode = 3; // mode for evaluation of the envelope - standard is set to 3

    Int_t dPtmax = 3;
    Int_t dhadPtmax = 3;
    Int_t dPtMin = 1;
    if(collsyst==2 && centbin==0) {dPtmax = 4; dhadPtmax = 7; dPtMin = 0;}
    if(collsyst==2 && centbin!=0) {dPtmax = 4; dhadPtmax = 3; dPtMin = 0;}    

    for(Int_t iDpt = dPtMin; iDpt <dPtmax; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <dhadPtmax; ihadpt++){


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
            if(collsyst==1) {
              if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 5;}
              if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 8;}
              if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 16;}            
            } else if(collsyst==2) {
              if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 4;}
              if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 7;}
              if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 10;}
              if(iDpt == 3) {Dptbin[0] = 11; Dptbin[1] = 11;}   
            } else {
              printf("Wrong collision system!\n");
              return;
            }
           
            // set correct paths
            inputcorrelation = Form("%s/%s%.0f_%.0f_%.1f_%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),Dptbin[0],Dptbin[1],hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            if(!oldnames) inputcorrelation = Form("%s/%s%dto%d_PoolInt_thr%.1fto%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have

            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            
            // set correct paths
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            // be careful... for this check, do not end the outputfilename with .root... the functions that are called will add a prefix based on the value of the v2 to be modulated
            
            TString binDIn = -1, binDOut = -1;
            TString binhIn = -1, binhOut = -1;
            if(iDpt == 0) {binDIn = "3"; binDOut = "5";}
            if(iDpt == 1) {binDIn = "5"; binDOut = "8";}
            if(iDpt == 2) {binDIn = "8"; binDOut = "16";}
            if(iDpt == 3) {binDIn = "16"; binDOut = "24";} 
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
            if(!purityOpt) SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode,oldnames,subtrMCclos,centbin,filepurity);
            if(!collsyst) continue; // if pp, stop here
                    
            //________________________________________________________________
            // v2 combination 1
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            if(!purityOpt) SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmax,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmax,systmode,oldnames,subtrMCclos,centbin,filepurity);
           
           
            //________________________________________________________________
            // v2 combination 2
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            if(!purityOpt) SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmax,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmax,systmode,oldnames,subtrMCclos,centbin,filepurity);
            
            
            //________________________________________________________________
            // v2 combination 3
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            if(!purityOpt) SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmin,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmin,systmode,oldnames,subtrMCclos,centbin,filepurity);
            
            //________________________________________________________________
            // v2 combination 4
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            if(!purityOpt) SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmin,systmode,oldnames,subtrMCclos,centbin);
            else SubtractFDexploitingClassDplus(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmin,systmode,oldnames,subtrMCclos,centbin,filepurity);
        }
        
    }
    
} // end function
//___________________________________________________________
void GetTotalEnvelopeFromV2(Int_t collsyst, Bool_t subtrMCclos, Bool_t oldnames, Int_t centbin, Int_t purityOpt){
    // note - this functions has to be runned only on p-Pb... setting it by mistake to pp should abort the process
    //Int_t collsyst = 2; // 0 is pp, 1 is p-Pb 2013, 2 is pPb 2016 (note that if you run on pp, it will perform only the v2=0 feeddown
    

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

    Double_t Dpt[] = {3,5,8,16,24}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1.0,2.0,3.0,1.0,2.0}; // set associated tracks pt bins (lower) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMaxInput[] = {99.0,0.99,99.0,99.0,99.0,1.99,2.99}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMax[] = {99.0,1.0,99.0,99.0,99.0,2.0,3.0}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t purities[7];

    if(collsyst==0){
     purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973;  //03-99, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; 
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");
    }
    else if(collsyst==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; //03-99, 03-1, 1-99
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst==2 && centbin==0){ //0-100
        purities[0] = 0.952; purities[1] = 0.949; purities[2] = 0.962; purities[3] = 0.970; purities[4] = 0.977; purities[5] = 0.958; purities[6] = 0.967; //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
        //CORRECTION FOR RATIO OF EFFICIENCIES WITH ONLY pi,K,p,e,mu (CORRECT) OVER ALL TRACKS (WRONG AND USED FOR RESULTS)
 /* DEACTIVATED */ //        purities[0]/=1.015;  purities[1]/= 1.010; purities[2]/= 1.042; purities[3]/= 1.063; purities[4]/= 1.069; purities[5]/= 1.038; purities[6]/= 1.063; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst==2 && centbin==1){ //0-20
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.961; purities[3] = 0.970; purities[4] = 0.977; purities[5] = 0.958; purities[6] = 0.967; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else if(collsyst==2 && centbin==2){ //20-60
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.962; purities[3] = 0.971; purities[4] = 0.979; purities[5] = 0.959; purities[6] = 0.968; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else if(collsyst==2 && centbin==3){ //60-100
        purities[0] = 0.953; purities[1] = 0.949; purities[2] = 0.964; purities[3] = 0.973; purities[4] = 0.979; purities[5] = 0.961; purities[6] = 0.971; //03-99, 03-1, 1-99, 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else cout << " ! Purity Value is not system wise " << endl;
    
    Int_t systmode = 3;

    Int_t dPtmax = 3;
    Int_t dhadPtmax = 3;
    Int_t dPtMin = 1;
    if(collsyst==2 && centbin==0) {dPtmax = 4; dhadPtmax = 7; dPtMin = 0;}
    if(collsyst==2 && centbin!=0) {dPtmax = 4; dhadPtmax = 3; dPtMin = 0;}    

    for(Int_t iDpt = dPtMin; iDpt <dPtmax; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <dhadPtmax; ihadpt++){
            
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
            if(collsyst==1) {
              if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 5;}
              if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 8;}
              if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 16;}           
            } else if(collsyst==2) {
              if(iDpt == 0) {Dptbin[0] = 3; Dptbin[1] = 4;}
              if(iDpt == 1) {Dptbin[0] = 5; Dptbin[1] = 7;}
              if(iDpt == 2) {Dptbin[0] = 8; Dptbin[1] = 10;}
              if(iDpt == 3) {Dptbin[0] = 11; Dptbin[1] = 11;}   
            } else {
              printf("Wrong collision system!\n");
              return;
            }         
            
            // set correct paths
            inputcorrelation = Form("%s/%s%.0f_%.0f_%.1f_%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),Dptbin[0],Dptbin[1],hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            if(!oldnames) inputcorrelation = Form("%s/%s%dto%d_PoolInt_thr%.1fto%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            // set correct paths
            outputfilename = "./Final_Plots_pPb/TotalEnvelope/pPb_FDsubDplus";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f.root",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            TString binDIn = -1, binDOut = -1;
            TString binhIn = -1, binhOut = -1;
            if(iDpt == 0) {binDIn = "3"; binDOut = "5";}
            if(iDpt == 1) {binDIn = "5"; binDOut = "8";}
            if(iDpt == 2) {binDIn = "8"; binDOut = "16";}
            if(iDpt == 3) {binDIn = "16"; binDOut = "24";} 
            if(ihadpt == 0) {binhIn = "03"; binhOut = "99";}
            if(ihadpt == 1) {binhIn = "03"; binhOut = "1";}
            if(ihadpt == 2) {binhIn = "1"; binhOut = "99";}
            if(ihadpt == 3) {binhIn = "2"; binhOut = "3";} //here I put 2-3 since this is what I have for templates (we don't use this bin) - it shall be 2-99
            if(ihadpt == 4) {binhIn = "3"; binhOut = "99";} 
            if(ihadpt == 5) {binhIn = "1"; binhOut = "2";}
            if(ihadpt == 6) {binhIn = "2"; binhOut = "3";}
            TString filepurity = Form("%s/DeltaPhi_%sto%s_%sto%s_RatioPrimOverAll.root",purityTempldir.Data(),binDIn.Data(),binDOut.Data(),binhIn.Data(),binhOut.Data());
            cout << "Purity option: " << purityOpt << " - FILE INPUT FOR PURITY TO BE LOADED: " << filepurity << endl;

            // running everything
            if(!purityOpt) SubtractFDexploitingClassDplusv2Modulations(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,systmode,oldnames,subtrMCclos,centbin); // collsyst is always set to 1!
            else SubtractFDexploitingClassDplusv2Modulations(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purityOpt,1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,systmode,oldnames,subtrMCclos,centbin,filepurity); // collsyst is always set to 1!
            
        }
        
    }
    
    
}


