//*******************************************************************************/
//*******************************************************************************/
//  MACRO TO PERFORM THE FD SUBTRACTION FOR D0-hadron correlations in pPb
//   Jitendra.Kumar (jitendra.kumar@cern.ch)  
//
//*******************************************************************************/

TString inputfc = "./Dzero/HFPtSpectrum_pp.root"; // input fprompt
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


void RunFeedown_pPb_Dzero(Int_t collsyst, Bool_t subtrMCclos){
  GetEnvelopeForEachV2(collsyst,subtrMCclos);
  GetTotalEnvelopeFromV2(collsyst,subtrMCclos);
}

void SetFDtemplateSystemString(TString str){
  strSystemFDtempl=str;
}

//_____________________________________________________________
void GetEnvelopeForEachV2(Int_t collsyst, Bool_t subtrMCclos){
    
    //**********************************
    // This function loops on all the templates, creating 5 envelopes for different v2 values.
    // The output file contains also the modulated templates
    //***********************************

    //Int_t collsyst = 2; // 0 is pp, 1 is p-Pb 2013, 2 is pPb 2016(note that if you run on pp, it will perform only the v2=0 feeddown
    if(!isLoadedFDsubtract){
      gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
      isLoadedFDsubtract=kTRUE;
    }


    SetSystemStringForTemplateFDnames(strSystemFDtempl.Data());
    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    TString inputcorrelation = ""; // input data file (not needed here)
    //    TString inputfc = "./Dzero/HFPtSpectrum_DrawFpromptVsRaaElossHypoCombined.root"; // input fprompt
    //    TString templatedir = "./Templates_pPb/"; // template path
  
    TString outputfilename = ""; //  (not needed here)

    Int_t oldnames=1; if(collsyst!=0 && collsyst!=1) oldnames=0;

    Double_t Dpt[] = {3,5,8,16,24}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1.0,2.0,3.0,1.0,2.0}; // set associated tracks pt bins (lower) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMaxInput[] = {99.0,0.99,99.0,99.0,99.0,1.99,2.99}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMax[] = {99.0,1.0,99.0,99.0,99.0,2.0,3.0}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t purities[7];
    if(collsyst ==0){
     purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973;  //03-99, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; 
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");
    }
    else if(collsyst ==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; //03-99, 03-1, 1-99
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst ==2){
        purities[0] = 0.958; purities[1] = 0.953; purities[2] = 0.973; purities[3] = 0.985; purities[4] = 0.990; purities[5] = 0.969; purities[6] = 0.982; //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
        //CORRECTION FOR RATIO OF EFFICIENCIES WITH ONLY pi,K,p,e,mu (CORRECT) OVER ALL TRACKS (WRONG AND USED FOR RESULTS)
  /*DEACTIVATED*/ //	purities[0]/=1.015;  purities[1]/= 1.010; purities[2]/= 1.042; purities[3]/= 1.063; purities[4]/= 1.069; purities[5]/= 1.038; purities[6]/= 1.063; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else cout << " ! Purity Value is not system wise " << endl;

    
    Int_t systmode = 3; // mode for evaluation of the envelope - standard is set to 3

    Int_t dPtmax = 3;
    Int_t dhadPtmax = 3;
    Int_t dPtMin = 1;
    if(collsyst==2) {dPtmax = 4; dhadPtmax = 7; dPtMin = 0;}

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
 
            //For D0 input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
	          if(collsyst==1) {
              if(iDpt == 0) {Dptbin[0] = 5; Dptbin[1] = 6;}
              if(iDpt == 1) {Dptbin[0] = 7; Dptbin[1] = 9;}
              if(iDpt == 2) {Dptbin[0] = 10; Dptbin[1] = 11;}                
            } else if(collsyst==2) {
              if(iDpt == 0) {Dptbin[0] = 4; Dptbin[1] = 5;}
              if(iDpt == 1) {Dptbin[0] = 6; Dptbin[1] = 8;}
              if(iDpt == 2) {Dptbin[0] = 9; Dptbin[1] = 11;}
              if(iDpt == 3) {Dptbin[0] = 12; Dptbin[1] = 12;}
            } else {
              printf("Wrong collision system!\n");
              return;
            }

            // set correct paths
            inputcorrelation = Form("%s/%s%dto%d_Limits_2_4_TreshPt_%.1f_to_%.2f_Data.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            if(collsyst>=2) inputcorrelation = Form("%s/%s%dto%d_PoolInt_thr%.1fto%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            
            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;
            
            
            // set correct paths
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);

	    //            outputfilename = "./Dzero/Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
	    //            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            // be careful... for this check, do not end the outputfilename with .root... the functions that are called will add a prefix based on the value of the v2 to be modulated
            
      
            // first one - no v2
            SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,0,0,systmode,oldnames,subtrMCclos);
            if(!collsyst) continue; // if pp, stop here
            
            //________________________________________________________________
            // v2 combination 1
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmax,systmode,oldnames,subtrMCclos);
           
           
            //________________________________________________________________
            // v2 combination 2
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmax,systmode,oldnames,subtrMCclos);
            
            
            //________________________________________________________________
            // v2 combination 3
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmax,v2hadmin,systmode,oldnames,subtrMCclos);
            
            //________________________________________________________________
            // v2 combination 4
            outputfilename = "./Final_Plots_pPb/Singlev2Envelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            SubtractFDexploitingClassDzero(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,v2Dmin,v2hadmin,systmode,oldnames,subtrMCclos);
        }
        
    }
    
} // end function
//___________________________________________________________
void GetTotalEnvelopeFromV2(Int_t collsyst, Bool_t subtrMCclos){
    // note - this functions has to be runned only on p-Pb... setting it by mistake to pp should abort the process
    //Int_t collsyst = 2; // 0 is pp, 1 is p-Pb 2013, 2 is pPb 2016 (note that if you run on pp, it will perform only the v2=0 feeddown
 
    gSystem->Exec("mkdir -p ./Final_Plots_pPb/TotalEnvelope/");
    gSystem->Sleep(100);
    if(!isLoadedFDsubtract){
      gROOT->LoadMacro(Form("%s/SubtractFD.C",fdsubtrmacrodir.Data()));
      isLoadedFDsubtract=kTRUE;
    }
    //exit(1);
    Double_t v2hadmin, v2hadmax, v2Dmin, v2Dmax;
    v2Dmin = 0.05; v2Dmax = 0.13;
    
    TString inputcorrelation = ""; // input data file
    //    TString inputfc = "HFPtSpectrum_DrawFpromptVsRaaElossHypoCombined.root"; // input fprompt
    //TString templatedir = "./Templates_pPb"; // template path
    
    TString outputfilename = "";
   
    Int_t oldnames=1; if(collsyst!=0 && collsyst!=1) oldnames=0;   
   
    Double_t Dpt[] = {3,5,8,16,24}; // set D meson pt bins
    Double_t hadpt[] = {0.3,0.3,1.0,2.0,3.0,1.0,2.0}; // set associated tracks pt bins (lower) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMaxInput[] = {99.0,0.99,99.0,99.0,99.0,1.99,2.99}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t hadptMax[] = {99.0,1.0,99.0,99.0,99.0,2.0,3.0}; // set associated tracks pt bins (upper) //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
    Double_t purities[7];

    if(collsyst ==0){
     purities[0] = 0.967; purities[1] = 0.963; purities[2] = 0.973;  //03-99, 03-1, 1-99
      //      purities[0] = 0.963; purities[1] = 0.967; purities[2] = 0.963; purities[3] = 0.973; 
      gSystem->Exec("mkdir -p Final_Plots_pp/Singlev2Envelope/");
    }
    else if(collsyst ==1){
        purities[0] = 0.965; purities[1] = 0.965; purities[2] = 0.965; //03-99, 03-1, 1-99
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }
    else if(collsyst ==2){
        purities[0] = 0.958; purities[1] = 0.953; purities[2] = 0.973; purities[3] = 0.985; purities[4] = 0.990; purities[5] = 0.969; purities[6] = 0.982; //03-99, 03-1, 1-99, 2-99, 3-99, 1-2, 2-3
        //CORRECTION FOR RATIO OF EFFICIENCIES WITH ONLY pi,K,p,e,mu (CORRECT) OVER ALL TRACKS (WRONG AND USED FOR RESULTS)
  /* DEACTIVATED*/ //	purities[0]/=1.015;  purities[1]/= 1.010; purities[2]/= 1.042; purities[3]/= 1.063; purities[4]/= 1.069; purities[5]/= 1.038; purities[6]/= 1.063; 
        gSystem->Exec("mkdir -p Final_Plots_pPb/Singlev2Envelope/");
    }    
    else cout << " ! Purity Value is not system wise " << endl;
    
    Int_t systmode = 3;

    Int_t dPtmax = 3;
    Int_t dhadPtmax = 3;
    Int_t dPtMin = 1;
    if(collsyst==2) {dPtmax = 4; dhadPtmax = 7; dPtMin = 0;}

    for(Int_t iDpt = dPtMin; iDpt <dPtmax; iDpt++){
        for(Int_t ihadpt = 0; ihadpt <dhadPtmax; ihadpt++){
            
            /*
            cout << " " << endl;
            cout << "====================================" << endl;
            cout << "BFeed down subtraction              " << endl;
            cout << "Dmeson pt " << Dpt[iDpt] << "-" << Dpt[iDpt+1] << endl;
            cout << "had pt = " << hadpt[ihadpt] << " to " << hadptMax[ihadpt] << " - applied purity " << purities[ihadpt] << endl;
            cout << "====================================" << endl;
            cout << " " << endl;
            */
            
            if(hadpt[ihadpt] < 0.6){v2hadmin = 0.05; v2hadmax = 0.07;}
            if(hadpt[ihadpt] > 0.6){v2hadmin = 0.08; v2hadmax = 0.13;}
            
            //For D0 input filenames! (bins and not pt)
            Double_t Dptbin[2] = {0,0};
            if(collsyst==1) {
              if(iDpt == 0) {Dptbin[0] = 5; Dptbin[1] = 6;}
              if(iDpt == 1) {Dptbin[0] = 7; Dptbin[1] = 9;}
              if(iDpt == 2) {Dptbin[0] = 10; Dptbin[1] = 11;}                
            } else if(collsyst==2) {
              if(iDpt == 0) {Dptbin[0] = 4; Dptbin[1] = 5;}
              if(iDpt == 1) {Dptbin[0] = 6; Dptbin[1] = 8;}
              if(iDpt == 2) {Dptbin[0] = 9; Dptbin[1] = 11;}
              if(iDpt == 3) {Dptbin[0] = 12; Dptbin[1] = 12;}
            } else {
              printf("Wrong collision system!\n");
              return;
            }          
            
             //1D_pPb_DplusCorr_3_5_0.3_1.0
            // set correct paths
            inputcorrelation = Form("%s/%s%dto%d_Limits_2_4_TreshPt_%.1f_to_%.2f_Data.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            if(collsyst>=2) inputcorrelation = Form("%s/%s%dto%d_PoolInt_thr%.1fto%.1f.root",inputcorrelationDir.Data(),inputfileroot.Data(),(int)Dptbin[0], (int)Dptbin[1], hadpt[ihadpt],hadptMaxInput[ihadpt]); // I guess all your input data files have
            cout << " inputcorrelation = " << inputcorrelation.Data() << endl;

            // set correct paths
            outputfilename = "./Final_Plots_pPb/TotalEnvelope/pPb_FDsubDzero";
            outputfilename += Form("Pt%dto%dassoc%1.1fto%1.1f.root",(int)Dpt[iDpt], (int)Dpt[iDpt+1], hadpt[ihadpt],hadptMax[ihadpt]);
            
            // be careful... for this check, END the outputfilename with .root... the functions that are called will NOT change the name of your output anymore
            //        GetTemplateFromFit(hFDtemplFile[0],hFDtempl[0],"cFitFD",0,v2D,v2Had);
printf("TOTAL ENVELOPE NOW!\n");
            // running everything
            SubtractFDexploitingClassDzerov2Modulations(Dpt[iDpt],Dpt[iDpt+1],hadpt[ihadpt],hadptMax[ihadpt],outputfilename.Data(),2,purities[ihadpt],1,inputcorrelation.Data(),inputfc.Data(),templatedir.Data(),collsyst,systmode,oldnames,subtrMCclos); // collsyst is always set to 1!
            
        }
        
    }
    
    
}


