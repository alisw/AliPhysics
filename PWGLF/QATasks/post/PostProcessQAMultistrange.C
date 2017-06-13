//////////////////////////////////////////////////
//
//  This macro was written by Domenico Colella (domenico.colella@cern.ch).
//  12 November 2013
//   - Updated for Run2 (24 October 2016)
//   - Centrality check removed (21 March 2017) 
//  11 May 2016
//   - Strong change in the output of the Task required a strange change in the 
//     structure of this macro. The object read from the task output are now TH2D.
//       
//  
//   ------------------------
//   ------ Arguments -------
//   ------------------------
//   --  collidingsystem   =  0) PbPb  1) pp/pPb
//   --  isMC              =  kTRUE if running on MC production 
//   --  mcass             =  kTRUE if want association for MC reconstructed particles
//   --  fileDir           =  "Input file directory"
//
//
//   -------------------------------------
//   ------ QATask output content --------
//   -------------------------------------
//   The output file produced by the AliAnalysisTaskQAMultistrange is a list of mono- end bi- dimensional histograms
//   that contain the main distribution of the variables used to select the cascade candidates. Here the full list of
//   histograms (given for XiMinus, but available also for XiPlus, OmegaMinus and OmegaPlus): 
//    -- PLOTS FOR RECONSTRUCTED CANDIDATES (ALSO ASSOCIATED IN CASE OF MC) 
//     -  fHistVarDcaCascDaughtXiMinus                [Max DCA Cascade Daughters]
//     -  fHistVarDcaBachToPrimVertexXiMinus          [Min DCA Bach To PV]
//     -  fHistVarCascCosineOfPointingAngleXiMinus    [Min Cascade Cosine Of PA]
//     -  fHistVarCascRadiusXiMinus                   [Min Cascade Radius Fid. Vol.]
//     -  fHistVarInvMassLambdaAsCascDghterXiMinus    [Window Invariant Mass Lambda]
//     -  fHistVarDcaV0DaughtersXiMinus               [Max DCA V0 Daughters]
//     -  fHistVarV0CosineOfPAToCascVertexXiMinus     [Min V0 Cosine Of PA To Cascade vertex]
//     -  fHistVarV0RadiusXiMinus                     [Min V0 Radius Fid. Vol.]
//     -  fHistVarDcaV0ToPrimVertexXiMinus            [Min DCA V0 To PV]
//     -  fHistVarDcaPosToPrimVertexXiMinus           [Min DCA Pos To PV]
//     -  fHistVarDcaNegToPrimVertexXiMinus           [Min DCA Neg To PV]
//     -  fHistMassXiMinus                            [Invariant Mass distribution]
//     -  fHistVarTransvMomentumXiMinus               [Transverse Momentum distribution]
//     -  fHistVarRapidityXiMinus                     [Rapidity distribution]
//     -  fHistVarCascProperLengthXiMinus             [Proper length distribution for the cascade]
//     -  fHistVarV0ProperLengthXiMinus               [Proper length distribution for the V0]
//    -- PLOTS FOR MC GENERATED PARTICLES (FILLED ONLY IN CASE OF MC PRODUCTION)
//     -  fHistCascadeMultiplicityXiMinus             [Distribution of the number of Xi minus per event]
//     -  fHistGenVarTotMomXiMinus                    [Total momentum]          
//     -  fHistGenVarTransvMomXiMinus                 [Transverse momentum]        
//     -  fHistGenVarYXiMinus                         [Rapidity]       
//     -  fHistGenVarEtaXiMinus                       [Pseudo-rapidity]      
//     -  fHistGenVarThetaXiMinus                     [Theta angle]      
//     -  fHistGenVarPhiXiMinus                       [Phi angle]      
//   Most of these histograms are TH2D with the present variable vs trasnverse momentum.
//
//
//   -----------------------------------
//   ------ Present Macro Checks -------
//   -----------------------------------
//   Using this macro many checks on the cascade topological reconstruction procedure
//   can be performed. In particular, the shape and the limit for the topological 
//   variable distributions as well as other kinematical variable distributions. The
//   reconstruction of the cascades are performed using two classes AliCascadeVertexer.cxx 
//   and AliV0vertexer.cxx contained in /STEER/ESD/ folder in Aliroot.
//   In the following are listed the contents of each page of the produced pdf:
//   
//   -- [Page 1] Distributions for the variables: 
//                DCA cascade daughters,  Bachelor IP to PV, 
//                Cascade cosine of PA,   Cascade radius of fiducial volume, 
//                Invariant mass Lambda,  DCA V0 daughters.
//   -- [Page 2] Distributions for the variables:
//                V0 cosine of PA to Casc Vtx,  Min V0 Radius fiducial volume, 
//                Min DCA V0 To PV,             Min DCA positive To PV, 
//                Min DCA negative To PV
//   -- [Page 3] Distributions for the variables;
//                Transverse momentum,    Rapidity
//                Cascade proper length,  V0 proper length.
//   -- [Page 4] Check on the invariant mass distribution fit
//   -- [Page 5] Generated cascade multiplicity distribution
//   -- [Page 6] Only in case of MC production, distributions for the MC generated
//               particles, of the variables:
//                Total momentum,         Transverse momentum,
//                Rapidity,               Pseudo-rapidity,
//                Theta angle,            Phi angle,
//
//////////////////////////////////////////////////////




class AliCFContainer;



//=================
// - MAIN FUNCTION
//=================
void PostProcessQAMultistrange(Int_t   collidingsystem = 1,                             // 0) PbPb  1) pp/pPb
                               Bool_t  isMC            = kTRUE,                        // kTRUE-->MC and kFALSE-->Exp.
                               Bool_t  mcass           = kTRUE,                         // if kTRUE use container with MC association
                               Char_t *fileDir         = ".",                           // Input file directory
                               Char_t *output          = "pdf"                          // "eps", "png" or "pdf"
                              ) {




     //______________________
     // - DEFINE DRAW OPTIONS
      myOptions();
      gROOT->ForceStyle();

     //__________________________
     // - SOURCE USEFUL LIBRARIES
     gSystem->Load("libANALYSIS.so");
     gSystem->Load("libANALYSISalice.so");
     gSystem->Load("libCORRFW.so");

     //______________
     // - DEFINE TEXT
     TLatex* t1 = new TLatex(0.6,0.7,"#color[3]{OK!!}");         myLatexMakeUp(t1,42,0.2,1);
     TLatex* t2 = new TLatex(0.6,0.7,"#color[2]{NOT OK!!}");     myLatexMakeUp(t2,42,0.2,2);
     TLatex* t31 = new TLatex(0.15,0.7,"#color[2]{TIGHTER CUT (WRT EXPECTED) IS USED!!}");               myLatexMakeUp(t31,42,0.1,2);
     TLatex* t32 = new TLatex(0.15,0.6,"#color[2]{PROBLEM FOR ANALYSIS. CHECK}"); myLatexMakeUp(t32,42,0.1,2);
     TLatex* t41 = new TLatex(0.15,0.7,"#color[42]{LOOSER CUT (WRT EXPECTED) IS USED!!}");                      myLatexMakeUp(t41,42,0.1,2);
     TLatex* t42 = new TLatex(0.15,0.6,"#color[42]{NOT AN ISSUE FOR ANALYSIS. CHECK}");  myLatexMakeUp(t42,42,0.1,2);
     Char_t *pname[4] = {"#color[1]{#Xi^{-}}", "#color[1]{#bar{#Xi}^{+}}", "#color[1]{#Omega^{-}}", "#color[1]{#bar{#Omega}^{+}}"};
     TLatex** tcas = new TLatex*[4];
     for (Int_t icas = 0; icas < 4; icas++) { tcas[icas] = new TLatex(0.8,0.3,Form("%s",pname[icas]));  myLatexMakeUp(tcas[icas],42,0.25,2); }
     const Int_t color[4] = {kRed+1,kOrange+1,kAzure+2,kViolet-4};
     Char_t *pdgmass[4] = {"PDG mass: 1.32171 GeV/c^{2}","PDG mass: 1.32171 GeV/c^{2}","PDG mass: 1.67245 GeV/c^{2}","PDG mass: 1.67245 GeV/c^{2}"};
     Char_t *cascade[4] = {"XiMinus","XiPlus","OmegaMinus","OmegaPlus"};

     //_____________________________________
     // - SOURCE THE FILE AND THE HISTOGRAMS
     TFile *f    = new TFile(Form("%s/AnalysisResults.root",fileDir));
     TList *list = (TList*) f->Get("PWGLFStrangeness.outputCheckCascade/fListHistMultistrangeQA");
     TH2F** h2dvar0  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar0[icas]  = list->FindObject(Form("fHistVarDcaCascDaught%s",cascade[icas])); 
       TH1F** hvar0  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar0[icas]  = h2dvar0[icas]->ProjectionY(Form("fHistVarDcaCascDaught[%i]",icas),0,-1);              
     TH2F** h2dvar1  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar1[icas]  = list->FindObject(Form("fHistVarDcaBachToPrimVertex%s",cascade[icas]));
       TH1F** hvar1  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar1[icas]  = h2dvar1[icas]->ProjectionY(Form("fHistVarDcaBachToPrimVertex[%i]",icas),0,-1);        
     TH2F** h2dvar2  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar2[icas]  = list->FindObject(Form("fHistVarCascCosineOfPointingAngle%s",cascade[icas]));
       TH1F** hvar2  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar2[icas]  = h2dvar2[icas]->ProjectionY(Form("fHistVarCascCosineOfPointingAngle[%i]",icas),0,-1);   
     TH2F** h2dvar3  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar3[icas]  = list->FindObject(Form("fHistVarCascRadius%s",cascade[icas]));       
       TH1F** hvar3  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar3[icas]  = h2dvar3[icas]->ProjectionY(Form("fHistVarCascRadius[%i]",icas),0,-1);                         
     TH2F** h2dvar4  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar4[icas]  = list->FindObject(Form("fHistVarInvMassLambdaAsCascDghter%s",cascade[icas])); 
       TH1F** hvar4  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar4[icas]  = h2dvar4[icas]->ProjectionY(Form("fHistVarInvMassLambdaAsCascDghter[%i]",icas),0,-1);  
     TH2F** h2dvar5  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar5[icas]  = list->FindObject(Form("fHistVarDcaV0Daughters%s",cascade[icas]));
       TH1F** hvar5  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar5[icas]  = h2dvar5[icas]->ProjectionY(Form("fHistVarDcaV0Daughters[%i]",icas),0,-1);             
     TH2F** h2dvar6  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar6[icas]  = list->FindObject(Form("fHistVarV0CosineOfPAToCascVertex%s",cascade[icas]));  
       TH1F** hvar6  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar6[icas]  = h2dvar6[icas]->ProjectionY(Form("fHistVarV0CosineOfPAToCascVertex[%i]",icas),0,-1);   
     TH2F** h2dvar7  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar7[icas]  = list->FindObject(Form("fHistVarV0Radius%s",cascade[icas]));                  
       TH1F** hvar7  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar7[icas]  = h2dvar7[icas]->ProjectionY(Form("fHistVarV0Radius[%i]",icas),0,-1);                   
     TH2F** h2dvar8  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar8[icas]  = list->FindObject(Form("fHistVarDcaV0ToPrimVertex%s",cascade[icas]));         
       TH1F** hvar8  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar8[icas]  = h2dvar8[icas]->ProjectionY(Form("fHistVarDcaV0ToPrimVertex[%i]",icas),0,-1);          
     TH2F** h2dvar9  = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar9[icas]  = list->FindObject(Form("fHistVarDcaPosToPrimVertex%s",cascade[icas]));        
       TH1F** hvar9  = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar9[icas]  = h2dvar9[icas]->ProjectionY(Form("fHistVarDcaPosToPrimVertex[%i]",icas),0,-1);         
     TH2F** h2dvar10 = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar10[icas] = list->FindObject(Form("fHistVarDcaNegToPrimVertex%s",cascade[icas]));        
       TH1F** hvar10 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar10[icas] = h2dvar10[icas]->ProjectionY(Form("fHistVarDcaNegToPrimVertex[%i]",icas),0,-1);        
       TH1F** hvar12 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar12[icas] = list->FindObject(Form("fHistVarTransvMomentum%s",cascade[icas]));                     
       TH1F** hvar13 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar13[icas] = list->FindObject(Form("fHistVarRapidity%s",cascade[icas]));                                         
     TH2F** h2dvar14 = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar14[icas] = list->FindObject(Form("fHistVarCascProperLength%s",cascade[icas]));          
       TH1F** hvar14 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar14[icas] = h2dvar14[icas]->ProjectionY(Form("fHistVarCascProperLength[%i]",icas),0,-1);          
     TH2F** h2dvar15 = new TH2F*[4];   for (Int_t icas = 0; icas < 4; icas++) h2dvar15[icas] = list->FindObject(Form("fHistVarV0ProperLength%s",cascade[icas]));            
       TH1F** hvar15 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar15[icas] = h2dvar15[icas]->ProjectionY(Form("fHistVarV0ProperLength[%i]",icas),0,-1);            
       TH1F** hvar16 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar16[icas] = list->FindObject(Form("fHistMass%s",cascade[icas]));                                  
       TH1F** hvar17 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar17[icas] = list->FindObject(Form("fHistCascadeMultiplicityMC%s",cascade[icas]));                 
       TH1F** hvar18 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar18[icas] = list->FindObject(Form("fHistGenVarTotMom%s",cascade[icas]));                          
       TH1F** hvar19 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar19[icas] = list->FindObject(Form("fHistGenVarTransvMom%s",cascade[icas]));                       
       TH1F** hvar20 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar20[icas] = list->FindObject(Form("fHistGenVarY%s",cascade[icas]));                               
       TH1F** hvar21 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar21[icas] = list->FindObject(Form("fHistGenVarEta%s",cascade[icas]));                             
       TH1F** hvar22 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar22[icas] = list->FindObject(Form("fHistGenVarTheta%s",cascade[icas]));                           
       TH1F** hvar23 = new TH1F*[4];   for (Int_t icas = 0; icas < 4; icas++) hvar23[icas] = list->FindObject(Form("fHistGenVarPhi%s",cascade[icas]));                             
    // - MAKEUP
    for (Int_t icas = 0; icas < 4; icas++) {

        myHistoMakeUp(hvar0[icas], color[icas], 0.06, 0.06, 0.06, 0.06);
        myHistoMakeUp(hvar1[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar2[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar3[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar4[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar5[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar6[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar7[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar8[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar9[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar10[icas], color[icas], 0.06, 0.06, 0.06, 0.06);          
        myHistoMakeUp(hvar12[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar13[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 

        myHistoMakeUp(hvar14[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar15[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar16[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar17[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar18[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar19[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar20[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar21[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar22[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 
        myHistoMakeUp(hvar23[icas], color[icas], 0.06, 0.06, 0.06, 0.06); 

    }

     //_____________________________________
     // - DEFINE TOPOLOGICAL VARIABLE VALUES                     [0]           [1]  
     //   --  0   = Max DCA Cascade Daughters                 PbPb: 0.3     pp/pPb: 2.0     
     //   --  1   = Min DCA Bach To PV                        PbPb: 0.03    pp/pPb: 0.01    
     //   --  2   = Min Cascade Cosine Of PA                  PbPb: 0.999   pp/pPb: 0.98   
     //   --  3   = Min Cascade Radius Fid. Vol.              PbPb: 0.9     pp/pPb: 0.2     
     //   --  4   = Window Invariant Mass Lambda              PbPb: 0.0008  pp/pPb: 0.008   
     //   --  5   = Max DCA V0 Daughters                      PbPb: 1.0     pp/pPb: 1.5     
     //   --  6   = Min V0 Cosine Of PA To PV                 PbPb: 0.98    pp/pPb: 0.9     
     //   --  7   = Min V0 Radius Fid. Vol.                   PbPb: 0.9     pp/pPb: 0.2     
     //   --  8   = Min DCA V0 To PV                          PbPb: 0.05    pp/pPb: 0.01    
     //   --  9   = Min DCA Pos To PV                         PbPb: 0.1     pp/pPb: 0.05    
     //   --  10  = Min DCA Neg To PV                         PbPb: 0.1     pp/pPb: 0.05    
     Float_t x0 = 0, x1 = 0, x2 = 0, x3 = 0, x41 = 0, x42 = 0, x5 = 0, x7 = 0, x8 = 0, x9 = 0, x10 = 0;
     if      (collidingsystem == 0) { x0 = 0.3; x1 = 0.03; x2 = 0.999; x3 = 0.9; x41 = 1.116 + 0.008; x42 = 1.115 - 0.008; x5 = 1.0; x7 = 0.9; x8 = 0.05; x9 = 0.10; x10 = 0.10;}
     else if (collidingsystem == 1) { x0 = 2; x1 = 0.01; x2 = 0.980; x3 = 0.2; x41 = 1.116 + 0.008; x42 = 1.115 - 0.008; x5 = 1.5; x7 = 0.2; x8 = 0.01; x9 = 0.05; x10 = 0.05;}
     TLine** line0  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line0[icas]  = new TLine(x0,0.,x0,hvar0[icas]->GetBinContent(hvar0[icas]->GetMaximumBin()));      myLineMakeUp(line0[icas],14,9,2.0);  }
     TLine** line1  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line1[icas]  = new TLine(x1,0.,x1,hvar1[icas]->GetBinContent(hvar1[icas]->GetMaximumBin()));      myLineMakeUp(line1[icas],14,9,2.0);  } 
     TLine** line2  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line2[icas]  = new TLine(x2,0.,x2,hvar2[icas]->GetBinContent(hvar2[icas]->GetMaximumBin()));      myLineMakeUp(line2[icas],14,9,2.0);  }
     TLine** line3  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line3[icas]  = new TLine(x3,0.,x3,hvar3[icas]->GetBinContent(hvar3[icas]->GetMaximumBin()));      myLineMakeUp(line3[icas],14,9,2.0);  }
     TLine** line41 = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line41[icas] = new TLine(x41,0.,x41,hvar4[icas]->GetBinContent(hvar4[icas]->GetMaximumBin()));    myLineMakeUp(line41[icas],14,9,2.0); }
     TLine** line42 = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line42[icas] = new TLine(x42,0.,x42,hvar4[icas]->GetBinContent(hvar4[icas]->GetMaximumBin()));    myLineMakeUp(line42[icas],14,9,2.0); }
     TLine** line5  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line5[icas]  = new TLine(x5,0.,x5,hvar5[icas]->GetBinContent(hvar5[icas]->GetMaximumBin()));      myLineMakeUp(line5[icas],14,9,2.0);  }
     TLine** line7  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line7[icas]  = new TLine(x7,0.,x7,hvar7[icas]->GetBinContent(hvar7[icas]->GetMaximumBin()));      myLineMakeUp(line7[icas],14,9,2.0);  }
     TLine** line8  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line8[icas]  = new TLine(x8,0.,x8,hvar8[icas]->GetBinContent(hvar8[icas]->GetMaximumBin()));      myLineMakeUp(line8[icas],14,9,2.0);  }
     TLine** line9  = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line9[icas]  = new TLine(x9,0.,x9,hvar9[icas]->GetBinContent(hvar9[icas]->GetMaximumBin()));      myLineMakeUp(line9[icas],14,9,2.0);  }
     TLine** line10 = new TLine*[4];  for (Int_t icas = 0; icas < 4; icas++) { line10[icas] = new TLine(x10,0.,x10,hvar10[icas]->GetBinContent(hvar10[icas]->GetMaximumBin()));  myLineMakeUp(line10[icas],14,9,2.0); }



     //_________________________________ 
     // - PAGE 1 : topological variables
     cout<<"\n--- BUILD THE FIRST PAGE: topological variable distributions for reco. candidate ---"<<endl;
     TCanvas** c1 = new TCanvas*[4];
     for (Int_t icas = 0; icas < 4; icas++)  {
           if      (icas == 0) cout<<"CASCADE: xi minus"<<endl;
           else if (icas == 1) cout<<"CASCADE: xi plus"<<endl;
           else if (icas == 2) cout<<"CASCADE: omega minus"<<endl;
           else if (icas == 3) cout<<"CASCADE: omega plus"<<endl;
           c1[icas] = new TCanvas(Form("c1_%i",icas),"",1200,800);
           c1[icas]->Divide(2,3); 
           // -- Pad 1: DCA cascade daughters
           c1[icas]->cd(1);
           myPadSetUp(gPad); 
            hvar0[icas]->Draw("histo");
            line0[icas]->Draw("same");
            if      (checkExactMaxLimit(hvar0[icas],x0) == 1) { cout<<"The cut on 'DCA cascade daughters' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMaxLimit(hvar0[icas],x0) == 2) { cout<<"The cut on 'DCA cascade daughters' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMaxLimit(hvar0[icas],x0) == 3) { cout<<"The cut on 'DCA cascade daughters' for cascade "<<icas<<" is FINE!!"<<endl;           t1->Draw();                }
            tcas[icas]->Draw();
           // -- Pad 2: Bachelor IP to PV
           c1[icas]->cd(2);
            myPadSetUp(gPad); 
            hvar1[icas]->GetXaxis()->SetRangeUser(0.,0.24);
            hvar1[icas]->Draw("histo");
            line1[icas]->Draw("same");
            if      (checkExactMinLimit(hvar1[icas],x1) == 1) { cout<<"The cut on 'Bachelor IP to PV' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw(); t42->Draw(); }
            else if (checkExactMinLimit(hvar1[icas],x1) == 2) { cout<<"The cut on 'Bachelor IP to PV' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw(); t32->Draw(); }
            else if (checkExactMinLimit(hvar1[icas],x1) == 3) { cout<<"The cut on 'Bachelor IP to PV' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();               }
           // -- Pad 3: Cascade cosine of Pointing Angle
           c1[icas]->cd(3);
            myPadSetUp(gPad); 
            hvar2[icas]->GetYaxis()->SetRangeUser(0.01,hvar2[icas]->GetMaximum()*1.5);
            hvar2[icas]->Draw("histo");
            line2[icas]->Draw("same");
            if      (checkExactMinLimit(hvar2[icas],x2) == 1) { cout<<"The cut on 'Cascade cosine of Pointing Angle' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw(); t42->Draw(); }
            else if (checkExactMinLimit(hvar2[icas],x2) == 2) { cout<<"The cut on 'Cascade cosine of Pointing Angle' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw(); t32->Draw(); }
            else if (checkExactMinLimit(hvar2[icas],x2) == 3) { cout<<"The cut on 'Cascade cosine of Pointing Angle' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();               }
           // -- Pad 4: Cascade radius of fiducial volume
           c1[icas]->cd(4);
            myPadSetUp(gPad); 
            hvar3[icas]->GetXaxis()->SetRangeUser(0.,3.8);
            hvar3[icas]->Draw("histo");
            line3[icas]->Draw("same");
            if      (checkExactMinLimit(hvar3[icas],x3) == 1) { cout<<"The cut on 'Cascade radius of fiducial volume' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMinLimit(hvar3[icas],x3) == 2) { cout<<"The cut on 'Cascade radius of fiducial volume' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMinLimit(hvar3[icas],x3) == 3) { cout<<"The cut on 'Cascade radius of fiducial volume' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad 5: Invariant mass Lambda
           c1[icas]->cd(5);
            myPadSetUp(gPad,kFALSE); 
            hvar4[icas]->Draw("histo");
            line41[icas]->Draw("same");
            line42[icas]->Draw("same");
            if      (checkExactMaxLimit(hvar4[icas],x41) == 1 && checkExactMinLimit(hvar4[icas],x42) == 1) { cout<<"The cut on 'Invariant mass Lambda' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMaxLimit(hvar4[icas],x41) == 2 && checkExactMinLimit(hvar4[icas],x42) == 2) { cout<<"The cut on 'Invariant mass Lambda' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMaxLimit(hvar4[icas],x41) == 3 && checkExactMinLimit(hvar4[icas],x42) == 3) { cout<<"The cut on 'Invariant mass Lambda' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad 6: DCA V0 daughters
           c1[icas]->cd(6);
            myPadSetUp(gPad); 
            hvar5[icas]->Draw("histo");
            line5[icas]->Draw("same");
            if      (checkExactMaxLimit(hvar5[icas],x5) == 1) { cout<<"The cut on 'DCA V0 daughters' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMaxLimit(hvar5[icas],x5) == 2) { cout<<"The cut on 'DCA V0 daughters' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMaxLimit(hvar5[icas],x5) == 3) { cout<<"The cut on 'DCA V0 daughters' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
          // -- OUTPUT
          if      (output == "png") c1[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page1_%i.png",icas));
          else if (output == "eps") c1[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page1_%i.eps",icas));
          else if (output == "pdf") c1[icas]->SaveAs("LF_QAanalysis_Multistrange.pdf(");
     }


     //_________________________________
     // - PAGE 2 : topological variables
     cout<<"\n--- BUILD THE SECOND PAGE: topological variable distributions for reco. candidate ---"<<endl;
     TCanvas** c2 = new TCanvas*[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           if      (icas == 0) cout<<"CASCADE: xi minus"<<endl;
           else if (icas == 1) cout<<"CASCADE: xi plus"<<endl;
           else if (icas == 2) cout<<"CASCADE: omega minus"<<endl;
           else if (icas == 3) cout<<"CASCADE: omega plus"<<endl;
           c2[icas] = new TCanvas(Form("c2_%i",icas),"",1200,800);
           c2[icas]->Divide(2,3);
           // -- Pad 1: V0 cosine of Pointing Angle to cascade vertex
           c2[icas]->cd(1);
            myPadSetUp(gPad);
            hvar6[icas]->GetYaxis()->SetRangeUser(0.01,(hvar6[icas]->GetBinContent(hvar6[icas]->GetMaximumBin()))*1.5);
            hvar6[icas]->Draw("histo");
            tcas[icas]->Draw();
           // -- Pad 2: Min V0 Radius Fid. Vol.  
           c2[icas]->cd(2);
            myPadSetUp(gPad); 
            hvar7[icas]->GetXaxis()->SetRangeUser(0.,3.0);
            hvar7[icas]->Draw("histo");
            line7[icas]->Draw("same");
            if      (checkExactMinLimit(hvar7[icas],x7) == 1) { cout<<"The cut on 'Min V0 Radius Fid. Vol.' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMinLimit(hvar7[icas],x7) == 2) { cout<<"The cut on 'Min V0 Radius Fid. Vol.' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMinLimit(hvar7[icas],x7) == 3) { cout<<"The cut on 'Min V0 Radius Fid. Vol.' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad3: Min DCA V0 To PV
           c2[icas]->cd(3);
            myPadSetUp(gPad); 
            hvar8[icas]->GetXaxis()->SetRangeUser(0.,0.3);
            hvar8[icas]->Draw("histo");
            line8[icas]->Draw("same");
            if      (checkExactMinLimit(hvar8[icas],x8) == 1) { cout<<"The cut on 'Min DCA V0 To PV' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMinLimit(hvar8[icas],x8) == 2) { cout<<"The cut on 'Min DCA V0 To PV' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMinLimit(hvar8[icas],x8) == 3) { cout<<"The cut on 'Min DCA V0 To PV' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad 4: Min DCA Pos To PV
           c2[icas]->cd(4);
            myPadSetUp(gPad); 
            hvar9[icas]->GetXaxis()->SetRangeUser(0.,0.2);
            hvar9[icas]->Draw("histo");
            line9[icas]->Draw("same");
            if      (checkExactMinLimit(hvar9[icas],x9) == 1) { cout<<"The cut on 'Min DCA Pos To PV' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMinLimit(hvar9[icas],x9) == 2) { cout<<"The cut on 'Min DCA Pos To PV' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMinLimit(hvar9[icas],x9) == 3) { cout<<"The cut on 'Min DCA Pos To PV' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad 5: Min DCA Neg To PV
           c2[icas]->cd(5);
            myPadSetUp(gPad); 
            hvar10[icas]->GetXaxis()->SetRangeUser(0.,0.2);
            hvar10[icas]->Draw("histo");
            line10[icas]->Draw("same");
            if      (checkExactMinLimit(hvar10[icas],x10) == 1) { cout<<"The cut on 'Min DCA Neg To PV' for cascade "<<icas<<" is +++LOOSER+++!!"<<endl;   t41->Draw();  t42->Draw(); }
            else if (checkExactMinLimit(hvar10[icas],x10) == 2) { cout<<"The cut on 'Min DCA Neg To PV' for cascade "<<icas<<" is ***TIGHTER***!!"<<endl;  t31->Draw();  t32->Draw(); }
            else if (checkExactMinLimit(hvar10[icas],x10) == 3) { cout<<"The cut on 'Min DCA Neg To PV' for cascade "<<icas<<" is FINE!!"<<endl;     t1->Draw();                }
           // -- Pad 6: Empty
           // -- OUTPUT
           if      (output == "png") c2[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page2_%i.png",icas));
           else if (output == "eps") c2[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page2_%i.eps",icas));
           else if (output == "pdf") c2[icas]->SaveAs("LF_QAanalysis_Multistrange.pdf");
     }


     //_____________________________
     // - PAGE 3 : further variables
     cout<<"\n--- BUILD THE THIRD PAGE: further variables for reco. candidate ---"<<endl;
     TCanvas** c3 = new TCanvas*[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           c3[icas] = new TCanvas(Form("c3_%i",icas),"",1200,800);
           c3[icas]->Divide(2,2);
           // -- Pad 1: Transverse momentum
           c3[icas]->cd(1);
            myPadSetUp(gPad,kFALSE);
            hvar12[icas]->Draw("histo");
            tcas[icas]->Draw();
           // -- Pad 2: Y
           c3[icas]->cd(2);
            myPadSetUp(gPad,kFALSE);
            hvar13[icas]->Draw("histo");
           // -- Pad 3: Cascade proper length
           c3[icas]->cd(3);
            myPadSetUp(gPad,kFALSE);
            hvar14[icas]->GetXaxis()->SetRangeUser(0.,90.);
            hvar14[icas]->Draw("histo");
           // -- Pad 4: V0 proper length 
           c3[icas]->cd(4);
            myPadSetUp(gPad,kFALSE);
            hvar15[icas]->GetXaxis()->SetRangeUser(0.,90.);
            hvar15[icas]->Draw("histo");
           // -- OUTPUT
           if      (output == "png") c3[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page3_%i.png",icas));
           else if (output == "eps") c3[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page3_%i.eps",icas));
           else if (output == "pdf") c3[icas]->SaveAs("LF_QAanalysis_Multistrange.pdf");
     }


     //____________________________________ 
     // - PAGE 4 : invariant mass fit check
     cout<<"\n--- BUILD THE FOURTH PAGE: invariant mass fit for reco. candidate ---"<<endl;
     TCanvas *c4 = new TCanvas("c4","",1200,800);
     c4->Divide(2,2);
     Float_t lowlimmass[4] = {1.300, 1.300, 1.645, 1.645};
     Float_t uplimmass[4]  = {1.340, 1.340, 1.700, 1.700};
     Double_t param1[4]    = {1.32171, 1.32171, 1.67245, 1.67245};
     Double_t param2min[4] = {1.310, 1.310, 1.664, 1.664};
     Double_t param2max[4] = {1.330, 1.330, 1.680, 1.680};
     Double_t param3min[4] = {0.001, 0.001, 0.0008, 0.0008};
     Double_t param3max[4] = {0.005, 0.005, 0.0060, 0.0060};     
     Float_t refwidth[4]   = {0.002, 0.002, 0.0025, 0.0025};
     TF1** fitinvmass = new TF1*[4];
     TPaveText** pave = new TPaveText[4];
     Double_t meanGauss[4], sigmaGauss[4];

     for (Int_t icas = 0; icas < 4; icas++)  {
           cout<<"\nAnalysing cascade: "<<icas<<endl;
           c4->cd(icas+1); 
            myPadSetUp(gPad,kFALSE);
            hvar16[icas]->Draw("histo");
            fitinvmass[icas] = new TF1(Form("fitinvmass[%i]",icas),"gaus(0)+pol2(3)",lowlimmass[icas],uplimmass[icas]);
            fitinvmass[icas]->SetParName(0, "cnstntG");
            fitinvmass[icas]->SetParName(1, "meanG");
            fitinvmass[icas]->SetParName(2, "sigmaG");
            fitinvmass[icas]->SetParLimits(0,0.,500000.);
            fitinvmass[icas]->SetParameter(1, param1[icas]);
            fitinvmass[icas]->SetParLimits(1, param2min[icas],param2max[icas]);
            fitinvmass[icas]->SetParLimits(2, param3min[icas],param3max[icas]);
            hvar16[icas]->Fit(Form("fitinvmass[%i]",icas),"rimeN");
            fitinvmass[icas]->SetLineColor(kRed);
            fitinvmass[icas]->Draw("same");
            meanGauss[icas]  = fitinvmass[icas]->GetParameter(1);
            sigmaGauss[icas] = fitinvmass[icas]->GetParameter(2);
            pave[icas] = new TPaveText(0.55,0.75,0.95,0.95,"NDC");
            myPaveMakeUp(pave[icas],Form("%s",pdgmass[icas]),0,22,0.04,0,0);
            pave[icas]->AddText(Form("#color[1]{Mass form Fit: %.5f #pm %.5f GeV/c^{2}}",meanGauss[icas],sigmaGauss[icas]));
            if (sigmaGauss[icas] > refwidth[icas] - 0.001 && sigmaGauss[icas] < refwidth[icas] + 0.001) pave[icas]->AddText("#color[3]{OK!! The width is compatible with standard.}");
            else                                                                                        pave[icas]->AddText("#color[2]{NOT OK!! Problem.}"); 
            pave[icas]->Draw("same");
            tcas[icas]->Draw();
     }
     // -- OUTPUT
     if      (output == "png") c4->SaveAs("LF_QAanalysis_Multistrange_page4.png");
     else if (output == "eps") c4->SaveAs("LF_QAanalysis_Multistrange_page4.eps");
     else if (output == "pdf" && isMC)  c4->SaveAs("LF_QAanalysis_Multistrange.pdf");
     else if (output == "pdf" && !isMC) c4->SaveAs("LF_QAanalysis_Multistrange.pdf");


     //_____________________________________________
     // - PAGE 5 : cascade multiplicity distribution
     cout<<"\n--- BUILD THE FIFTH PAGE: generated cascade multiplicity distribution ---"<<endl;
     if (isMC) {
        TList *l = (TList*) f->Get("PWGLFStrangeness.outputCheckCascade/fListHistMultistrangeQA");
        TCanvas *c5 = new TCanvas("c5","",1200,800);
        c5->Divide(2,2);
        Double_t mean[4] = {};
        Double_t integral[4] = {};
        Int_t    entries[4] = {};
        TPaveText** pave1 = new TPaveText[4];

        for (Int_t icas = 0; icas < 4; icas++)  {
              c5->cd(icas+1);
               myPadSetUp(gPad,kFALSE); 
               hvar17[icas]->Draw();
               for (Int_t j = 0; j < hvar17[icas]->GetNbinsX(); j++) {
                     mean[icas]     = mean[icas] + hvar17[icas]->GetBinCenter(j) * hvar17[icas]->GetBinWidth(j) * hvar17[icas]->GetBinContent(j);
                     integral[icas] = integral[icas] + hvar17[icas]->GetBinContent(j) * hvar17[icas]->GetBinWidth(j);
                     entries[icas] = entries[icas] + hvar17[icas]->GetBinContent(j);
               }
               mean[icas] = mean[icas] / integral[icas];
               tcas[icas]->Draw();  
               pave1[icas] = new TPaveText(0.55,0.75,0.95,0.95,"NDC");
               myPaveMakeUp(pave1[icas],Form("Mean : %.2f",mean[icas]),0,22,0.1,0,0);
               pave1[icas]->Draw("same");
        }
        // -- OUTPUT
        if      (output == "png") c5->SaveAs("LF_QAanalysis_Multistrange_page5.png");
        else if (output == "eps") c5->SaveAs("LF_QAanalysis_Multistrange_page5.eps");
        else if (output == "pdf") c5->SaveAs("LF_QAanalysis_Multistrange.pdf");
     }


     //_______________________________________________
     // - PAGE 6 : MC generated particles check
     cout<<"\n--- BUILD THE SEVENTH PAGE: general variables for MC generated particles ---"<<endl;
     if (isMC) { 
          AliCFContainer *cfMC = (AliCFContainer*) (f->Get("PWGLFStrangeness.outputCheckCascade/fCFContCascadeMCgen"));
          TCanvas** c6 = new TCanvas*[4];

          for (Int_t icas = 0; icas < 4; icas++)  {
                c6[icas] = new TCanvas(Form("c6_%i",icas),"",1200,800);
                c6[icas]->Divide(2,3);
                // -- Pad 1: Total Momentum
                c6[icas]->cd(1);
                 myPadSetUp(gPad,kFALSE);
                 hvar18[icas]->Draw("histo");
                 tcas[icas]->Draw();
                // -- Pad 2: Transverse Momentum
                c6[icas]->cd(2);
                 myPadSetUp(gPad,kFALSE);
                 hvar19[icas]->Draw("histo");
                // -- Pad 3: Rapidity (y)
                c6[icas]->cd(3);
                 myPadSetUp(gPad,kFALSE);
                 hvar20[icas]->Draw("histo");
                // -- Pad 4: Pseudo-rapidity (eta)
                c6[icas]->cd(4);
                 myPadSetUp(gPad,kFALSE);
                 hvar21[icas]->Draw("histo");
                // -- Pad 5: Theta
                c6[icas]->cd(5);
                 myPadSetUp(gPad,kFALSE);
                 hvar22[icas]->Draw("histo");
                // -- Pad 6: Phi
                c6[icas]->cd(6);
                 myPadSetUp(gPad,kFALSE);
                 hvar23[icas]->Draw("histo");
                // -- OUTPUT
                if      (output == "png") c6[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page6_%i.png",icas));
                else if (output == "eps") c6[icas]->SaveAs(Form("LF_QAanalysis_Multistrange_page6_%i.eps",icas));
                else if (output == "pdf") {
                    if (icas < 3)  c6[icas]->SaveAs("LF_QAanalysis_Multistrange.pdf");
                    else           c6[icas]->SaveAs("LF_QAanalysis_Multistrange.pdf)");
                }
          }
     }

}




//====================
// - USEFUL FUNCTIONS
//====================
//______________________
Int_t checkExactMinLimit(TH1F *lHist, Float_t limit) {
         Int_t checkOk = -1;
         Int_t lastbin  = lHist->GetNbinsX();
         Int_t content = 0;
         for (Int_t i = 1; i < lastbin; i++) {
              content = 0;
              content = lHist->GetBinContent(i);
              if (content != 0) break; 
         }
         Float_t binEdge = lHist->GetBinLowEdge(i);
         if      (binEdge < limit) checkOk = 1;
         else if (binEdge > limit) checkOk = 2; 
         else                      checkOk = 3;
         return checkOk;
}
//_____________________
Int_t checkExactMaxLimit(TH1F *lHist, Float_t limit) {
         Int_t checkOk = -1;
         Int_t lastbin  = lHist->GetNbinsX();
         Int_t content = 0;
         for (Int_t i = lastbin; i >= 0; i--) {
              content = 0;
              content = lHist->GetBinContent(i);
              if (content != 0) break;  
         } 
         Float_t binEdge = lHist->GetBinLowEdge(i+1);
         if      (binEdge > limit) checkOk = 1;
         else if (binEdge < limit) checkOk = 2;
         else                      checkOk = 3;
         return checkOk;
}
//-----------------------------
void myOptions(Int_t lStat=0) {
        cout << "Set my personal style!" << endl;
        gStyle->SetOptStat("ie");
        gStyle->SetOptStat(kFALSE);
        gStyle->SetOptTitle(kFALSE);
        gStyle->SetOptLogy(kFALSE);
        gStyle->SetFrameLineWidth(2.5);
        gStyle->SetHistLineWidth(2.5);
        // Set gStyle
        int font = 42;
        // From plain
        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameFillColor(0);
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(0);
        gStyle->SetCanvasColor(0);
        gStyle->SetTitleFillColor(0);
        gStyle->SetTitleBorderSize(1);
        gStyle->SetStatColor(0);
        gStyle->SetStatBorderSize(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetDrawBorder(0);
        gStyle->SetTextFont(font);
        gStyle->SetStatFont(font);
        gStyle->SetStatFontSize(0.05);
        gStyle->SetStatX(0.97);
        gStyle->SetStatY(0.98);
        gStyle->SetStatH(0.03);
        gStyle->SetStatW(0.3);
        gStyle->SetTickLength(0.02,"y");
        gStyle->SetEndErrorSize(3);
        gStyle->SetLabelSize(0.05,"xyz");
        gStyle->SetLabelFont(font,"xyz");
        gStyle->SetLabelOffset(0.01,"xyz");
        gStyle->SetTitleFont(font,"xyz");
        gStyle->SetTitleOffset(1.0,"xyz");
        gStyle->SetTitleSize(0.06,"xyz");
        gStyle->SetMarkerSize(1);
        gStyle->SetGridColor(kGray);
        gStyle->SetPalette(1,0);
        if (lStat){
              gStyle->SetOptTitle(1);
              gStyle->SetOptStat(1111);
              gStyle->SetOptFit(1111);
        } else {
              gStyle->SetOptTitle(0);
              gStyle->SetOptStat(0);
              gStyle->SetOptFit(0);
        }
}
//----------------------------------------------------------
void myLatexMakeUp(TLatex *currentLatex, Int_t textfont, Double_t textsize, Int_t textcolor){
       currentLatex->SetNDC();
       currentLatex->SetTextFont(textfont);
       currentLatex->SetTextSize(textsize);
       currentLatex->SetTextColor(textcolor);
       return;
}
//----------------------------------------------------------
void myLineMakeUp(TLine *currentLine, Int_t linecolor, Int_t linestyle, Double_t linewidth){
       currentLine->SetLineColor(linecolor);
       currentLine->SetLineStyle(linestyle);
       currentLine->SetLineWidth(linewidth);
       return;
}
//-----------------------------------------------------------
void myPadSetUp(TPad *currentPad, Bool_t setlogy=kTRUE, Float_t currentLeft=0.11, Float_t currentTop=0.04, Float_t currentRight=0.04, Float_t currentBottom=0.15){
        if (setlogy) currentPad->SetLogy();
        currentPad->SetLeftMargin(currentLeft);
        currentPad->SetTopMargin(currentTop);
        currentPad->SetRightMargin(currentRight);
        currentPad->SetBottomMargin(currentBottom);
        return;
}
//----------------------------------------------------------
void myPaveMakeUp(TPaveText *currentPave, Char_t *texture, Double_t textangle, Int_t textalign, Double_t textsize, Int_t pavefillcolor, Double_t pavebordersize){
       TText *text = currentPave->AddText(Form("%s",texture));
       text->SetTextAngle(textangle);
       text->SetTextAlign(textalign);
       text->SetTextSize(textsize);
       currentPave->SetFillColor(pavefillcolor);
       currentPave->SetFillStyle(0);
       currentPave->SetBorderSize(pavebordersize);
       return;
}
//----------------------------------------------------------
void myHistoMakeUp(TH1F *currentTH1, Int_t linecolor, Double_t labelsizex, Double_t titlesizex, Double_t labelsizey, Double_t titlesizey){
//void myHistoMakeUp(TH1F *currentTH1, Int_t linestyle=1, Int_t linecolor, Int_t markersymbol, Int_t markercolor, Float_t markersize){

//       currentTH1->SetLineStyle(linestyle);
       currentTH1->SetLineColor(linecolor);
       currentTH1->SetLineWidth(1.0);

//       currentTH1->SetMarkerStyle(markersymbol);
//       currentTH1->SetMarkerColor(markercolor);
//       currentTH1->SetMarkerSize(markersize);

       currentTH1->GetXaxis()->SetLabelSize(labelsizex);
//       currentTH1->GetXaxis()->SetLabelOffset(labeloffsetx);
//       currentTH1->GetXaxis()->SetTitle(title);
       currentTH1->GetXaxis()->SetTitleSize(titlesizex);

       currentTH1->GetYaxis()->SetLabelSize(labelsizey);
//       currentTH1->GetYaxis()->SetLabelOffset(labeloffsety);
//       currentTH1->GetYaxis()->SetTitle(title);
       currentTH1->GetYaxis()->SetTitleSize(titlesizey);
       return;
}
