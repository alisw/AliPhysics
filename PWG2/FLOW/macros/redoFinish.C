// Macro redoFinish.C is typically used after the merging macros (mergeOutput.C or
// mergeOutputOnGrid.C) have been used to produce the merged, large statistics
// file of flow analysis. Results stored in merged file are WRONG because after
// merging the results from small statistics files are trivially summed up in all
// histograms. This is taken into account and corrected for with macro redoFinish.C.
// Another typical use of the macro redoFinish.C is to repeat the call to Finish()
// in all classes, but with different values of some settings which might modify
// the final results (Example: redo the Finish() and apply correction for detector
// effects in QC code because by default this correction is switched off).

// Name of the merged, large statistics file obtained with the merging macros:
TString mergedFileName = "output.root";
// Final output file name holding correct final results for large statistics sample:
TString outputFileName = "AnalysisResults.root";

Bool_t bApplyCorrectionForNUA = kFALSE; // apply correction for non-uniform acceptance
Bool_t bApplyCorrectionForNUAVsM = kFALSE; // apply correction for non-uniform acceptance in each multiplicity bin independently
Bool_t bPropagateErrorAlsoFromNIT = kFALSE; // propagate error also from non-isotropic terms
Bool_t bMinimumBiasReferenceFlow = kTRUE; // store in CRH for reference flow the result obtained wihout rebinning in multiplicity (kTRUE)

void redoFinish()
{
  LoadLibraries();

  TString mergedFileFullPathName(gSystem->pwd());
  mergedFileFullPathName+="/";
  mergedFileFullPathName+=mergedFileName.Data();
  TFile *mergedFile = NULL;
  if(gSystem->AccessPathName(mergedFileFullPathName.Data(),kFileExists))
  {
    cout<<endl;
    cout<<" WARNING: Couldn't find a file: "<<mergedFileName.Data()<<endl;
    cout<<"          in directory "<<gSystem->pwd()<<" !!!!"<<endl;
    cout<<endl;
    exit(0);
  }
  else
  {
    // Create temporarily copy of <mergedFileName> if neccessary:
    if(!(mergedFileName == outputFileName))
    {
      TSystemFile *fileTemp = new TSystemFile(mergedFileFullPathName.Data(),".");
      fileTemp->Copy("mergedAnalysisResultsTemp.root");
      delete fileTemp;
    }
    // Access <mergedFileName>:
    mergedFile = TFile::Open(mergedFileFullPathName.Data(),"UPDATE");
  }

  // Access from <mergedFileName> the merged TDirectoryFile's for each method and from them the lists holding histograms:
  
  TList* mergedFileKeys = mergedFile->GetListOfKeys();
  for(Int_t i=0; i<mergedFileKeys->GetEntries(); i++)
  {
    TDirectory* directory = dynamic_cast<TDirectory*>(mergedFile->Get(mergedFileKeys->At(i)->GetName()));
    if (!directory) continue;

    TList* listTemp = directory->GetListOfKeys();
    for (Int_t icent=0; icent<listTemp->GetEntries(); icent++)
    {
      TList* list = dynamic_cast<TList*>(directory->Get(listTemp->At(icent)->GetName()));
      if (!list) continue;

      ////////////////////
      if(TString(list->GetName()).Contains("MCEP"))
      {
        AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
        mcep->GetOutputHistograms(list);
        mcep->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // SP:
      else if(TString(list->GetName()).Contains("SP"))
      {
        AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
        sp->GetOutputHistograms(list);
        sp->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // GFC:
      else if(TString(list->GetName()).Contains("GFC"))
      {
        AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
        gfc->GetOutputHistograms(list);
        gfc->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // QC:
      else if(TString(list->GetName()).Contains("QC"))
      {
        AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
        qc->GetOutputHistograms(list);
        qc->GetIntFlowFlags()->SetBinContent(3,(Int_t)bApplyCorrectionForNUA);
        qc->GetIntFlowFlags()->SetBinContent(8,(Int_t)bApplyCorrectionForNUAVsM);
        qc->GetIntFlowFlags()->SetBinContent(9,(Int_t)bPropagateErrorAlsoFromNIT);
        qc->GetIntFlowFlags()->SetBinContent(11,(Int_t)bMinimumBiasReferenceFlow);
        qc->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // FQD:
      else if(TString(list->GetName()).Contains("FQD"))
      {
        AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
        fqd->GetOutputHistograms(list);
        fqd->Finish(kTRUE);
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // LYZ1SUM:
      else if(TString(list->GetName()).Contains("LYZ1SUM"))
      {
        AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
        lyz1sum->SetFirstRun(kTRUE);
        lyz1sum->SetUseSum(kTRUE);
        lyz1sum->GetOutputHistograms(list);
        lyz1sum->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // LYZ2SUM:
      else if(TString(list->GetName()).Contains("LYZ2SUM"))
      {
        AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
        lyz2sum->SetFirstRun(kFALSE);
        lyz2sum->SetUseSum(kTRUE);
        lyz2sum->GetOutputHistograms(list);
        lyz2sum->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // LYZ1PROD:
      else if(TString(list->GetName()).Contains("LYZ1PROD"))
      {
        AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
        lyz1prod->SetFirstRun(kTRUE);
        lyz1prod->SetUseSum(kFALSE);
        lyz1prod->GetOutputHistograms(list);
        lyz1prod->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // LYZ2PROD:
      else if(TString(list->GetName()).Contains("LYZ2PROD"))
      {
        AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
        lyz2prod->SetFirstRun(kFALSE);
        lyz2prod->SetUseSum(kFALSE);
        lyz2prod->GetOutputHistograms(list);
        lyz2prod->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // LYZEP:
      else if(TString(list->GetName()).Contains("LYZEP"))
      {
        AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
        lyzep->GetOutputHistograms(list);
        lyzep->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // MH:
      else if(TString(list->GetName()).Contains("MH"))
      {
        AliFlowAnalysisWithMixedHarmonics* mh = new AliFlowAnalysisWithMixedHarmonics();
        mh->GetOutputHistograms(list);
        mh->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
      // NL:
      else if(TString(list->GetName()).Contains("NL"))
      {
        AliFlowAnalysisWithNestedLoops* nl = new AliFlowAnalysisWithNestedLoops();
        nl->GetOutputHistograms(list);
        nl->Finish();
        directory->Add(list,kTRUE);
        directory->Write(directory->GetName(),TObject::kSingleKey+TObject::kWriteDelete);
      }
    }//for icent
  }//

    // Close the final output file:
    delete mergedFile;

    // Giving the final names if neccessary:
    if(!(mergedFileName == outputFileName))
    {
      TSystemFile *outputFileFinal = new TSystemFile(mergedFileName.Data(),".");
      outputFileFinal->Rename(outputFileName.Data());
      delete outputFileFinal;
      TSystemFile *mergedFileFinal = new TSystemFile("mergedAnalysisResultsTemp.root",".");
      mergedFileFinal->Rename(mergedFileName.Data());
      delete mergedFileFinal;
    } // end of if(!(mergedFileName == outputFileName))

    cout<<endl;

} // end of void redoFinish(Int_t mode=mLocal)

void LoadLibraries()
{
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");

  // for AliRoot
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG2flowCommon");
  gSystem->Load("libPWG2flowTasks");
} // end of void LoadLibrariesRF(const libModes mode)


