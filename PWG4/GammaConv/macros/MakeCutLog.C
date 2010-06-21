// provided by Gamma Conversion Group, PWG4, Kathrin Koch, kkoch@physi.uni-heidelberg.de

#include <Riostream>
#include <fstream>
using namespace std;

void MakeCutLog(const char *inputRootFile = "AnalysisResults",const char *path = "./", const char* outputDir="./Output/"){

  fstream outputFile(Form("%sCutSelection.log",outputDir),ios::out);
  if(!outputFile.is_open()){
    cout<<"Problem opening file"<<endl;
    return;
  }

  //  Char_t filename_input1[200] = (Form("%s%s",path,input1));	
  TString filename = Form("%s%s.root",path,inputRootFile);	
  TFile f(filename.Data());  

  TList *directories = f.GetListOfKeys(); // get the list of directories in the file
  
  for(Int_t entFile=0;entFile<directories->GetEntries();entFile++){

    TObject * o = f.Get(directories->At(entFile)->GetName()); // get the object in the base directory

    if(TString(o->IsA()->GetName())=="TDirectoryFile"){ // means that this is a directory (PWG4......)
      
      TDirectory *pwg4dir =(TDirectory*)o;
 
      TString baseDirName = pwg4dir->GetName();
      
      TString reconstructionFlagString = ""; // this is for new scheme where also the flags are coded in numbers in the PWG4.... name

      if(baseDirName.Length()>31){
	reconstructionFlagString = baseDirName(baseDirName.Index("GammaConversion_")+16,8);
      }
      
      TList *pwg4list = pwg4dir->GetListOfKeys(); // list of the yeys inside the base directory

      for(Int_t entHist=0;entHist<pwg4list->GetEntries();entHist++){
	TString name = pwg4list->At(entHist)->GetName();

	if(name.Contains("container")==0){ // does not try to read the container (get errors if tried)
	  TObject * oHist = pwg4dir->Get(pwg4list->At(entHist)->GetName()); // get the object 
	  
	  if(TString(oHist->IsA()->GetName())=="TList"){ // check if the object is a TList
	    
	    TString listname = oHist->GetName();
	    cout<<"Reading: "<<listname.Data()<<endl;
	    
	    TString cutString = listname(listname.Index("_")+1,listname.Length()) + "\n";// get the Cut string from the name


	    outputFile << cutString.Data();
	  }
	}
      }
    }
  }
  outputFile.close();
}
