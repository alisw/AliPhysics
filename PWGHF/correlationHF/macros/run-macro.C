//-*- Mode: C++ -*-
// $Id$

/// @file   run-macro.C
/// @author Matthias.Richter@ift.uib.no
/// @date   2012-04-06
/// @brief  Compile and run a macro
///
/// Dependency header files and classes are automatically searched in the
/// current directory and in all AliRoot libraries

#if defined(__CINT__) && !defined(__MAKECINT__)
TString GetIncludeHeaders(const char* filename, TString& headers, TString& libs, bool loadClass=true);

void run_macro(const char* macro, bool bExecute=true, const char* includePath=NULL, const char* libraries=NULL)
{
  gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include ");
  if (includePath) gSystem->AddIncludePath(includePath);

  TString strLibraries(libraries);
  TObjArray* tokens=strLibraries.Tokenize(" ");
  if (tokens) {
    TIter next(tokens);
    TObject* object=NULL;
    while ((object=next())!=NULL) {
      gSystem->Load(object->GetName());
    }
    delete tokens;
  }

  TString dependencyHeader;
  TString dependencyLibraries;
  GetIncludeHeaders(macro, dependencyHeader, dependencyLibraries, true);

  TString command(macro);
  Int_t error=0;
  command+="+g";
  if (bExecute) {
    gROOT->Macro(command); 
  } else {
    gROOT->LoadMacro(command);
  }
}

TString GetIncludeHeaders(const char* filename, TString& headers, TString& libs, bool loadClass)
{
  // scan the file and add all include headers found by path
  // to the parameter headers
  ifstream input(filename);
  
  if (input.bad()) {
    cerr << "failed to open file " << filename << endl;
    return headers;
  }
  TString line; 
  while (!line.ReadLine(input).eof()) {
    if (!line.Contains("#include") || !line.Contains(".h")) continue;
    line=line(0, line.Index(".h"));line+=".h";
    line.Replace(0, line.Index("#include"), "");
    line.ReplaceAll("#include", "");
    line.ReplaceAll(" ", "");
    line.ReplaceAll("\"", "");
    if (!line.BeginsWith("Ali") && !line.BeginsWith("T")) continue;
    if (gSystem->AccessPathName(line)!=0) {
      // not an include file in the current directory, check if class
      // is available or find library
      line.ReplaceAll(".h","");
      //cout << "checking class " << line << endl;
      if (TClass::GetClass(line)==NULL) {
	TString command;
	TString resfilename(gSystem->TempDirectory()); resfilename+="/findlib.txt";
	command.Form("for lib in $ALICE_ROOT/lib/*/lib*.so; do (nm $lib | grep %s | grep ' T ' | grep Class_Name > /dev/null) && echo $lib > %s; done", line.Data(), resfilename.Data());
	gSystem->Exec(command);
	ifstream resfile(resfilename.Data());
	if (resfile.good()) {
	  TString result;
	  if (!result.ReadLine(resfile).eof()) {
	    Ssiz_t haveSlash=-1;
	    while ((haveSlash=result.First('/'))>=0) result.Replace(0, haveSlash+1, "");
	    if (!libs.Contains(result)) {
	      cout << "loading dependency library '" << result << "' for class '" << line << "'" << endl;
	      gSystem->Load(result);
	      if (!libs.IsNull()) libs+=" ";
	      libs+=result;
	    }
	  }
	  command="rm "; command+=resfilename;
	  gSystem->Exec(command);
	}
      }
    } else {
      if (headers.Contains(line)) {
        if (!headers.BeginsWith(line)) {
          headers.ReplaceAll(line, "");
          if (!headers.IsNull()) headers.Insert(0, " ");
          headers.Insert(0, line);
        }
        continue;
      }
      if (!headers.IsNull()) headers.Insert(0, " ");
      headers.Insert(0, line);
      TString source=line; source.ReplaceAll(".h", ".cxx");
      if (gSystem->AccessPathName(source)==0) {
	GetIncludeHeaders(source, headers, libs);
      }
      GetIncludeHeaders(line, headers, libs);
      if (loadClass && gSystem->AccessPathName(source)==0) {
	line.ReplaceAll(".h", "");
	if (TClass::GetClass(line)==NULL) {
	  source+="+g";
	  gROOT->LoadMacro(source);
	}
      }
    }
  }
  return headers;
}

void run_macro()
{
  cout << endl;
  cout << "run-macro.C" << endl;
  cout << "Compile and run a macro, include header dependencies are searched in all" << endl;
  cout << "AliRoot libraries and the local directory." << endl;
  cout << endl;
  cout << "Usage:" << endl;
  cout << " aliroot -l run-macro.C'(\"macro.C\, execute, \"include\", \"libraries\")'" << endl;
  cout << endl;
  cout << "Parameters:" << endl;
  cout << " macro path (mandatory)" << endl;
  cout << " execute (optional, default 'true')" << endl;
  cout << " include (optional, e.g. -I$ALICE_ROOT/TPC)" << endl;
  cout << " libraries (optional, e.g. libTPCbase.so)" << endl;
  cout << endl;
}

#elif
{
  cerr << "this macro can not be compiled, remove option '+'" << endl;
}
#endif
