#!/bin/bash

Usage() {
    echo Usage: mrgfl.sh [-f] outFile inpFile1 inpFile2 ....
    echo Option -f will force overwriting existing output file
    exit 1
}
#
if [ $# -lt 2 ] ; then Usage ;fi
#
force=0
if [ "$1" == "-f" ] ; then
    force=1
    shift 1 ;
    if [ $# -lt 2 ] ; then Usage ;fi
fi
#
mrglst="_mrglst.txt"
if [ -e $mrglst ] ; then rm $mrglst ;fi
#
outf=$1
shift 1 ;

while [ $# -gt 0 ] ; do
    echo $1
    shift 1 ;
done > $mrglst

macroName="TmpMergeMacro"
if [ -e ${macroName}.C ] ; then rm ${macroName}.C ;fi
#
#cat > ${macroName}.C << "EOF"
echo \
void ${macroName}'(const char* outFile,const char* mrgLst,Bool_t force) {
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  ifstream in;
  in.open(mrgLst);
  TFileMerger m(1);
  TString line; 
  while (in.good()) {
    in >> line;
    if (line.Length() == 0) continue;
    TString fileName;
    fileName.Form("%s", line.Data());
    Printf("%s", fileName.Data());
    m.AddFile(fileName);
  }
  m.SetFastMethod();
  m.OutputFile(outFile,force);
  m.Merge();
}' > ${macroName}.C

aliroot -b -q ${macroName}.C\(\"${outf}\",\"${mrglst}\",${force}\)

rm ${mrglst}
rm ${macroName}.C
