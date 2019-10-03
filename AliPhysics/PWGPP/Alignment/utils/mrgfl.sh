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
ver=`date +%s`
macroName="TmpMergeMacro_${ver}"
mrglst="_mrglst_${ver}"
cnt=0
macroNameS="${macroName}"
mrglstS="${mrglst}"
while [ -f "${mrglstS}.txt" ] || [ -f "${macroNameS}.C" ] ; do
    mrglstS="${mrglst}"_"$cnt"
    macroNameS="${macroName}"_"$cnt"
    cnt=$((cnt+1))
done
#
macroName="$macroNameS"
mrglst="${mrglstS}.txt"
#
outf=$1
shift 1 ;
#
while [ $# -gt 0 ] ; do
    echo $1
    shift 1 ;
done > $mrglst

#
#cat > ${macroName}.C << "EOF"
echo \
void ${macroName}'(const char* outFile,const char* mrgLst,Bool_t force) {
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libAlg");
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
    if (line.Contains("alien") && !gGrid) TGrid::Connect("alien://");
    m.AddFile(fileName);
  }
  m.SetFastMethod();
  m.OutputFile(outFile,force);
  m.Merge();
}' > ${macroName}.C

aliroot -b -q ${macroName}.C\(\"${outf}\",\"${mrglst}\",${force}\)

rm ${mrglst}
rm ${macroName}.C
