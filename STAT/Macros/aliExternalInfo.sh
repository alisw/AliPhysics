#!/usr/bin/env bash
# macro to cache production infromation and to create an "materrialized view
# Input:
#    MonALISA production tables - aceeesed using AliExternalInfo interface
# Output:
#    mcRunQATrending.root       - MC detectorQA+production per run
#    mcPeriodQATrending.root    - MC detectorQA+production per period
#    realDataRunTrending.root   - Real data QA+production per run
#    realDataPeriodTrending.root- Real data QA+production per period
# ( source $AliRoot_SRC/STAT/Macros/aliExternalInfo.sh; cacheMCInfo; cacheRealDataInfo )

# cache MC trending trees and create summary trees with run statistic
cacheMCInfo(){
    aliroot -n -b -l <<\EOF 2>&1 | tee cacheMC.log
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C+
     CacheTestMCProductions("QA.TPC",0);
     CacheTestMCProductions("QA.ITS",0);
     CacheTestMCProductions("QA.TRD",0);
    .q
EOF
    # 2.) parse log to get production table info
    echo 'period/C:type/C:nRuns/D:nRunsProd/D:runList/C'  > mcPeriodQATrending.mif
    cat  cacheMC.log |grep "CacheTestMCProductionsPeriod::"  | sed s_".*CacheTestMCProductionsPeriod::\ *"__  | sed s_"\t"_"\""_g >> mcPeriodQATrending.mif
    echo 'period/C:type/C:nRuns/D:nRunsProd/D/:index/D:run/D'  > mcRunQATrending.mif
    cat  cacheMC.log |grep "CacheTestMCProductionsRun::"  | sed s_".*CacheTestMCProductionsRun::\ *"__  | sed s_"\t"_"\""_g >> mcRunQATrending.mif
    #
    # 3.) convert ascii format to root file
    aliroot -n -b -l <<\EOF 2>&1 | tee -a cacheMC.log
    TTree treeRun; treeRun.ReadFile("mcRunQATrending.mif","",'\"');
    TTree treePeriod; treePeriod.ReadFile("mcPeriodQATrending.mif","",'\"');
    TFile *f0 = TFile::Open("mcRunQATrending.root","recreate"); treeRun->Write("trending"); delete f0;
    TFile *f1 = TFile::Open("mcPeriodQATrending.root","recreate"); treePeriod->Write("trending"); delete f1;
EOF
}

# cache Real trending trees and create summary trees with run statistic
# 1.) cache data
# 2.) parse log to get production table info
# 3.) convert ascii format to root file
cacheRealDataInfo(){
# 1.) cache data
    aliroot -n -b -l <<\EOF 2>&1 | tee cacheRealData.log
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C+
        CacheTrendingProductions("QA.TPC");
    .q
EOF
aliroot -n -b -l <<\EOF 2>&1 | tee -a cacheRealData.log
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C+
    CacheTrendingProductions("QA.TRD");
EOF
aliroot -n -b -l <<\EOF 2>&1 | tee -a cacheRealData.log
    .L $AliRoot_SRC/STAT/Macros/aliExternalInfo.C+
    CacheTrendingProductions("QA.ITS");
    .q
EOF
# 2.) parse log to get production table info
   echo 'period/C:pass/C:type/C:nRuns/D:nRunsProd/D:runList/C'  > realDataPeriodTrending.mif
   cat  cacheRealData.log |grep "CacheTrendingProductionsPeriod::"  | sed s_".*CacheTrendingProductionsPeriod::\ *"__ | sed s_"\t"_"\""_g >> realDataPeriodTrending.mif
   echo 'period/C:pass/C:type/C:nRuns/D:nRunsProd/D:run/C'  > realDataRunTrending.mif
   cat  cacheRealData.log |grep "CacheTrendingProductionsRun::"  | sed s_".*CacheTrendingProductionsRun::\ *"__  | sed s_"\t"_"\""_g >> realDataRunTrending.mif
# 3.) convert ascii format to root file
    aliroot -n -b -l <<\EOF 2>&1 | tee -a cacheRealData.log
    TTree treeRun; treeRun.ReadFile("realDataRunTrending.mif","",'\"');
    TTree treePeriod; treePeriod.ReadFile("realDataPeriodTrending.mif","",'\"');
    TFile *f0 = TFile::Open("realDataRunTrending.root","recreate"); treeRun->Write("trending"); delete f0;
    TFile *f1 = TFile::Open("realDataPeriodTrending.root","recreate"); treePeriod->Write("trending"); delete f1;
EOF
}
