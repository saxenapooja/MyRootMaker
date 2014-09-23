#!/bin/bash
date

################################################
# Sunday, April 13, 2014: 12:48 hr             #                                                                                                                          
################################################


usage () {
    echo "@@@@@@ Higgs_tau_tau "
    echo "              ==Possible arguments=="
    echo "  makelib                                 - to build lib of MainAnalyzer"
    echo "  skim                                    - to execute skim"
    echo "  main                                    - to execute main"
    echo "  help                                    - shows this help"
}


SetupEnv() {
#    source $SETUP
    source $LIB
}

MakeLib() {
    cd /nfs/dust/cms/user/anayak/CMS/CMSSW_5_3_14_h2tau/src/MyRootMaker/MyRootMaker/test/ANAL/AnalysisTool
    source AnalysisToolUseThis
    ./configure --prefix=/nfs/dust/cms/user/anayak/CMS/CMSSW_5_3_14_h2tau/src/MyRootMaker/MyRootMaker/test/ANAL/lib_AnalysisTool
    make
    make install
    echo '@@@@@@@@@@@@@ LIBERARIES ARE BUILD @@@@@@@@@@@@@@@@@@@@@@'
}

Skim() {
    make -B main_skim
# test prupose (single sample)
# ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root  /nfs/dust/cms/user/pooja/scratch/plot-macro/armin-file/TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6_*.root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg

# test purpose (huge dataset)
#    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root /nfs/dust/cms/user/pooja/samples/h2tautau-8TeV/TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6_28*root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
    
    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV*root /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
#    ./main_skim ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root  TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6_29.root  /nfs/dust/cms/user/pooja/scratch/plot-macro/tau-hadronic/h2tautau-analysis/Config/skim/common.cfg
}


Main() {
    make -B main
#    ./main skimmed.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main finalskim.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_patch2/src/MyRootMaker/MyRootMaker/test/GluGluToHToTauTau_M-125.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_H2Tau/src/MyRootMaker/MyRootMaker/test/VBF_HToTauTau_M-125.root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6_10.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/vbf_tautau_reproduce/TTFH_MC_VBF_HToTauTau_M-125_*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_9_patch2/src/MyRootMaker/MyRootMaker/test/GluGluToHToTauTau_M-125_full.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/VBF_HToTauTau_M-125_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/GluGluToHToTauTau_M-125_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/SUSYGluGluToHToTauTau_M-120_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg  
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/DYJetsToLL_M-50_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_DYJetsToLL_M-50_TuneZ2star_8TeV_madgraph.root ../Config/analysis/2012.cfg
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/VBF_HToTauTau_M-125_tauPolarOff_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg
#    ./main /nfs/dust/cms/user/pooja/samples/higgs_cp_study/GluGluToHToTauTau_M-125_tauPolarOff_MC_v8_vtxWithBS/GluGluToHToTauTau_M*root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
#    ./main /nfs/dust/cms/user/pooja/scratch/h2tautau/CMSSW_5_3_14/src/MyRootMaker/MyRootMaker/test/GluGluToHToTauTau_M-125_EvntAll.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg
    ./main /nfs/dust/cms/user/pooja/samples/cpStudy_Aug2014/GluGluToHToTauTau_M-125_MC_v8_vtxWithBS_test/GluGluToHToTauTau_M-125*.root ../../Lumi/new/LUMI_INFO_TTFH_MC_GluGluToHToTauTau_M-90_8TeV-powheg-pythia6.root ../Config/analysis/2012.cfg 
}




#if [ $# -lt 1 ]; then
#    usage
#    exit 1
#fi

if [ "$1" = "help" ]; then
    usage

elif [ "$1" = "makelib" ]; then
    SetupEnv
    MakeLib

#main_skim
elif [ "$1" = "skim" ]; then
    SetupEnv
    Skim

    
#main
elif [ "$1" = "main" ]; then
    SetupEnv
    Main
fi

