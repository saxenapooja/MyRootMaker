// Headers files
#include "AnalysisTool/Analyse.h"
#include "Common/Config.h"
#include "Common/PileupWeight.h"
//#include "Common/FakeRateFunctions.h"
#include "Common/HistoCall.h"
//#include "Common/MvaTree.h"
#include "Common/CutFlowTree.h"
#include "Common/MvaEvaluator.h"
#include "Common/TauWP.h"
#include "Common/diTauCand.h"
#include "Common/GenTools.h"
#include "Common/ElectronId.h"
#include "Common/Cutflow.h"
#include "Common/IrredMVA.h"
#include "Common/CombinedMVA.h"
#include "Common/METWeight.h"
#include "Common/MuonId.h"

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <Math/VectorUtil.h>
#include <Math/Vector3D.h>
#include <TROOT.h>

// C & C++
#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include <map> 
using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LVec;
typedef ROOT::Math::XYZVectorD PVec;
typedef ROOT::Math::XYZPointD Point3D;
#define MaxCount 500

// global variables
// doDebug    = false;
bool doFillTree = false;

// class declaration
class MyAnalysis : public Analyse {
public:
  TFile *file;
  TTree *tree;

  Int_t   count;
  Float_t angle;
  Float_t phiStar_gen;//[MaxCount];
  Float_t phiLab_gen; //[MaxCount];
  Float_t CPPhiPrime_gen; //[MaxCount];
  Float_t CPPsiPrime_gen; //[MaxCount];
  Float_t phiLab;
  Float_t CPPhiStar_gen;
  Float_t CPPhiStar_gen_CP;

  Float_t genPionPos_pt;
  Float_t genPionNeg_pt;
  Float_t genPionNeg_Eta;
  Float_t genPionPos_Eta;
  Float_t genTauPos_pt;
  Float_t genTauNeg_pt;
  Float_t genTauNeg_Eta;
  Float_t genTauPos_Eta;

  double tau1_Pt;
  double tau1_Eta;
  double tau1_Phi;
  double tau1_Iso;
  double tau1_RelIso;
  double tau2_Pt;
  double tau2_Eta;
  double tau2_Phi;
  double tau2_Iso;
  double tau2_RelIso;
  double nPV;
  //double met; 
  double tautau_Mass;
  double tautau_Pt;
  double tautau_Eta;
  double tautau_Phi;
  double tautau_dEta;
  double tautau_dR;
  double PVndof;
  double PVz;
  double PVy;
  double PVx;
  double PVr;


private:
  const bool ENABLE_PU_REWEIGHTING;
  const bool ENABLE_SKIM;
  const bool ENABLE_TRIGGER;
  const bool ENABLE_MATCHING; //CP Study at GenLevel
  const bool VERTEX_STUDY;    //CP Study at RecoLevel
  const bool EARLY_MATCHING;
  const bool SAME_SIGN_TAU_PAIRS;
  const bool HIGHPT_SELECTION;
  const bool ONLY_BEST_TRIPLET;
  const bool ONLY_BEST_DOUBLET;
  const bool INCLUSIVE_FAKERATES;
  const bool CAPPED_FAKERATES;
  const std::string DATA_ERA;

  const TauWP OSTauWP;
  const TauWP SSTauWP;

  const diTauCand *ditau;

  const double LEADTAU_PT_CUT;
  const double SUBTAU_PT_CUT;
  const double MET_CUT;
  const double MASS_CUT;

  const double TAU_SCALE;
  const double MET_SCALE;

  const PileupWeight puWeight;
  const METWeight    metWeight;
  const IrredMVA     irredMVA;
  const CombinedMVA  combinedMVA;
  //  const FakeRateFunctions fakeRateFuncsNormal;

private:
  UInt_t currun;
  UInt_t curlumi;

  double cPVNdof;
  double cPVz;
  double cPVr ;

  // define outputfile 
  TFile* histfile;
  std::ofstream evfile;
  std::ofstream outfile;
  std::set<std::string> processedFiles;

  //Event before
  TH1D* hPVN;
  TH1D* hPVN_u;
  TH1D* hPVNdof;
  TH1D* hPVz;
  TH1D* hPVr;
  TH1D *vSubTauPt_;  

  // Cutflow
  CutflowRegions cutflow;
  FillInfo::RegionDef def;
  const CutflowRegions::Cut* cutSkimmed;
  const CutflowRegions::Cut* cutTrigger;
  const CutflowRegions::Cut* cutVtxNDof;
  const CutflowRegions::Cut* cutVtxZ;
  const CutflowRegions::Cut* cutVtxR;
  const CutflowRegions::Cut* cutMuPt;
  const CutflowRegions::Cut* cutMuEta;
  const CutflowRegions::Cut* cutMuID;
  const CutflowRegions::Cut* cutMuNumStations;
  const CutflowRegions::Cut* cutMuNumChambers;
  const CutflowRegions::Cut* cutMuNumLayers;
  const CutflowRegions::Cut* cutMuNumPixel;
  const CutflowRegions::Cut* cutMuChi2;
  const CutflowRegions::Cut* cutMuDxy;
  const CutflowRegions::Cut* cutMuDz;
  const CutflowRegions::Cut* cutMuIso;
  const CutflowRegions::Cut* cutOSTauMatching;
  const CutflowRegions::Cut* cutOSTauPt;
  const CutflowRegions::Cut* cutOSTauEta;
  const CutflowRegions::Cut* cutOSTauCharge;
  const CutflowRegions::Cut* cutOSTauChHadCandPt;
  const CutflowRegions::Cut* cutOSTauDecayModeFinding;
  const CutflowRegions::Cut* cutOSTauIsolation;
  const CutflowRegions::Cut* cutOSTauAntiMu;
  const CutflowRegions::Cut* cutOSTauAntiE;
  const CutflowRegions::Cut* cutOSTauDz;
  const CutflowRegions::Cut* cutOSTauDxy;
  const CutflowRegions::Cut* cutSSTauMatching;
  const CutflowRegions::Cut* cutSSTauPt;
  const CutflowRegions::Cut* cutSSTauEta;
  const CutflowRegions::Cut* cutSSTauCharge;
  const CutflowRegions::Cut* cutSSTauChHadCandPt;  
  const CutflowRegions::Cut* cutSSTauDecayModeFinding;
  const CutflowRegions::Cut* cutSSTauIsolation;
  const CutflowRegions::Cut* cutSSTauAntiMu;
  const CutflowRegions::Cut* cutSSTauAntiE;
  const CutflowRegions::Cut* cutSSTauDz;
  const CutflowRegions::Cut* cutSSTauDxy;
  const CutflowRegions::Cut* cutDrOSTauMu;
  const CutflowRegions::Cut* cutDrSSTauMu;
  const CutflowRegions::Cut* cutDrOSTauSSTau;
  const CutflowRegions::Cut* cutMuonVeto;
  const CutflowRegions::Cut* cutElVeto;
  const CutflowRegions::Cut* cutZVeto;
  const CutflowRegions::Cut* cutBJetVeto;
  const CutflowRegions::Cut* cutMET;
  const CutflowRegions::Cut* cutMt;
  const CutflowRegions::Cut* cutMVA;
  const CutflowRegions::Cut* cutCombinatorics;


  // Event after
  VariableBySelection<TH1D> vNPV;
  VariableBySelection<TH1D> vMET;
  VariableBySelection<TH1D> vMETTYPE1;
  VariableBySelection<TH1D> vNTriplets;
  // Tau
  VariableBySelection<TH1D> vLeadTauPt;
  VariableBySelection<TH1D> vLeadTauEta;
  VariableBySelection<TH1D> vSubTauPt;
  VariableBySelection<TH1D> vSubTauEta;
  VariableBySelection<TH1D> vOSTauPt;
  VariableBySelection<TH1D> vSSTauPt;
  VariableBySelection<TH1D> vOSTauEta;
  VariableBySelection<TH1D> vSSTauEta;
  VariableBySelection<TH1D> vOSTauIso;
  VariableBySelection<TH1D> vSSTauIso;
  VariableBySelection<TH1D> vOSTauRelIso;
  VariableBySelection<TH1D> vSSTauRelIso;
  // Tau Tau
  VariableBySelection<TH1D> vTauTauVisMass;
  VariableBySelection<TH1D> vTauTauPt;
  VariableBySelection<TH1D> vTauTauEta;
  VariableBySelection<TH1D> vTauTauDPhi;
  VariableBySelection<TH1D> vTauTauDEta;
  VariableBySelection<TH1D> vTauTauDR;
  VariableBySelection<TH1D> vTauTauPtRatio;
  VariableBySelection<TH1D> vTauTauAngle;
  VariableBySelection<TH1D> vTauTauSVfitMass;
  // MVA
  VariableBySelection<TH1D> vBDT8_v2;
  VariableBySelection<TH1D> vIrredMVA;
  VariableBySelection<TH1D> vCombinedMVA;
  // Profile
  VariableBySelection<TProfile> pDiTauPtOSTauMuMass;
  VariableBySelection<TProfile> pLeadTauPtSubTauPt;
  VariableBySelection<TProfile> pLeadTauPtMET;
  VariableBySelection<TProfile> pLeadTauPtDeltaR;
  VariableBySelection<TProfile> pLeadTauPtPtRatio;
  VariableBySelection<TProfile> pSubTauPtMET;
  VariableBySelection<TProfile> pSubTauPtDeltaR;
  VariableBySelection<TProfile> pSubTauPtPtRatio;
  VariableBySelection<TProfile> pMETDeltaR;
  VariableBySelection<TProfile> pMETPtRatio;
  VariableBySelection<TProfile> pDeltaRPtRatio;

  // number of processed events
  TH1D* nEvents;
  // pile-up
  Double_t weight;
  // MVA tree
  CutFlowTree mvaTree;
  const char *sample;
  //  const diTau * DITAU;
  
public:
  struct triplet {
    int region;
    unsigned int lep;
    const TauCand* OStau;
    const TauCand* SStau;
    bool lepJetFake;
    bool OSTauJetFake;
    bool SSTauJetFake;
    bool wjetsVetoMT;
    bool OSTauTrackFake;
    bool SSTauTrackFake;
    CutFlowTree::Variables mvaVars;
    double mvaValue;
  };
  
  struct DiTau {
    int region;
    const TauCand* OStau;
    const TauCand* SStau;
  };
  
  MyAnalysis(const Config& config);
  virtual ~MyAnalysis();
  
  // AnalyzeEvent is a virtual function which is called for each event.
  virtual Int_t AnalyzeEvent();
  int  Cutflow_PV(int regions);
  int  Cutflow_Muon(int mu,  int regions);
  int  Cutflow_OSTau(int mu, const TauCand& OStau, int regions);
  int  Cutflow_OSTau(const TauCand& OStau, int regions);
  int  Cutflow_SSTau(int mu, const TauCand& OStau, const TauCand& SStau, int regions);
  int  Cutflow_SSTau(const TauCand& OStau, const TauCand& SStau, int regions);
  int  Cutflow_Combined(int mu, const TauCand& OStau, const TauCand& SStau, int regions);
  int  Cutflow_FinalSelection(int mu, const TauCand& OStau, const TauCand& SStau, int regions);
  int  Cutflow_FinalSelection(const TauCand& OStau, const TauCand& SStau, int regions);
  int  Cutflow_recoTauP(const TauCand& tauP, int regions);
  int  Cutflow_recoTauN(const TauCand& tauN, int regions);
  bool hasMuon(const TauCand& OStau, const Muon& muon, int regions);
  bool hasEle(const TauCand& OStau, const Electron& ele, int regions);
  void SetSample(std::string flag);
  void FillHistogramsFromTriplets(const std::vector<triplet>& triplets, bool mvaCut);
  void FillHistogramsFromDoublets(const std::vector<DiTau>& doublets);
  void FillHistograms(const DiTau& ditau);
  void FillHistograms(const DiTau& ditau,bool mvaCut);
  TLorentzVector MET() const;
  std::map<std::string, TH1F*>  CreateHistograms();

  TTree* syncTree;
  unsigned int run;
  unsigned int lumi;
  unsigned int event;

  float muonPt;
  float muonEta;
  float muonPhi;
  float muonIso;
  int muonCharge;
  bool muonGlobal;
  bool muonTracker;
  bool muonPF;
  unsigned int muonPixelHits;
  unsigned int muonNumStations;
  unsigned int muonNumChamberHits;
  unsigned int muonPixelLayers;
  unsigned int muonStripLayers;
  float muonChi2OverNdof;

  float tauOSPt;
  float tauOSEta;
  float tauOSPhi;
  float tauOSIso3Hits;
  bool tauOSIso3HitsLoose;
  bool tauOSIso3HitsMedium;
  bool tauOSIso3HitsTight;
  bool tauOSAgainstElectronLoose;
  bool tauOSAgainstElectronLooseMVA3;
  bool tauOSAgainstElectronMediumMVA3;
  bool tauOSAgainstElectronTightMVA3;
  bool tauOSAgainstMuonLoose;
  bool tauOSAgainstMuonMedium;
  bool tauOSAgainstMuonTight;

  float tauSSPt;
  float tauSSEta;
  float tauSSPhi;
  float tauSSIso3Hits;
  bool tauSSIso3HitsLoose;
  bool tauSSIso3HitsMedium;
  bool tauSSIso3HitsTight;
  bool tauSSAgainstElectronLoose;
  bool tauSSAgainstElectronLooseMVA3;
  bool tauSSAgainstElectronMediumMVA3;
  bool tauSSAgainstElectronTightMVA3;
  bool tauSSAgainstMuonLoose;
  bool tauSSAgainstMuonMedium;
  bool tauSSAgainstMuonTight;

  float met;
  float metPhi;
  float mt;

};


static double calcMt(double px1, double py1, double px2, double py2)
{
  double pt1 = TMath::Sqrt(px1*px1 + py1*py1);
  double pt2 = TMath::Sqrt(px2*px2 + py2*py2);
  
  double p1Dotp2 = px1*px2 + py1*py2;
  double cosAlpha = p1Dotp2/(pt1*pt2);
  
  return TMath::Sqrt(2*pt1*pt2*(1 - cosAlpha));
}


static const int HISTOGRAMS = Variable<TH1D>::SIGNAL
	| Variable<TH1D>::FAKE_12_WZ
	| Variable<TH1D>::FAKE_12_UNWEIGHTED
	| Variable<TH1D>::FAKE_SS_WZ
	| Variable<TH1D>::FAKE_SS_ERR
	| Variable<TH1D>::FAKE_SS_UNWEIGHTED;

MyAnalysis::MyAnalysis(const Config& config) : Analyse(), currun(0), curlumi(0),
  ENABLE_PU_REWEIGHTING(config.get<bool>("enable_pu_reweighting")),
  ENABLE_SKIM(config.get<bool>("enable_skim", true)),
  ENABLE_TRIGGER(config.get<bool>("enable_trigger", false)),
  ENABLE_MATCHING(config.get<bool>("enable_matching", false)), 
  VERTEX_STUDY(config.get<bool>("vertex_study", true)), 
  EARLY_MATCHING(true),
  SAME_SIGN_TAU_PAIRS(config.get<bool>("same_sign_tau_pairs", false)),
  HIGHPT_SELECTION(config.get<bool>("highpt_selection", false)),
  ONLY_BEST_TRIPLET(config.get<bool>("only_best_triplet", true)),
  ONLY_BEST_DOUBLET(config.get<bool>("only_best_triplet", true)),
  INCLUSIVE_FAKERATES(config.get<bool>("inclusive_fakerates", false)),
  CAPPED_FAKERATES(config.get<bool>("capped_fakerates", false)),
  DATA_ERA(config.get<std::string>("data_era")),

  OSTauWP(config, "tauOS", "loose_3hits", "loose", "loose"), //(config, prefix, default_isolation, antiE, antiMu)
  SSTauWP(config, "tauSS", "loose_3hits", "loose", "loose"),

  LEADTAU_PT_CUT(HIGHPT_SELECTION ? FillInfo::HIGHPT_TAU1_CUT : 40),
  SUBTAU_PT_CUT(HIGHPT_SELECTION ? FillInfo::HIGHPT_TAU2_CUT : 40),
  MET_CUT(HIGHPT_SELECTION ? FillInfo::HIGHPT_MET_CUT : 0),
  MASS_CUT(config.get<double>("mass_cut", 85.0)),

  TAU_SCALE(config.get<double>("tau_scale", 1.0)),
  MET_SCALE(config.get<double>("met_scale", 1.0)),

  puWeight(DATA_ERA.c_str()),
  metWeight(DATA_ERA.c_str(), "mu"),
					       //  fakeRateFuncsNormal(config.get<std::string>("fakerate_directory", "plots").c_str(), DATA_ERA.c_str(), (std::string("wmuJetsSS") + (CAPPED_FAKERATES ? "Capped" : "") + "AgainstElecMu").c_str(), OSTauWP.antie_frfunc_name(), OSTauWP.isolation_frfunc_name(), SSTauWP.antie_frfunc_name(), SSTauWP.isolation_frfunc_name(), true, INCLUSIVE_FAKERATES, CAPPED_FAKERATES),
  //  fakeRateFuncsInclusive(config.get<std::string>("fakerate_directory", "plots").c_str(), DATA_ERA.c_str(), (std::string("wmuJetsSS") + (CAPPED_FAKERATES ? "Capped" : "") + "AgainstElecMu").c_str(), OSTauWP.antie_frfunc_name(), OSTauWP.isolation_frfunc_name(), SSTauWP.antie_frfunc_name(), SSTauWP.isolation_frfunc_name(), true, true, CAPPED_FAKERATES),
  //  mvaEvaluatorBDT8_v1("v1-std", "BDT8", "sane"),
  //  mvaEvaluatorBDT8_v2("v2-std", "BDT8", "sane"),
  //  mvaEvaluatorBDT8_v3("v3-std", "BDT8", "all"),
  //  mvaEvaluatorBDTG_v2("v2-std", "BDTG", "sane"),
  //  mvaEvaluatorBDT8Mt_v2("v2-std", "BDT8", "sane-plus-mt"),
  //  mvaEvaluatorBDTGMt_v2("v2-std", "BDTG", "sane-plus-mt"),
  histfile(new TFile("plots/ditau.root", "RECREATE")),
  def(cutflow),

  vNPV(HISTOGRAMS, "vNPV", "Number of Primary Vertices", "N_{PV}", NULL, "Candidates", 101, -0.5, 100.5),
  vMET(HISTOGRAMS, "vMET", "Missing transverse energy", "MET", "GeV", "Candidates", 350, 0.0, 350.0),
  vMETTYPE1(HISTOGRAMS, "vMETTYPE1", "Missing transverse energy", "METTYPE1", "GeV", "Candidates", 350, 0.0, 350.0),
  vNTriplets(HISTOGRAMS, "vNTriplets", "Number of lep-tau-tau triplets in an event", "N_{Triplets}", "", "Candidates", 101, -0.5, 100.5),
  vLeadTauPt(HISTOGRAMS, "vLeadTauPt", "Leading tau transverse momentum", "p_{T}^{#tau1}", "GeV/c", "Candidates", 400, 0.0, 400.0),
  vLeadTauEta(HISTOGRAMS, "vLeadTauEta", "Leading tau pseudorapidity", "#eta^{#tau1}", NULL, "Candidates", 60, -3.0, 3.0),
  vSubTauPt(HISTOGRAMS, "vSubTauPt", "Subleading tau transverse momentum", "p_{T}^{#tau2}", "GeV/c", "Candidates", 400, 0.0, 400.0),
  vSubTauEta(HISTOGRAMS, "vSubTauEta", "Subleading tau pseudorapidity", "#eta^{#tau2}", NULL, "Candidates", 60, -3.0, 3.0),
  vOSTauPt(HISTOGRAMS, "vOSTauPt", "OS tau transverse momentum", "p_{T}^{#tauOS}", "GeV/c", "Candidates", 400, 0.0, 400.0),
  vSSTauPt(HISTOGRAMS, "vSSTauPt", "SS tau transverse momentum", "p_{T}^{#tauSS}", "GeV/c", "Candidates", 400, 0.0, 400.0),
  vOSTauEta(HISTOGRAMS, "vOSTauEta", "OS tau pseudorapidity", "#eta^{#tauOS}", NULL, "Candidates", 60, -3.0, 3.0),
  vSSTauEta(HISTOGRAMS, "vSSTauEta", "SS tau pseudorapidity", "#eta^{#tauSS}", NULL, "Candidates", 60, -3.0, 3.0),
  vOSTauIso(HISTOGRAMS, "vOSTauIso", "OS tau isolation", "I_{PF}", "GeV/c", "Candidates", 500, 0.0, 50.0),
  vSSTauIso(HISTOGRAMS, "vSSTauIso", "SS tau isolation", "I_{PF}", "GeV/c", "Candidates", 500, 0.0, 50.0),
  vOSTauRelIso(HISTOGRAMS, "vOSTauRelIso", "OS tau relative isolation", "I_{PF}^{rel}", NULL, "Candidates", 500, 0.0, 5.0),
  vSSTauRelIso(HISTOGRAMS, "vSSTauRelIso", "SS tau relative isolation", "I_{PF}^{rel}", NULL, "Candidates", 500, 0.0, 5.0),
  vTauTauVisMass(HISTOGRAMS | Variable<TH1D>::MVA_OPTIMIZATION, "vTauTauVisMass", "Di-tau visible mass", "m^{#tau#tau}", "GeV/c^{2}", "Candidates", 400, 0.0, 400.0),
  vTauTauPt(HISTOGRAMS, "vTauTauPt", "Di-tau transverse momentum", "p_{T}^{#tau#tau}", "GeV/c", "Candidates", 400, 0.0, 400),
  vTauTauEta(HISTOGRAMS, "vTauTauEta", "Di-tau pseudorapidity", "#eta^{#tau#tau}", NULL, "Candidates", 60, -3.0, 3.0),
  vTauTauDPhi(HISTOGRAMS, "vTauTauDPhi", "#Delta #phi between the taus", "#Delta #phi", NULL, "Candidates", 100, -M_PI, M_PI),
  vTauTauDEta(HISTOGRAMS, "vTauTauDEta", "#Delta #eta between the taus", "#Delta #eta", NULL, "Candidates", 100, -5.0, 5.0),
  vTauTauDR(HISTOGRAMS, "vTauTauDR", "#Delta R between the taus", "#Delta R", NULL, "Candidates", 100, 0.0, 10.0),
  vTauTauPtRatio(HISTOGRAMS, "vTauTauPtRatio", "#Di-tau p_{T} over scalar p_{T} sum", "p_{T} ratio", NULL, "Candidates", 100, 0.0, 1.0),
  vTauTauAngle(HISTOGRAMS, "vTauTauAngle", "Strange angle I have no idea what it is", "#alpha", "rad", "Candidates", 100, -M_PI, M_PI),
  vTauTauSVfitMass(HISTOGRAMS, "vTauTauSVfitMass", "Di-tau SVfit mass", "m_{svfit}^{#tau#tau}", "GeV/c^{2}", "Candidates", 400, 0.0, 400.0),
  //vBDT8_v1(HISTOGRAMS, "vBDT8_v1", "BDT8_v1 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  vBDT8_v2(HISTOGRAMS, "vBDT8_v2", "BDT8_v2 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  //vBDT8_v3(HISTOGRAMS, "vBDT8_v3", "BDT8_v3 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  //vBDTG_v2(HISTOGRAMS, "vBDTG_v2", "BDTG_v2 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  //vBDT8Mt_v2(HISTOGRAMS, "vBDT8Mt_v2", "BDT8Mt_v2 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  //vBDTGMt_v2(HISTOGRAMS, "vBDTGMt_v2", "BDTGMt_v2 Output", "BDT Output", NULL, "Candidates", 200, -1.0, 1.0),
  vIrredMVA(HISTOGRAMS, "vIrredMVA", "Irreducible MVA Output", "MVA Output", NULL, "Candidates", 100, 0.0, 1.0),
  vCombinedMVA(HISTOGRAMS, "vCombinedMVA", "Combined MVA Output", "MVA Output", NULL, "Candidates", 200, -1.0, 1.0),
  pDiTauPtOSTauMuMass(HISTOGRAMS, "pDiTauPtOSTauMuMass", "Ditau pT vs. OSTauMu mass", "p_{T}^{#tau#tau}", "GeV/c", "M_{#tauOS,#mu}", 100, 0.0, 100.0),
  pLeadTauPtSubTauPt(HISTOGRAMS, "pLeadTauPtSubTauPt", "Leadtau pT vs. Subtau pT", "p_{T}^{#tau1}", "GeV/c", "p_{T}^{#tau2}", 200, 0.0, 200.0),
  pLeadTauPtMET(HISTOGRAMS, "pLeadTauPtMET", "Leadtau pT vs. MET", "p_{T}^{#tau1}", "GeV/c", "MET", 200, 0.0, 200.0),
  pLeadTauPtDeltaR(HISTOGRAMS, "pLeadTauPtDeltaR", "Leadtau pT vs. DeltaR", "p_{T}^{#tau1}", "GeV/c", "#DeltaR", 200, 0.0, 200.0),
  pLeadTauPtPtRatio(HISTOGRAMS, "pLeadTauPtPtRatio", "Leadtau pT vs. Pt Ratio", "p_{T}^{#tau1}", "GeV/c", "Pt Ratio", 200, 0.0, 200.0),
  pSubTauPtMET(HISTOGRAMS, "pSubTauPtMET", "Subtau pT vs. MET", "p_{T}^{#tau2}", "GeV/c", "MET", 200, 0.0, 200.0),
  pSubTauPtDeltaR(HISTOGRAMS, "pSubTauPtDeltaR", "Subtau pT vs DeltaR", "p_{T}^{#tau2}", "GeV/c", "DeltaR", 200, 0.0, 200.0),
  pSubTauPtPtRatio(HISTOGRAMS, "pSubTauPtPtRatio", "Subtau pT vs. PtRatio", "p_{T}^{#tau2}", "GeV/c", "PtRatio", 200, 0.0, 200.0),
  pMETDeltaR(HISTOGRAMS, "pMETDeltaR", "MET vs. #DeltaR", "MET", "GeV/c", "DeltaR", 200, 0.0, 200.0),
  pMETPtRatio(HISTOGRAMS, "pMETPtRatio", "MET vs. PtRatio", "MET", "GeV/c", "PtRatio", 200, 0.0, 200.0),
  pDeltaRPtRatio(HISTOGRAMS, "pDeltaRPtRatio", "#DeltaR vs. PtRatio", "#DeltaR", NULL, "PtRatio", 100, 0.0, 10.0),
  mvaTree("mvaTree", "MVA Tree")
{
  FillInfo::MASS_CUT = MASS_CUT;
  
  // Set conditional cuts to null
  cutTrigger = cutOSTauMatching    = cutSSTauMatching    = NULL;
  
  // Initialize cutflow
  if(ENABLE_SKIM)                           cutSkimmed       = cutflow.AddCut("Skimmed");
  if(ENABLE_TRIGGER)                        cutTrigger       = cutflow.AddCut("Trigger");
  if(ENABLE_MATCHING && EARLY_MATCHING)     cutOSTauMatching = cutflow.AddCut("OSTau Matching");
  if(ENABLE_MATCHING && EARLY_MATCHING)     cutSSTauMatching = cutflow.AddCut("SSTau Matching");
  cutVtxNDof                                                 = cutflow.AddCut("Primary Vertex NDof");
  cutVtxZ                                                    = cutflow.AddCut("Primary Vertex Z");
  cutVtxR                                                    = cutflow.AddCut("Primary Vertex R");
  
  if(ENABLE_MATCHING && !EARLY_MATCHING)   cutOSTauMatching = cutflow.AddCut("OSTau Matching");
  cutOSTauPt                                                 = cutflow.AddCut("OSTau pT");
  cutOSTauEta                                                = cutflow.AddCut("OSTau eta");
  cutOSTauCharge                                             = cutflow.AddCut("OSTau charge"); 
  cutOSTauChHadCandPt                                        = cutflow.AddCut("OSTau leadChHadCand pt");
  cutOSTauDecayModeFinding                                   = cutflow.AddCut("OSTau Decay Mode Finding");
  cutOSTauIsolation                                          = cutflow.AddCut(TString::Format("OSTau Isolation%s", OSTauWP.isolation_hr_name()));
  cutOSTauDz                                                 = cutflow.AddCut("OSTau Dz");
  cutOSTauDxy                                                = cutflow.AddCut("OSTau Dxy");
  cutOSTauAntiMu                                             = cutflow.AddCut(TString::Format("OSTau AntiMu%s", OSTauWP.antimu_hr_name()));
  cutOSTauAntiE                                              = cutflow.AddCut(TString::Format("OSTau AntiE%s", OSTauWP.antie_hr_name()));

  if(ENABLE_MATCHING && !EARLY_MATCHING)    cutSSTauMatching = cutflow.AddCut("SSTau Matching");
  cutSSTauPt                                                 = cutflow.AddCut("SSTau pT");
  cutSSTauEta                                                = cutflow.AddCut("SSTau eta");
  cutSSTauCharge                                             = cutflow.AddCut("SSTau charge");
  cutSSTauChHadCandPt                                        = cutflow.AddCut("SSTau leadChHadCand pt");
  cutSSTauDecayModeFinding                                   = cutflow.AddCut("SSTau Decay Mode Finding");
  cutSSTauIsolation                                          = cutflow.AddCut(TString::Format("SSTau Isolation%s", SSTauWP.isolation_hr_name()));
  cutSSTauDz                                                 = cutflow.AddCut("SSTau Dz");
  cutSSTauDxy                                                = cutflow.AddCut("OSTau Dxy");
  cutSSTauAntiMu                                             = cutflow.AddCut(TString::Format("SSTau AntiMu%s", SSTauWP.antimu_hr_name()));
  cutSSTauAntiE                                              = cutflow.AddCut(TString::Format("SSTau AntiE%s", SSTauWP.antie_hr_name()));

  //cutDrSSTauMu               = cutflow.AddCut("DeltaR SSTau Muon");
  //cutDrOSTauSSTau            = cutflow.AddCut("DeltaR OSTau SSTau");
  //cutMuonVeto                = cutflow.AddCut("Muon Veto");
  //cutElVeto                  = cutflow.AddCut("Electron Veto");
  // cutZVeto                   = cutflow.AddCut("Z Veto");
  // cutBJetVeto                = cutflow.AddCut("BJet Veto");
  // cutMET                     = cutflow.AddCut("MET");
  // cutMt                      = cutflow.AddCut("Mt");
  // cutMVA                     = cutflow.AddCut("MVA Cut");
  // cutCombinatorics           = cutflow.AddCut("Combinatorics");


  //_________PV cuts___________
  cPVNdof           = 7;
  cPVz              = 24;
  cPVr              = 2;//2;
  
  // load needed informations
  LoadTrigger();
  LoadBeamSpot();
  LoadPrimVertices();
  LoadElectrons();
  LoadMuons();
  LoadTaus();
  LoadTracks();
  LoadAK5PFJets();
  LoadMET();
  LoadAllGenParticles();
  LoadGenParticles();
  LoadGenInfo();

  // output file
  outfile.open("plots/taumu.log");
  evfile.open("missingEvents.log");

  // Event histograms
  hPVN      = new TH1D("hPVN", "Number of PV;#PV;Candidates/1.0", 40, 0., 40.);
  hPVN_u    = new TH1D("hPVN_u", "Number of PV unweighted;#PV;Candidates/1.0", 40, 0., 40.);
  hPVNdof   = new TH1D("hPVNdof", "PV Ndof;Ndf_{PV};Candidates/1.0", 250, 0., 250.);
  hPVz      = new TH1D("hPVz", "PV z;z_{PV} [cm];Candidates/1.0 [cm]", 60, -30., 30.);
  hPVr      = new TH1D("hPVr", "PV r;r_{PV} [cm];Candidates/0.01 [cm]", 300, 0., 3.);
  nEvents   = new TH1D("nEvents", "Number of processed events", 2, -0.5, 1.5);

  syncTree  = new TTree("syncTree", "syncTree");
  syncTree->Branch("run", &run, "run/i");
  syncTree->Branch("lumi", &lumi, "lumi/i");
  syncTree->Branch("event", &event, "event/i");

  syncTree->Branch("muonPt", &muonPt, "muonPt/F");
  syncTree->Branch("muonEta", &muonEta, "muonEta/F");
  syncTree->Branch("muonPhi", &muonPhi, "muonPhi/F");
  syncTree->Branch("muonIso", &muonIso, "muonIso/F");
  syncTree->Branch("muonCharge", &muonCharge, "muonCharge/I");
  syncTree->Branch("muonGlobal", &muonGlobal, "muonGlobal/O");
  syncTree->Branch("muonTracker", &muonTracker, "muonTracker/O");
  syncTree->Branch("muonPF", &muonPF, "muonPF/O");
  syncTree->Branch("muonPixelHits", &muonPixelHits, "muonPixelHits/i");
  syncTree->Branch("muonNumStations", &muonNumStations, "muonNumStations/i");
  syncTree->Branch("muonNumChamberHits", &muonNumChamberHits, "muonNumChamberHits/i");
  syncTree->Branch("muonPixelLayers", &muonPixelLayers, "muonPixelLayers/i");
  syncTree->Branch("muonStripLayers", &muonStripLayers, "muonStripLayers/i");
  syncTree->Branch("muonChi2OverNdof", &muonChi2OverNdof, "muonChi2OverNdof/F");

  syncTree->Branch("tauOSPt", &tauOSPt, "tauOSPt/F");
  syncTree->Branch("tauOSEta", &tauOSEta, "tauOSEta/F");
  syncTree->Branch("tauOSPhi", &tauOSPhi, "tauOSPhi/F");
  syncTree->Branch("tauOSIso3Hits", &tauOSIso3Hits, "tauOSIso3Hits/F");
  syncTree->Branch("tauOSIso3HitsLoose", &tauOSIso3HitsLoose, "tauOSIso3HitsLoose/O");
  syncTree->Branch("tauOSIso3HitsMedium", &tauOSIso3HitsMedium, "tauOSIso3HitsMedium/O");
  syncTree->Branch("tauOSIso3HitsTight", &tauOSIso3HitsTight, "tauOSIso3HitsTight/O");
  syncTree->Branch("tauOSAgainstElectronLoose", &tauOSAgainstElectronLoose, "tauOSAgainstElectronLoose/O");
  syncTree->Branch("tauOSAgainstElectronLooseMVA3", &tauOSAgainstElectronLooseMVA3, "tauOSAgainstElectronLooseMVA3/O");
  syncTree->Branch("tauOSAgainstElectronMediumMVA3", &tauOSAgainstElectronMediumMVA3, "tauOSAgainstElectronMediumMVA3/O");
  syncTree->Branch("tauOSAgainstElectronTightMVA3", &tauOSAgainstElectronTightMVA3, "tauOSAgainstElectronTightMVA3/O");
  syncTree->Branch("tauOSAgainstMuonLoose", &tauOSAgainstMuonLoose, "tauOSAgainstMuonLoose/O");
  syncTree->Branch("tauOSAgainstMuonMedium", &tauOSAgainstMuonMedium, "tauOSAgainstMuonMedium/O");
  syncTree->Branch("tauOSAgainstMuonTight", &tauOSAgainstMuonTight, "tauOSAgainstMuonTight/O");

  syncTree->Branch("tauSSPt", &tauSSPt, "tauSSPt/F");
  syncTree->Branch("tauSSEta", &tauSSEta, "tauSSEta/F");
  syncTree->Branch("tauSSPhi", &tauSSPhi, "tauSSPhi/F");
  syncTree->Branch("tauSSIso3Hits", &tauSSIso3Hits, "tauSSIso3Hits/F");
  syncTree->Branch("tauSSIso3HitsLoose", &tauSSIso3HitsLoose, "tauSSIso3HitsLoose/O");
  syncTree->Branch("tauSSIso3HitsMedium", &tauSSIso3HitsMedium, "tauSSIso3HitsMedium/O");
  syncTree->Branch("tauSSIso3HitsTight", &tauSSIso3HitsTight, "tauSSIso3HitsTight/O");
  syncTree->Branch("tauSSAgainstElectronLoose", &tauSSAgainstElectronLoose, "tauSSAgainstElectronLoose/O");
  syncTree->Branch("tauSSAgainstElectronLooseMVA3", &tauSSAgainstElectronLooseMVA3, "tauSSAgainstElectronLooseMVA3/O");
  syncTree->Branch("tauSSAgainstElectronMediumMVA3", &tauSSAgainstElectronMediumMVA3, "tauSSAgainstElectronMediumMVA3/O");
  syncTree->Branch("tauSSAgainstElectronTightMVA3", &tauSSAgainstElectronTightMVA3, "tauSSAgainstElectronTightMVA3/O");
  syncTree->Branch("tauSSAgainstMuonLoose", &tauSSAgainstMuonLoose, "tauSSAgainstMuonLoose/O");
  syncTree->Branch("tauSSAgainstMuonMedium", &tauSSAgainstMuonMedium, "tauSSAgainstMuonMedium/O");
  syncTree->Branch("tauSSAgainstMuonTight", &tauSSAgainstMuonTight, "tauSSAgainstMuonTight/O");

  syncTree->Branch("met", &met, "met/F");
  syncTree->Branch("metPhi", &metPhi, "metPhi/F");
  syncTree->Branch("mt", &mt, "mt/F");


//   // output file
//   outfile.open("plots/ditau.log");

//   // Event histograms
//   vSubTauPt_ = new TH1D("vSubTauPt_","Subleading tau transverse momentum", 100, 0.0, 400.0);
//   hPVN      = new TH1D("hPVN", "Number of PV;#PV;Candidates/1.0", 40, 0., 40.);
//   hPVN_u    = new TH1D("hPVN_u", "Number of PV unweighted;#PV;Candidates/1.0", 40, 0., 40.);
//   hPVNdof   = new TH1D("hPVNdof", "PV Ndof;Ndf_{PV};Candidates/1.0", 250, 0., 250.);
//   hPVz      = new TH1D("hPVz", "PV z;z_{PV} [cm];Candidates/1.0 [cm]", 60, -30., 30.);
//   hPVr      = new TH1D("hPVr", "PV r;r_{PV} [cm];Candidates/0.01 [cm]", 300, 0., 3.);
//   nEvents   = new TH1D("nEvents", "Number of processed events", 2, -0.5, 1.5);

//   file = new TFile("ditauplots.root","RECREATE");
//   tree = new TTree("tree","A Tree with Events");

//   //branches

//   if(ENABLE_MATCHING) {
//     tree->Branch("count" ,           &count ,  "count/i");
//     tree->Branch("angle" ,           &angle ,  "angle/f");
//     tree->Branch("phiStar_gen",&phiStar_gen,"phiStar_gen/f");
//     tree->Branch("phiLab_gen",&phiLab_gen,"phiLab_gen/f");
//     tree->Branch("CPPhiStar_gen",&CPPhiStar_gen,"CPPhiStar_gen/f");
//     tree->Branch("CPPhiStar_gen_CP",&CPPhiStar_gen_CP,"CPPhiStar_gen_CP/f");
//     tree->Branch("CPPhiPrime_gen",&CPPhiPrime_gen,"CPPhiPrime_gen/f");
//     tree->Branch("CPPsiPrime_gen",&CPPsiPrime_gen,"CPPsiPrime_gen/f");
//     tree->Branch("CPPsiPrime_gen",&CPPsiPrime_gen,"CPPsiPrime_gen/f");
//     tree->Branch("phiLab",&phiLab,"phiLab/f");
//     tree->Branch("genPionNeg_Eta",&genPionNeg_Eta,"genPionNeg_Eta/f");
//     tree->Branch("genPionPos_Eta",&genPionPos_Eta,"genPionPos_Eta/f");
//     tree->Branch("genPionPos_pt",&genPionPos_pt,"genPionPos_pt/f");
//     tree->Branch("genPionNeg_pt",&genPionNeg_pt,"genPionNeg_pt/f");
//     tree->Branch("genTauPos_pt",&genTauPos_pt,"genTauPos_pt/f");
//     tree->Branch("genTauNeg_pt",&genTauNeg_pt,"genTauNeg_pt/f");
//     tree->Branch("genTauPos_Eta",&genTauPos_Eta,"genTauPos_Eta/f");
//     tree->Branch("genTauNeg_Eta",&genTauNeg_Eta,"genTauNeg_Eta/f");
//   }
//   tree->Branch("tau1_Pt",&tau1_Pt,"tau1_Pt/D");
//   tree->Branch("tau1_Eta",&tau1_Eta,"tau1_Eta/D");
//   tree->Branch("tau1_Phi",&tau1_Phi,"tau1_Phi/D");
//   tree->Branch("tau1_Iso",&tau1_Iso,"tau1_Iso/D");
//   tree->Branch("tau1_RelIso",&tau1_RelIso,"tau1_RelIso/D");
//   tree->Branch("tau2_Pt",&tau2_Pt,"tau2_Pt/D");
//   tree->Branch("tau2_Eta",&tau2_Eta,"tau2_Eta/D");
//   tree->Branch("tau2_Phi",&tau2_Phi,"tau2_Phi/D");
//   tree->Branch("tau2_Iso",&tau2_Iso,"tau2_Iso/D");
//   tree->Branch("tau2_RelIso",&tau2_RelIso,"tau2_RelIso/D");

//   tree->Branch("nPV",&nPV,"nPV/D");
//   tree->Branch("met",&met,"met/D");
//   tree->Branch("tautau_Mass",&tautau_Mass,"tautau_Mass/D");
//   tree->Branch("tautau_Pt",&tautau_Pt,"tautau_Pt/D");
//   tree->Branch("tautau_Eta",&tautau_Eta,"tautau_Eta/D");
//   tree->Branch("tautau_Phi",&tautau_Phi,"tautau_Phi/D");
//   tree->Branch("tautau_dEta",&tautau_dEta,"tautau_dEta/D");
//   tree->Branch("tautau_dR",&tautau_dR,"tautau_dR/D");
//   tree->Branch("PVndof",&PVndof,"PVndof/D");
//   tree->Branch("PVz",&PVz,"PVz/D");
//   tree->Branch("PVy",&PVy,"PVy/D");
//   tree->Branch("PVx",&PVx,"PVx/D");
//   tree->Branch("PVr",&PVr,"PVr/D");
}


// Destructor:
MyAnalysis::~MyAnalysis()
{
  TH1D* hCutflowSignal            = cutflow.GetCutflowHistogram("Signal");
  TH1D* hCutflowBackground        = cutflow.GetCutflowHistogram("Background");
  TH1D* hCutflowControl           = cutflow.GetCutflowHistogram("Control");
  TH1D* hCutflowControlBackground = cutflow.GetCutflowHistogram("ControlBackground");
  
  hCutflowSignal->SetBinContent(1, nEvents->GetBinContent(1));
  cout<<"sig->events :"<< (nEvents->GetBinContent(1)) << " getNbins : "<< (hCutflowSignal->GetNbinsX())<< endl;
  hCutflowBackground->SetBinContent(1, nEvents->GetBinContent(1));
  hCutflowControl->SetBinContent(1, nEvents->GetBinContent(1));
  hCutflowControlBackground->SetBinContent(1, nEvents->GetBinContent(1));

  std::cout << std::setw(30) << "Selection" << std::setw(25) << std::setprecision(3) << "Events" << std::setw(25) << "Efficiency [%]" << std::endl;
  for(unsigned int i = 2; i <= hCutflowSignal->GetNbinsX(); ++i)
    {
      if(!hCutflowSignal->GetXaxis()->GetBinLabel(i) || !*hCutflowSignal->GetXaxis()->GetBinLabel(i)) break;
      std::cout << std::setw(30)
		<< hCutflowSignal->GetXaxis()->GetBinLabel(i)
		<< std::setw(25)
		<< (unsigned int)hCutflowSignal->GetBinContent(i)
		<< std::setw(25)
		<< 100. * hCutflowSignal->GetBinContent(i) / hCutflowSignal->GetBinContent(i-1)
		<< std::endl;
    }
  
  gROOT->GetListOfFiles()->Remove(histfile); // Trick
  histfile->Write();
  histfile->Close();
  //   file->Write();
  //   file->Close();
}


// Analysis
Int_t MyAnalysis::AnalyzeEvent()
{
  if(doDebug) cout<<"-----------------------------------------------------------------------------------------------"<< event << endl;
  count = 0;
  
  // We start in all regions
  const int ALL_REGIONS = cutflow.AllRegions(); //15
  
  // Event counting
  if(processedFiles.find(GetCurrentFileName()) == processedFiles.end())
    {
      nEvents->SetBinContent(1, nEvents->GetBinContent(1) + GetOrigEvents());
      processedFiles.insert(GetCurrentFileName());
    }
  
  nEvents->Fill(1.0);
  if(doDebug)  cout<<"CutSkimmed cut :"<< ALL_REGIONS<< endl;
  if(ENABLE_SKIM)  cutflow.Pass(cutSkimmed, ALL_REGIONS);
  
  // Trigger
  bool OK_Trigger = false;
  Int_t mtrigger = GetTriggerSelection("TauTrigger")->Result();
  // cout<< "name :"<< (GetTriggerSelection("TauTrigger")->GetTriggerName()) << endl;
  if(mtrigger == 0)
    std::cout << "WARNING!!! No trigger in run=" << Run() << ", lumi=" << LumiBlock() << std::endl;
  if(abs(mtrigger) > 1)
    std::cout << "WARNING!!! The trigger is prescaled!" << std::endl;
  if(!ENABLE_TRIGGER || mtrigger == 1) OK_Trigger = true;
  if(doDebug)  cout<<"Trigger condition is passed with : "<< OK_Trigger << endl;

  std::vector<TauCand>   OStaus, SStaus ;
  std::vector<TauCand> L1taus, L2taus; 
  std::vector<Muon>      muons;
  std::vector<triplet>   triplets;
  std::vector<triplet>   triplets_mva;
  std::vector<DiTau>     doublet, finaldoublet;
  
  // Trigger pass
  if(OK_Trigger)
    {
      if(doDebug)  cout<<"Ok_trigger with region :"<< ALL_REGIONS << endl;
      if(ENABLE_TRIGGER) cutflow.Pass(cutTrigger, ALL_REGIONS); //15, all_cuts displayed are 43
      weight = 1.0f;
      if(doDebug)    cout<<"Trigger pass, will proceed further"<< endl;
      
      if(ENABLE_PU_REWEIGHTING && strcmp(sample, "Data") != 0)
	weight = puWeight.getWeight(NumTruePileUpInteractions());
      
      hPVN->Fill(NumGoodPrimVertices(), weight);
      hPVN_u->Fill(NumGoodPrimVertices());
      
      bool haveEarlyMatch = false;
      int OSTauIndex = -1, SSTauIndex = -1;
      int tauPIndex = -1, tauNIndex = -1;      
      
      if(ENABLE_MATCHING)
	{
	  // First, find generator level particles for the matching
	  GenParticle GenTau, GenTauPos, GenTauNeg;
	  vector<GenParticle> GenTausFromHiggs;
	  bool TauFound = false, PionPosFound = false, PionNegFound = false, TauPosFound = false, TauNegFound = false;
	  TVector3 GenHiggs_Vtx = (0, 0, 0);
	  
	  if(doDebug) cout<<"All genParticle size() : " << NumAllGenParticles() << endl;  

	  // loop over allgenparticles
	  for(unsigned int i = 0 ; i < NumAllGenParticles() ; i++)
	    {
	      doFillTree = false;
	      if(doDebug)	      cout<<i <<"  ---th genParticle "<< endl;
	      GenParticle Gen = AllGenParticles(i);
	      if(Gen.Status() != 3) continue;
	      
	      // Select SM/BSM higgs
	      if(TMath::Abs(Gen.PDGId()) == 23 ||  TMath::Abs(Gen.PDGId()) == 35 || TMath::Abs(Gen.PDGId()) == 36 || TMath::Abs(Gen.PDGId()) == 25)
		{
		  if(doDebug)  cout<<"Higgs found -> " <<Gen.PDGId()<< " No of daugthers are :" <<  Gen.NumDaughters()<< "with pdgid : "<< Gen.PDGId()  << endl;
		  
		  // Taus from Higgs
		  std::vector<GenParticle> GenTausFromHiggs(gen::GetGenTausFromAllHiggs(Gen));
		  if(doDebug) cout<<" GenTausFromHiggs->size : "<< GenTausFromHiggs.size() << endl;
		  if(GenTausFromHiggs.size() == 0) continue;
		  if(GenTausFromHiggs.size() != 2) cout<<"MORE HIGGS CHILDREN MODE "<< GenTausFromHiggs.size() << endl;

		  assert( GenTausFromHiggs.size() == 2);
		  
		  // rest frame of higgs, provided the four-momenta of lab-frame
		  ROOT::Math::Boost Boost_RestFrame_Higgs(Gen.P4().BoostToCM());
		  GenHiggs_Vtx =  Gen.Vertex();		  

		  // initialiazation
		  LVec genTauP_LabFrame(0, 0, 0, 0);  LVec genTauN_LabFrame(0, 0, 0, 0);
		  LVec genPionP_LabFrame(0, 0, 0, 0); LVec genPionN_LabFrame(0, 0, 0, 0);
		  //LVec genNuP_LabFrame(0, 0, 0, 0); LVec genNuN_LabFrame(0, 0, 0, 0);
		  TVector3 genPionP_Vtx(0,0,0);  TVector3 genPionN_Vtx(0,0,0);
		  TVector3 GenTauP_Vtx(0,0,0);  TVector3 GenTauN_Vtx(0,0,0);
		  GenParticle GenNPion, GenPPion, GenNuN, GenNuP;
		  
		  for(std::vector<GenParticle>::iterator genTau=GenTausFromHiggs.begin(); genTau != GenTausFromHiggs.end(); genTau++)
		    {
		      // Lab-frame of taus provided the four-momenta
		      if((*genTau).Charge() == -1) {
			//cout<<"charge of -tau is :"<< (*genTau).Charge() <<" pdgid : "<< (*genTau).PDGId() << endl;
			assert((*genTau).PDGId() == 15);
			genTauN_LabFrame = (*genTau).P4();
			GenTauN_Vtx = (*genTau).Vertex();
			GenNPion = gen::GetTauPionDecayProduct(*genTau);
			GenNuN = gen::GetTauInvisibleDecayProduct(*genTau);
			genTauNeg_pt = (*genTau).Pt();
			genTauNeg_Eta = (*genTau).Eta();
		      }
		      
		      if((*genTau).Charge() == +1)   {
			assert((*genTau).PDGId() == -15);
			genTauP_LabFrame = (*genTau).P4();
			GenTauP_Vtx = (*genTau).Vertex();
			GenPPion = gen::GetTauPionDecayProduct(*genTau);
			GenNuP = gen::GetTauInvisibleDecayProduct(*genTau);
			genTauPos_pt = (*genTau).Pt();
			genTauPos_Eta = (*genTau).Eta();
		      }
		    } // for(std::vector<GenParticle>::iterator genTau=GenTausFromHiggs.begin(); genTau != GenTausFromHiggs.end(); genTau++)
		  
		  // check point
		  if(GenNPion.PDGId()  != -211 || GenPPion.PDGId()  != 211) continue;
		  assert( GenNPion.PDGId() == -211 && GenPPion.PDGId() == 211);
		  PionPosFound = true; PionNegFound = true;		  
		  if(doDebug) cout<<"Pair of pi+pi- found"<< endl;

		  genPionPos_pt  = GenPPion.Pt();
		  genPionNeg_pt  = GenNPion.Pt();
		  genPionPos_Eta = GenPPion.Eta();
		  genPionNeg_Eta = GenNPion.Eta();

		  // pion vertex
		  genPionP_Vtx = GenPPion.Vertex();
		  genPionN_Vtx = GenNPion.Vertex();
		  
		  // lab-frame of pions provided the four momenta
		  genPionP_LabFrame = GenPPion.P4();
		  genPionN_LabFrame = GenNPion.P4();
		  
		  // rest-frame of tau+/ tau- system
		  ROOT::Math::Boost Boost_RestFrame_TauP(genTauP_LabFrame.BoostToCM());
		  ROOT::Math::Boost Boost_RestFrame_TauN(genTauN_LabFrame.BoostToCM());
		  
		  // four momentum of taus in the higgs rest frame 
		  LVec genTauP_RestFrame = genTauP_LabFrame.pt() > 0 ? Boost_RestFrame_Higgs(genTauP_LabFrame) : LVec(0, 0, 0, 0);
		  LVec genTauN_RestFrame = genTauN_LabFrame.pt() > 0 ? Boost_RestFrame_Higgs(genTauN_LabFrame) : LVec(0, 0, 0, 0);
		  
		  // four momentum of pions in higgs rest frame
		  LVec genPionP_RestFrame  = genPionP_LabFrame.pt() > 0 ? Boost_RestFrame_TauP(genPionP_LabFrame) : LVec(0, 0, 0, 0);
		  LVec genPionN_RestFrame  = genPionN_LabFrame.pt() > 0 ? Boost_RestFrame_TauN(genPionN_LabFrame) : LVec(0, 0, 0, 0);
		  
		  if( genTauP_RestFrame.pt()  > 0 && genTauN_RestFrame.pt()  > 0 && 
		      genPionP_RestFrame.pt() > 0 && genPionN_RestFrame.pt() > 0 ) 
		    { 
		      
		      // 3-momentum of taus/pions in higgs rest-frame
		      PVec genTauP_RestFrame3d  =  genTauP_RestFrame.Vect();
		      PVec genTauN_RestFrame3d  =  genTauN_RestFrame.Vect();
		      PVec genPionP_RestFrame3d =  genPionP_RestFrame.Vect();
		      PVec genPionN_RestFrame3d =  genPionN_RestFrame.Vect();
		      
		      // Get the normal to the decay plane
 		      PVec dPlaneNormal_1 = genPionP_RestFrame3d.Cross(genTauP_RestFrame3d);
 		      PVec dPlaneNormal_2 = genTauN_RestFrame3d.Cross(genPionN_RestFrame3d);
		      PVec dPlaneNormal_unitVec1 = dPlaneNormal_1/sqrt(dPlaneNormal_1.mag2());
		      PVec dPlaneNormal_unitVec2 = dPlaneNormal_2/sqrt(dPlaneNormal_2.mag2());
		      
		      // Get CP angle in ZMF
		      double Phi_GenCP = TMath::ACos(dPlaneNormal_unitVec1.Dot(dPlaneNormal_unitVec2));
		      if(doDebug) cout<<"angle is :"<< Phi_GenCP << endl;
		      angle = Phi_GenCP;
		      doFillTree = true;
		      if(doDebug) cout<<"----------------------------"<< endl;
		    }
		  
		  // phi_Star in pi-pi+ rest-frame
		  if(genPionP_LabFrame.pt() > 0 && genPionN_LabFrame.pt() > 0)
		    {
		      // creating a boost in rest frame of pi+pi-
		      if(doDebug) cout<<" Calculating the phistar in dipion rest frame"<< endl;
		      LVec gen_diPion_LabFrame  = genPionP_LabFrame + genPionN_LabFrame;
		      ROOT::Math::Boost Boost_RestFrame_diPion(gen_diPion_LabFrame.BoostToCM());
    
		      // calculating impact parameter vectors n1, n2
		      PVec pionP = genPionP_LabFrame.Vect();
		      PVec pionN = genPionN_LabFrame.Vect();
		      PVec TauP = genTauP_LabFrame.Vect();
		      PVec TauN = genTauN_LabFrame.Vect();
		      
		      PVec n1 = TauP - ((TauP.Dot(pionP)) / ( pionP.Dot(pionP))) * pionP;
		      PVec n2 = TauN - ((TauN.Dot(pionN)) / ( pionN.Dot(pionN))) * pionN;
		      n1 = n1.Unit();
		      n2 = n2.Unit();
		      phiLab = TMath::ACos(n1.Dot(n2));

 		      // Boosting 4-vectors (n1,0), (n2,0), p1, p2 with rest_frame_boost
		      LVec n1_mu(n1.X(), n1.Y(), n1.Z(), 0); //tau+
		      LVec n2_mu(n2.X(), n2.Y(), n2.Z(), 0); //tau-
		      n1_mu = Boost_RestFrame_diPion(n1_mu); 
		      n2_mu = Boost_RestFrame_diPion(n2_mu);

		      // pions in dipion rest-frame
		      LVec p1 = Boost_RestFrame_diPion(genPionP_LabFrame); //pi+
		      LVec p2 = Boost_RestFrame_diPion(genPionN_LabFrame); //pi-

 		      // Calculation of the transverse component of n1, n2 to p1, p2 (after Boosting)                                 
		      PVec n1_rf(n1_mu.Px(), n1_mu.Py(), n1_mu.Pz());
		      PVec n2_rf(n2_mu.Px(), n2_mu.Py(), n2_mu.Pz());
		      PVec p1_rf(p1.Px(),p1.Py(), p1.Pz());
		      PVec p2_rf(p2.Px(),p2.Py(), p2.Pz());

 		      PVec n1t = n1_rf - ((n1_rf.Dot(p1_rf)) / (p1_rf.Dot(p1_rf))) * p1_rf;
 		      PVec n2t = n2_rf - ((n2_rf.Dot(p2_rf)) / (p2_rf.Dot(p2_rf))) * p2_rf;

 		      n1t = n1t.Unit();
 		      n2t = n2t.Unit();
 		      PVec p1n = p1_rf.Unit();

		      if(doDebug)     cout <<  acos(n1t.Dot(p1)) << "             " << acos(n2t.Dot(p2)) << std::endl;
		      if(doDebug)     cout<< acos(p1n.Dot(n1t.Cross(n2t))) << endl;

 		      // Calculating Phi* and Psi*CP
		      phiLab_gen  = acos(n1t.Dot(n2t));
 		      phiStar_gen = acos(p1n.Dot(n2t.Cross(n1t)));


		      //Alternative method used by Arun
		      LVec gen_pion_pair_lab = genPionP_LabFrame + genPionN_LabFrame;
		      ROOT::Math::Boost boost_to_rf_gen(gen_pion_pair_lab.BoostToCM());
		      
		      if(doDebug) cout<<"Higgs Vertex ("<< GenHiggs_Vtx.x()<<","<< GenHiggs_Vtx.y()<<","<< GenHiggs_Vtx.z()<<")"<<std::endl;
		      if(doDebug) cout<<"lab frame : pionP pt = "<<genPionP_LabFrame.pt()<<" pionN pt "<<genPionN_LabFrame.pt()<<std::endl;
		      
		      //boost to rf of two pion system
		      LVec genpionP_rf_2p = boost_to_rf_gen(genPionP_LabFrame);
		      LVec genpionN_rf_2p = boost_to_rf_gen(genPionN_LabFrame);
		      if(doDebug) cout<<"p1+p2 rest frame : pionP pt = "<<genpionP_rf_2p.pt()<<" pionN pt "<<genpionN_rf_2p.pt()<<std::endl;
		      
		      //get unit momentum vectors
		      PVec genpionP_lab_3d = genPionP_LabFrame.Vect();
		      PVec genpionN_lab_3d = genPionN_LabFrame.Vect();
		      PVec genpionP_lab_3d_u = genpionP_lab_3d/sqrt(genpionP_lab_3d.mag2());
		      PVec genpionN_lab_3d_u = genpionN_lab_3d/sqrt(genpionN_lab_3d.mag2());
		      
		      if(doDebug) cout<<"genTau+ Vertex ("<<GenTauP_Vtx.x()<<","<<GenTauP_Vtx.y()<<","<<GenTauP_Vtx.z()<<")"<<std::endl;
		      if(doDebug) cout<<"genTau- Vertex ("<<GenTauN_Vtx.x()<<","<<GenTauN_Vtx.y()<<","<<GenTauN_Vtx.z()<<")"<<std::endl;
		      
		      if(doDebug) cout<<"genPion+ Vertex ("<<genPionP_Vtx.x()<<","<<genPionP_Vtx.y()<<","<<genPionP_Vtx.z()<<")"<<std::endl;
		      if(doDebug) cout<<"genPion- Vertex ("<<genPionN_Vtx.x()<<","<<genPionN_Vtx.y()<<","<<genPionN_Vtx.z()<<")"<<std::endl;
		      
		      //Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
		      double tP = ( genpionP_lab_3d_u.x()*(GenHiggs_Vtx.x() - genPionP_Vtx.x()) + 
				    genpionP_lab_3d_u.y()*(GenHiggs_Vtx.y() - genPionP_Vtx.y()) + 
				    genpionP_lab_3d_u.z()*(GenHiggs_Vtx.z() - genPionP_Vtx.z()) );
		      double tN = ( genpionN_lab_3d_u.x()*(GenHiggs_Vtx.x() - genPionN_Vtx.x()) +  
				    genpionN_lab_3d_u.y()*(GenHiggs_Vtx.y() - genPionN_Vtx.y()) + 	
				    genpionN_lab_3d_u.z()*(GenHiggs_Vtx.z() - genPionN_Vtx.z()) );
		      
		      Point3D pcaP(genPionP_Vtx.x() + genpionP_lab_3d_u.x()*tP, genPionP_Vtx.y() + genpionP_lab_3d_u.y()*tP, genPionP_Vtx.z() + genpionP_lab_3d_u.z()*tP);
		      Point3D pcaN(genPionN_Vtx.x() + genpionN_lab_3d_u.x()*tN, genPionN_Vtx.y() + genpionN_lab_3d_u.y()*tN, genPionN_Vtx.z() + genpionN_lab_3d_u.z()*tN);
		      
		      //if(doDebug) cout<<"PCA Pion+ ("<<pcaP.x()<<","<<pcaP.y()<<","<<pcaP.z()<<")"<<std::endl;
		      //if(doDebug) cout<<"PCA Pion- ("<<pcaN.x()<<","<<pcaN.y()<<","<<pcaN.z()<<")"<<std::endl;
		      
		      //Get the normalized IP vectors
		      PVec ipvP_lab_gen_3d(pcaP.x() - GenHiggs_Vtx.x(), pcaP.y() - GenHiggs_Vtx.y(), pcaP.z() - GenHiggs_Vtx.z());
		      PVec ipvN_lab_gen_3d(pcaN.x() - GenHiggs_Vtx.x(), pcaN.y() - GenHiggs_Vtx.y(), pcaN.z() - GenHiggs_Vtx.z());
		      
		      if(doDebug) cout<<"IPV lab Pion+ ("<<ipvP_lab_gen_3d.x()<<","<<ipvP_lab_gen_3d.y()<<","<<ipvP_lab_gen_3d.z()<<")"<<std::endl;
		      if(doDebug) cout<<"IPV lab Pion- ("<<ipvN_lab_gen_3d.x()<<","<<ipvN_lab_gen_3d.y()<<","<<ipvN_lab_gen_3d.z()<<")"<<std::endl;
		      
		      LVec ipvP_lab_gen(ipvP_lab_gen_3d.x()/sqrt(ipvP_lab_gen_3d.mag2()), 
					ipvP_lab_gen_3d.y()/sqrt(ipvP_lab_gen_3d.mag2()), 
					ipvP_lab_gen_3d.z()/sqrt(ipvP_lab_gen_3d.mag2()), 0);
		      LVec ipvN_lab_gen(ipvN_lab_gen_3d.x()/sqrt(ipvN_lab_gen_3d.mag2()), 
					ipvN_lab_gen_3d.y()/sqrt(ipvN_lab_gen_3d.mag2()), 
					ipvN_lab_gen_3d.z()/sqrt(ipvN_lab_gen_3d.mag2()), 0);
		      
		      //boost to ZMF
		      LVec ipvP_rf_gen = boost_to_rf_gen(ipvP_lab_gen);
		      LVec ipvN_rf_gen = boost_to_rf_gen(ipvN_lab_gen);
		      
		      //Get only the position component of the vector
		      PVec ipvP_rf_gen_3d = ipvP_rf_gen.Vect();
		      PVec ipvN_rf_gen_3d = ipvN_rf_gen.Vect();
		      
		      //if(doDebug) cout<<"IPV Pion+ ("<<ipvP_rf_gen_3d.x()<<","<<ipvP_rf_gen_3d.y()<<","<<ipvP_rf_gen_3d.z()<<")"<<std::endl;
		      //if(doDebug) cout<<"IPV Pion- ("<<ipvN_rf_gen_3d.x()<<","<<ipvN_rf_gen_3d.y()<<","<<ipvN_rf_gen_3d.z()<<")"<<std::endl;
		      
		      PVec genpionP_rf_2p_3d = genpionP_rf_2p.Vect();
		      PVec genpionN_rf_2p_3d = genpionN_rf_2p.Vect();
		      
		      //Get the unit (normalized) vector along pion momentum
		      PVec genpionP_rf_2p_3d_u = genpionP_rf_2p_3d/TMath::Sqrt(genpionP_rf_2p_3d.mag2());
		      PVec genpionN_rf_2p_3d_u = genpionN_rf_2p_3d/TMath::Sqrt(genpionN_rf_2p_3d.mag2());
		      
		      //Get the longitudinal component of IP vector parallel to pion momenta
		      PVec ipvP_rf_gen_3d_l = ipvP_rf_gen_3d.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
		      PVec ipvN_rf_gen_3d_l = ipvN_rf_gen_3d.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;
		      
		      //if(doDebug) cout<<"LIPV Pion+ ("<<ipvP_rf_gen_3d_l.x()<<","<<ipvP_rf_gen_3d_l.y()<<","<<ipvP_rf_gen_3d_l.z()<<")"<<std::endl;
		      //if(doDebug) cout<<"LIPV Pion- ("<<ipvN_rf_gen_3d_l.x()<<","<<ipvN_rf_gen_3d_l.y()<<","<<ipvN_rf_gen_3d_l.z()<<")"<<std::endl;
		      
		      //Get IP vector normal to pion momenta
		      PVec ipvP_rf_gen_3d_t = ipvP_rf_gen_3d - ipvP_rf_gen_3d_l;
		      PVec ipvN_rf_gen_3d_t = ipvN_rf_gen_3d - ipvN_rf_gen_3d_l;
		      
		      //if(doDebug) cout<<"TIPV Pion+ ("<<ipvP_rf_gen_3d_t.x()<<","<<ipvP_rf_gen_3d_t.y()<<","<<ipvP_rf_gen_3d_t.z()<<")"<<std::endl;
		      //if(doDebug) cout<<"TIPV Pion- ("<<ipvN_rf_gen_3d_t.x()<<","<<ipvN_rf_gen_3d_t.y()<<","<<ipvN_rf_gen_3d_t.z()<<")"<<std::endl;
		      
		      //Get normalized normal IP vector
		      PVec ipvP_rf_gen_3d_t_u = ipvP_rf_gen_3d_t/TMath::Sqrt(ipvP_rf_gen_3d_t.mag2());
		      PVec ipvN_rf_gen_3d_t_u = ipvN_rf_gen_3d_t/TMath::Sqrt(ipvN_rf_gen_3d_t.mag2());
		      
		      //Get CP angle in ZMF
		      double phi_star_gen = TMath::ACos(ipvP_rf_gen_3d_t_u.Dot(ipvN_rf_gen_3d_t_u));
		      CPPhiStar_gen = phi_star_gen;

		      //Get CP angle in lab frame
		      double phi_lab_gen = TMath::ACos(ipvP_lab_gen_3d.Dot(ipvN_lab_gen_3d)/(sqrt(ipvP_lab_gen_3d.mag2())*sqrt(ipvN_lab_gen_3d.mag2())));
		      double OstarCP_gen = genpionN_rf_2p_3d_u.Dot(ipvP_rf_gen_3d_t_u.Cross(ipvN_rf_gen_3d_t_u));
		      double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_gen : (TMath::TwoPi() - phi_star_gen);
		      CPPhiStar_gen_CP = phi_star_gen_cp;
		      
		      //Compute the angle between pi+-pi- momentum in respective tau+/- rf.  check arxiv:1108.0670
		      //It should also be same as the angle bewteen the decay planes. 
		      //ROOT::Math::Boost Boost_RestFrame_TauP(gentauP_lab.BoostToCM());
		      //ROOT::Math::Boost Boost_RestFrame_TauN(gentauN_lab.BoostToCM());
		      
		      //boost to rf of respective taus
		      LVec genpionP_rf_tau = Boost_RestFrame_TauP(genPionP_LabFrame);
		      LVec genpionN_rf_tau = Boost_RestFrame_TauN(genPionN_LabFrame);
		      
		      PVec genpionP_rf_tau_3d = genpionP_rf_tau.Vect();
		      PVec genpionN_rf_tau_3d = genpionN_rf_tau.Vect();
		      
		      PVec genpionP_rf_tau_3d_u = genpionP_rf_tau_3d/sqrt(genpionP_rf_tau_3d.mag2());
		      PVec genpionN_rf_tau_3d_u = genpionN_rf_tau_3d/sqrt(genpionN_rf_tau_3d.mag2());
		      
		      double phi_prime_gen = TMath::ACos(genpionP_rf_tau_3d_u.Dot(genpionN_rf_tau_3d_u));
		      
		      PVec gentauN_rf_3d = genTauN_RestFrame.Vect();
		      PVec gentauN_rf_3d_u = gentauN_rf_3d/sqrt(gentauN_rf_3d.mag2());
		      double psi_prime_gen = TMath::ACos(gentauN_rf_3d_u.Dot(genpionP_rf_tau_3d_u.Cross(genpionN_rf_tau_3d_u)));
		      
		      CPPhiPrime_gen = phi_prime_gen;
		      CPPsiPrime_gen =  psi_prime_gen;

		      //cross check of IP vs Momenta
		      //check if tau vertex direction is parallel to tau momenta
		      PVec tauPDecayVect(genPionP_Vtx.x() - GenHiggs_Vtx.x(), genPionP_Vtx.y() - GenHiggs_Vtx.y(), genPionP_Vtx.z() - GenHiggs_Vtx.z());
		      PVec tauNDecayVect(genPionN_Vtx.x() - GenHiggs_Vtx.x(), genPionN_Vtx.y() - GenHiggs_Vtx.y(), genPionN_Vtx.z() - GenHiggs_Vtx.z());
// 		      histos_["xcheck_angle1"]->Fill(TMath::ACos(tauPDecayVect.Dot(gentauP_lab_3d)/(sqrt(tauPDecayVect.mag2())*sqrt(gentauP_lab_3d.mag2()))));
// 		      histos_["xcheck_angle1"]->Fill(TMath::ACos(tauNDecayVect.Dot(gentauN_lab_3d)/(sqrt(tauNDecayVect.mag2())*sqrt(gentauN_lab_3d.mag2()))));
//uncomment it later		      
// 		      if(doDebug) cout<<"tau+ Dir. Vect ("<<tauPDecayVect.Unit().x()<<", "<<tauPDecayVect.Unit().y()<<", "<<tauPDecayVect.Unit().z()<<")"<<std::endl;
// 		      if(doDebug) cout<<"tau+ Mom. Vect ("<<gentauP_lab_3d.Unit().x()<<", "<<gentauP_lab_3d.Unit().y()<<", "<<gentauP_lab_3d.Unit().z()<<")"<<std::endl;
// 		      if(doDebug) cout<<"tau- Dir. Vect ("<<tauNDecayVect.Unit().x()<<", "<<tauNDecayVect.Unit().y()<<", "<<tauNDecayVect.Unit().z()<<")"<<std::endl;
// 		      if(doDebug) cout<<"tau- Mom. Vect ("<<gentauN_lab_3d.Unit().x()<<", "<<gentauN_lab_3d.Unit().y()<<", "<<gentauN_lab_3d.Unit().z()<<")"<<std::endl;
		      
// 		      //check if IP of pion is parallel to normal componennt of tau vector
// 		      histos_["xcheck_angle2"]->Fill(TMath::ACos(ipvP_lab_gen_3d.Dot(gentauP_lab_3d_n_u)/sqrt(ipvP_lab_gen_3d.mag2())));
// 		      histos_["xcheck_angle2"]->Fill(TMath::ACos(ipvN_lab_gen_3d.Dot(gentauN_lab_3d_n_u)/sqrt(ipvN_lab_gen_3d.mag2())));
		      
// 		      //check if IP of pion is parallel to normal componennt of tau vector, in pi-pi+ RF
// 		      histos_["xcheck_angle3"]->Fill(TMath::ACos(ipvP_rf_gen_3d.Dot(gentauP_rf_2p_n_3d)/(sqrt(ipvP_rf_gen_3d.mag2())*sqrt(gentauP_rf_2p_n_3d.mag2()))));
// 		      histos_["xcheck_angle3"]->Fill(TMath::ACos(ipvN_rf_gen_3d.Dot(gentauN_rf_2p_n_3d)/(sqrt(ipvN_rf_gen_3d.mag2())*sqrt(gentauN_rf_2p_n_3d.mag2()))));
		    } // if(genPionP_LabFrame.pt() > 0 && genPionN_LabFrame.pt() > 0)
		} // TMath::Abs(Gen.PDGId()) == 35

	      count++;
	      if(doDebug)	      cout<<"doFillTree while filling tree: "<< doFillTree << " event " << count << endl;
	      if(doFillTree) tree->Fill();	      
 	    } // for(unsigned int i = 0 ; i < NumAllGenParticles() ; i++)	      
	} // if(ENABLE_MATCHING)
      // end of generator level studies
      
      // ditau selection
      if(!ENABLE_MATCHING || !EARLY_MATCHING || haveEarlyMatch) 
	{

	  if(doDebug) cout<<"NumDiTaus() :"<< NumDiTaus() << endl;
	  if(VERTEX_STUDY) {
	    OStaus.reserve(NumDiTaus());
	    SStaus.reserve(NumDiTaus());
	    
	    for(unsigned int k = 0; k < NumDiTaus(); ++k)
	      {
		if(doDebug)	cout<<"Looping over NumDiTaus "<< k << " of :"<< NumDiTaus() << endl;
		const diTau ditau = diTaus(k);
		const Tau tauP = Taus(ditau.diTau_Leg1Index());
		const Tau tauN = Taus(ditau.diTau_Leg2Index());
		if(doDebug)	cout<<"(leg1, leg2) :"<< ditau.diTau_Leg1Index() <<", "<< ditau.diTau_Leg2Index() << endl;
		OStaus.push_back(TauCand(OSTauWP, tauP, TAU_SCALE));
		SStaus.push_back(TauCand(SSTauWP, tauN, TAU_SCALE));
	      } 
	  } // VERTEX_STUDY
	  
	  //making sure primery vertex selection
	  const int VERTEX_REGIONS = Cutflow_PV(ALL_REGIONS); //15
	  if(doDebug) 	cout<<"VERTEX_REGIONS :"<< VERTEX_REGIONS << endl;
	  
	  if(VERTEX_REGIONS != 0)
	    {
	      //first tau
	      for(unsigned int j = 0; j < NumDiTaus(); ++j)
		{
		  if(doDebug) cout<< "j ------> "<< j << endl;
		  if( ENABLE_MATCHING && j != OSTauIndex ) continue;
		  if( !EARLY_MATCHING ) cutflow.Pass(cutOSTauMatching, VERTEX_REGIONS);
		  if(doDebug)  cout<<"CutOSMatching cut :"<< VERTEX_REGIONS<< endl;
		  
		  const TauCand& OSTau = OStaus[j];
		  if (! OSTau.L1trigger_match ) continue;
		  
		  const int OSTAU_REGIONS = Cutflow_OSTau(OSTau, VERTEX_REGIONS);
		  if(doDebug) 		cout<<"OSTAU_REGIONS :"<< OSTAU_REGIONS <<endl;
		  //	  if( OSTAU_REGIONS == 0) continue;
		  if(!OSTAU_REGIONS) continue;

		  if(doDebug)    cout<<"@@@@   FIRST TAU SELECTED "<< endl;
		  
		  // second tau
		  for(unsigned int k = j+1; k < NumDiTaus(); ++k) if(j != k)
		    {
		      if(doDebug)  cout<<"k ------>" <<  k << endl;
		      const TauCand& SSTau = OStaus[k];
		      if (! SSTau.L1trigger_match ) continue;
		      
		      // Avoid double counting in same sign case
		      // if(SAME_SIGN_TAU_PAIRS && k< j && Cutflow_OSTau(i, SSTau, 0) && Cutflow_SSTau(i, SSTau, OSTau, 0) && Cutflow_Combined(i, SSTau, OSTau, 0)) continue;
		      if(ENABLE_MATCHING && k != SSTauIndex) continue;                    // enable_matching == false
		      if(!EARLY_MATCHING) cutflow.Pass(cutSSTauMatching, OSTAU_REGIONS);  // early_matching  == true
		      
		      const int SSTAU_REGIONS = Cutflow_SSTau(OSTau, SSTau, OSTAU_REGIONS);
		      if(doDebug)    cout<<"SSTAU_REGIONS :"<< SSTAU_REGIONS <<endl;
		      if(!SSTAU_REGIONS) continue;
		      
		      if(doDebug)    cout<<"@@@@   SECOND TAU SELECTED "<< endl;
		      
		      DiTau ditau = {SSTAU_REGIONS, &OSTau, &SSTau};		    
		      doublet.push_back(ditau); 
		      const TauCand *OSTau_ = ditau.OStau;
		      
		      if(doublet.size() < 2 ) continue;
		      if(doDebug) cout<<"          good doublet is selected"<< endl;
		      
		      const DiTau* best_doublet = NULL;
		      //cout<<"size :"<< doublet.size() << endl;
		      for(unsigned int i = 0; i < doublet.size(); ++i)
			{
			  if(!best_doublet ||  doublet[i].OStau->Pt() * doublet[i].SStau->Pt() > best_doublet->OStau->Pt() * best_doublet->SStau->Pt())
			    best_doublet = &doublet[i];
			  // cout<<" doublet[i].OStau->Pt() * doublet[i].SStau->Pt() :"<< (doublet[i].OStau->Pt() * doublet[i].SStau->Pt() )<< endl;
			  // cout<<"best_doublet->OStau->Pt() * best_doublet->SStau->Pt() :"<< (best_doublet->OStau->Pt() * best_doublet->SStau->Pt()) << endl;		
			}
		    
		      if(best_doublet == NULL) continue;
		      
		      if(doDebug) cout<<" BEST_DOUBLET is about to be SELECTED "<< endl;
		      //const int FINAL_REGIONS = Cutflow_FinalSelection(OSTau, SSTau, SSTAU_REGIONS);
		      const int FINAL_REGIONS = SSTAU_REGIONS; //Cutflow_FinalSelection(*best_doublet->OStau, *best_doublet->SStau, SSTAU_REGIONS);
		      if(doDebug)		    cout<<"FINAL_RGIONS "<< FINAL_REGIONS << endl;
		      
		      if(!FINAL_REGIONS) continue;
		      if( (FINAL_REGIONS - 1) != 0 ) continue;
		      
		      // One doublet can only end up in one region
		      assert((FINAL_REGIONS & (FINAL_REGIONS - 1)) == 0);
		      
		      // keep selected ditau events
		      DiTau selected_ditau = {FINAL_REGIONS, &OSTau, &SSTau};
		      finaldoublet.push_back(selected_ditau);  //vector filling the ditau event
		      if(doDebug)    cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ TREE ENTRY IS BEING FILLED"<< endl;
		      
		      //     vSubTauPt_->Fill(selected_ditau.OStau->Pt());//, infos);
		      
		      FillHistogramsFromDoublets(finaldoublet);
		    } //  for(unsigned int k = 0; k < NumTaus(); ++k) if(j != k)
		} // for(unsigned int j = 0; j < NumTaus(); ++j)
	    } // if(VERTEX_REGIONS != 0)
	} // if(!ENABLE_MATCHING || !EARLY_MATCHING || haveEarlyMatch)
    } //  if(OK_Trigger)
  
    //// fill histograms
    // FillHistogramsFromTriplets(triplets    , false);
    // FillHistogramsFromTriplets(triplets_mva, true );
    // FillHistograms(finaldoublet);
    // cutflow.GetCutflow("Signal")->Print();
  
#if 0
  // WZ if(Number() == 711216|| Number() == 1412003)
  //if(Number() == 1121870|| Number() == 1937053|| Number() == 1882145|| Number() == 1982077|| Number() == 148464|| Number() == 316326|| Number() == 1083104)
  //if(Number() == 1692418 || Number() == 1554907 || Number() == 2126397 || Number() == 459336 || Number() == 638573 || Number() == 869633 || Number() == 1155256 || Number() == 1953408 || Number() == 1399060)
  //if(Number() == 5913|| Number() == 34236|| Number() == 197625|| Number() == 51992|| Number() == 37623|| Number() == 89391|| Number() == 51073|| Number() == 35199|| Number() == 2524|| Number() == 49537|| Number() == 183987|| Number() == 110767|| Number() == 180225|| Number() == 9940|| Number() == 32660|| Number() == 36341|| Number() == 128089)
  if( Number() == 73416 || Number() == 75218 || Number() == 75510 || Number() == 76686 || Number() == 76879 || Number() == 76953 || Number() == 80117 || Number() == 80250)
  {
    const std::string trigname = GetTriggerSelection("TauTrigger")->GetTriggerName();

    std::cout << "Event " << Run() << ":" << LumiBlock() << ":" << Number() << ":\n  ";
    cutflow.GetCutflow("Signal")->Print();
    std::cout << "  Trigger=" << trigname << std::endl;
    std::cout << "  MET=" << MET().Et() << ", phi=" << MET().Phi() << std::endl;

    std::cout << "  Triplets: " << std::endl;
    for(unsigned int i = 0; i < triplets.size(); ++i)
      std::cout << "    Triplet " << i << ": muPt=" << Muons(triplets[i].lep).Pt() << ", tauOSPt=" << triplets[i].OStau->Pt() << ", tauSSPt=" << triplets[i].SStau->Pt() << ", region=" << triplets[i].region << ", mvaValue=" << triplets[i].mvaValue << ", mass=" << (triplets[i].OStau->p4 + triplets[i].SStau->p4).M() << std::endl;

    std::cout << "  Muons: " << std::endl;
    for(unsigned int i = 0; i < NumMuons(); ++i)
      std::cout << "    Muon " << i << ": pT=" << Muons(i).Pt() << ", eta=" << Muons(i).Eta() << ", phi=" << Muons(i).Phi() << ", isolation=" << Muons(i).PFIsoR4RelDB() << " (charged=" << Muons(i).PFIsoR4ChargedHad() << ", neutral=" << Muons(i).PFIsoR4NeutralHad() << ", gamma=" << Muons(i).PFIsoR4Photon() << ", PU=" << Muons(i).PFIsoR4SumPUPt() << "), Dxy=" << Muons(i).InnerTrack().Dxy() << ", Dz=" << Muons(i).InnerTrack().Dz() << ", NumStations=" << Muons(i).NumStations() << ", NumChambers=" << Muons(i).NumChamberHits() << ", Chi2OverNdof=" << Muons(i).Chi2OverNdof() << ", NPixelHits=" << Muons(i).InnerTrack().NPixelHits() << std::endl;
    std::cout << "  Electrons: " << std::endl;
    for(unsigned int i = 0; i < NumElectrons(); ++i)
      std::cout << "    Electron " << i << ": pT=" << Electrons(i).Pt() << ", eta=" << Electrons(i).Eta() << ", phi=" << Electrons(i).Phi() << ", charge=" << Electrons(i).Charge() << ", Isolation=" << Electrons(i).PFIsoR4RelDB() << ", Mt=" << calcMt(Electrons(i).Px(), Electrons(i).Py(), MET().Px(), MET().Py()) << ", MVAIdNonTrig=" << Electrons(i).MVAIdNonTrig() << ", Dz=" << Electrons(i).Dz() << ", Dxy=" << Electrons(i).Dxy() << std::endl;
    std::cout << "  Taus: " << std::endl;
    for(unsigned int i = 0; i < NumTaus(); ++i)
      if(Taus(i).TauDiscriminator("decayModeFinding") && Taus(i).NumTracks() > 0 && Taus(i).Pt() > 20.)
        std::cout << "    Tau " << i << ": pT=" << Taus(i).Pt() << ", eta=" << Taus(i).Eta() << ", phi=" << Taus(i).Phi() << ", charge=" << Taus(i).Charge() << ", trigger=" << Taus(i).Trigger(trigname) << ", AntiELoose=" << Taus(i).TauDiscriminator("againstElectronLoose") << ", AntiEMedium=" << Taus(i).TauDiscriminator("againstElectronMedium") << ", AntiETight=" << Taus(i).TauDiscriminator("againstElectronTight") << ", isolation=" << Taus(i).Iso3Hits() << ", Dz=" << Taus(i).LeadingTrack().Dz() << ", Dxy=" << Taus(i).LeadingTrack().Dxy() << ", loose=" << Taus(i).TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr") << ", medium=" << Taus(i).TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr") << ", tight=" << Taus(i).TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr") << ", loose3Hits=" << Taus(i).TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr3Hits") << ", medium3Hits=" << Taus(i).TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr3Hits") << ", tight3Hits=" << Taus(i).TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr3Hits") << ", antiMuTight=" << Taus(i).TauDiscriminator("againstMuonTight") << std::endl;
    std::cout << "  Jets: " << std::endl;
    for(unsigned int i = 0; i < NumAK5PFJets(); ++i)
      if(AK5PFJets(i).Pt() > 20)
        std::cout << "    Jet " << i << ": pT=" << AK5PFJets(i).Pt() << ", eta=" << AK5PFJets(i).Eta() << ", phi=" << AK5PFJets(i).Phi() << ", btag=" << AK5PFJets(i).combinedSecondaryVertexBJetTags() << std::endl;
  }
#endif

  cutflow.FinishEvent();
 
#if 0
  if(false) //Number() == 149833) //false) //Number() == 172415)
    for(unsigned int i = 0; i < NumTaus(); ++i)
      std::cout << "Tau " << i << ": " << Taus(i).Pt() << ", " << Taus(i).Eta() << ", " << Taus(i).Phi() << ", decay mode: " << Taus(i).TauDiscriminator("decayModeFinding") << std::endl;

  if(Number() == 42580 || Number() == 168440)
  {
    std::cout << Number() << ", in file " << GetCurrentFileName() << ":" << std::endl;
    std::cout << "MET: " << MET().Pt() << std::endl;
    std::cout << "RAW PF MET: " << PFMET().Pt() << std::endl;
    std::cout << "TYPE-1 PF MET: " << PFMETTYPE1().Pt() << std::endl;
    for(unsigned int i = 0; i < NumMuons(); ++i)
    {
      std::cout << "Mt " << i << ": " << calcMt(Muons(i).Px(), Muons(i).Py(), MET().Px(), MET().Py()) << std::endl;
      std::cout << "Mt PF " << i << ": " << calcMt(Muons(i).Px(), Muons(i).Py(), PFMET().Px(), PFMET().Py()) << std::endl;
      std::cout << "Mt T1 " << i << ": " << calcMt(Muons(i).Px(), Muons(i).Py(), PFMETTYPE1().Px(), PFMETTYPE1().Py()) << std::endl;
    }
  }

  if(false) //Number() == 140157) //false) //Number() == 7575008 || Number() == 8985493 || Number() == 3218007)
  {
    std::cout << Number() << ", in file " << GetCurrentFileName() << ":" << std::endl;
    for(unsigned int i = 0; i < NumElectrons(); ++i)
    {
      std::cout << "Electron " << i << ": pT=" << Electrons(i).Pt() << ", Eta=" << Electrons(i).Eta() << ", Phi=" << Electrons(i).Phi() << std::endl;
      std::cout << "Charge: " << Electrons(i).Charge() << std::endl;
      std::cout << "ElectronID Loose: " << ElectronId::loose(Electrons(i)) << ", MVA Value=" << Electrons(i).MVAIdNonTrig() << std::endl;
      std::cout << "Isolation: " << Electrons(i).PFIsoR4Rel() << std::endl;
      std::cout << "Isolation (DB-corrected): " << Electrons(i).PFIsoR4RelDB() << std::endl;
      std::cout << "Mt(el,MET): " << calcMt(Electrons(i).Px(), Electrons(i).Py(), MET().Px(), MET().Py()) << std::endl;

      for(unsigned int j = 0; j < NumElectrons(); ++j) if(i!=j)
        std::cout << "\tDi-El mass to " << j << ": " << (Electrons(i)+Electrons(j)).M() << std::endl;
      for(unsigned int j = 0; j < NumTaus(); ++j)
        std::cout << "\tDz to tau " << j << ": " << Taus(j).LeadingTrack().Dz() - Electrons(i).Dz() << std::endl;
    }
  }

  if(false) //Number() == 4080839)
  {
    std::cout << Number() << ", in file " << GetCurrentFileName() << ":" << std::endl;
    for( unsigned int i = 0; i < NumAK5PFJets(); ++i)
    {
      std::cout << "Jet " << i << ": pT=" << AK5PFJets(i).Pt() << ", Eta=" << AK5PFJets(i).Eta() << ", Phi=" << AK5PFJets(i).Phi() << std::endl;
      std::cout << "PUJetFullLoose: " << AK5PFJets(i).puJetFullLoose() << std::endl;
      std::cout << "BTag: " << AK5PFJets(i).combinedSecondaryVertexBJetTags() << std::endl;
    }
  }

  if(false) //Number() == 4218245 || Number() == 1609728) //Number() == 149833) //false)//Number() == 283078600)//163089721 || Number() == 2188659530)
  {
    std::cout << Number() << ", in file " << GetCurrentFileName() << ":" << std::endl;
    std::cout << "  PtMu: " << OK_PtMu << std::endl;
    std::cout << "  EtaMu: " << OK_EtaMu << std::endl;
    std::cout << "  Glotr: " << OK_glotr << std::endl;
    std::cout << "  NStationsMu: " << OK_NStationsMu << std::endl;
    std::cout << "  NChamberMu: " << OK_NChamberMu << std::endl;
    std::cout << "  TrackerLayersMu: " << OK_TrackerLayersMu << std::endl;
    std::cout << "  HitsPixelMu: " << OK_HitsPixelMu << std::endl;
    std::cout << "  Chi2NdfMu: " << OK_Chi2NdfMu << std::endl;
    std::cout << "  DxyMu: " << OK_DxyMu << std::endl;
    std::cout << "  DzMu: " << OK_DzMu << std::endl;
    std::cout << "  IsoMu: " << OK_IsoMu << std::endl;
    std::cout << "  LeadPtTau: " << OK_LeadPtTau << std::endl;
    std::cout << "  LeadEtaTau: " << OK_LeadEtaTau << std::endl;
    std::cout << "  LeadDecayModeDisc: " << OK_LeadDecayModeDisc << std::endl;
    std::cout << "  LeadTightIsoDisc: " << OK_LeadTightIsoDisc << std::endl;
    std::cout << "  LeadTightMuonDisc: " << OK_LeadTightMuonDisc << std::endl;
    std::cout << "  LeadLooseElDisc: " << OK_LeadLooseElDisc << std::endl;
    std::cout << "  LeadDz: " << OK_LeadDz << std::endl;
    std::cout << "  SubPtTau: " << OK_SubPtTau << std::endl;
    std::cout << "  SubEtaTau: " << OK_SubEtaTau << std::endl;
    std::cout << "  SubDecayModeDisc: " << OK_SubDecayModeDisc << std::endl;
    std::cout << "  SubMediumIsoDisc: " << OK_SubMediumIsoDisc << std::endl;
    std::cout << "  SubTightMuonDisc: " << OK_SubTightMuonDisc << std::endl;
    std::cout << "  SubMediumElDisc: " << OK_SubMediumElDisc << std::endl;
    std::cout << "  SubDz: " << OK_SubDz << std::endl;
    std::cout << "  Rleadtaumu: " << OK_Rleadtaumu << std::endl;
    std::cout << "  Rsubtaumu: " << OK_Rsubtaumu << std::endl;
    std::cout << "  Rleadtausubtau: " << OK_Rleadtausubtau << std::endl;
    std::cout << "  OppositeCharge: " << OK_oppCharge << std::endl;
    std::cout << "  MuonVeto: " << OK_MuonVeto << std::endl;
    std::cout << "  ElVeto: " << OK_ElVeto << std::endl;
    std::cout << "  ZVeto: " << OK_ZVeto << std::endl;
    std::cout << "  BJetVeto: " << OK_BJetVeto << std::endl;
    std::cout << "  MET: " << OK_MET << std::endl;
    std::cout << "  MtMuMET: " << OK_MtMuMET << std::endl;
    std::cout << "----------------------------" << OK_MtMuMET << std::endl;
  }
#endif

  return(1);
} // Int_t MyAnalysis::AnalyzeEvent()



int main(int argc, char* argv[])
{
  string namebuf;
  string filename;

  std::auto_ptr<Config> config;
  std::vector<std::string> lumiFiles;
  std::vector<std::string> files;

  for(int i = 1; i < argc; ++i)
  {
    namebuf = argv[i];

    UInt_t slashpos = namebuf.find_last_of("/");
    if(slashpos == namebuf.size())
      filename = namebuf;
    else
      filename = namebuf.substr(slashpos+1);

    if(filename.find("LUMI_") == 0)
      lumiFiles.push_back(namebuf);
    else if(filename.find(".cfg") != std::string::npos)
      if(config.get() == NULL)
        config.reset(new Config(namebuf.c_str()));
      else
        config->merge(Config(namebuf.c_str()));
    else
      files.push_back(namebuf);
  }
  
  if(!config.get())
    throw std::runtime_error("No configuration file provided!");

  FillInfo::HIGHPT_TAU1_CUT = 45.;
  FillInfo::HIGHPT_TAU2_CUT = 30.;
  FillInfo::HIGHPT_MET_CUT = 15.;
  
  TH1::SetDefaultSumw2(true);
  MyAnalysis ana(*config);
  
  for(std::vector<std::string>::const_iterator iter = lumiFiles.begin(); iter != lumiFiles.end(); ++iter)
    ana.AddLumiFile(*iter);
  for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
    ana.AddFile(*iter);
  
  if(filename.find("Data")== 0  || filename.find("TTFH_Data")== 0 ||filename.find("LUMI_INFO_TTFH_Data")== 0 ){
    ana.SetSample("Data");
    cout << "Data" << endl;
  }
  else if(filename.find("QCD")== 0 || filename.find("MC_QCD")== 0 || filename.find("TTFH_QCD")== 0 || filename.find("LUMI_INFO_TTFH_QCD")== 0 || filename.find("LUMI_INFO_TTFH_MC_QCD")== 0) {
    ana.SetSample("QCD");
    cout << "QCD" << endl;
  }
  else if(filename.find("MC_TTJets")== 0 || filename.find("TTFH_MC_TTJets")== 0 || filename.find("LUMI_INFO_TTFH_MC_TTJets")== 0 ) {
    ana.SetSample("TTJets");
    cout << "TTJets" << endl;
  }
  else if(filename.find("MC_DYJets")== 0 || filename.find("TTFH_MC_DYJets")== 0 || filename.find("LUMI_INFO_TTFH_MC_DYJets")== 0 ) {
    ana.SetSample("DYJets");
    cout << "DYJets" << endl;
  }
  else if(filename.find("MC_WJets")== 0 || filename.find("TTFH_MC_WJets")== 0 || filename.find("LUMI_INFO_TTFH_MC_WJets")== 0  || filename.find("MC_W1Jet")== 0 || filename.find("TTFH_MC_W1Jet")== 0 || filename.find("LUMI_INFO_TTFH_MC_W1Jet" || filename.find("MC_W2Jets")== 0 || filename.find("TTFH_MC_W2Jets")== 0 || filename.find("LUMI_INFO_TTFH_MC_W2Jets") || filename.find("MC_W3Jets")== 0 || filename.find("TTFH_MC_W3Jets")== 0 || filename.find("LUMI_INFO_TTFH_MC_W3Jets") || filename.find("MC_W4Jets")== 0 || filename.find("TTFH_MC_W4Jets")== 0 || filename.find("LUMI_INFO_TTFH_MC_W4Jets"))) {
    ana.SetSample("WJets");
    cout << "WJets" << endl;
  }
  else if(filename.find("MC_WW")== 0 || filename.find("TTFH_MC_WW")== 0 || filename.find("LUMI_INFO_TTFH_MC_WW")== 0 ) {
    ana.SetSample("WW");
    cout << "WW" << endl;
  }
  else if(filename.find("MC_WZ")== 0 || filename.find("TTFH_MC_WZ")== 0 || filename.find("LUMI_INFO_TTFH_MC_WZ")== 0 ) {
    ana.SetSample("WZ");
    cout << "WZ" << endl;
  }
  else if(filename.find("MC_ZZ")== 0 || filename.find("TTFH_MC_ZZ")== 0 || filename.find("LUMI_INFO_TTFH_MC_ZZ")== 0 ) {
    ana.SetSample("ZZ");
    cout << "ZZ" << endl;
  }
  else if(filename.find("MC_WH")== 0 || filename.find("TTFH_MC_WH")== 0 || filename.find("LUMI_INFO_TTFH_MC_WH")== 0 ) {
    ana.SetSample("WH");
    cout << "WH" << endl;
  }
  else{
    cout << "failed to set sample:" << endl;
    return(1);
  }

  // ana.SetSample("WH");
  // Loop will start to run the analysis on the specified range or on 
  // all events if no range is given. 
  ana.SetPrintInfo(10000);
  ana.EnableDuplicateCheck();
  ana.PrintLumiOfRuns();

  if(doDebug)  cout<<"trying to trigger the event"<< endl;
  vector<string> vimtrigger;
  vimtrigger.push_back("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v.*");
  vimtrigger.push_back("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v.*");
  vimtrigger.push_back("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v.*");
  vimtrigger.push_back("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v.*");
  vimtrigger.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_v*");
  ana.AddTriggerSelection("TauTrigger", vimtrigger);
  ana.GetTriggerSelection("TauTrigger")->PrintInfo();
  ana.Loop();
} //int main(int argc, char* argv[])


void MyAnalysis::SetSample( std::string flag)
{
  sample = flag.c_str();
}

int MyAnalysis::Cutflow_PV(int regions)
{
#if 0
  double cPVNdof = -1; /* TEMP! */
  double cPVz = 100;
  double cPVr = 100;
#endif
  
  const double rPV = sqrt(pow(PrimVertex().x(),2)+pow(PrimVertex().y(),2));
  
  if( !( PrimVertex().Ndof()               > cPVNdof)) return 0; cutflow.Pass(cutVtxNDof, regions);
  if(doDebug) cout<<"PV.ndof :"<< regions<<endl;
  if( !( TMath::Abs(PrimVertex().z())      < cPVz))    return 0; cutflow.Pass(cutVtxZ, regions);
  if(doDebug) cout<<"PV.z :"<< regions<<endl;
  if( !( TMath::Abs(rPV)                   < cPVr))    return 0; cutflow.Pass(cutVtxR, regions);
  if(doDebug) cout<<"PV.R :"<< regions<<endl;
  // PV
  hPVNdof->Fill(PrimVertex().Ndof(), weight);
  hPVz->Fill(PrimVertex().z(), weight);
  hPVr->Fill(rPV, weight);
  if(doDebug) cout<<"PV return region:"<< regions<<endl;
  return regions;
}

int MyAnalysis::Cutflow_Muon(int mu, int regions)
{
//   if( !( Muons(mu).Pt()                                            > 24))                                return 0; cutflow.Pass(cutMuPt, regions);
//   if( !( TMath::Abs(Muons(mu).Eta())                               < 2.1))                               return 0; cutflow.Pass(cutMuEta, regions);
//   if( !( (Muons(mu).IsGlobal() && Muons(mu).IsPF())))                                                    return 0; cutflow.Pass(cutMuID, regions);
//   if( !( Muons(mu).NumStations()                                   >= 2))                                return 0; cutflow.Pass(cutMuNumStations, regions);
//   if( !( Muons(mu).NumChamberHits()                                >= 1))                                return 0; cutflow.Pass(cutMuNumChambers, regions);
//   if( !( Muons(mu).InnerTrack().NPixelLayers() + Muons(mu).InnerTrack().NStripLayers() > 5))             return 0; cutflow.Pass(cutMuNumLayers, regions);
//   if( !( Muons(mu).InnerTrack().NPixelHits()                       >= 1))                                return 0; cutflow.Pass(cutMuNumPixel, regions);
//   if( !( Muons(mu).Chi2OverNdof()                                  < 10.0))                              return 0; cutflow.Pass(cutMuChi2, regions);
//   if( !( TMath::Abs(Muons(mu).InnerTrack().Dxy())                  < 0.045))                              return 0; cutflow.Pass(cutMuDxy, regions);
//   if( !( TMath::Abs(Muons(mu).InnerTrack().Dz())                   < 0.2))                               return 0; cutflow.Pass(cutMuDz, regions);
//   if( !( Muons(mu).PFIsoR4RelDB()                                  < 0.1))                               return 0; cutflow.Pass(cutMuIso, regions);
  return regions;
}

bool MyAnalysis::hasEle(const TauCand& tau, const Electron& ele, int regions)
{
  bool hasele = false;
  hasele = ( ele.Pt() > 10 &&
	     TMath::Abs(ele.Eta()) < 2.5 &&
	     ElectronId::loose(ele) &&
	     ( ROOT::Math::VectorUtil::DeltaR(tau.p4, ele) < 0.3 ) );
  return hasele;
}

bool MyAnalysis::hasMuon(const TauCand& tau, const Muon& muon, int regions)
{
  bool hasMu = false;
  hasMu = ( muon.Pt() > 10 &&
	    TMath::Abs(muon.Eta()) < 2.4 &&
	    MuonId::loose(muon) &&
	    ( ROOT::Math::VectorUtil::DeltaR(tau.p4, muon) < 0.3 ) );
  return hasMu;
}


int MyAnalysis::Cutflow_OSTau(const TauCand& OStau, int regions)
{
  if(doDebug)  cout<<"inside CutFlow_OSTAU"<< endl;

  if( !( OStau.Pt()                                > SUBTAU_PT_CUT))   return 0; cutflow.Pass(cutOSTauPt              , regions);

  if(doDebug)    cout<<" OStau pt pass - check! " << OStau.Pt() <<", " << regions << endl;


  if( !( TMath::Abs(OStau.Eta())                   < 2.1          ))   return 0; cutflow.Pass(cutOSTauEta             , regions);
  if(doDebug)    cout<<" OStau Eta pass - check! " << OStau.Eta() <<", " << regions << endl;

  if( !( SAME_SIGN_TAU_PAIRS || OStau.Charge       < 0            ))   return 0; cutflow.Pass(cutOSTauCharge          , regions);
  if(doDebug)    cout<<" OStau Charge pass - check! " << OStau.Charge <<", "<< regions << endl;

  if( !(OStau.leadpfchargedhadrcandpt              > 0            ))   return 0; cutflow.Pass( cutOSTauChHadCandPt    , regions);

  if( !( OStau.PFChargedHadrCands_size ==1 && OStau.PFGammaCands_size == 0 ))  return 0; cutflow.Pass(cutOSTauDecayModeFinding, regions);
  //if( !( OStau.decayModeFinding                                   ))   return 0; cutflow.Pass(cutOSTauDecayModeFinding, regions);
  if(doDebug)    cout<<" OStau decaymodeFinding pass - check! " << regions << endl;

  if( CAPPED_FAKERATES && OStau.Iso3Hits            > 10.          )   return 0;

  if( !( OStau.*OSTauWP.ISOLATION_FLAG                            ))   return 0; cutflow.Pass(cutOSTauIsolation        , regions);
  if(doDebug)    cout<<" OStau isolation ("<< OSTauWP.ISOLATION_FLAG <<") pass - check! "<< regions << endl;

  if( !( OStau.*OSTauWP.ISOLATION_FLAG                            ))   regions &= ~(def.SIGNAL  | def.BACKGROUND         );   
  else                                                                 regions &= ~(def.CONTROL | def.CONTROL_BACKGROUND );
  if(doDebug)  cout<<"OStau --> (sig, bkg, cont, cont_bkg) : ("<< def.SIGNAL <<", "<< def.BACKGROUND<< ", "<< def.CONTROL <<", "<< def.CONTROL_BACKGROUND <<")"<< endl;

  if(!(TMath::Abs(OStau.Dz)                         < 0.2         ))   return 0; cutflow.Pass(cutOSTauDz              , regions);
  if(doDebug)    cout<<" OStau Dz pass - check! " << regions << endl;

  if(!(TMath::Abs(OStau.Dxy)                        < 0.2         ))   return 0; cutflow.Pass(cutOSTauDxy             , regions);
  if(doDebug)    cout<<" OStau Dxy pass - check! " << regions << endl;
  if( !( OStau.*OSTauWP.ANTIMU_FLAG                                                ))   return 0; cutflow.Pass(cutOSTauAntiMu, regions);
  if( !( OStau.*OSTauWP.ANTIE_FLAG                                                 ))   return 0; cutflow.Pass(cutOSTauAntiE, regions);
  
  cutflow.Pass(cutOSTauIsolation, regions);
  return regions;
}


int MyAnalysis::Cutflow_SSTau(const TauCand& OStau, const TauCand& SStau, int regions)
{
  if(doDebug)  cout<<"inside CutFlow_SSTAU"<< endl;
  const double PT_CUT = (OStau.Pt() > LEADTAU_PT_CUT ? SUBTAU_PT_CUT : LEADTAU_PT_CUT);
  if( !( SStau.Pt()                                > PT_CUT       ))   return 0; cutflow.Pass(cutSSTauPt              , regions);
  if(doDebug)    cout<<" SStau pt pass - check! " << regions << endl;
  
  if( !( TMath::Abs(SStau.Eta())                   < 2.1          ))   return 0; cutflow.Pass(cutSSTauEta             , regions);
  if(doDebug)    cout<<" SStau Eta pass - check! " << regions << endl;

  if( !( (!SAME_SIGN_TAU_PAIRS && SStau.Charge     > 0)     ||
         (SAME_SIGN_TAU_PAIRS && SStau.Charge * OStau.Charge > 0) ))   return 0; cutflow.Pass(cutSSTauCharge          , regions);
  if(doDebug)    cout<<" SStau Charge pass - check! " << regions << endl;

  if( !(SStau.leadpfchargedhadrcandpt              > 0            ))   return 0; cutflow.Pass( cutSSTauChHadCandPt    , regions);

  if( !( SStau.PFChargedHadrCands_size ==1 && SStau.PFGammaCands_size == 0 ))  return 0; cutflow.Pass(cutSSTauDecayModeFinding, regions);
  //if( !( SStau.decayModeFinding                                   ))   return 0; cutflow.Pass(cutSSTauDecayModeFinding, regions);
  if(doDebug) cout<<" SSTau decayModeFinding pass - check! " << regions << endl;

  if( CAPPED_FAKERATES && SStau.Iso3Hits            > 10.          )   return 0;
  if( !( SStau.*OSTauWP.ISOLATION_FLAG                            ))   return 0; cutflow.Pass(cutSSTauIsolation        , regions);
  if(doDebug)    cout<<" SStau isolation passs - check!" << regions << endl;

  if(!(SStau.*SSTauWP.ISOLATION_FLAG))                                 regions &= ~(def.SIGNAL     | def.CONTROL               );
  else                                                                 regions &= ~(def.BACKGROUND | def.CONTROL_BACKGROUND    );
  if(doDebug)  cout<<"SStau --> (sig, bkg, cont, cont_bkg) : ("<< def.SIGNAL <<", "<< def.BACKGROUND<< ", "<< def.CONTROL <<", "<< def.CONTROL_BACKGROUND <<")"<< endl;

  if(!(TMath::Abs(SStau.Dz)                         < 0.2         ))   return 0; cutflow.Pass(cutSSTauDz              , regions);
  if(doDebug)    cout<<" SStau Dz pass - check! " << regions << endl;

  if(!(TMath::Abs(SStau.Dxy)                        < 0.2         ))   return 0; cutflow.Pass(cutSSTauDxy             , regions);
  if(doDebug)    cout<<" SStau Dxy pass - check! " << regions << endl;

  if( !( SStau.*OSTauWP.ANTIMU_FLAG                                                ))   return 0; cutflow.Pass(cutSSTauAntiMu, regions);
  if( !( SStau.*OSTauWP.ANTIE_FLAG                                                 ))   return 0; cutflow.Pass(cutSSTauAntiE, regions);
  cutflow.Pass(cutSSTauIsolation, regions);
  return regions;
}


#if 0
int MyAnalysis::Cutflow_recoTauP(const TauCand& tauP, int regions)
{
  if(doDebug)  cout<<"inside CutFlow_recoTauP"<< endl;
  if( !( tauP.Pt()                                                                 > 20 ))  return 0; cutflow.Pass(cutRecoTauPpt, regions);
  if(doDebug) cout<<" tauP_pt :"<< regions << endl;

  if( !( TMath::Abs(tauP.Eta())                                                    < 2.1))   return 0; cutflow.Pass(cutRecoTauPeta, regions);
  if(doDebug) cout<<" tauP_eta :"<< regions << endl;
  
  if( tauP.Charge < 0)                                                                      return 0; cutflow.Pass(cutRecoTauPcharge, regions);
  
  if( !( tauP.byLooseCombinedIsolationDeltaBetaCorr3Hits                                ))   return 0; cutflow.Pass(cutRecoTauPdecayModeFinding, regions);
  if(doDebug) cout<<" tauP_decay :"<< regions << endl;
  
  if( CAPPED_FAKERATES && tauP.leadpfchargedhadrcandpt < 0.)                                return 0;
  
  if(doDebug)  cout<<"sig :"<< def.SIGNAL <<" bkg :"<< def.BACKGROUND<< " cont : "<< def.CONTROL <<" cont_bkg :"<< def.CONTROL_BACKGROUND << endl;
  
//   if( !( tauP.*TauPWP.ISOLATION_FLAG  ))     regions &= ~(def.SIGNAL | def.BACKGROUND);   
//   else                                        regions &= ~(def.CONTROL | def.CONTROL_BACKGROUND);
//   if(doDebug)  cout<<"OS_iso :"<< regions << endl;
  
  //cutflow.Pass(cutRecoTauPisolation, regions);
  if(doDebug)  cout<<"tauP return region :"<< regions << endl;

  return regions;
}


int MyAnalysis::Cutflow_recoTauN(const TauCand& tauN, int regions)
{
  if(doDebug)  cout<<"inside CutFlow_recoTauN"<< endl;
  if( !( tauN.Pt()                                                                 > 20 ))  return 0; cutflow.Pass(cutRecoTauNpt, regions);
  if(doDebug) cout<<" tauN_pt :"<< regions << endl;

  if( !( TMath::Abs(tauN.Eta())                                                    < 2.1))   return 0; cutflow.Pass(cutRecoTauNeta, regions);
  if(doDebug) cout<<" tauN_eta :"<< regions << endl;
  
  if( tauN.Charge > 0)                                                                      return 0; cutflow.Pass(cutRecoTauNcharge, regions);
  
  if( !( tauN.byLooseCombinedIsolationDeltaBetaCorr3Hits                                ))   return 0; cutflow.Pass(cutRecoTauNdecayModeFinding, regions);
  if(doDebug) cout<<" tauN_decay :"<< regions << endl;
  
  if( CAPPED_FAKERATES && tauN.leadpfchargedhadrcandpt < 0.)                                return 0;
    if(doDebug)  cout<<"sig :"<< def.SIGNAL <<" bkg :"<< def.BACKGROUND<< " cont : "<< def.CONTROL <<" cont_bkg :"<< def.CONTROL_BACKGROUND << endl;

//   if( !( tauN.*TauPWP.ISOLATION_FLAG  ))      regions &= ~(def.SIGNAL | def.BACKGROUND);   
//   else                                        regions &= ~(def.CONTROL | def.CONTROL_BACKGROUND);
//   if(doDebug)  cout<<"tauN_iso :"<< regions << endl;
  
    // cutflow.Pass(cutRecoTauNisolation, regions);
  if(doDebug)  cout<<"tauN return region :"<< regions << endl;

  return regions;
}

#endif



int MyAnalysis::Cutflow_SSTau(int mu, const TauCand& OStau, const TauCand& SStau, int regions)
{
  const double PT_CUT = (OStau.Pt() > LEADTAU_PT_CUT ? SUBTAU_PT_CUT : LEADTAU_PT_CUT);
  if( !( SStau.Pt()                                                            > PT_CUT)) return 0; cutflow.Pass(cutSSTauPt, regions);
  if( !( TMath::Abs(SStau.Eta())                                               < 2.3)) return 0; cutflow.Pass(cutSSTauEta, regions);
  if( !( (!SAME_SIGN_TAU_PAIRS && Muons(mu).Charge() * SStau.Charge > 0) ||
         (SAME_SIGN_TAU_PAIRS && SStau.Charge * OStau.Charge > 0)                   )) return 0; cutflow.Pass(cutSSTauCharge, regions);
  if( !( SStau.decayModeFinding                                                    ))  return 0; cutflow.Pass(cutSSTauDecayModeFinding, regions);
  if(CAPPED_FAKERATES && SStau.Iso3Hits > 10.) return 0;

  if(!(SStau.*SSTauWP.ISOLATION_FLAG))
    regions &= ~(def.SIGNAL | def.CONTROL); //bitwise OR (|) 
  else
    regions &= ~(def.BACKGROUND | def.CONTROL_BACKGROUND);
  cutflow.Pass(cutSSTauIsolation, regions);

  if( !( SStau.*SSTauWP.ANTIMU_FLAG                                               )) return 0; cutflow.Pass(cutSSTauAntiMu, regions);
  if( !( SStau.*SSTauWP.ANTIE_FLAG                                                )) return 0; cutflow.Pass(cutSSTauAntiE, regions);
  //  if( !(TMath::Abs(SStau.Dz)                                                 < 0.2))  return 0; cutflow.Pass(cutSSTauDz, regions);
  return regions;
}


int MyAnalysis::Cutflow_OSTau(int mu, const TauCand& OStau, int regions)
{
  if( !( OStau.Pt()                                                             > SUBTAU_PT_CUT)) return 0; cutflow.Pass(cutOSTauPt, regions);
  if( !( TMath::Abs(OStau.Eta())                                                < 2.3))   return 0; cutflow.Pass(cutOSTauEta, regions);
  //if( !( SAME_SIGN_TAU_PAIRS || Muons(mu).Charge() * OStau.Charge < 0))                   return 0; cutflow.Pass(cutOSTauCharge, regions);
  //  if( !( SAME_SIGN_TAU_PAIRS ))                                       return 0; cutflow.Pass(cutOSTauCharge, regions);
  if( !( OStau.decayModeFinding                                                      ))   return 0; cutflow.Pass(cutOSTauDecayModeFinding, regions);
  //if(CAPPED_FAKERATES && OStau.Iso3Hits > 10.) return 0;
  if( !( OStau.*OSTauWP.ISOLATION_FLAG                                             ))
    regions &= ~(def.SIGNAL | def.BACKGROUND);
  else
    regions &= ~(def.CONTROL | def.CONTROL_BACKGROUND);
  cutflow.Pass(cutOSTauIsolation, regions);
  
  if( !( OStau.*OSTauWP.ANTIMU_FLAG                                                ))   return 0; cutflow.Pass(cutOSTauAntiMu, regions);
  if( !( OStau.*OSTauWP.ANTIE_FLAG                                                 ))   return 0; cutflow.Pass(cutOSTauAntiE, regions);
  // if(!(TMath::Abs(OStau.Dz)                                                     < 0.2))   return 0; cutflow.Pass(cutOSTauDz, regions);
  return regions;
}


int MyAnalysis::Cutflow_FinalSelection(const TauCand& OStau, const TauCand& SStau, int regions)
{ 

  //////////////////////////////////////////////////////////////////////  
  // Muon Veto
  bool OK_MuonVetoEvent=false;
  for( unsigned int i2 = 0; i2 < NumMuons(); ++i2){
    if(Muons(i2).Pt()                       > 10  && 
       TMath::Abs(Muons(i2).Eta())          < 2.4 &&
       MuonId::loose(Muons(i2))                   &&
       ( (ROOT::Math::VectorUtil::DeltaR(OStau.p4, Muons(i2))  < 0.3) || (ROOT::Math::VectorUtil::DeltaR(SStau.p4, Muons(i2))  < 0.3) )
       ) OK_MuonVetoEvent=true;
  }
  if(OK_MuonVetoEvent==true) return 0;
  cutflow.Pass(cutMuonVeto, regions);
  

  //////////////////////////////////////////////////////////////////////
  // Electron Veto
  bool OK_ElVetoEvent=false;
  for( unsigned int j2 = 0; j2 < NumElectrons(); ++j2){
    if(Electrons(j2).Pt()                       > 10  &&
       TMath::Abs(Electrons(j2).Eta())          < 2.5 &&
       ElectronId::loose(Electrons(j2))               &&
       ( (ROOT::Math::VectorUtil::DeltaR(SStau.p4, Electrons(j2)) < 0.3) || (ROOT::Math::VectorUtil::DeltaR(SStau.p4, Electrons(j2))  < 0.3) )
       ) OK_ElVetoEvent=true;
  }
  if(OK_ElVetoEvent==true) return 0;
  cutflow.Pass(cutElVeto, regions);
  
  if(doDebug)  cout<<"regions return by Cutflow_finalSelection :"<< regions<< endl;
  return regions;
}


int MyAnalysis::Cutflow_FinalSelection(int mu, const TauCand& OStau, const TauCand& SStau, int regions)
{ 
   if( !( ROOT::Math::VectorUtil::DeltaR(OStau.p4, Muons(mu))                   > 0.5)) return 0; cutflow.Pass(cutDrOSTauMu, regions);
   if( !( ROOT::Math::VectorUtil::DeltaR(SStau.p4, Muons(mu))                   > 0.5)) return 0; cutflow.Pass(cutDrSSTauMu, regions);
   // if( !( ROOT::Math::VectorUtil::DeltaR(OStau.p4, SStau.p4)                    > 0.5)) return 0; cutflow.Pass(cutDrOSTauSSTau, regions);

  //////////////////////////////////////////////////////////////////////  
  // Muon Veto
  bool OK_MuonVetoEvent=false;
  for( unsigned int i2 = 0; i2 < NumMuons(); ++i2){
    if(i2==mu) continue;
    if(Muons(i2).Pt()                       > 10  && 
       TMath::Abs(Muons(i2).Eta())          < 2.4 &&
       MuonId::loose(Muons(i2))) OK_MuonVetoEvent=true;
  }
  if(OK_MuonVetoEvent==true) return 0;
  cutflow.Pass(cutMuonVeto, regions);

  //////////////////////////////////////////////////////////////////////
  // Electron Veto
  bool OK_ElVetoEvent=false;
  for( unsigned int j2 = 0; j2 < NumElectrons(); ++j2){
    if(Electrons(j2).Pt()                       > 10 &&
       TMath::Abs(Electrons(j2).Eta())          < 2.5 &&
       ElectronId::loose(Electrons(j2))) OK_ElVetoEvent=true;
  }
  if(OK_ElVetoEvent==true) return 0;
  cutflow.Pass(cutElVeto, regions);
  
  return regions;
}



int MyAnalysis::Cutflow_Combined(int mu, const TauCand& OStau, const TauCand& SStau, int regions)
{ 
  if( !( ROOT::Math::VectorUtil::DeltaR(OStau.p4, Muons(mu))                   > 0.5)) return 0; cutflow.Pass(cutDrOSTauMu, regions);
  if( !( ROOT::Math::VectorUtil::DeltaR(SStau.p4, Muons(mu))                   > 0.5)) return 0; cutflow.Pass(cutDrSSTauMu, regions);
  if( !( ROOT::Math::VectorUtil::DeltaR(OStau.p4, SStau.p4)                    > 0.5)) return 0; cutflow.Pass(cutDrOSTauSSTau, regions);

  //////////////////////////////////////////////////////////////////////  
  // Muon Veto
  bool OK_MuonVetoEvent=false;
  for( unsigned int i2 = 0; i2 < NumMuons(); ++i2){
    if(i2==mu)continue;
    if(Muons(i2).Pt()                       > 10 && 
       TMath::Abs(Muons(i2).Eta())          < 2.1 &&
       Muons(i2).PFIsoR4RelDB()             < 0.3 &&
       Muons(i2).IsGlobal() && Muons(i2).IsPF() &&
       TMath::Abs(Muons(i2).InnerTrack().Dz()) < 0.2) OK_MuonVetoEvent=true;
  }
  if(OK_MuonVetoEvent==true) return 0;
  cutflow.Pass(cutMuonVeto, regions);

  //////////////////////////////////////////////////////////////////////
  // Electron Veto
  bool OK_ElVetoEvent=false;
  for( unsigned int j2 = 0; j2 < NumElectrons(); ++j2){
    if(Electrons(j2).Pt()                      >10 &&
       TMath::Abs(Electrons(j2).Eta())         < 2.5 &&
       Electrons(j2).PFIsoR4RelDB()             < 0.3 &&
       ElectronId::loose(Electrons(j2)) && 
       TMath::Abs(Electrons(j2).Dz()) < 0.2) OK_ElVetoEvent=true;
  }
  if(OK_ElVetoEvent==true) return 0;
  cutflow.Pass(cutElVeto, regions);

  //////////////////////////////////////////////////////////////////////
  // Z Veto
  cutflow.Pass(cutZVeto, regions);

  //////////////////////////////////////////////////////////////////////
  // BTag Veto
  bool OK_BJetVetoEvent=false;
  for( unsigned int k2 = 0; k2 < NumAK5PFJets(); ++k2){
    if(AK5PFJets(k2).Pt()                      >20 && 
       TMath::Abs(AK5PFJets(k2).Eta())         <2.4 &&
       AK5PFJets(k2).puJetFullLoose()               &&
       AK5PFJets(k2).combinedSecondaryVertexBJetTags() > 0.898) OK_BJetVetoEvent=true;
  }
  if(OK_BJetVetoEvent==true) return 0;
  cutflow.Pass(cutBJetVeto, regions);

  //////////////////////////////////////////////////////////////////////
  // MET Cut
  if(!(MET().Et() > MET_CUT)) return 0;
  cutflow.Pass(cutMET, regions);

  //////////////////////////////////////////////////////////////////////
  // Mt Cut
  if(!(calcMt(Muons(mu).Px(), Muons(mu).Py(), MET().Px(), MET().Py()) > 30)) return 0;
  cutflow.Pass(cutMt, regions);

  return regions;
}

void MyAnalysis::FillHistogramsFromDoublets(const std::vector<DiTau>& doublets)
{
  if(doDebug)  cout<<"Inside filling loop for ehceking doublets" << endl;
  const DiTau* best_doublet = NULL;
  if(doDebug)  cout<<"doublets.size() :"<< doublets.size() << endl;

  for(unsigned int i = 0; i < doublets.size(); ++i) 
    {
      best_doublet = &doublets[i];
      if(doDebug) cout<<" inside doublet_pt2 :"<<  (best_doublet->OStau->Pt() * best_doublet->SStau->Pt()) <<" ttpt2 :"<< (doublets[i].OStau->Pt() * doublets[i].SStau->Pt() ) << endl; 
    }
  if(doDebug)  cout<<" outside doublet_pt2 :"<<  (best_doublet->OStau->Pt() * best_doublet->SStau->Pt()) <<" ttpt2 :"<< (doublets[0].OStau->Pt() * doublets[0].SStau->Pt() ) << endl;

  if(ONLY_BEST_DOUBLET && best_doublet) FillHistograms(*best_doublet, false);//, doublets.size());
}


//void MyAnalysis::FillHistograms(const triplet& trip, unsigned int n_triplets, bool mvaCut)
void MyAnalysis::FillHistograms(const DiTau& ditau, bool mvaCut)
{
  if(doDebug) cout<<"inside the FillHistograms(const DiTau& ditau, bool mvaCut) "<< endl;
  //  const Muon mu = Muons(trip.lep);
  const TauCand& OSTau        = *ditau.OStau;
  const TauCand& SSTau        = *ditau.SStau;
 
  const TLorentzVector tautau = OSTau.p4 + SSTau.p4;
  const double Lt             = OSTau.Pt() + SSTau.Pt();
  double SVfitMass = -1.0;


//   for(unsigned int dt = 0; dt < NumMuTauTauPairs(); ++dt)
//   {
//     const MuTauTauPair taupair = MuTauTauPairs(dt);
//     if(ROOT::Math::VectorUtil::DeltaR(taupair.Muon(), mu) > 0.001) continue;
//     if(ROOT::Math::VectorUtil::DeltaR(taupair.Leg1(), OSTau.p4) > 0.001 && ROOT::Math::VectorUtil::DeltaR(taupair.Leg2(), OSTau.p4) > 0.001) continue;
//     if(ROOT::Math::VectorUtil::DeltaR(taupair.Leg1(), SSTau.p4) > 0.001 && ROOT::Math::VectorUtil::DeltaR(taupair.Leg2(), SSTau.p4) > 0.001) continue;

//     SVfitMass = taupair.SVFitMassInt();
//     break;
//   }

  if(mvaCut) cutflow.Pass(cutCombinatorics, ditau.region);

//  if(SVfitMass < 0.0)
//    std::cerr << "WARNING!!! Event with no SVfit mass, run=" << Run() << ", lumi=" << LumiBlock() << ", event=" << Number() << ": " << NumMuTauTauPairs() << " di tau pair candidates. Fakes: " << mu_fake << ", " << leadtau_fake << ", " << subtau_fake << std::endl;

  if(mvaCut)
  {
    if(ditau.region & def.SIGNAL)
      outfile << "SIGNAL Run_" << Run() << "_LS_" << LumiBlock() << "_Event_" << Number() << std::endl;
    if(ditau.region & def.BACKGROUND)
      outfile << "BACKGROUND Run_" << Run() << "_LS_" << LumiBlock() << "_Event_" << Number() << std::endl;
    if(ditau.region & def.CONTROL)
      outfile << "CONTROL Run_" << Run() << "_LS_" << LumiBlock() << "_Event_" << Number() << std::endl;
    if(ditau.region & def.CONTROL_BACKGROUND)
      outfile << "CTRLBKG Run_" << Run() << "_LS_" << LumiBlock() << "_Event_" << Number() << std::endl;
  }

  FillInfo infoNormal(def, OSTau, SSTau);
  //   const FillInfo infoInclusive(def, fakeRateFuncsInclusive, mvaEvaluatorBDT8_v2, irredMVA, NULL, trip.mvaVars, mvaCut, mu, OSTau, SSTau);
  //   const FillInfo infoMETWeighted(def, fakeRateFuncsNormal, mvaEvaluatorBDT8_v2, irredMVA, &metWeight, trip.mvaVars, mvaCut, mu, OSTau, SSTau);
  
  std::vector<const FillInfo*> infos;
  infos.push_back(&infoNormal);
  // infos.push_back(&infoInclusive);
  // infos.push_back(&infoMETWeighted);


//   // Log the event.
//   if(!trip.lepJetFake && !trip.OSTauJetFake && !trip.SSTauJetFake && (trip.wjetsVetoMT || !trip.OSTauTrackFake || !trip.SSTauTrackFake))
//   if(mvaVars.mvaVars.leadPt > 45.0 && mvaVars.mvaVars.subPt > 30.0 && mvaVars.mvaVars.Met > 20.0)
//   outfile << "Run_" << Run() << "_LS_" << LumiBlock() << "_Event_" << Number() << std::endl;
  
  if(doDebug)  cout<<"filling the Event /variables"<< endl;
  // Event
  if(doDebug)  cout<<"first round"<< endl;
  vNPV.fill(NumGoodPrimVertices(), weight, infos);
  vMET.fill(MET().Et(), weight, infos);
  vMETTYPE1.fill(PFMETTYPE1().Et(), weight, infos);
  //  vNTriplets.fill(n_triplets, weight, infos);
  // Tau
  const double LeadTauPt = OSTau.Pt() > SSTau.Pt() ? OSTau.Pt() : SSTau.Pt();
  const double SubTauPt = OSTau.Pt() < SSTau.Pt() ? OSTau.Pt() : SSTau.Pt();
  const double LeadTauEta = OSTau.Pt() > SSTau.Pt() ? OSTau.Eta() : SSTau.Eta();
  const double SubTauEta = OSTau.Pt() < SSTau.Pt() ? OSTau.Eta() : SSTau.Eta();

  if(doDebug)  cout<<"second round"<< endl;
  vLeadTauPt.fill(LeadTauPt, weight, infos);
  vLeadTauEta.fill(LeadTauEta, weight, infos);
  vSubTauPt.fill(SubTauPt, weight, infos);
  vSubTauEta.fill(SubTauEta, weight, infos);
  vOSTauPt.fill(OSTau.Pt(), weight, infos);
  vSSTauPt.fill(SSTau.Pt(), weight, infos);
  vOSTauEta.fill(OSTau.Eta(), weight, infos);
  vSSTauEta.fill(SSTau.Eta(), weight, infos);
  vOSTauIso.fill(OSTau.IsoNeutralsPt + OSTau.IsoChargedPt + OSTau.IsoGammaPt, weight, infos);
  vSSTauIso.fill(SSTau.IsoNeutralsPt + SSTau.IsoChargedPt + SSTau.IsoGammaPt, weight, infos);
  vOSTauRelIso.fill( (OSTau.IsoNeutralsPt + OSTau.IsoChargedPt + OSTau.IsoGammaPt) / OSTau.Pt(), weight, infos);
  vSSTauRelIso.fill( (SSTau.IsoNeutralsPt + SSTau.IsoChargedPt + SSTau.IsoGammaPt) / SSTau.Pt(), weight, infos);

  // Tau Tau
  if(doDebug)  cout<<"third round"<< endl;
  vTauTauVisMass.fill(tautau.M(), weight, infos);
  vTauTauPt.fill(tautau.Pt(), weight, infos);
  vTauTauEta.fill(tautau.Eta(), weight, infos);
//   vTauTauDPhi.fill(ditau.mvaVars.mvaVars.deltaPhiDiTau, weight, infos);
//   vTauTauDEta.fill(ditau.mvaVars.mvaVars.deltaEtaDiTau, weight, infos);
//   vTauTauDR.fill(ditau.mvaVars.mvaVars.deltaRDiTau, weight, infos);
//   vTauTauPtRatio.fill(ditau.mvaVars.mvaVars.ptRatio, weight, infos);
//   vTauTauAngle.fill(ditau.mvaVars.mvaVars.angle, weight, infos);
  vTauTauSVfitMass.fill(SVfitMass, weight, infos);
  // Tau Tau Muon
  //vBDT8_v1.fill(mvaEvaluatorBDT8_v1.mvaValue(mvaVars.mvaVars), weight, infos);
  //const double mvaValue = mvaEvaluatorBDT8_v2.mvaValue(trip.mvaVars.mvaVars);
  //vBDT8_v2.fill(mvaValue, weight, infos);
  //vBDT8_v3.fill(mvaEvaluatorBDT8_v3.mvaValue(mvaVars.mvaVars), weight, infos);
  //vBDTG_v2.fill(mvaEvaluatorBDTG_v2.mvaValue(mvaVars.mvaVars), weight, infos);
  //vBDT8Mt_v2.fill(mvaEvaluatorBDT8Mt_v2.mvaValue(mvaVars.mvaVars), weight, infos);
  //vBDTGMt_v2.fill(mvaEvaluatorBDTGMt_v2.mvaValue(mvaVars.mvaVars), weight, infos);
  //vIrredMVA.fill(infoNormal.irredMVAValue, weight, infos);
  //vCombinedMVA.fill(combinedMVA.mvaValue(trip.mvaValue, infoNormal.irredMVAValue), weight, infos);
  // Profiles
  //  pDiTauPtOSTauMuMass.fill(Variable<TProfile>::Filler(tautau.Pt(), (OSTau.p4 + mu).M()), weight, infos);
  pLeadTauPtSubTauPt.fill(Variable<TProfile>::Filler(LeadTauPt, SubTauPt), weight, infos);
  pLeadTauPtMET.fill(Variable<TProfile>::Filler(LeadTauPt, MET().Pt()), weight, infos);
  pSubTauPtMET.fill(Variable<TProfile>::Filler(SubTauPt, MET().Pt()), weight, infos);
//   pLeadTauPtDeltaR.fill(Variable<TProfile>::Filler(LeadTauPt, ditau.mvaVars.mvaVars.deltaRDiTau), weight, infos);
//   pLeadTauPtPtRatio.fill(Variable<TProfile>::Filler(LeadTauPt, ditau.mvaVars.mvaVars.ptRatio), weight, infos);
//   pSubTauPtDeltaR.fill(Variable<TProfile>::Filler(SubTauPt, ditau.mvaVars.mvaVars.deltaRDiTau), weight, infos);
//   pSubTauPtPtRatio.fill(Variable<TProfile>::Filler(SubTauPt, ditau.mvaVars.mvaVars.ptRatio), weight, infos);
//   pMETDeltaR.fill(Variable<TProfile>::Filler(MET().Pt(), ditau.mvaVars.mvaVars.deltaRDiTau), weight, infos);
//   pMETPtRatio.fill(Variable<TProfile>::Filler(MET().Pt(), ditau.mvaVars.mvaVars.ptRatio), weight, infos);
  //pDeltaRPtRatio.fill(Variable<TProfile>::Filler(trip.mvaVars.mvaVars.deltaRDiTau, trip.mvaVars.mvaVars.ptRatio), weight, infos);

  // MVA Tree
  // if(!mvaCut) mvaTree.fill(trip.mvaVars, mvaValue, weight);
  
  if(!mvaCut && (ditau.region & def.SIGNAL))
  {
    if(doDebug)  cout<<"inside !mcaCut loop fill"<< endl; 
    run = Run();
    lumi = LumiBlock();
    event = Number();

//     muonPt = mu.Pt();
//     muonEta = mu.Eta();
//     muonPhi = mu.Phi();
//     muonIso = mu.PFIsoR4RelDB();
//     muonCharge = mu.Charge();
//     muonGlobal = mu.IsGlobal();
//     muonTracker = mu.IsTracker();
//     muonPF = mu.IsPF();
//     muonPixelHits = mu.InnerTrack().NPixelHits();
//     muonNumStations = mu.NumStations();
//     muonNumChamberHits = mu.NumChamberHits();
//     muonPixelLayers = mu.InnerTrack().NPixelLayers();
//     muonStripLayers = mu.InnerTrack().NStripLayers();
//     muonChi2OverNdof = mu.Chi2OverNdof();

    tauOSPt = OSTau.Pt();
    tauOSEta = OSTau.Eta();
    tauOSPhi = OSTau.Phi();
    tauOSIso3Hits = OSTau.Iso3Hits;
    tauOSIso3HitsLoose = OSTau.byLooseCombinedIsolationDeltaBetaCorr3Hits;
    tauOSIso3HitsMedium = OSTau.byMediumCombinedIsolationDeltaBetaCorr3Hits;
    tauOSIso3HitsTight = OSTau.byTightCombinedIsolationDeltaBetaCorr3Hits;
    tauOSAgainstElectronLoose = OSTau.againstElectronLoose;
    tauOSAgainstElectronLooseMVA3 = OSTau.againstElectronLooseMVA3;
    tauOSAgainstElectronMediumMVA3 = OSTau.againstElectronMediumMVA3;
    tauOSAgainstElectronTightMVA3 = OSTau.againstElectronTightMVA3;
    tauOSAgainstMuonLoose = OSTau.againstMuonLoose;
    tauOSAgainstMuonMedium = OSTau.againstMuonMedium;
    tauOSAgainstMuonTight = OSTau.againstMuonTight;

    tauSSPt = OSTau.Pt();
    tauSSEta = OSTau.Eta();
    tauSSPhi = OSTau.Phi();
    tauSSIso3Hits = OSTau.Iso3Hits;
    tauSSIso3HitsLoose = OSTau.byLooseCombinedIsolationDeltaBetaCorr3Hits;
    tauSSIso3HitsMedium = OSTau.byMediumCombinedIsolationDeltaBetaCorr3Hits;
    tauSSIso3HitsTight = OSTau.byTightCombinedIsolationDeltaBetaCorr3Hits;
    tauSSAgainstElectronLoose = OSTau.againstElectronLoose;
    tauSSAgainstElectronLooseMVA3 = OSTau.againstElectronLooseMVA3;
    tauSSAgainstElectronMediumMVA3 = OSTau.againstElectronMediumMVA3;
    tauSSAgainstElectronTightMVA3 = OSTau.againstElectronTightMVA3;
    tauSSAgainstMuonLoose = OSTau.againstMuonLoose;
    tauSSAgainstMuonMedium = OSTau.againstMuonMedium;
    tauSSAgainstMuonTight = OSTau.againstMuonTight;

    met = MET().Pt();
    metPhi = MET().Phi();
    //  mt = calcMt(MET().Px(), MET().Py());//, mu.Px(), mu.Py());
    if(doDebug)    cout<<"suncTree is being filled"<< endl;
    syncTree->Fill();

  }
  if(doDebug) cout<<"exiting the Loop  for filling the ditau information synTree"<< endl;
}


void MyAnalysis::FillHistograms(const DiTau& ditau)
 {
   if(doDebug) cout<<"Finally filling the entry----"<< endl;
  const TauCand& OSTau = *ditau.OStau;
  const TauCand& SSTau = *ditau.SStau;

  const TLorentzVector tautau = OSTau.p4   + SSTau.p4;
  const double Lt             = OSTau.Pt() + SSTau.Pt();
  const double LeadTauPt      = OSTau.Pt() > SSTau.Pt() ? OSTau.Pt()  : SSTau.Pt();
  const double SubTauPt       = OSTau.Pt() < SSTau.Pt() ? OSTau.Pt()  : SSTau.Pt();
  const double LeadTauEta     = OSTau.Pt() > SSTau.Pt() ? OSTau.Eta() : SSTau.Eta();
  const double SubTauEta      = OSTau.Pt() < SSTau.Pt() ? OSTau.Eta() : SSTau.Eta();
  
  nPV = NumGoodPrimVertices();
  met = MET().Et(); 
  tau1_Pt = LeadTauPt;
  tau1_Eta = LeadTauEta;
  tau1_Phi = OSTau.Phi();
  tau1_Iso = (OSTau.IsoNeutralsPt + OSTau.IsoChargedPt + OSTau.IsoGammaPt);
  tau1_RelIso = ((OSTau.IsoNeutralsPt + OSTau.IsoChargedPt + OSTau.IsoGammaPt) / OSTau.Pt());
  tau2_Pt = SubTauPt;
  tau2_Eta = SubTauEta;
  tau2_Phi = SSTau.Phi();
  tau2_Iso = (SSTau.IsoNeutralsPt + SSTau.IsoChargedPt + SSTau.IsoGammaPt);
  tau2_RelIso = ((SSTau.IsoNeutralsPt + SSTau.IsoChargedPt + SSTau.IsoGammaPt) / SSTau.Pt());
  tautau_Mass = tautau.M();
  tautau_Pt = tautau.Pt();
  tautau_Eta = tautau.Eta();
  //   tautau_Phi = ditau.deltaPhiDiTau;
  //   tautau_dEta = ditau.deltaEtaDiTau;
  //   tautau_dR =ditau.deltaRDiTau;
  PVndof = PrimVertex().Ndof(); 
  PVz = PrimVertex().z();  
  PVy = PrimVertex().y(); 
  PVx = PrimVertex().x();
  PVr = sqrt(pow(PrimVertex().x(),2)+pow(PrimVertex().y(),2));
  
  tree->Fill();
} // FillHistograms(const triplet& trip, unsigned int n_triplets)

#if 0
void MyAnalysis::FillHistogramsFromTriplets(const std::vector<triplet>& triplets, bool mvaCut)
{
  const triplet* best_triplet = NULL;
  for(unsigned int i = 0; i < triplets.size(); ++i)
    {
      if(!best_triplet || 
	 triplets[i].OStau->Pt() * triplets[i].SStau->Pt() > best_triplet->OStau->Pt() * best_triplet->SStau->Pt()
	 )
	best_triplet = &triplets[i];
      
      if(!ONLY_BEST_TRIPLET)
	FillHistograms(triplets[i], triplets.size(), mvaCut);
    }
  
  if(ONLY_BEST_TRIPLET && best_triplet)
    FillHistograms(*best_triplet, triplets.size(), mvaCut);
}
#endif


TLorentzVector MyAnalysis::MET() const
{
  return PFMETMVA() * MET_SCALE;
}


std::map<std::string, TH1F*> MyAnalysis::CreateHistograms() {
  std::map<std::string, TH1F*> histos_;

  histos_["DeltaPVX"] = new TH1F("DeltaPVX", "", 100, -0.1, 0.1);  
  histos_["DeltaPVY"] = new TH1F("DeltaPVY", "", 100, -0.1, 0.1);
  histos_["DeltaPVZ"] = new TH1F("DeltaPVZ", "", 100, -0.1, 0.1);

  histos_["DeltaPVwrtGenX"] = new TH1F("DeltaPVwrtGenX", "", 100, -0.1, 0.1);
  histos_["DeltaPVwrtGenY"] = new TH1F("DeltaPVwrtGenY", "", 100, -0.1, 0.1);
  histos_["DeltaPVwrtGenZ"] = new TH1F("DeltaPVwrtGenZ", "", 100, -0.1, 0.1);

  histos_["DeltaRfPVwrtGenX"] = new TH1F("DeltaRfPVwrtGenX", "", 100, -0.1, 0.1);
  histos_["DeltaRfPVwrtGenY"] = new TH1F("DeltaRfPVwrtGenY", "", 100, -0.1, 0.1);
  histos_["DeltaRfPVwrtGenZ"] = new TH1F("DeltaRfPVwrtGenZ", "", 100, -0.1, 0.1);

  histos_["ResPVwrtGenX"] = new TH1F("ResPVwrtGenX", "", 100, -10.0, 10.0);
  histos_["ResPVwrtGenY"] = new TH1F("ResPVwrtGenY", "", 100, -10.0, 10.0);
  histos_["ResPVwrtGenZ"] = new TH1F("ResPVwrtGenZ", "", 100, -10.0, 10.0);

  histos_["ResRfPVwrtGenX"] = new TH1F("ResRfPVwrtGenX", "", 100, -10.0, 10.0);
  histos_["ResRfPVwrtGenY"] = new TH1F("ResRfPVwrtGenY", "", 100, -10.0, 10.0);
  histos_["ResRfPVwrtGenZ"] = new TH1F("ResRfPVwrtGenZ", "", 100, -10.0, 10.0);

  histos_["VtxXRes"] = new TH1F("VtxXRes",  "", 100, -0.1, 0.1);
  histos_["VtxYRes"] = new TH1F("VtxYRes",  "", 100, -0.1, 0.1);
  histos_["VtxZRes"] = new TH1F("VtxZRes",  "", 100, -0.1, 0.1);

  histos_["ReFitVtxXRes"] = new TH1F("ReFitVtxXRes",  "", 100, -0.1, 0.1);
  histos_["ReFitVtxYRes"] = new TH1F("ReFitVtxYRes",  "", 100, -0.1, 0.1);
  histos_["ReFitVtxZRes"] = new TH1F("ReFitVtxZRes",  "", 100, -0.1, 0.1);
  
  histos_["VtxXErr"] = new TH1F("VtxXErr",  "", 100, -0.05, 0.05);
  histos_["VtxYErr"] = new TH1F("VtxYErr",  "", 100, -0.05, 0.05);
  histos_["VtxZErr"] = new TH1F("VtxZErr",  "", 100, -0.05, 0.05);

  histos_["ReFitVtxXErr"] = new TH1F("ReFitVtxXErr",  "", 100, -0.05, 0.05);
  histos_["ReFitVtxYErr"] = new TH1F("ReFitVtxYErr",  "", 100, -0.05, 0.05);
  histos_["ReFitVtxZErr"] = new TH1F("ReFitVtxZErr",  "", 100, -0.05, 0.05);

  histos_["CPPhiStar"] = new TH1F("CPPhiStar",  "", 20, 0.0, 4.0);
  histos_["CPPhiLab"] = new TH1F("CPPhiLab",  "", 20, 0.0, 4.0);

  histos_["CPPhiStar_opv"] = new TH1F("CPPhiStar_opv",  "", 20, 0.0, 4.0);
  histos_["CPPhiStar_gen"] = new TH1F("CPPhiStar_gen",  "", 20, 0.0, 4.0);
  histos_["CPPhiLab_gen"] = new TH1F("CPPhiLab_gen",  "", 20, 0.0, 4.0);
  histos_["CPPhiCP_gen"] = new TH1F("CPPhiCP_gen",  "", 20, 0.0, 4.0);
  
  return histos_;
}
  
										      
