#ifndef _TTCPANALYSIS_VARIABLE_H_
#define _TTCPANALYSIS_VARIABLE_H_
//
#include "Common/FakeRateFunctions.h"
#include "Common/TauCand.h"
#include "Common/MvaEvaluator.h"
#include "Common/IrredMVA.h"
#include "Common/METWeight.h"
#include "Common/Cutflow.h"
#include "AnalysisTool/Analyse.h"
#include <TLorentzVector.h>
#include <TH1D.h>

struct FillInfo
{
  struct Lepton: TLorentzVector
  {
    Lepton(const Muon&     muon): TLorentzVector(muon), Charge(muon.Charge()) {}
    Lepton(const Electron& elec): TLorentzVector(elec), Charge(elec.Charge()) {}
    const int Charge;
  };
  
  struct RegionDef
  {
    RegionDef(CutflowRegions & regions);
    const int SIGNAL;
    const int BACKGROUND;
    const int CONTROL;
    const int CONTROL_BACKGROUND;
  };
  
  FillInfo(const RegionDef& regionDef,  const TauCand& ostau,  const TauCand& sstau);
  const RegionDef& def;
  const TauCand& ostau;
  const TauCand& sstau;


/*   const double wosw; */
/*   const double wssw; */
/*   const double wosz; */
/*   const double wssz; */

/*   const double wosw_up; */
/*   const double wosw_down; */
/*   const double wssw_up; */
/*   const double wssw_down; */
/*   const double wosz_up; */
/*   const double wosz_down; */
/*   const double wssz_up; */
/*   const double wssz_down; */

  static double HIGHPT_TAU1_CUT;
  static double HIGHPT_TAU2_CUT;
  static double HIGHPT_MET_CUT;
  static double MASS_CUT;
}; //FillInfo


//type-1
template<typename HistoType>
class Variable
{
 public:
  static const int SIGNAL                      = 1 << 0; //1
  static const int FAKE_SS_WZ                  = 1 << 2; //4
  static const int FAKE_SS_ERR                 = 1 << 3; //8
  static const int FAKE_SS_UNWEIGHTED          = 1 << 4; //16
  static const int FAKE_12_WZ                  = 1 << 7; //128
  static const int FAKE_12_ERR                 = 1 << 8; //256
  static const int FAKE_12_UNWEIGHTED          = 1 << 9; //512
  static const int MVA_OPTIMIZATION            = 1 << 16; //65536
  
  Variable(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);
  
  struct Filler {
    Filler(double x, double y = 0);
    void operator()(HistoType* h, double w) const;
    double x, y;
  };
  
  void fill(const Filler& filler, double weight, const FillInfo& info) const;

private:
  HistoType* HistoSignal;
  HistoType* HistoSignalNoPU;
  HistoType* HistoTauSSFakeWJets;
  HistoType* HistoTauSSFakeZJets;
  HistoType* HistoTauSSFakeErrUpWJets;
  HistoType* HistoTauSSFakeErrDownWJets;
  HistoType* HistoTauSSFakeErrUpZJets;
  HistoType* HistoTauSSFakeErrDownZJets;
  HistoType* HistoTauSSFakeUnweighted;
  HistoType* HistoTauSSFakeUnweightedNoPU;
  HistoType* HistoTauOSFakeWJets;
  HistoType* HistoTauOSFakeZJets;
  HistoType* HistoTauOSFakeErrUpWJets;
  HistoType* HistoTauOSFakeErrDownWJets;
  HistoType* HistoTauOSFakeErrUpZJets;
  HistoType* HistoTauOSFakeErrDownZJets;
  HistoType* HistoTauOSFakeUnweighted;
  HistoType* HistoTauOSFakeUnweightedNoPU;
  HistoType* HistoTau12FakeWJets;
  HistoType* HistoTau12FakeZJets;
  HistoType* HistoTau12FakeErrUpWJets;
  HistoType* HistoTau12FakeErrDownWJets;
  HistoType* HistoTau12FakeErrUpZJets;
  HistoType* HistoTau12FakeErrDownZJets;
  HistoType* HistoTau12FakeUnweighted;
  HistoType* HistoTau12FakeUnweightedNoPU;
  HistoType* HistoTau12FakeTauOSUnweightedWJets;
  HistoType* HistoTau12FakeTauOSUnweightedErrUpWJets;
  HistoType* HistoTau12FakeTauOSUnweightedErrDownWJets;
  HistoType* HistoTau12FakeTauSSUnweightedWJets;
  HistoType* HistoTau12FakeTauOSUnweightedZJets;
  HistoType* HistoTau12FakeTauSSUnweightedZJets;
};

//type-2
template<typename HistoType>
class VariableByTrackIso
{
 public:
  VariableByTrackIso(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);
  
  void fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& infos) const;
 private:
  std::auto_ptr<Variable<HistoType> > normal;
  std::auto_ptr<Variable<HistoType> > inclusive;
  std::auto_ptr<Variable<HistoType> > metWeighted;
  std::auto_ptr<Variable<HistoType> > mu;
  std::auto_ptr<Variable<HistoType> > muInclusive;
};


//type-3
template<typename HistoType>
class VariableByCategory
{
 public:
  VariableByCategory(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);
  
  void fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& infos) const;
 private:
  
  std::auto_ptr<VariableByTrackIso<HistoType> > inclusive;
  
  std::auto_ptr<VariableByTrackIso<HistoType> > low_mass;
  std::auto_ptr<VariableByTrackIso<HistoType> > high_mass;
  
  std::auto_ptr<VariableByTrackIso<HistoType> > low_irred;
  std::auto_ptr<VariableByTrackIso<HistoType> > high_irred;
  
  std::auto_ptr<VariableByTrackIso<HistoType> > mass3_low;
  std::auto_ptr<VariableByTrackIso<HistoType> > mass3_medium;
  std::auto_ptr<VariableByTrackIso<HistoType> > mass3_high;
};


//type-4
template<typename HistoType>
class VariableBySelection
{
public:
  VariableBySelection(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);

  void fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& infos) const;
  //  void fill(const typename Variable<HistoType>::Filler& filler, double weight) const;

private:
  std::auto_ptr<VariableByCategory<HistoType> > uncut;
  std::auto_ptr<VariableByCategory<HistoType> > cut;
  std::auto_ptr<VariableByCategory<HistoType> > mva;
  std::auto_ptr<VariableByCategory<HistoType> > mva_control;
};

#endif
