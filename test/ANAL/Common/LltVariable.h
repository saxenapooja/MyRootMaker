
#ifndef _WHANALYSIS_LLTVARIABLE_H_
#define _WHANALYSIS_LLTVARIABLE_H_

#include "Common/TauCand.h"
#include "Common/Cutflow.h"
#include "Common/kNNEvaluator.h"
#include "Common/MuonId.h"
#include "Common/FakeRateFunctions.h"

#include "AnalysisTool/Analyse.h"

#include <TLorentzVector.h>
#include <TH1D.h>

struct LltFillInfo
{
  struct Lepton: TLorentzVector
  {
    Lepton(const Muon& muon, bool genMatched, const Analyse& ana): TLorentzVector(muon), Charge(muon.Charge()), Isolation(muon.PFIsoR4DB()), JetPt(MuonId::jetPt(ana, muon)), GenMatched(genMatched) {}
    Lepton(const Electron& elec, bool genMatched, const Analyse& ana): TLorentzVector(elec), Charge(elec.Charge()), Isolation(elec.PFIsoR4DB()), JetPt(-1.), GenMatched(genMatched) {}

    const int Charge;
    const double Isolation;
    const double JetPt;
    const bool GenMatched;
  };

  struct RegionDef
  {
    RegionDef(CutflowRegions& regions, const char* data_era);

    const int SIGNAL;
    const int SIGNAL_LEADMU_FAKE;
    const int SIGNAL_SUBMU_FAKE;
    const int SIGNAL_BOTHMU_FAKE;

    const int F3;
    const int F3_LEADMU_FAKE;
    const int F3_SUBMU_FAKE;
    const int F3_BOTHMU_FAKE;

    double wzLumiScale;
    double zzLumiScale;
    double wjetsDataLead;
    double wjetsDataSub;
    double wjetsWZLead;
    double wjetsWZSub;
    double wjetsZZLead;
    double wjetsZZSub;
  };

  LltFillInfo(const RegionDef& regionDef, const FakeRateFunctions& fr1Jet, const FakeRateFunctions& fr2Jet, const Basic_kNNEvaluator& kNNData, const Basic_kNNEvaluator* kNNWZ, const Basic_kNNEvaluator* kNNZZ, const Lepton& leadLep, const Lepton& subLep, const TauCand& tau, double met, unsigned int nJets20, int regions);

  const RegionDef& def;

  const Lepton leadlep;
  const Lepton sublep;
  const TauCand& tau;
  const int regions;

  const double lt;
  const double met;
  unsigned int njets20;

  const double wfr1jet;
  const double wfr2jet;

  const double wkNNleadData;
  const double wkNNsubData;

  const double wkNNleadWZ;
  const double wkNNsubWZ;

  const double wkNNleadZZ;
  const double wkNNsubZZ;

  const double wkNNlead;
  const double wkNNsub;
};

template<typename HistoType>
class LltVariable
{
public:
  static const int SIGNAL = 1 << 0;
  static const int F3 = 1 << 1;

  LltVariable(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);

  struct Filler {
    Filler(double x, double y = 0);
    void operator()(HistoType* h, double w) const;

    double x, y;
  };

  void fill(const Filler& filler, double weight, const LltFillInfo& info) const;

private:
  HistoType* HistoSignal;
  HistoType* HistoSignalNoPU;

  HistoType* HistoSignalLeadMuFakeUnweighted;
  HistoType* HistoSignalLeadMuFakekNNWeighted;
  HistoType* HistoSignalLeadMuFakekNNWeightedData;

  HistoType* HistoSignalSubMuFakeUnweighted;
  HistoType* HistoSignalSubMuFakekNNWeighted;
  HistoType* HistoSignalSubMuFakekNNWeightedData;

  HistoType* HistoSignalBothMuFakeUnweighted;
  HistoType* HistoSignalBothMuFakekNNWeighted;
  HistoType* HistoSignalBothMuFakekNNWeightedData;

  HistoType* HistoF3;
  HistoType* HistoF3NoPU;
  HistoType* HistoF3W1JetWeighted;
  HistoType* HistoF3W2JetWeighted;

  HistoType* HistoF3LeadMuFakeUnweighted;
  HistoType* HistoF3LeadMuFakekNNWeighted;
  HistoType* HistoF3LeadMuFakekNNWeightedData;

  HistoType* HistoF3SubMuFakeUnweighted;
  HistoType* HistoF3SubMuFakekNNWeighted;
  HistoType* HistoF3SubMuFakekNNWeightedData;

  HistoType* HistoF3BothMuFakeUnweighted;
  HistoType* HistoF3BothMuFakekNNWeighted;
  HistoType* HistoF3BothMuFakekNNWeightedData;
};

template<typename HistoType>
class LltVariableByTrackIso
{
public:
  LltVariableByTrackIso(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);

  void fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& infos) const;
private:
  std::auto_ptr<LltVariable<HistoType> > frOrderedkNN;
  std::auto_ptr<LltVariable<HistoType> > frUnorderedkNN;
  std::auto_ptr<LltVariable<HistoType> > frDY;
  std::auto_ptr<LltVariable<HistoType> > frTT;
};

template<typename HistoType>
class LltVariableByCategory
{
public:
  LltVariableByCategory(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);

  void fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& infos) const;
private:

  std::auto_ptr<LltVariableByTrackIso<HistoType> > inclusive;

  std::auto_ptr<LltVariableByTrackIso<HistoType> > low_lt;
  std::auto_ptr<LltVariableByTrackIso<HistoType> > high_lt;

  std::auto_ptr<LltVariableByTrackIso<HistoType> > lt80;
};

template<typename HistoType>
class LltVariableBySelection
{
public:
  LltVariableBySelection(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax);

  void fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& infos) const;

private:
  std::auto_ptr<LltVariableByCategory<HistoType> > uncut;
};

#endif // _WHANALYSIS_LLTVARIABLE_H_
