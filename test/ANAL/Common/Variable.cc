#include "Common/Variable.h"
#include "MariaFunctions.h"
#include <TProfile.h>

namespace
{
  double w(double f) { 
    assert(f >= 0.0); 
    if(f >= 1.0) return 0.0; 
    assert(f < 1.0); 
    return f/(1.0 - f); 
  }
}

double FillInfo::MASS_CUT = -1.0;
double FillInfo::HIGHPT_TAU1_CUT = -1.0;
double FillInfo::HIGHPT_TAU2_CUT = -1.0;
double FillInfo::HIGHPT_MET_CUT = -1.0;

FillInfo::RegionDef::RegionDef(CutflowRegions& regions):
  SIGNAL(regions.AddRegion("Signal")),
  BACKGROUND(regions.AddRegion("Background")),
  CONTROL(regions.AddRegion("Control")),
  CONTROL_BACKGROUND(regions.AddRegion("ControlBackground"))
{
}

FillInfo::FillInfo(const RegionDef& regionDef, const FakeRateFunctions& fr, const MvaEvaluator& mvaEvaluator, const IrredMVA& irredMVA, const METWeight* metWeight, const MvaTree::Variables& vars, bool mvaCut, const Lepton& lepton, const TauCand& ostau, const TauCand& sstau):

  def(regionDef),

  lepton(lepton), ostau(ostau), sstau(sstau),
  regions(vars.regions), case2(vars.case2),

  met(vars.mvaVars.Met),
  mvaCut(mvaCut),
  irredMVAValue(irredMVA.mvaValue(lepton.Pt(), ostau.Pt(), sstau.Pt(), vars.mvaVars.Met)),

  ssMassWeight(metWeight == NULL ? 1.0 : (ostau.p4 + sstau.p4).M() > 60 && (ostau.p4 + sstau.p4).M() < 80 ? 1.20 : 1.),
  ssFakeMetWeight(1.0), //metWeight == NULL ? 1.0 : metWeight->getWeight(met)),

  wosw(fr.TauOS->getPtEta(FakeRateFunctions::WJETS_FUNC, ostau.Pt(), ostau.Eta(), case2)),
  wssw(fr.TauSS->getPtEta(FakeRateFunctions::WJETS_FUNC, sstau.Pt(), sstau.Eta(), case2)),
  wosz(fr.TauOS->getPtEta(FakeRateFunctions::ZJETS_FUNC, ostau.Pt(), ostau.Eta(), case2)),
  wssz(fr.TauSS->getPtEta(FakeRateFunctions::ZJETS_FUNC, sstau.Pt(), sstau.Eta(), case2)),

  wosw_up(fr.TauOS->getPtEta(FakeRateFunctions::WJETS_ERR_UP, ostau.Pt(), ostau.Eta(), case2)),
  wosw_down(fr.TauOS->getPtEta(FakeRateFunctions::WJETS_ERR_DOWN, ostau.Pt(), ostau.Eta(), case2)),
  wssw_up(fr.TauSS->getPtEta(FakeRateFunctions::WJETS_ERR_UP, sstau.Pt(), sstau.Eta(), case2)),
  wssw_down(fr.TauSS->getPtEta(FakeRateFunctions::WJETS_ERR_DOWN, sstau.Pt(), sstau.Eta(), case2)),
  wosz_up(fr.TauOS->getPtEta(FakeRateFunctions::ZJETS_ERR_UP, ostau.Pt(), ostau.Eta(), case2)),
  wosz_down(fr.TauOS->getPtEta(FakeRateFunctions::ZJETS_ERR_DOWN, ostau.Pt(), ostau.Eta(), case2)),
  wssz_up(fr.TauSS->getPtEta(FakeRateFunctions::ZJETS_ERR_UP, sstau.Pt(), sstau.Eta(), case2)),
  wssz_down(fr.TauSS->getPtEta(FakeRateFunctions::ZJETS_ERR_DOWN, sstau.Pt(), sstau.Eta(), case2))
{
}


// FillInfo::FillInfo(const RegionDef& regionDef, const TauCand& ostau, const TauCand& sstau):
//   def(regionDef),
//   ostau(ostau), 
//   sstau(sstau),
//   lepton(lepton), 
//   regions(NULL), case2(NULL),

//   met(NULL),
//   mvaCut(NULL),
//   irredMVAValue(NULL),

//   ssMassWeight(0),
//   ssFakeMetWeight(0),
//   wosw(NULL),
//   wssw(NULL),
//   wosz(NULL),
//   wssz(NULL),

//   wosw_up(NULL),
//   wosw_down(NULL),
//   wssw_up(NULL),
//   wssw_down(NULL),
//   wosz_up(NULL),
//   wosz_down(NULL),
//   wssz_up(NULL),
//   wssz_down(NULL)
// {
// }


template<typename HistoType>
Variable<HistoType>::Filler::Filler(double x, double y):
  x(x), y(y)
{
}

template<>
void Variable<TH1D>::Filler::operator()(TH1D* histo, double weight) const
{
  histo->Fill(x, weight);
}

template<>
void Variable<TProfile>::Filler::operator()(TProfile* prof, double weight) const
{
  prof->Fill(x, y, weight);
}

template<typename HistoType>
Variable<HistoType>::Variable(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax):
  HistoSignal(NULL), HistoSignalNoPU(NULL),
  HistoTauSSFakeWJets(NULL), HistoTauSSFakeZJets(NULL),
  HistoTauSSFakeErrUpWJets(NULL), HistoTauSSFakeErrDownWJets(NULL), HistoTauSSFakeErrUpZJets(NULL), HistoTauSSFakeErrDownZJets(NULL),
  HistoTauSSFakeUnweighted(NULL), HistoTauSSFakeUnweightedNoPU(NULL),
  HistoTauOSFakeWJets(NULL), HistoTauOSFakeZJets(NULL),
  HistoTauOSFakeErrUpWJets(NULL), HistoTauOSFakeErrDownWJets(NULL), HistoTauOSFakeErrUpZJets(NULL), HistoTauOSFakeErrDownZJets(NULL),
  HistoTauOSFakeUnweighted(NULL), HistoTauOSFakeUnweightedNoPU(NULL),
  HistoTau12FakeWJets(NULL), HistoTau12FakeZJets(NULL),
  HistoTau12FakeErrUpWJets(NULL), HistoTau12FakeErrDownWJets(NULL), HistoTau12FakeErrUpZJets(NULL), HistoTau12FakeErrDownZJets(NULL),
  HistoTau12FakeUnweighted(NULL), HistoTau12FakeUnweightedNoPU(NULL),
  HistoTau12FakeTauOSUnweightedWJets(NULL), HistoTau12FakeTauSSUnweightedWJets(NULL),
  HistoTau12FakeTauOSUnweightedZJets(NULL), HistoTau12FakeTauSSUnweightedZJets(NULL),
  HistoTau12FakeTauOSUnweightedErrUpWJets(NULL), HistoTau12FakeTauOSUnweightedErrDownWJets(NULL)
{
  if(strchr(desc, ';') != NULL || strchr(xtitle, ';') != NULL || strchr(ytitle, ';') != NULL)
    throw std::runtime_error("Histogram titles must not contain ';' characters");

  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  std::string xtitle_str = std::string(xtitle);
  std::string ytitle_str = std::string(ytitle);
  if(xunit) xtitle_str = xtitle_str + " [" + xunit + "]";
  //if(yunit) ytitle_str = ytitle_str + " [" + yunit + "]";

  const std::string axes = ";" + xtitle_str + ";" + ytitle_str;

  if(histograms & SIGNAL)
  {
    HistoSignal = new HistoType("signal", (desc + axes).c_str(), nbins, xmin, xmax);
    HistoSignalNoPU = new HistoType("signalNoPU", (std::string(desc) + ", no PU weight" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & FAKE_SS_WZ)
  {
    HistoTauSSFakeWJets = new HistoType("tauSSfakeWJets", (std::string(desc) + " of same sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauSSFakeZJets = new HistoType("tauSSfakeZJets", (std::string(desc) + " of same sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeWJets = new HistoType("tauOSfakeWJets", (std::string(desc) + " of opposite sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeZJets = new HistoType("tauOSfakeZJets", (std::string(desc) + " of opposite sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & (FAKE_SS_WZ | FAKE_SS_ERR) == (FAKE_SS_WZ | FAKE_SS_ERR))
  {
    HistoTauSSFakeErrUpWJets = new HistoType("tauSSfakeErrUpWJets", (std::string(desc) + " upper error of same sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauSSFakeErrDownWJets = new HistoType("tauSSfakeErrDownWJets", (std::string(desc) + " lower error of same sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauSSFakeErrUpZJets = new HistoType("tauSSfakeErrUpZJets", (std::string(desc) + " upper error of same sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauSSFakeErrDownZJets = new HistoType("tauSSfakeErrDownZJets", (std::string(desc) + " lower error of same sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeErrUpWJets = new HistoType("tauOSfakeErrUpWJets", (std::string(desc) + " upper error of oposite sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeErrDownWJets = new HistoType("tauOSfakeErrDownWJets", (std::string(desc) + " lower error of opposite sign #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeErrUpZJets = new HistoType("tauOSfakeErrUpZJets", (std::string(desc) + " upper error of opposite sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeErrDownZJets = new HistoType("tauOSfakeErrDownZJets", (std::string(desc) + " lower error of opposite sign #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & FAKE_SS_UNWEIGHTED)
  {
    HistoTauSSFakeUnweighted = new HistoType("tauSSfakeUnweighted", (std::string(desc) + " of unweighted same sign #tau fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauSSFakeUnweightedNoPU = new HistoType("tauSSfakeUnweightedNoPU", (std::string(desc) + " of unweighted same sign #tau fakes, no PU weight" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeUnweighted = new HistoType("tauOSfakeUnweighted", (std::string(desc) + " of unweighted opposite sign #tau fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTauOSFakeUnweightedNoPU = new HistoType("tauOSfakeUnweightedNoPU", (std::string(desc) + " of unweighted opposite sign #tau fakes, no PU weight" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & FAKE_12_WZ)
  {
    HistoTau12FakeWJets = new HistoType("tau12fakeWJets", (std::string(desc) + " of both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTau12FakeZJets = new HistoType("tau12fakeZJets", (std::string(desc) + " of both #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & (FAKE_12_WZ | FAKE_12_ERR) == (FAKE_12_WZ | FAKE_12_ERR))
  {
    HistoTau12FakeErrUpWJets = new HistoType("tau12fakeErrUpWJets", (std::string(desc) + " upper error of both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTau12FakeErrDownWJets = new HistoType("tau12fakeDownUpWJets", (std::string(desc) + " lower error of both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTau12FakeErrUpZJets = new HistoType("tau12fakeErrUpZJets", (std::string(desc) + " upper error of both #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTau12FakeErrDownZJets = new HistoType("tau12fakeErrDownZJets", (std::string(desc) + " lower error of both #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & FAKE_12_UNWEIGHTED)
  {
    HistoTau12FakeUnweighted = new HistoType("tau12fakeUnweighted", (std::string(desc) + " of unweighted both #tau fakes" + axes).c_str(), nbins, xmin, xmax);
    HistoTau12FakeUnweightedNoPU = new HistoType("tau12fakeUnweightedNoPU", (std::string(desc) + " of unweighted both #tau fakes, no PU weight" + axes).c_str(), nbins, xmin, xmax);

    if(histograms & FAKE_12_WZ)
    {
      HistoTau12FakeTauSSUnweightedWJets = new HistoType("tau12fakeTauSSUnweightedWJets",  (std::string(desc) + " of SS tau unweighted both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
      HistoTau12FakeTauOSUnweightedWJets = new HistoType("tau12fakeTauOSUnweightedWJets", (std::string(desc) + " of OS tau unweighted both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
      HistoTau12FakeTauSSUnweightedZJets = new HistoType("tau12fakeTauSSUnweightedZJets",  (std::string(desc) + " of SS tau unweighted both #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
      HistoTau12FakeTauOSUnweightedZJets = new HistoType("tau12fakeTauOSUnweightedZJets", (std::string(desc) + " of OS tau unweighted both #tau Z fakes" + axes).c_str(), nbins, xmin, xmax);
    }

    if(histograms & (FAKE_12_WZ | FAKE_12_ERR) == (FAKE_12_WZ | FAKE_12_ERR))
    {
      HistoTau12FakeTauOSUnweightedErrUpWJets = new HistoType("tau12fakeTauOSUnweightedErrUpWJets",  (std::string(desc) + " upper error of OS tau unweighted both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
      HistoTau12FakeTauOSUnweightedErrDownWJets = new HistoType("tau12fakeTauOSUnweightedErrDownWJets",  (std::string(desc) + " lower error of OS tau unweighted both #tau W fakes" + axes).c_str(), nbins, xmin, xmax);
    }
  }

  prevDirectory->cd();
}

template<typename HistoType>
void Variable<HistoType>::fill(const Filler& filler, double weight, const FillInfo& info) const
{
  const int regions = info.regions;
  const double ssFakeMetWeight = info.ssMassWeight * info.ssFakeMetWeight;

  if(regions & info.def.SIGNAL)
  {
    if(HistoSignal) filler(HistoSignal, weight);
    if(HistoSignalNoPU) filler(HistoSignalNoPU, 1.0);
  }
  if(regions & info.def.BACKGROUND)
  {
    if(HistoTauSSFakeWJets) filler(HistoTauSSFakeWJets, weight * ssFakeMetWeight * w(info.wssw));
    if(HistoTauSSFakeZJets) filler(HistoTauSSFakeZJets, weight * ssFakeMetWeight * w(info.wssz));
    if(HistoTauSSFakeErrUpWJets) filler(HistoTauSSFakeErrUpWJets, weight * ssFakeMetWeight * w(info.wssw_up));
    if(HistoTauSSFakeErrDownWJets) filler(HistoTauSSFakeErrDownWJets, weight * ssFakeMetWeight * w(info.wssw_down));
    if(HistoTauSSFakeErrUpZJets) filler(HistoTauSSFakeErrUpZJets, weight * ssFakeMetWeight * w(info.wssz_up));
    if(HistoTauSSFakeErrDownZJets) filler(HistoTauSSFakeErrDownZJets, weight * ssFakeMetWeight * w(info.wssz_down));
    if(HistoTauSSFakeUnweighted) filler(HistoTauSSFakeUnweighted, weight * ssFakeMetWeight);
    if(HistoTauSSFakeUnweightedNoPU) filler(HistoTauSSFakeUnweightedNoPU, ssFakeMetWeight);
  }
  if(regions & info.def.CONTROL)
  {
    if(HistoTauOSFakeWJets) filler(HistoTauOSFakeWJets, weight * w(info.wosw));
    if(HistoTauOSFakeZJets) filler(HistoTauOSFakeZJets, weight * w(info.wosz));
    if(HistoTauOSFakeErrUpWJets) filler(HistoTauOSFakeErrUpWJets, weight * w(info.wosw_up));
    if(HistoTauOSFakeErrDownWJets) filler(HistoTauOSFakeErrDownWJets, weight * w(info.wosw_down));
    if(HistoTauOSFakeErrUpZJets) filler(HistoTauOSFakeErrUpZJets, weight * w(info.wosz_up));
    if(HistoTauOSFakeErrDownZJets) filler(HistoTauOSFakeErrDownZJets, weight * w(info.wosz_down));
    if(HistoTauOSFakeUnweighted) filler(HistoTauOSFakeUnweighted, weight);
    if(HistoTauOSFakeUnweightedNoPU) filler(HistoTauOSFakeUnweightedNoPU, 1.0);
  }
  if(regions & info.def.CONTROL_BACKGROUND)
  {
    if(HistoTau12FakeWJets) filler(HistoTau12FakeWJets, weight * ssFakeMetWeight * w(info.wosw) * w(info.wssw));
    if(HistoTau12FakeZJets) filler(HistoTau12FakeZJets, weight * ssFakeMetWeight * w(info.wosz) * w(info.wssz));
    if(HistoTau12FakeErrUpWJets) filler(HistoTau12FakeErrUpWJets, weight * ssFakeMetWeight * w(info.wosw_up) * w(info.wssw_up));
    if(HistoTau12FakeErrDownWJets) filler(HistoTau12FakeErrDownWJets, weight * ssFakeMetWeight * w(info.wosw_down) * w(info.wssw_down));
    if(HistoTau12FakeErrUpZJets) filler(HistoTau12FakeErrUpZJets, weight * ssFakeMetWeight * w(info.wosz_up) * w(info.wssz_up));
    if(HistoTau12FakeErrDownZJets) filler(HistoTau12FakeErrDownZJets, weight * ssFakeMetWeight * w(info.wosz_down) * w(info.wssz_down));
    if(HistoTau12FakeUnweighted) filler(HistoTau12FakeUnweighted, weight * ssFakeMetWeight);
    if(HistoTau12FakeUnweightedNoPU) filler(HistoTau12FakeUnweightedNoPU, ssFakeMetWeight);
    if(HistoTau12FakeTauSSUnweightedWJets) filler(HistoTau12FakeTauSSUnweightedWJets, weight * ssFakeMetWeight * w(info.wosw));
    if(HistoTau12FakeTauOSUnweightedWJets) filler(HistoTau12FakeTauOSUnweightedWJets, weight * ssFakeMetWeight * w(info.wssw));
    if(HistoTau12FakeTauOSUnweightedErrUpWJets) filler(HistoTau12FakeTauOSUnweightedErrUpWJets, weight * ssFakeMetWeight * w(info.wssw_up));
    if(HistoTau12FakeTauOSUnweightedErrDownWJets) filler(HistoTau12FakeTauOSUnweightedErrDownWJets, weight * ssFakeMetWeight * w(info.wssw_down));
    if(HistoTau12FakeTauSSUnweightedZJets) filler(HistoTau12FakeTauSSUnweightedZJets, weight * ssFakeMetWeight * w(info.wosz));
    if(HistoTau12FakeTauOSUnweightedZJets) filler(HistoTau12FakeTauOSUnweightedZJets, weight * ssFakeMetWeight * w(info.wssz));
  }
}

template<typename HistoType>
VariableByTrackIso<HistoType>::VariableByTrackIso(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  normal.reset(new Variable<HistoType>(histograms, "normal", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  inclusive.reset(new Variable<HistoType>(histograms, "inclusive", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  metWeighted.reset(new Variable<HistoType>(histograms, "metWeighted", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  mu.reset(new Variable<HistoType>(histograms, "mu", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  muInclusive.reset(new Variable<HistoType>(histograms, "muInclusive", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void VariableByTrackIso<HistoType>::fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& infos) const
{
  if(infos.size() > 0) normal->fill(filler, weight, *infos[0]);
  if(infos.size() > 1) inclusive->fill(filler, weight, *infos[1]);
  if(infos.size() > 2) metWeighted->fill(filler, weight, *infos[2]);
  if(infos.size() > 3) mu->fill(filler, weight, *infos[3]);
  if(infos.size() > 4) muInclusive->fill(filler, weight, *infos[4]);
}

template<typename HistoType>
VariableByCategory<HistoType>::VariableByCategory(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  inclusive.reset(new VariableByTrackIso<HistoType>(histograms, "inclusive", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  low_mass.reset(new VariableByTrackIso<HistoType>(histograms, "low_mass", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  high_mass.reset(new VariableByTrackIso<HistoType>(histograms, "high_mass", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  //low_irred.reset(new VariableByTrackIso(histograms, "low_irred", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  //high_irred.reset(new VariableByTrackIso(histograms, "high_irred", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  //mass3_low.reset(new VariableByTrackIso(histograms, "mass3_low", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  //mass3_medium.reset(new VariableByTrackIso(histograms, "mass3_medium", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  //mass3_high.reset(new VariableByTrackIso(histograms, "mass3_high", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void VariableByCategory<HistoType>::fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& info) const
{
  assert(FillInfo::MASS_CUT >= 0.0);
  const double M = (info[0]->ostau.p4 + info[0]->sstau.p4).M();
  const double Irred = info[0]->irredMVAValue;

  inclusive->fill(filler, weight, info);
  if(M < FillInfo::MASS_CUT) low_mass->fill(filler, weight, info);
  else high_mass->fill(filler, weight, info);

  //if(M < 70.) mass3_low->fill(value, weight, info);
  //else if(M < 90.) mass3_medium->fill(value, weight, info);
  //else mass3_high->fill(value, weight, info);

  //if(Irred < 0.775) low_irred->fill(value, weight, info);
  //else high_irred->fill(value, weight, info);
}


template<typename HistoType>
VariableBySelection<HistoType>::VariableBySelection(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  uncut.reset(new VariableByCategory<HistoType>(histograms, "uncut", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  cut.reset(new VariableByCategory<HistoType>(histograms, "cut", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  mva.reset(new VariableByCategory<HistoType>(histograms, "mva", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  //mva_control.reset(new VariableByCategory<HistoType>(histograms, "mva_control", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void VariableBySelection<HistoType>::fill(const typename Variable<HistoType>::Filler& filler, double weight, const std::vector<const FillInfo*>& info) const
{
  assert(FillInfo::HIGHPT_TAU1_CUT >= 0.0 && FillInfo::HIGHPT_TAU2_CUT >= 0.0 && FillInfo::HIGHPT_MET_CUT >= 0.0);

  if(!info[0]->mvaCut)
  {
    uncut->fill(filler, weight, info);
    if(cut.get() && (TMath::Max(info[0]->ostau.Pt(), info[0]->sstau.Pt()) > FillInfo::HIGHPT_TAU1_CUT && TMath::Min(info[0]->ostau.Pt(), info[0]->sstau.Pt()) > FillInfo::HIGHPT_TAU2_CUT && info[0]->met > FillInfo::HIGHPT_MET_CUT))
      cut->fill(filler, weight, info);
  }
  else
  {
    if(mva.get()) mva->fill(filler, weight, info);
  }
}

template class Variable<TH1D>::Filler;
template class VariableByCategory<TH1D>;
template class VariableBySelection<TH1D>;

template class Variable<TProfile>::Filler;
template class VariableByCategory<TProfile>;
template class VariableBySelection<TProfile>;
