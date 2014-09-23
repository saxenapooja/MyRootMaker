
#include "Common/Cutflow.h"
#include "Common/LltVariable.h"

#include "MariaFunctions.h"

#include <TProfile.h>

namespace
{
  double w(double f) { if(f < 0.) return 0.; assert(f >= 0.0); if(f >= 1.0) return 0.; assert(f < 1.0); return f/(1.0 - f); }
}

LltFillInfo::RegionDef::RegionDef(CutflowRegions& regions, const char* data_era):
  SIGNAL(regions.AddRegion("Signal")),
  SIGNAL_LEADMU_FAKE(regions.AddRegion("SignalLeadMuFake")),
  SIGNAL_SUBMU_FAKE(regions.AddRegion("SignalSubMuFake")),
  SIGNAL_BOTHMU_FAKE(regions.AddRegion("SignalBothMuFake")),

  F3(regions.AddRegion("F3")),
  F3_LEADMU_FAKE(regions.AddRegion("F3LeadMuFake")),
  F3_SUBMU_FAKE(regions.AddRegion("F3SubMuFake")),
  F3_BOTHMU_FAKE(regions.AddRegion("F3BothMuFake"))
{
  if(strcmp(data_era, "2011") == 0 || strcmp(data_era, "fall11") == 0)
  {
    wzLumiScale = 0.00354118714899;
    zzLumiScale = 0.000445565029346;
    wjetsDataLead = 830.;
    wjetsDataSub = 3894.;
    wjetsWZLead = 4340. * wzLumiScale;
    wjetsWZSub = 7541. * wzLumiScale;
    wjetsZZLead = 6131. * zzLumiScale;
    wjetsZZSub = 9036. * zzLumiScale;
  }
  else if(strcmp(data_era, "2012") == 0 || strcmp(data_era, "summer12") == 0)
  {
    wzLumiScale = 0.010324537742;
    zzLumiScale = 0.000532967626359;
    wjetsDataLead = 5665.;
    wjetsDataSub = 25400.;
    wjetsWZLead = 7780. * wzLumiScale;
    wjetsWZSub = 13333. * wzLumiScale;
    wjetsZZLead = 17683. * zzLumiScale;
    wjetsZZSub = 29438. * zzLumiScale;
  }
  else
    throw std::runtime_error("Unknown data era: " + std::string(data_era));
}

LltFillInfo::LltFillInfo(const RegionDef& regionDef, const FakeRateFunctions& fr1Jet, const FakeRateFunctions& fr2Jet, const Basic_kNNEvaluator& kNNData, const Basic_kNNEvaluator* kNNWZ, const Basic_kNNEvaluator* kNNZZ, const Lepton& leadLep, const Lepton& subLep, const TauCand& tau, double met, unsigned int nJets20, int regions):
  def(regionDef), leadlep(leadLep), sublep(subLep), tau(tau),
  regions(regions),

  lt(leadlep.Pt() + subLep.Pt() + tau.Pt()),
  met(met),
  njets20(nJets20),

  wfr1jet(fr1Jet.TauOS->getPtEta(FakeRateFunctions::WJETS_FUNC, tau.Pt(), tau.Eta(), false)),
  wfr2jet(fr2Jet.TauOS->getPtEta(FakeRateFunctions::WJETS_FUNC, tau.Pt(), tau.Eta(), false)),

  wkNNleadData(kNNData.getLeadingMuon(leadLep.Pt(), leadLep.JetPt, nJets20)),
  wkNNsubData(kNNData.getSubleadingMuon(subLep.Pt(), subLep.JetPt, nJets20)),

  wkNNleadWZ(kNNWZ ? kNNWZ->getLeadingMuon(leadLep.Pt(), leadLep.JetPt, nJets20) : 0.),
  wkNNsubWZ(kNNWZ ? kNNWZ->getSubleadingMuon(subLep.Pt(), subLep.JetPt, nJets20) : 0.),

  wkNNleadZZ(kNNZZ ? kNNZZ->getLeadingMuon(leadLep.Pt(), leadLep.JetPt, nJets20) : 0.),
  wkNNsubZZ(kNNZZ ? kNNZZ->getSubleadingMuon(subLep.Pt(), subLep.JetPt, nJets20) : 0.),

  wkNNlead(wkNNleadData - (def.wjetsWZLead / def.wjetsDataLead) * wkNNleadWZ - (def.wjetsZZLead / def.wjetsDataLead) * wkNNleadZZ),
  wkNNsub(wkNNsubData - (def.wjetsWZSub / def.wjetsDataSub) * wkNNsubWZ - (def.wjetsZZSub / def.wjetsDataSub) * wkNNsubZZ)
{
}

template<typename HistoType>
LltVariable<HistoType>::Filler::Filler(double x, double y):
  x(x), y(y)
{
}

template<>
void LltVariable<TH1D>::Filler::operator()(TH1D* histo, double weight) const
{
  histo->Fill(x, weight);
}

template<>
void LltVariable<TProfile>::Filler::operator()(TProfile* prof, double weight) const
{
  prof->Fill(x, y, weight);
}

template<typename HistoType>
LltVariable<HistoType>::LltVariable(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax):
  HistoSignal(NULL), HistoSignalNoPU(NULL)
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

    HistoSignalLeadMuFakeUnweighted = new HistoType("signalLeadMuFakeUnweighted", (std::string(desc) + ", lead mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalLeadMuFakekNNWeighted = new HistoType("signalLeadMuFakekNNWeighted", (std::string(desc) + ", lead mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalLeadMuFakekNNWeightedData = new HistoType("signalLeadMuFakekNNWeightedData", (std::string(desc) + ", lead mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);

    HistoSignalSubMuFakeUnweighted = new HistoType("signalSubMuFakeUnweighted", (std::string(desc) + ", sub mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalSubMuFakekNNWeighted = new HistoType("signalSubMuFakekNNWeighted", (std::string(desc) + ", sub mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalSubMuFakekNNWeightedData = new HistoType("signalSubMuFakekNNWeightedData", (std::string(desc) + ", sub mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);

    HistoSignalBothMuFakeUnweighted = new HistoType("signalBothMuFakeUnweighted", (std::string(desc) + ", both mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalBothMuFakekNNWeighted = new HistoType("signalBothMuFakekNNWeighted", (std::string(desc) + ", both mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoSignalBothMuFakekNNWeightedData = new HistoType("signalBothMuFakekNNWeightedData", (std::string(desc) + ", both mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
  }

  if(histograms & F3)
  {
    HistoF3 = new HistoType("f3", (std::string(desc) + ", F3" + axes).c_str(), nbins, xmin, xmax);
    HistoF3NoPU = new HistoType("f3NoPU", (std::string(desc) + ", F3, no PU weight" + axes).c_str(), nbins, xmin, xmax);
    HistoF3W1JetWeighted = new HistoType("f3W1JetWeighted", (std::string(desc) + ", F3, 1-jet weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3W2JetWeighted = new HistoType("f3W2JetWeighted", (std::string(desc) + ", F3, 2-jet weighted" + axes).c_str(), nbins, xmin, xmax);

    HistoF3LeadMuFakeUnweighted = new HistoType("f3LeadMuFakeUnweighted", (std::string(desc) + ", F3, lead mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3LeadMuFakekNNWeighted = new HistoType("f3LeadMuFakekNNWeighted", (std::string(desc) + ", F3, lead mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3LeadMuFakekNNWeightedData = new HistoType("f3LeadMuFakekNNWeightedData", (std::string(desc) + ", F3, lead mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);

    HistoF3SubMuFakeUnweighted = new HistoType("f3SubMuFakeUnweighted", (std::string(desc) + ", F3, sub mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3SubMuFakekNNWeighted = new HistoType("f3SubMuFakekNNWeighted", (std::string(desc) + ", F3, sub mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3SubMuFakekNNWeightedData = new HistoType("f3SubMuFakekNNWeightedData", (std::string(desc) + ", F3, sub mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);

    HistoF3BothMuFakeUnweighted = new HistoType("f3BothMuFakeUnweighted", (std::string(desc) + ", F3, both mu fake, unweighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3BothMuFakekNNWeighted = new HistoType("f3BothMuFakekNNWeighted", (std::string(desc) + ", F3, both mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
    HistoF3BothMuFakekNNWeightedData = new HistoType("f3BothMuFakekNNWeightedData", (std::string(desc) + ", F3, both mu fake, kNN weighted" + axes).c_str(), nbins, xmin, xmax);
  }

  prevDirectory->cd();
}

template<typename HistoType>
void LltVariable<HistoType>::fill(const Filler& filler, double weight, const LltFillInfo& info) const
{
  const int regions = info.regions;

  // Signal
  if(regions & info.def.SIGNAL)
  {
    if(HistoSignal) filler(HistoSignal, weight);
    if(HistoSignalNoPU) filler(HistoSignalNoPU, 1.0);
  }

  if(regions & info.def.SIGNAL_LEADMU_FAKE)
  {
    if(HistoSignalLeadMuFakeUnweighted) filler(HistoSignalLeadMuFakeUnweighted, weight);
    if(HistoSignalLeadMuFakekNNWeighted) filler(HistoSignalLeadMuFakekNNWeighted, weight * w(info.wkNNlead));
    if(HistoSignalLeadMuFakekNNWeightedData) filler(HistoSignalLeadMuFakekNNWeightedData, weight * w(info.wkNNleadData));
  }

  if(regions & info.def.SIGNAL_SUBMU_FAKE)
  {
    if(HistoSignalSubMuFakeUnweighted) filler(HistoSignalSubMuFakeUnweighted, weight);
    if(HistoSignalSubMuFakekNNWeighted) filler(HistoSignalSubMuFakekNNWeighted, weight * w(info.wkNNsub));
    if(HistoSignalSubMuFakekNNWeightedData) filler(HistoSignalSubMuFakekNNWeightedData, weight * w(info.wkNNsubData));
  }

  if(regions & info.def.SIGNAL_BOTHMU_FAKE)
  {
    if(HistoSignalBothMuFakeUnweighted) filler(HistoSignalBothMuFakeUnweighted, weight);
    if(HistoSignalBothMuFakekNNWeighted) filler(HistoSignalBothMuFakekNNWeighted, weight * w(info.wkNNlead) * w(info.wkNNsub));
    if(HistoSignalBothMuFakekNNWeightedData) filler(HistoSignalBothMuFakekNNWeightedData, weight * w(info.wkNNleadData) * w(info.wkNNsubData));
  }

  // F3
  if(regions & info.def.F3)
  {
    if(HistoF3) filler(HistoF3, weight);
    if(HistoF3NoPU) filler(HistoF3NoPU, 1.0);
    if(HistoF3W1JetWeighted) filler(HistoF3W1JetWeighted, weight * w(info.wfr1jet));
    if(HistoF3W2JetWeighted) filler(HistoF3W2JetWeighted, weight * w(info.wfr2jet));
  }

  if(regions & info.def.F3_LEADMU_FAKE)
  {
    if(HistoF3LeadMuFakeUnweighted) filler(HistoF3LeadMuFakeUnweighted, weight);
    if(HistoF3LeadMuFakekNNWeighted) filler(HistoF3LeadMuFakekNNWeighted, weight * w(info.wkNNlead));
    if(HistoF3LeadMuFakekNNWeightedData) filler(HistoF3LeadMuFakekNNWeightedData, weight * w(info.wkNNleadData));
  }

  if(regions & info.def.F3_SUBMU_FAKE)
  {
    if(HistoF3SubMuFakeUnweighted) filler(HistoF3SubMuFakeUnweighted, weight);
    if(HistoF3SubMuFakekNNWeighted) filler(HistoF3SubMuFakekNNWeighted, weight * w(info.wkNNsub));
    if(HistoF3SubMuFakekNNWeightedData) filler(HistoF3SubMuFakekNNWeightedData, weight * w(info.wkNNsubData));
  }

  if(regions & info.def.F3_BOTHMU_FAKE)
  {
    if(HistoF3BothMuFakeUnweighted) filler(HistoF3BothMuFakeUnweighted, weight);
    if(HistoF3BothMuFakekNNWeighted) filler(HistoF3BothMuFakekNNWeighted, weight * w(info.wkNNlead) * w(info.wkNNsub));
    if(HistoF3BothMuFakekNNWeightedData) filler(HistoF3BothMuFakekNNWeightedData, weight * w(info.wkNNleadData) * w(info.wkNNsubData));
  }
}

template<typename HistoType>
LltVariableByTrackIso<HistoType>::LltVariableByTrackIso(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  frOrderedkNN.reset(new LltVariable<HistoType>(histograms, "frOrderedkNN", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  frUnorderedkNN.reset(new LltVariable<HistoType>(histograms, "frUnorderedkNN", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  frDY.reset(new LltVariable<HistoType>(histograms, "frDY", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  frTT.reset(new LltVariable<HistoType>(histograms, "frTT", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void LltVariableByTrackIso<HistoType>::fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& infos) const
{
  if(infos.size() > 0) frOrderedkNN->fill(filler, weight, *infos[0]);
  if(infos.size() > 1) frUnorderedkNN->fill(filler, weight, *infos[1]);
  if(infos.size() > 2) frDY->fill(filler, weight, *infos[2]);
  if(infos.size() > 3) frTT->fill(filler, weight, *infos[3]);
}

template<typename HistoType>
LltVariableByCategory<HistoType>::LltVariableByCategory(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  inclusive.reset(new LltVariableByTrackIso<HistoType>(histograms, "inclusive", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  low_lt.reset(new LltVariableByTrackIso<HistoType>(histograms, "low_lt", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));
  high_lt.reset(new LltVariableByTrackIso<HistoType>(histograms, "high_lt", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  lt80.reset(new LltVariableByTrackIso<HistoType>(histograms, "lt80", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void LltVariableByCategory<HistoType>::fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& info) const
{
  inclusive->fill(filler, weight, info);

  if(info[0]->lt < 130.) low_lt->fill(filler, weight, info);
  else high_lt->fill(filler, weight, info);

  if(info[0]->lt > 80.) lt80->fill(filler, weight, info);
}

template<typename HistoType>
LltVariableBySelection<HistoType>::LltVariableBySelection(int histograms, const char* name, const char* desc, const char* xtitle, const char* xunit, const char* ytitle, unsigned int nbins, double xmin, double xmax)
{
  TDirectory* prevDirectory = gDirectory;

  gDirectory->mkdir(name);
  gDirectory->cd(name);

  uncut.reset(new LltVariableByCategory<HistoType>(histograms, "uncut", desc, xtitle, xunit, ytitle, nbins, xmin, xmax));

  prevDirectory->cd();
}

template<typename HistoType>
void LltVariableBySelection<HistoType>::fill(const typename LltVariable<HistoType>::Filler& filler, double weight, const std::vector<const LltFillInfo*>& info) const
{
  uncut->fill(filler, weight, info);
}

template class LltVariable<TH1D>::Filler;
template class LltVariableByCategory<TH1D>;
template class LltVariableBySelection<TH1D>;

template class LltVariable<TProfile>::Filler;
template class LltVariableByCategory<TProfile>;
template class LltVariableBySelection<TProfile>;
