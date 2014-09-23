#include "Common/FakeRateFunctions.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <cassert>

namespace
{
  struct WZFuncs {
    FakeRateFunctions::FakeRateFunc& W;
    FakeRateFunctions::FakeRateFunc& Z;
  };
  
  double evalFunc(const WZFuncs& funcs, FakeRateFunctions::Function func, double pt)
  {
    pt = std::min(100., pt);
    
    switch(func)
      {
      case FakeRateFunctions::WJETS_FUNC:     return funcs.W.func->Eval(pt);
      case FakeRateFunctions::WJETS_ERR_UP:   return funcs.W.func->Eval(pt) + funcs.W.err->Eval(pt);
      case FakeRateFunctions::WJETS_ERR_DOWN: return funcs.W.func->Eval(pt) - funcs.W.err->Eval(pt);
      case FakeRateFunctions::ZJETS_FUNC:     return funcs.Z.func->Eval(pt);
      case FakeRateFunctions::ZJETS_ERR_UP:   return funcs.Z.func->Eval(pt) + funcs.Z.err->Eval(pt);
      case FakeRateFunctions::ZJETS_ERR_DOWN: return funcs.Z.func->Eval(pt) - funcs.Z.err->Eval(pt);
      default: assert(false);                 return -1.;
      }
  }
  
  const char* cappedString(bool capped)
  {
    if(capped) return "Capped";
    else return "";
  }
  
  const char* lepTauString(bool lepTauSS)
  {
    if(lepTauSS) return "lepTauSS";
    else return "lepTauOS";
  }
}

FakeRateFunctions::FakeRateFunc::FakeRateFunc(TFile* file, const char* label)
{
  std::cout << label << "... ";
  func = dynamic_cast<TF1*>(file->Get( (std::string(label) + "/landau").c_str()));
  err = dynamic_cast<TF1*>(file->Get( (std::string(label) + "/landauErr").c_str()));
  if(!func || !err)
    throw std::runtime_error("Not all Fake Rate functions could be loaded!");
  std::cout << "done!" << std::endl;
}

FakeRateFunctions::FakeRateSet::FakeRateSet(TFile* file, const char* data_era, const char* channel, const char* antie, const char* isolation, bool lepTauSS, bool capped):
  WJetsInclusive(file, (std::string(data_era) + "_" + std::string(channel) + "_" + std::string(antie) + "_" + std::string(lepTauString(lepTauSS)) + "_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_inclusive").c_str()),
  WJetsBarrel(file, (std::string(data_era) + "_" + std::string(channel) + "_" + std::string(antie) + "_" + std::string(lepTauString(lepTauSS)) + "_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_barrel").c_str()),
  WJetsTransition(file, (std::string(data_era) + "_" + std::string(channel) + "_" + std::string(antie) + "_" + std::string(lepTauString(lepTauSS)) + "_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_transition").c_str()),
  WJetsForward(file, (std::string(data_era) + "_" + std::string(channel) + "_" + std::string(antie) + "_" + std::string(lepTauString(lepTauSS)) + "_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_forward").c_str()),

  ZJetsInclusive(file, (std::string(data_era) + "_zJetsAll_" + std::string(antie) + "_all_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_inclusive").c_str()),
  ZJetsBarrel(file, (std::string(data_era) + "_zJetsAll_" + std::string(antie) + "_all_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_barrel").c_str()),
  ZJetsTransition(file, (std::string(data_era) + "_zJetsAll_" + std::string(antie) + "_all_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_transition").c_str()),
  ZJetsForward(file, (std::string(data_era) + "_zJetsAll_" + std::string(antie) + "_all_fakeable" + cappedString(capped) + "_" + std::string(isolation) + "_byEta_projX_forward").c_str())
{
}

FakeRateFunctions::FakeRateIsolations::FakeRateIsolations(TFile* file, const char* data_era, const char* channel, const char* antie, const char* isolation, bool lepTauSS, bool inclusive, bool capped):
  INCLUSIVE(inclusive),
  RegularIso(file, data_era, channel, antie, isolation, lepTauSS, capped)
{
  //TrackIso.reset(new FakeRateSet(file, data_era, channel, antie, "trackIsolated", capped));
  //TrackIsoIfRegularIso.reset(new FakeRateSet(file, data_era, channel, antie, (std::string(isolation) + "TrackIsolated", capped).c_str()));
}

double FakeRateFunctions::FakeRateIsolations::getPtOnly(Function func, double pt, bool case2)
{
  double p;
  assert(!case2 || TrackIsoIfRegularIso.get() != NULL);

  const WZFuncs isoFuncs = { RegularIso.WJetsInclusive, RegularIso.ZJetsInclusive };

  if(case2)
  {
    const WZFuncs trackIsoIfRegularIsoFuncs = { TrackIsoIfRegularIso->WJetsInclusive, TrackIsoIfRegularIso->ZJetsInclusive };
    p = evalFunc(isoFuncs, func, pt) * evalFunc(trackIsoIfRegularIsoFuncs, func, pt);
  }
  else
  {
    p = evalFunc(isoFuncs, func, pt);
  }

  return std::min(std::max(0.0, p), 0.999);
}

double FakeRateFunctions::FakeRateIsolations::getPtEta(Function func, double pt, double eta, bool case2)
{
  double p;
  assert(!case2); // TODO: Fix for case2

  WZFuncs isoFuncsInclusive = { RegularIso.WJetsInclusive, RegularIso.ZJetsInclusive };
  WZFuncs trackIsoIfRegularIsoFuncsInclusive = { RegularIso.WJetsInclusive, RegularIso.ZJetsInclusive };
  WZFuncs isoFuncsBarrel = { RegularIso.WJetsBarrel, RegularIso.ZJetsBarrel };
  WZFuncs trackIsoIfRegularIsoFuncsBarrel = { RegularIso.WJetsBarrel, RegularIso.ZJetsBarrel };
  WZFuncs isoFuncsTransition = { RegularIso.WJetsTransition, RegularIso.ZJetsTransition };
  WZFuncs trackIsoIfRegularIsoFuncsTransition = { RegularIso.WJetsTransition, RegularIso.ZJetsTransition };
  WZFuncs isoFuncsForward = { RegularIso.WJetsForward, RegularIso.ZJetsForward };
  WZFuncs trackIsoIfRegularIsoFuncsForward = { RegularIso.WJetsForward, RegularIso.ZJetsForward };

  const WZFuncs* isoFuncs = NULL;
  const WZFuncs* trackIsoIfRegularIsoFuncs = NULL;
  if(INCLUSIVE) { isoFuncs = &isoFuncsInclusive; trackIsoIfRegularIsoFuncs = &trackIsoIfRegularIsoFuncsInclusive; }
  else if(fabs(eta) < 0.8) { isoFuncs = &isoFuncsBarrel; trackIsoIfRegularIsoFuncs = &trackIsoIfRegularIsoFuncsBarrel; }
  else if(fabs(eta) < 1.6) { isoFuncs = &isoFuncsTransition; trackIsoIfRegularIsoFuncs = &trackIsoIfRegularIsoFuncsTransition; }
  else { isoFuncs = &isoFuncsForward; trackIsoIfRegularIsoFuncs = &trackIsoIfRegularIsoFuncsForward; }

  if(case2)
    p = evalFunc(*isoFuncs, func, pt) * evalFunc(*trackIsoIfRegularIsoFuncs, func, pt);
  else
    p = evalFunc(*isoFuncs, func, pt);

  return std::min(std::max(0.0, p), 1.0);
}

FakeRateFunctions::FakeRateFunctions(const char* directory, const char* data_era, const char* channel, const char* tauOS_antie, const char* tauOS_isolation, const char* tauSS_antie, const char* tauSS_isolation, bool lepTauSS, bool inclusive, bool capped)
{
  std::cout << "Loading " << data_era << " Fake Rates... " << std::endl;

  std::string filename = "/nfs/dust/cms/user/aburgm/tautau-fullyhadronic/WH_Analysis/WH_Analysis/FakeRate/" + std::string(directory) + "/" + std::string(data_era) + ".root";
  fileFakeRate = new TFile(filename.c_str(), "READ");
  if(fileFakeRate->IsZombie())
    throw std::runtime_error("Fake Rate function file could not be opened!");

  if(tauOS_antie != NULL && tauOS_isolation != NULL)
    TauOS.reset(new FakeRateIsolations(fileFakeRate, data_era, channel, tauOS_antie, tauOS_isolation, lepTauSS, inclusive, capped));
  if(tauSS_antie != NULL && tauSS_isolation != NULL)
    TauSS.reset(new FakeRateIsolations(fileFakeRate, data_era, channel, tauSS_antie, tauSS_isolation, lepTauSS, inclusive, capped));

  std::cout << "Done!" << std::endl;
}
