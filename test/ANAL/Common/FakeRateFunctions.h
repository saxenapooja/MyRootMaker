#ifndef _WHANALYSIS_FAKERATEFUNCTIONS_H_
#define _WHANALYSIS_FAKERATEFUNCTIONS_H_

#include <TFile.h>
#include <TF1.h>
#include <memory>

class FakeRateFunctions
{
public:
  FakeRateFunctions(const char* directory, const char* data_era, const char* channel, const char* tauOS_antie, const char* tauOS_isolation, const char* tauSS_antie, const char* tauSS_isolation, bool lepTauSS, bool inclusive, bool capped);
  
  enum Function {
    WJETS_FUNC,
    WJETS_ERR_UP,
    WJETS_ERR_DOWN,
    ZJETS_FUNC,
    ZJETS_ERR_UP,
    ZJETS_ERR_DOWN
  };

private:
  TFile* fileFakeRate;

public:
  struct FakeRateFunc {
    FakeRateFunc(TFile* file, const char* label);
    TF1* func;
    TF1* err;
  };

  struct FakeRateSet {
    FakeRateSet(TFile* file, const char* data_era, const char* channel, const char* antie, const char* isolation, bool lepTauSS, bool capped);
    FakeRateFunc WJetsInclusive;
    FakeRateFunc WJetsBarrel;
    FakeRateFunc WJetsTransition;
    FakeRateFunc WJetsForward;
    
    FakeRateFunc ZJetsInclusive;
    FakeRateFunc ZJetsBarrel;
    FakeRateFunc ZJetsTransition;
    FakeRateFunc ZJetsForward;
  };
  
  struct FakeRateIsolations {
    FakeRateIsolations(TFile* file, const char* data_era, const char* channel, const char* antie, const char* isolation, bool lepTauSS, bool inclusive, bool capped);
    double getPtOnly(Function func, double pt, bool case2);
    double getPtEta(Function func, double pt, double eta, bool case2);
    
    FakeRateSet RegularIso;
    std::auto_ptr<FakeRateSet> TrackIso; // optional
    std::auto_ptr<FakeRateSet> TrackIsoIfRegularIso; // optional

    const bool INCLUSIVE;
  };

  std::auto_ptr<FakeRateIsolations> TauOS;
  std::auto_ptr<FakeRateIsolations> TauSS;
};

#endif // _WHANALYSIS_FAKERATEFUNCTIONS_H_
