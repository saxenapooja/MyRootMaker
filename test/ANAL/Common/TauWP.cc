#include "Common/TauWP.h"
#include <cassert>
#include <stdexcept>

TauWP::Isolation string_to_isolation(const std::string& str)
{
  if(str == "loose")             return TauWP::ISO_LOOSE;
  else if(str == "loose_3hits")  return TauWP::ISO_LOOSE_3HITS;
  else if(str == "medium")       return TauWP::ISO_MEDIUM;
  else if(str == "medium_3hits") return TauWP::ISO_MEDIUM_3HITS;
  else if(str == "tight")        return TauWP::ISO_TIGHT;
  else if(str == "tight_3hits")  return TauWP::ISO_TIGHT_3HITS;
  else throw std::runtime_error("No such tau isolation working point: " + str);
}

TauWP::AntiE string_to_antie(const std::string& str)
{
  if(str == "loose")            return TauWP::ANTIE_LOOSE;
  else if(str == "loose_mva3")  return TauWP::ANTIE_LOOSE_MVA3;
  else if(str == "medium")      return TauWP::ANTIE_MEDIUM;
  else if(str == "medium_mva3") return TauWP::ANTIE_MEDIUM_MVA3;
  else if(str == "tight")       return TauWP::ANTIE_TIGHT;
  else if(str == "tight_mva3")  return TauWP::ANTIE_TIGHT_MVA3;
  else throw std::runtime_error("No such electron rejection working point: " + str);
}

TauWP::AntiMu string_to_antimu(const std::string& str)
{
  if(str == "loose")        return TauWP::ANTIMU_LOOSE;
  else if(str == "medium")  return TauWP::ANTIMU_MEDIUM;
  else if(str == "tight")   return TauWP::ANTIMU_LOOSE;
  else throw std::runtime_error("No such muon rejection working point: " + str);
}

const bool TauCand::*get_isolation_flag(TauWP::Isolation isolation)
{
  switch(isolation)
  {
  case TauWP::ISO_LOOSE:        return &TauCand::byLooseCombinedIsolationDeltaBetaCorr;
  case TauWP::ISO_LOOSE_3HITS:  return &TauCand::byLooseCombinedIsolationDeltaBetaCorr3Hits;
  case TauWP::ISO_MEDIUM:       return &TauCand::byMediumCombinedIsolationDeltaBetaCorr;
  case TauWP::ISO_MEDIUM_3HITS: return &TauCand::byMediumCombinedIsolationDeltaBetaCorr3Hits;
  case TauWP::ISO_TIGHT:        return &TauCand::byTightCombinedIsolationDeltaBetaCorr;
  case TauWP::ISO_TIGHT_3HITS:  return &TauCand::byTightCombinedIsolationDeltaBetaCorr3Hits;
  default: assert(false);       return NULL;
  }
}

const bool TauCand::*get_antie_flag(TauWP::AntiE antie)
{
  switch(antie)
    {
    case TauWP::ANTIE_LOOSE: return &TauCand::againstElectronLoose;
    case TauWP::ANTIE_LOOSE_MVA3: return &TauCand::againstElectronLooseMVA3;
    case TauWP::ANTIE_MEDIUM: return &TauCand::againstElectronMedium;
    case TauWP::ANTIE_MEDIUM_MVA3: return &TauCand::againstElectronMediumMVA3;
    case TauWP::ANTIE_TIGHT: return &TauCand::againstElectronTight;
    case TauWP::ANTIE_TIGHT_MVA3: return &TauCand::againstElectronTightMVA3;
    default: assert(false); return NULL;
    }
}

const bool TauCand::*get_antimu_flag(TauWP::AntiMu antimu)
{
  switch(antimu)
    {
    case TauWP::ANTIMU_LOOSE: return &TauCand::againstMuonLoose;
    case TauWP::ANTIMU_MEDIUM: return &TauCand::againstMuonMedium;
      //case TauWP::ANTIMU_TIGHT: return &TauCand::againstMuonTight;
    case TauWP::ANTIMU_TIGHT: return &TauCand::againstMuonLoose;
    default: assert(false); return NULL;
    }
}

TauWP::TauWP(const Config& cfg,  const char* prefix,  const char* defaultIsolation,    const char* defaultAntiE,   const char* defaultAntiMu):
  ISOLATION(string_to_isolation(cfg.get<std::string>(std::string(prefix) + "_isolation", defaultIsolation))),
  ANTIE(string_to_antie(cfg.get<std::string>(std::string(prefix) + "_antie", defaultAntiE))),
  ANTIMU(string_to_antimu(cfg.get<std::string>(std::string(prefix) + "_antimu", defaultAntiMu))),
  ISOLATION_FLAG(get_isolation_flag(ISOLATION)),
  ANTIE_FLAG(get_antie_flag(ANTIE)),
  ANTIMU_FLAG(get_antimu_flag(ANTIMU))
{
}

// TauWP::TauWP(const Config& cfg,  const char* prefix,  const char* defaultIsolation):
//   ISOLATION(string_to_isolation(cfg.get<std::string>(std::string(prefix) + "_isolation", defaultIsolation))),
//   ISOLATION_FLAG(get_isolation_flag(ISOLATION))
// {
// }



const char* TauWP::isolation_hr_name() const
{
  switch(ISOLATION)
    {
    case ISO_LOOSE: return "Loose";
    case ISO_LOOSE_3HITS: return "Loose3Hits";
    case ISO_MEDIUM: return "Medium";
    case ISO_MEDIUM_3HITS: return "Medium3Hits";
    case ISO_TIGHT: return "Tight";
    case ISO_TIGHT_3HITS: return "Tight3Hits";
    default: assert(false); return "";
    }
}

const char* TauWP::antie_hr_name() const
{
  switch(ANTIE)
    {
    case ANTIE_LOOSE: return "Loose";
    case ANTIE_LOOSE_MVA3: return "LooseMVA3";
    case ANTIE_MEDIUM: return "Medium";
    case ANTIE_MEDIUM_MVA3: return "MediumMVA3";
    case ANTIE_TIGHT: return "Tight";
    case ANTIE_TIGHT_MVA3: return "TightMVA3";
    default: assert(false); return "";
    }
}

const char* TauWP::antimu_hr_name() const
{
  switch(ANTIMU)
    {
    case ANTIMU_LOOSE: return "Loose";
    case ANTIMU_MEDIUM: return "Medium";
    case ANTIMU_TIGHT: return "Tight";
  }
}

const char* TauWP::isolation_frfunc_name() const
{
  switch(ISOLATION)
    {
    case ISO_LOOSE: return "combinedIsolationLoose";
    case ISO_LOOSE_3HITS: return "combinedIsolationLoose3Hits";
    case ISO_MEDIUM: return "combinedIsolationMedium";
    case ISO_MEDIUM_3HITS: return "combinedIsolationMedium3Hits";
    case ISO_TIGHT: return "combinedIsolationTight";
    case ISO_TIGHT_3HITS: return "combinedIsolationTight3Hits";
    default: assert(false); return "";
    }
}

const char* TauWP::antie_frfunc_name() const
{
  switch(ANTIMU)
    {
    case ANTIMU_LOOSE:
      switch(ANTIE)
	{
	case ANTIE_LOOSE: return "electronLoose_muonLoose";
	case ANTIE_LOOSE_MVA3: return "electronLooseMVA3_muonLoose";
	case ANTIE_MEDIUM: return "electronMedium_muonLoose";
	case ANTIE_MEDIUM_MVA3: return "electronMediumMVA3_muonLoose";
	case ANTIE_TIGHT: return "electronTight_muonLoose";
	case ANTIE_TIGHT_MVA3: return "electronTightMVA3_muonLoose";
	default: assert(false); return "";
	}
      break;
    case ANTIMU_MEDIUM:
      switch(ANTIE)
	{
	case ANTIE_LOOSE: return "electronLoose_muonMedium";
	case ANTIE_LOOSE_MVA3: return "electronLooseMVA3_muonMedium";
	case ANTIE_MEDIUM: return "electronMedium_muonMedium";
	case ANTIE_MEDIUM_MVA3: return "electronMediumMVA3_muonMedium";
	case ANTIE_TIGHT: return "electronTight_muonMedium";
	case ANTIE_TIGHT_MVA3: return "electronTightMVA3_muonMedium";
	default: assert(false); return "";
	}
      break;
    case ANTIMU_TIGHT:
      switch(ANTIE)
	{
	case ANTIE_LOOSE: return "electronLoose_muonTight";
	case ANTIE_LOOSE_MVA3: return "electronLooseMVA3_muonTight";
	case ANTIE_MEDIUM: return "electronMedium_muonTight";
	case ANTIE_MEDIUM_MVA3: return "electronMediumMVA3_muonTight";
	case ANTIE_TIGHT: return "electronTight_muonTight";
	case ANTIE_TIGHT_MVA3: return "electronTightMVA3_muonTight";
	default: assert(false); return "";
	}
      break;
    }
}
