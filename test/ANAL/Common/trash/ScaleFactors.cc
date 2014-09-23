#include "Common/ScaleFactors.h"

#include <cmath>
#include <cassert>
#include <stdexcept>

#include <TH1.h>
#include "MuonPOG2011/efficiencyFunctions.C"
#include "MuonPOG2011/efficiencyFunctionsIDIso.C"

namespace
{

double Corr_Trg_Mu_2012_53X(double mupt, double mueta)
{
  if(mupt >= 10 && mupt < 15)
  {
    if (fabs(mueta) < 0.8) return 0.9818;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9713;
    else if (1.2 <= fabs(mueta)) return 0.9675;
  }
  else if(mupt >= 15 && mupt < 20)
  {
    if (fabs(mueta) < 0.8) return 0.9781;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9782;
    else if (1.2 <= fabs(mueta) ) return 0.9587;
  }
  else if(mupt >= 20 && mupt < 25)
  {
    if (fabs(mueta) < 0.8) return 0.9873;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9532;
    else if (1.2 <= fabs(mueta) ) return 0.9605;
  }
  else if(mupt >= 25 && mupt < 30)
  {
    if (fabs(mueta) < 0.8) return 0.9755;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9818;
    else if (1.2 <= fabs(mueta) ) return 0.9632;
  }
  else if(mupt >= 30)
  {
    if (fabs(mueta) < 0.8) return 0.9956;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9644;
    else if (1.2 <= fabs(mueta) ) return 0.9530;
  }

  assert(false);
}

double Corr_IdIso_Mu_2012_53X(double mupt, double mueta)
{
  if(mupt >= 10 && mupt < 15)
  {
    if (fabs(mueta) < 0.8) return 0.9771;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9746;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9644;
    else if (1.6 <= fabs(mueta)) return 0.9891;
  }
  else if(mupt >= 15 && mupt < 20)
  {
    if (fabs(mueta) < 0.8) return 0.9548;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9701;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9766;
    else if (1.6 <= fabs(mueta) ) return 0.9892;
  }
  else if(mupt >= 20 && mupt < 25)
  {
    if (fabs(mueta) < 0.8) return 0.9648;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9836;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9820;
    else if (1.6 <= fabs(mueta) ) return 0.9909;
  }
  else if(mupt >= 25 && mupt < 30)
  {
    if (fabs(mueta) < 0.8) return 0.9676;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9817;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9886;
    else if (1.6 <= fabs(mueta) ) return 0.9883;
  }
  else if(mupt >= 30 && mupt < 35)
  {
    if (fabs(mueta) < 0.8) return 0.9730;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9833;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9910;
    else if (1.6 <= fabs(mueta) ) return 0.9900;
  }
  else if(mupt >= 35)
  {
    if (fabs(mueta) < 0.8) return 0.9826;
    else if (0.8 <= fabs(mueta) && fabs(mueta) < 1.2) return 0.9841;
    else if (1.2 <= fabs(mueta) && fabs(mueta) < 1.6) return 0.9900;
    else if (1.6 <= fabs(mueta) ) return 0.9886;
  }

  assert(false);
}

}

ScaleFactors::ScaleFactors(const char* data_era)
{
	if(strcmp(data_era, "fall11") == 0 || strncmp(data_era, "fall11-", 7) == 0)
	{
		era_ = Era2011;
	}
	else if(strcmp(data_era, "summer12") == 0 || strncmp(data_era, "summer12-", 9) == 0)
	{
		era_ = Era2012;
	}
	else if(strcmp(data_era, "2011") == 0 || strcmp(data_era, "2012") == 0)
	{
		era_ = Data;
	}
	else
	{
		throw std::runtime_error("No such data era: " + std::string(data_era));
	}
}

float ScaleFactors::getDoubleMuTrigWeight(const char* trigName, float muPt1, float muEta1, float muPt2, float muEta2) const
{
	switch(era_)
	{
	case Era2011:
		if(strncmp(trigName, "HLT_Mu13_Mu8_v", 14) == 0)
			return Eff_HLT_Mu13_Mu8_2011_TPfit_RunAB_EtaEta_DATAoverMC(muEta1, muEta2);
		else if(strncmp(trigName, "HLT_Mu17_Mu8_v", 14) == 0)
			return Eff_HLT_Mu17_Mu8_2011_TPfit_RunAB_EtaEta_DATAoverMC(muEta1, muEta2);
		else if(strncmp(trigName, "HLT_DoubleMu7_v", 15) == 0)
			// We don't have dedicated DoubleMu7 scale factors:
			return Eff_HLT_Mu13_Mu8_2011_TPfit_RunAB_EtaEta_DATAoverMC(muEta1, muEta2);
		else
			assert(false);
		break;
	case Era2012:
		// This is derived as an OR of HLT_Mu17_Mu8 and HLT_Mu17_TkMu8:
		assert(strncmp(trigName, "HLT_Mu13_Mu8_v", 14) == 0 || strncmp(trigName, "HLT_Mu17_Mu8_v", 14) == 0);
		return Corr_Trg_Mu_2012_53X(muPt1, muEta1) * Corr_Trg_Mu_2012_53X(muPt2, muEta2);
	case Data:
		return 1.;
	default:
		assert(false);
		break;
	}
}

float ScaleFactors::getMuIdIsoWeight(float muPt, float muEta) const
{
	switch(era_)
	{
	case Era2011:
		{
			const double id2011A = PFTIGHT_nL8_2011A(muPt, muEta);
			const double id2011B = PFTIGHT_nL8_2011B(muPt, muEta);
			const double iso2011A = combRelPFISO12_2011A(muPt, muEta);
			const double iso2011B = combRelPFISO12_2011B(muPt, muEta);
			const double lumi2011A = 215.630 + 955.212 + 389.876 + 706.700;
			const double lumi2011B = 2714.453;

			const double id = (id2011A * lumi2011A + id2011B * lumi2011B) / (lumi2011A + lumi2011B);
			const double iso = (iso2011A * lumi2011A + iso2011B * lumi2011B) / (lumi2011A + lumi2011B);
			return id * iso;
		}
		break;
	case Era2012:
		return Corr_IdIso_Mu_2012_53X(muPt, muEta);
	case Data:
		return 1.;
	default:
		assert(false);
		break;
	}
}


