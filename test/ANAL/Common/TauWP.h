#ifndef _WHANALYSIS_TAUWP_H_
#define _WHANALYSIS_TAUWP_H_

#include "Common/Config.h"
#include "Common/TauCand.h"

class TauWP
{
public:
  enum Isolation {
    ISO_LOOSE,
    ISO_LOOSE_3HITS,
    ISO_MEDIUM,
    ISO_MEDIUM_3HITS,
    ISO_TIGHT,
    ISO_TIGHT_3HITS
  };

  enum AntiE {
    ANTIE_LOOSE,
    ANTIE_LOOSE_MVA3,
    ANTIE_MEDIUM,
    ANTIE_MEDIUM_MVA3,
    ANTIE_TIGHT,
    ANTIE_TIGHT_MVA3
  };

  enum AntiMu {
    ANTIMU_LOOSE,
    ANTIMU_MEDIUM,
    ANTIMU_TIGHT
  };

  TauWP(const Config& cfg, 
	const char* prefix, 
	const char* defaultIsolation, 
	const char* defaultAntiE, 
	const char* defaultAntiMu);


  TauWP(const Config& cfg, 
	const char* prefix, 
	const char* defaultIsolation);

  const Isolation ISOLATION;
  const AntiE ANTIE;
  const AntiMu ANTIMU;

  const bool TauCand::*ISOLATION_FLAG;
  const bool TauCand::*ANTIE_FLAG;
  const bool TauCand::*ANTIMU_FLAG;

  const char* isolation_hr_name() const;
  const char* antie_hr_name() const;
  const char* antimu_hr_name() const;

  const char* isolation_frfunc_name() const;
  const char* antie_frfunc_name() const;
};

#endif // _WHANALYSIS_TAUWP_H_
