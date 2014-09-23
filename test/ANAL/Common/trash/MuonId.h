
#ifndef _WHANALYSIS_MUONID_H_
#define _WHANALYSIS_MUONID_H_

#include "AnalysisTool/DataForm.h"

namespace MuonId
{

bool loose(const Muon& mu);
bool pfTight(const Muon& mu);
bool matchedJet(const Analyse& ana, const Muon& muon, Jet& j);
double jetPt(const Analyse& ana, const Muon& muon);
bool bjetVeto(const Analyse& ana, const Muon& muon);

}

#endif // _WHANALYSIS_MUONID_H_

