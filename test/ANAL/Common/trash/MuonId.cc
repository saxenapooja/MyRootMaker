#include "Common/MuonId.h"
#include "AnalysisTool/Analyse.h"

#include <Math/VectorUtil.h>

namespace MuonId
{

bool loose(const Muon& muon)
{
  if(!muon.IsGlobal() && !muon.IsTracker()) return false;
  if(muon.InnerTrack().NPixelHits() < 1) return false;
  return true;
}

bool pfTight(const Muon& muon)
{
  if(!loose(muon)) return false;
  if(!muon.IsPF()) return false;
  if(muon.NumStations() < 2) return false;
  if(muon.NumChamberHits() < 1) return false;
  if(muon.InnerTrack().NPixelLayers() + muon.InnerTrack().NStripLayers() <= 5) return false;
  if(muon.Chi2OverNdof() >= 10.0) return false;
  return true;
}

bool matchedJet(const Analyse& ana, const Muon& muon, Jet& j, double minPt)
{
	Jet closestJet;
	double closestDR = -1.;
	for(unsigned int i = 0; i < ana.NumAK5PFJets(); ++i)
	{
		const Jet jet = ana.AK5PFJets(i);
		if(jet.Pt() > minPt)
		{
			const double dR = ROOT::Math::VectorUtil::DeltaR(muon, jet);
			if(dR < 0.5 && (closestDR < 0. || dR < closestDR))
			{
				closestJet = jet;
				closestDR = dR;
			}
		}
	}

	if(closestDR < 0.) return false;

	j = closestJet;
	return true;
}

double jetPt(const Analyse& ana, const Muon& muon)
{
	Jet closestJet;
	if(!matchedJet(ana, muon, closestJet, 10.)) return muon.Pt();
	return std::max(closestJet.Pt(), muon.Pt());
}

bool bjetVeto(const Analyse& ana, const Muon& muon)
{
	Jet closestJet;
	if(!matchedJet(ana, muon, closestJet, 12.)) return true;
	return closestJet.combinedSecondaryVertexBJetTags() < 0.8;
}

} // namespace MuonId
