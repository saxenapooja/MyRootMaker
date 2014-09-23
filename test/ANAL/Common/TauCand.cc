
#include "AnalysisTool/DataForm.h"
#include "Common/TauCand.h"
#include "Common/TauWP.h"

static TLorentzVector correct(const TLorentzVector& p4, double corr)
{
	const double px = p4.Px() * corr;
	const double py = p4.Py() * corr;
	const double pz = p4.Pz() * corr;
	const double M = p4.M();
	return TLorentzVector(px, py, pz, TMath::Sqrt(M*M+px*px+py*py+pz*pz));
}

TauCand::TauCand(const TauWP& wp, const Tau& tau, double sf):
	trigObj(tau), Charge(tau.Charge()), NewEMFraction(tau.NewEMFraction()), NumTracks(tau.NumTracks()),
	Dxy(tau.Dxy()), Dz(tau.Dz()),
	IsoNeutralsPt(tau.IsoNeutralsPt()), IsoChargedPt(tau.IsoChargedPt()), IsoGammaPt(tau.IsoGammaPt()), Iso3Hits(tau.Iso3Hits()), TrackIso(0.0f),
	decayModeFinding(tau.TauDiscriminator("decayModeFinding")),
	byLooseCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr") == 1),
	byLooseCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr3Hits") == 1),
	byMediumCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr") == 1),
	byMediumCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr3Hits") == 1),
	byTightCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr") == 1),
	byTightCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr3Hits") == 1),
	againstElectronLoose(tau.TauDiscriminator("againstElectronLoose") == 1),
	againstElectronLooseMVA3(tau.TauDiscriminator("againstElectronLooseMVA3") == 1),
	againstElectronMedium(tau.TauDiscriminator("againstElectronMedium") == 1),
	againstElectronMediumMVA3(tau.TauDiscriminator("againstElectronMediumMVA3") == 1),
	againstElectronTight(tau.TauDiscriminator("againstElectronTight") == 1),
	againstElectronTightMVA3(tau.TauDiscriminator("againstElectronTightMVA3") == 1),
	againstMuonLoose(tau.TauDiscriminator("againstMuonLoose") == 1),
	againstMuonMedium(tau.TauDiscriminator("againstMuonMedium") == 1),
	againstMuonTight(tau.TauDiscriminator("againstMuonTight") == 1),
	corr(sf), p4(correct(tau, corr)),
	leadpfchargedhadrcandpt(tau.leadpfchargedhadrcandpt()),
	leadpfchargedhadrcandp4(tau.leadpfchargedhadrcandp4()),
	PFChargedHadrCands_size(tau.signalPFChargedHadrCands_size()),
	PFGammaCands_size(tau.signalPFGammaCands_size()),
	L1trigger_match(tau.TauL1trigger_match())
{
}


TauCand::TauCand(const Tau& tau, double sf):
	trigObj(tau), Charge(tau.Charge()), NewEMFraction(tau.NewEMFraction()), NumTracks(tau.NumTracks()),
	Dxy(tau.Dxy()), Dz(tau.Dz()),
	IsoNeutralsPt(tau.IsoNeutralsPt()), IsoChargedPt(tau.IsoChargedPt()), IsoGammaPt(tau.IsoGammaPt()), Iso3Hits(tau.Iso3Hits()), TrackIso(0.0f),
	decayModeFinding(tau.TauDiscriminator("decayModeFinding")),
	byLooseCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr") == 1),
	byLooseCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byLooseCombinedIsolationDeltaBetaCorr3Hits") == 1),
	byMediumCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr") == 1),
	byMediumCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byMediumCombinedIsolationDeltaBetaCorr3Hits") == 1),
	byTightCombinedIsolationDeltaBetaCorr(tau.TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr") == 1),
	byTightCombinedIsolationDeltaBetaCorr3Hits(tau.TauDiscriminator("byTightCombinedIsolationDeltaBetaCorr3Hits") == 1),
	againstElectronLoose(tau.TauDiscriminator("againstElectronLoose") == 1),
	againstElectronLooseMVA3(tau.TauDiscriminator("againstElectronLooseMVA3") == 1),
	againstElectronMedium(tau.TauDiscriminator("againstElectronMedium") == 1),
	againstElectronMediumMVA3(tau.TauDiscriminator("againstElectronMediumMVA3") == 1),
	againstElectronTight(tau.TauDiscriminator("againstElectronTight") == 1),
	againstElectronTightMVA3(tau.TauDiscriminator("againstElectronTightMVA3") == 1),
	againstMuonLoose(tau.TauDiscriminator("againstMuonLoose") == 1),
	againstMuonMedium(tau.TauDiscriminator("againstMuonMedium") == 1),
	againstMuonTight(tau.TauDiscriminator("againstMuonTight") == 1),
	corr(sf), p4(correct(tau, corr)),
	leadpfchargedhadrcandpt(tau.leadpfchargedhadrcandpt()),
	leadpfchargedhadrcandp4(tau.leadpfchargedhadrcandp4()),
	L1trigger_match(tau.TauL1trigger_match())
{
}
