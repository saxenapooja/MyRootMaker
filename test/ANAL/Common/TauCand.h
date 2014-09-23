#ifndef _WHANALYSIS_TAUCAND_H_
#define _WHANALYSIS_TAUCAND_H_

#include "AnalysisTool/DataForm.h"
#include "Common/Config.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class TauCand
{
public:
  TauCand(const class TauWP& wp, const Tau& tau, double sf);
  TauCand(const Tau& tau, double sf);

  float Pt()  const { return p4.Pt();  }
  float Eta() const { return p4.Eta(); }
  float Phi() const { return p4.Phi(); }
  float E()   const { return p4.E();   }

  TriggerObject trigObj;

  // These are not const so that they can be used with a std::vector<>.
  // You should only use const TauCands, though.
  int Charge;
  float NewEMFraction;
  unsigned int NumTracks;
  float Dxy;
  float Dz;

  int PFChargedHadrCands_size;
  int PFGammaCands_size;

  float IsoNeutralsPt;
  float IsoChargedPt;
  float IsoGammaPt;
  float Iso3Hits;
  float TrackIso;
  float leadpfchargedhadrcandpt;
  LorentzVector leadpfchargedhadrcandp4;

  bool L1trigger_match;
  bool decayModeFinding;
  bool byLooseCombinedIsolationDeltaBetaCorr;
  bool byLooseCombinedIsolationDeltaBetaCorr3Hits;
  bool byMediumCombinedIsolationDeltaBetaCorr;
  bool byMediumCombinedIsolationDeltaBetaCorr3Hits;
  bool byTightCombinedIsolationDeltaBetaCorr;
  bool byTightCombinedIsolationDeltaBetaCorr3Hits;

  bool againstElectronLoose;
  bool againstElectronLooseMVA3;
  bool againstElectronMedium;
  bool againstElectronMediumMVA3;
  bool againstElectronTight;
  bool againstElectronTightMVA3;

  bool againstMuonLoose;
  bool againstMuonMedium;
  bool againstMuonTight;

  double corr;
  TLorentzVector p4;

};

#endif // _WHANALYSIS_TAUCAND_H_
