
#ifndef _WHANALYSIS_MVATREE_H_
#define _WHANALYSIS_MVATREE_H_

#include "Common/FakeRateFunctions.h"
#include "Common/TauCand.h"

#include "AnalysisTool/Analyse.h"

#include <TLorentzVector.h>
#include <TTree.h>

class MvaTree
{
public:
  MvaTree(const FakeRateFunctions& fr, const char* name, const char* desc);

  class Variables
  {
  public:
    Variables(const Muon& muon, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30, int regions, bool case2); //bool lepton_fake, bool OStauJetFake, bool SStauJetFake, bool wjets_veto_mt, bool OStauTrackFake, bool SStauTrackFake);
    Variables(const Electron& electron, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30, int regions, bool case2); //bool lepton_fake, bool OStauJetFake, bool SStauJetFake, bool wjets_veto_mt, bool OStauTrackFake, bool SStauTrackFake);

    struct MvaVars
    {
      float leadPt;
      float subPt;
      float OSPt;
      float OSEta;
      float SSPt;
      float SSEta;
      //float Lt;
      float Met;
      float deltaRDiTau;
      float deltaPhiDiTau;
      float deltaEtaDiTau;
      //float etaDiTau;
      //float etaLepTauTau;
      //float deltaRLepOSTau;
      //float deltaRLepSSTau;
      float ptRatio;
      float angle;

      //unsigned int nJets30;
      //float etaLeadTau;
      //float etaSubTau;
      float transMassLepMET;
      //float transMassOSTauMET;
      //float transMassSSTauMET;
      //float deltaPhiLepMET;

      //float dxyOSTau;
      //float dxySSTau;
      //float dzOSTau;
      //float dzSSTau;
    };

/*    const bool leptonFake;
    const bool OSTauFake;
    const bool SSTauFake;
    const bool wJetsVetoMT;
    const bool OSTauTrackFake;
    const bool SSTauTrackFake;*/
/*    const int lepCharge;
    const int leadTauCharge;
    const int subTauCharge;*/
    int regions;
    bool case2;
    MvaVars mvaVars;
  private:
    void init(const TLorentzVector& lepton, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30);
  };

  void fill(const Variables& vars, double weight);

private:
  const FakeRateFunctions& fakeRateFuncs;

  struct TreeVars
  {
    Variables::MvaVars mvaVars;
    float weight;
    float Wweight;
    float Zweight;
  };

  TTree* mvaTree;
  TTree* mvaTreeTau12Fake;
  TTree* mvaTreeTauSSFake;

  TreeVars mvaVars;
  TreeVars mvaVarsTau12Fake;
  TreeVars mvaVarsTauSSFake;
};

#endif // _WHANALYSIS_MVATREE_H_
