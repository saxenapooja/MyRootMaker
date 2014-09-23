/*This is an example how to use the AC1B Root-Tree structure. We provide this interface to have an easy access to 
 * data stored in the new Root-Trees.*/

//Include Analyse.h to make the interface available. 
#include "AnalysisTool/Analyse.h"
#include "Common/Config.h"
#include "Common/LumiFilter.h"
#include <TH1D.h>
#include <iostream>
#include <vector>

using namespace std;

//You should derive your own class from Analyse.
class MyAnalysis : public Analyse
{
public:
  enum GenFilter { 
    NONE, 
    W_MU_H_TAUTAU,
    W_MU_H_TAUTAU_SKIMMED, 
    W_TAU_H_MUTAU, 
    Z_TAUTAU_H_MUTAU, 
    Z_MUTAU_H_TAUTAU, 
    WH_LEPDECAY, 
    WH_HADDECAY, 
    ZH, 
    TTH, 
    EWK_EE, 
    EWK_MM, 
    EWK_TT 
  };
  
private:
  std::set<std::string> processedFiles;
  const bool LOWPT_SKIM;
  const GenFilter GEN_FILTER;

public:
  MyAnalysis(const Config& cfg, const std::string& database);
  virtual ~MyAnalysis();
  virtual Int_t AnalyseEvent();
  
  double nEv_Skim;
  double nEV_lead;
  double nEV_sub;
  double nEV_muons;

  unsigned int nW;
  unsigned int nZ;
  unsigned int nT;
  unsigned int nO;

  TH1D* nEvents;
};//class description over


MyAnalysis::GenFilter string_to_filter(const std::string& str)
{
  std::cout << "--- Filter --> " << str << std::endl;
  if(str == "none") return MyAnalysis::NONE;
  else if(str == "w_lep_h_tautau")         return MyAnalysis::W_MU_H_TAUTAU;
  else if(str == "w_lep_h_tautau_skimmed") return MyAnalysis::W_MU_H_TAUTAU_SKIMMED;
  else if(str == "w_tau_h_leptau")         return MyAnalysis::W_TAU_H_MUTAU;
  else if(str == "z_tautau_h_leptau")      return MyAnalysis::Z_TAUTAU_H_MUTAU;
  else if(str == "z_leptau_h_tautau")      return MyAnalysis::Z_MUTAU_H_TAUTAU;
  else if(str == "wh_lepdecay")            return MyAnalysis::WH_LEPDECAY;
  else if(str == "wh_haddecay")            return MyAnalysis::WH_HADDECAY;
  else if(str == "zh")                     return MyAnalysis::ZH;
  else if(str == "tth")                    return MyAnalysis::TTH;
  else if(str == "ewk_ee")                 return MyAnalysis::EWK_EE;
  else if(str == "ewk_mm")                 return MyAnalysis::EWK_MM;
  else if(str == "ewk_tt")                 return MyAnalysis::EWK_TT;
  throw std::runtime_error("Invalid Filter: " + str);
}

//Constructor:
MyAnalysis::MyAnalysis(const Config& cfg,  const std::string& database): LOWPT_SKIM(cfg.get<bool>("lowpt_skim")),
									 GEN_FILTER(string_to_filter(cfg.get<std::string>("gen_filter", "none")))
{
  //cout<<"(lowpt_skim, gen_filter) : ("<< LOWPT_SKIM << ", "<< GEN_FILTER <<")"<< endl;
  LoadTrigger();
  LoadBeamSpot();
  LoadPrimVertices();
  LoadTaus();
  LoadMuons();
  LoadAK5PFJets();
  LoadElectrons();
  LoadMET();
  LoadGenParticles();
  LoadAllGenParticles();
  LoadGenInfo();
  
  PrepareSkimming("skimmed.root");
  nEvents = new TH1D("nEvents", "Number of events", 1, -1, 1);
  nEv_Skim = nEV_lead = nEV_sub = nEV_muons = 0.0;
  nW = nZ = nT = nO = 0;
}


//Destructor:
MyAnalysis::~MyAnalysis(){
  
  cout<<setw(20)<<"Selection"         <<setw(15)<<setprecision(5)<<" "            <<setw(15)<<"Efficiency"<<endl;
  cout<<setw(20)<<"nEv_Skim:"         <<setw(15)<<nEv_Skim      <<setw(15)<<1<<endl;                        //  
  cout<<"__________________________________________________"<<endl;
  //cout<<setw(20)<<"Muons"             <<setw(15)<<nEV_muons     <<setw(15)<<nEV_muons/nEv_Skim<<endl;               //  
  cout<<setw(20)<<"Tau1"              <<setw(15)<<nEV_sub       <<setw(15)<<nEV_sub/nEv_Skim<<endl;
  cout<<setw(20)<<"Tau2"              <<setw(15)<<nEV_lead      <<setw(15)<<nEV_lead/nEV_sub<<endl;
  cout<<"__________________________________________________"<<endl;
  cout<<setw(20)<<"Total"             <<setw(15)<< " "          <<setw(15)<<nEV_lead/nEv_Skim<<endl;
  cout << endl;
  
  cout << "W: " << nW << std::endl;
  cout << "Z: " << nZ << std::endl;
  cout << "T: " << nT << std::endl;
  cout << "O: " << nO << std::endl;
  
  nEvents->Write();
}//MyAnalysis::~MyAnalysis()



Int_t MyAnalysis::AnalyseEvent() {
  
  if(IsLumiAvailable() != 2){
    cout << "NoLumiInfo: " << Run() << " "<<LumiBlock()<<" "<<Number() << endl;
    return(1);
  }
  
  // Lumi Filter
  if(!LumiFilter::pass(Run(), LumiBlock()))  return 1;
  
  // If we don't apply a filter on the WH dataset, write out the number of original events
  // If we do apply such a filter, we store the number of events passing it further below.
  const bool nEventsFromSkim = (GEN_FILTER != W_MU_H_TAUTAU         && 
				GEN_FILTER != W_MU_H_TAUTAU_SKIMMED && 
				GEN_FILTER != W_TAU_H_MUTAU         && 
				GEN_FILTER != Z_TAUTAU_H_MUTAU      && 
				GEN_FILTER != Z_MUTAU_H_TAUTAU      && 
				GEN_FILTER != WH_LEPDECAY           && 
				GEN_FILTER != WH_HADDECAY           && 
				GEN_FILTER != ZH                    && 
				GEN_FILTER != TTH );
  //  cout<<" nEventsFromSkim "<< nEventsFromSkim << endl;
  
  if(nEventsFromSkim)
    {
      if(processedFiles.find(GetCurrentFileName()) == processedFiles.end())
	{
	  nEvents->SetBinContent(1, nEvents->GetBinContent(1) + GetOrigEvents());
	  processedFiles.insert(GetCurrentFileName());
	}
    } // if(nEventsFromSkim)
  

  // Apply generator filter
  if(GEN_FILTER != NONE)
    {
      unsigned int haveW = 0;
      unsigned int haveZ = 0;
      unsigned int haveT = 0;
      unsigned int haveH = 0;
      unsigned int haveEl = 0;
      unsigned int haveMu = 0;
      unsigned int haveTau = 0;
      unsigned int haveElFromW = 0;
      unsigned int haveMuFromW = 0;
      unsigned int haveTauFromW = 0;
      unsigned int haveMuFromZ = 0;
      unsigned int haveTauFromZ = 0;
      unsigned int haveElFromH = 0;
      unsigned int haveMuFromH = 0;
      unsigned int haveTauFromH = 0;

      if(NumAllGenParticles() > 0)
	{
	  for(unsigned int i = 0; i < NumAllGenParticles(); ++i) // returning filled tree variable: 'genallparticles_count'
	    {
	      const GenParticle part = AllGenParticles(i);
	      if(TMath::Abs(part.PDGId()) == 24 && part.Status() == 3 && part.NumMothers() == 2) ++haveW;
	      if(TMath::Abs(part.PDGId()) == 23 && part.Status() == 3) ++haveZ;
	      if(TMath::Abs(part.PDGId()) == 6  && part.Status() == 3) ++haveT;
	      if(TMath::Abs(part.PDGId()) == 25 && part.Status() == 3) ++haveH;

	      if(TMath::Abs(part.PDGId()) == 11 && part.Status() == 3) ++haveEl;
	      if(TMath::Abs(part.PDGId()) == 11 && part.Status() == 3)
		{
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 24 && part.GetMother(0).NumMothers() == 2)
		    ++haveElFromW;
		}

	      if(TMath::Abs(part.PDGId()) == 13 && part.Status() == 3) ++haveMu;
	      if(TMath::Abs(part.PDGId()) == 13 && part.Status() == 3)
		{
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 24 && part.GetMother(0).NumMothers() == 2)
		    ++haveMuFromW;
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 23)
		    ++haveMuFromZ;
		}

	      if(TMath::Abs(part.PDGId()) == 15 && part.Status() == 3) ++haveTau;
	      if(TMath::Abs(part.PDGId()) == 15 && part.Status() == 3)
		{
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 24 && part.GetMother(0).NumMothers() == 2)
		    ++haveTauFromW;
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 23)
		    ++haveTauFromZ;
		  if(part.NumMothers() == 1 && TMath::Abs(part.GetMother(0).PDGId()) == 25)
		    {
		      GenParticle particle = part;
		      while(particle.NumDaughters() == 1 && TMath::Abs(particle.GetDaughter(0).PDGId()) == 15)
			particle = particle.GetDaughter(0);
		      
		      bool haveEl = false, haveMu = false;
		      for(unsigned int j = 0; j < particle.NumDaughters(); ++j)
			{
			  if(TMath::Abs(particle.GetDaughter(j).PDGId()) == 11) haveEl = true;
			  if(TMath::Abs(particle.GetDaughter(j).PDGId()) == 13) haveMu = true;
			}
		      
		      if(haveEl)              ++haveElFromH;
		      if(haveMu)              ++haveMuFromH;
		      if(!haveEl && !haveMu)  ++haveTauFromH;
		    }
		}
	    }//for(unsigned int i = 0; i < NumAllGenParticles(); ++i)
	}//if(NumAllGenParticles() > 0)
      else if(NumGenParticles() > 0) //return filled variable 'genparticles_count'
	{
	  for(unsigned int i = 0; i < NumGenParticles(); ++i)
	    {
	      const GenLightParticle part = GenParticles(i);
	      if(TMath::Abs(part.PDGId()) == 24 && part.Status() == 3) ++haveW;
	      if(TMath::Abs(part.PDGId()) == 23 && part.Status() == 3) ++haveZ;
	      if(TMath::Abs(part.PDGId()) == 6  && part.Status() == 3) ++haveT;
	      if(TMath::Abs(part.PDGId()) == 25 && part.Status() == 3) ++haveH;
	      if(TMath::Abs(part.PDGId()) == 11 && part.Status() == 3) ++haveEl;
	      if(TMath::Abs(part.PDGId()) == 13 && part.Status() == 3) ++haveMu;
	      if(TMath::Abs(part.PDGId()) == 15 && part.Status() == 3) ++haveTau;
	      if(TMath::Abs(part.PDGId()) == 11 && part.Status() == 3 && part.MoreInfo(GenLightParticle::FromW)) ++haveElFromW;
	      if(TMath::Abs(part.PDGId()) == 13 && part.Status() == 3 && part.MoreInfo(GenLightParticle::FromW)) ++haveMuFromW;
	      if(TMath::Abs(part.PDGId()) == 15 && part.Status() == 3 && part.MoreInfo(GenLightParticle::FromW)) ++haveTauFromW;
	      //if(TMath::Abs(part.PDGId()) == 13 && part.Status() == 3 && part.MoreInfo(GenLightParticle::FromZ)) ++haveMuFromZ;
	      //if(TMath::Abs(part.PDGId()) == 15 && part.Status() == 3 && part.MoreInfo(GenLightParticle::FromZ)) ++haveTauFromZ;
	      haveMuFromZ = 0; // TODO
	      haveMuFromZ = 0; // TODO
	      haveTauFromH = 2; // TODO
	      haveElFromH = 0; // TODO
	      haveMuFromH = 0; // TODO
	    }// for(unsigned int i = 0; i < NumGenParticles(); ++i)
	}//else if(NumGenParticles() > 0) 
      else
	{
	  std::cerr << "No Generator information available!" << std::endl;
	  assert(false);
	  return 1;
	}
    
      if(haveW == 1) ++nW;
      else if(haveZ == 1) ++nZ;
      else if(haveT == 2) ++nT;
      else  std::cout << "Other: T=" << haveT << ", Z=" << haveZ << ", W=" << haveW << ", H=" << haveH << std::endl; 
      
      switch(GEN_FILTER)
	{
	case W_MU_H_TAUTAU:

	case W_MU_H_TAUTAU_SKIMMED:
	  if(haveT != 0 || haveZ != 0 || haveW != 1 || haveH != 1 || haveMuFromW != 1 || haveTauFromW != 0 || haveMuFromH != 0 || haveTauFromH != 2) return 1;
	  break;

	case W_TAU_H_MUTAU:
	  if(haveT != 0 || haveZ != 0 || haveW != 1 || haveH != 1 || haveMuFromW != 0 || haveTauFromW != 1 || haveMuFromH != 1 || haveTauFromH != 1) {
	    cout<<" @@@@@@ Filter : W_TAU_H_MUTA" << endl; 
	    return 1;
	  }
	  break;

	case Z_TAUTAU_H_MUTAU:
	  if(haveT != 0 || haveZ != 1 || haveW != 0 || haveH != 1 || haveMuFromZ != 0 || haveTauFromZ != 2 || haveMuFromH != 1 || haveTauFromH != 1) return 1;
	  break;

	case Z_MUTAU_H_TAUTAU:
	  // TODO: Doesn't work -- need to trace tau decay products!
	  if(haveT != 0 || haveZ != 1 || haveW != 0 || haveH != 1 || haveMuFromZ != 1 || haveTauFromZ != 1 || haveMuFromH != 0 || haveTauFromH != 2) return 1;
	  break;

	case WH_LEPDECAY:
	  if(haveT != 0 || haveZ != 0 || haveW != 1 || haveH != 1 || (haveElFromW + haveMuFromW + haveTauFromW) != 1) return 1;
	  break;

	case WH_HADDECAY:
	  if(haveT != 0 || haveZ != 0 || haveW != 1 || haveH != 1 || (haveElFromW + haveMuFromW + haveTauFromW) != 0) return 1;
	  break;

	case ZH:
	  if(haveT != 0 || haveZ != 1 || haveW != 0 || haveH != 1) return 1;
	  break;

	case TTH:
	  if(haveT != 2 || haveZ != 0 || haveW != 0 || haveH != 1) return 1;
	  break;

	case EWK_EE:
	  if(haveEl == 0) return 1;
	  break;

	case EWK_MM:
	  if(haveMu == 0) return 1;
	  break;

	case EWK_TT:
	  if(haveTau == 0) return 1;
	  break;
	default:
	  assert(false);
	  break;
	}//switch(GEN_FILTER)
    }//if(GEN_FILTER != NONE)

  nEv_Skim++;
  if(!nEventsFromSkim)   nEvents->Fill(0.0);
  
  // loop over the taus
  std::vector< Tau > taus1;
  std::vector< Tau > taus2;

  //first tau
  for( unsigned int i = 0; i < NumTaus(); i++){

    if( (LOWPT_SKIM
	 && Taus(i).Pt() > 20.
	 && Abs(Taus(i).Eta()) < 2.3
	 && Taus(i).TauDiscriminator("decayModeFinding")==1
	 && (Taus(i).TauDiscriminator("againstElectronLoose") == 1 || Taus(i).TauDiscriminator("againstElectronMVA") || Taus(i).TauDiscriminator("againstElectronVLooseMVA2") == 1 || Taus(i).TauDiscriminator("againstElectronLooseMVA3") == 1)
	 && (Taus(i).TauDiscriminator("againstMuonLoose") == 1 || Taus(i).TauDiscriminator("againstMuonLoose2") == 1)
	 ) || (!LOWPT_SKIM
	       && Taus(i).Pt() > 30.
	       && Abs(Taus(i).Eta()) < 2.3
	       && Taus(i).TauDiscriminator("decayModeFinding")==1
	       ))
      {
	taus1.push_back(Taus(i));

	//second tau
	for( unsigned int j = 0; j < NumTaus(); j++) if(i!=j) {
	  if( (LOWPT_SKIM
               && Taus(j).Pt() > 25.
               && Abs(Taus(j).Eta()) < 2.3
               && Taus(j).TauDiscriminator("decayModeFinding")==1
               && (Taus(j).TauDiscriminator("againstElectronLoose") == 1 || Taus(j).TauDiscriminator("againstElectronMVA") || Taus(j).TauDiscriminator("againstElectronVLooseMVA2") == 1 || Taus(j).TauDiscriminator("againstElectronLooseMVA3") == 1)
               && (Taus(j).TauDiscriminator("againstMuonLoose") == 1 || Taus(j).TauDiscriminator("againstMuonLoose2") == 1)
	       ) || (!LOWPT_SKIM
		     && Taus(j).Pt() > 40.
		     && Abs(Taus(j).Eta()) < 2.3
		     && Taus(j).TauDiscriminator("decayModeFinding")==1
		     ))
	    {
	      taus2.push_back(Taus(j));
	      i = NumTaus();
	      break;
	    }
	}// for( unsigned int j = 0; j < NumTaus(); j++) if(i!=j)
      }
  }// for( unsigned int i = 0; i < NumTaus(); i++)
  
  if(false)
    {
      std::cout << "Run=" << Run() << ", Lumi=" << LumiBlock() << ", Event=" << Number() << std::endl;
      std::cout << "LeadTau: " << (taus2.size()>=1) << std::endl;
      std::cout << "SubTau: " << (taus1.size()>=1) << std::endl;
    }
  
  const bool doFilter = GEN_FILTER != W_MU_H_TAUTAU    && 
                        GEN_FILTER != W_TAU_H_MUTAU    && 
                        GEN_FILTER != Z_TAUTAU_H_MUTAU &&  
                        GEN_FILTER != Z_MUTAU_H_TAUTAU;

  if(doFilter && taus1.size() < 1 )    return(1);
  nEV_sub++;
  if(doFilter && taus2.size() < 1 )    return(1);
  nEV_lead++;
  
  SkimEvent();
  return(1);
}

//main 
int main(int argc, char* argv[])
{
  //initialiazation
  TTree::SetMaxTreeSize(500000000);
  string namebuf;
  string filename;
  std::auto_ptr<Config> config;
  std::vector<std::string> lumiFiles;
  std::vector<std::string> files;
  std::string database;

  for(int i = 1; i < argc; ++i)
  {
    namebuf = argv[i];

    UInt_t slashpos = namebuf.find_last_of("/");
    if(slashpos == namebuf.size())
      filename = namebuf;
    else
      filename = namebuf.substr(slashpos+1);

    if(filename.find("LUMI_") == 0)
      lumiFiles.push_back(namebuf);
    else if(filename.find(".cfg") != std::string::npos)
      if(config.get() == NULL)
        config.reset(new Config(namebuf.c_str()));
      else
        config->merge(Config(namebuf.c_str()));
    else if(filename.find(".db") != std::string::npos)
      database = namebuf;
    else
      files.push_back(namebuf);
  }
  
  if(!config.get())  throw std::runtime_error("No configuration file provided!");
  
  TH1::SetDefaultSumw2(true);
  MyAnalysis ana(*config, database);
  
  for(std::vector<std::string>::const_iterator iter = lumiFiles.begin(); iter != lumiFiles.end(); ++iter)
    ana.AddLumiFile(*iter);
  for(std::vector<std::string>::const_iterator iter = files.begin(); iter != files.end(); ++iter)
    ana.AddFile(*iter);
  
  // Loop will start to run the analysis on the specified range or on all events if no range is given. 
  ana.SetPrintInfo(10000);
  ana.EnableDuplicateCheck();
  ana.PrintLumiOfRuns();
  
  ana.Loop(); //0,10); 
  ana.GetLumi(1);
}
