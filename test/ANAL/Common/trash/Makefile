Common.a: Config.o PileupWeight.o FakeRateFunctions.o Variable.o MvaTree.o MvaEvaluator.o TauWP.o TauCand.o GenTools.o ElectronId.o TrackIsoReader.o LumiFilter.o Cutflow.o IrredMVA.o CombinedMVA.o METWeight.o LltVariable.o MuonId.o kNNEvaluator.o ScaleFactors.o
	ar cru Common.a Config.o PileupWeight.o FakeRateFunctions.o Variable.o MvaTree.o MvaEvaluator.o TauWP.o TauCand.o GenTools.o ElectronId.o TrackIsoReader.o LumiFilter.o Cutflow.o IrredMVA.o CombinedMVA.o METWeight.o LltVariable.o MuonId.o kNNEvaluator.o ScaleFactors.o

Config.o: Config.cc Config.h
	g++ Config.cc -c -o Config.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

PileupWeight.o: PileupWeight.cc PileupWeight.h
	g++ PileupWeight.cc -c -o PileupWeight.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

FakeRateFunctions.o: FakeRateFunctions.cc FakeRateFunctions.h
	g++ FakeRateFunctions.cc -c -o FakeRateFunctions.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

Variable.o: Variable.cc Variable.h MvaEvaluator.h MariaFunctions.h MariaFunctions_ZMedium.h MariaFunctions_ZTight.h MariaFunctions_WMedium.h MariaFunctions_WTight.h
	g++ Variable.cc -c -o Variable.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

MvaTree.o: MvaTree.cc MvaTree.h
	g++ MvaTree.cc -c -o MvaTree.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

MvaEvaluator.o: MvaEvaluator.cc MvaEvaluator.h MvaTree.h
	g++ MvaEvaluator.cc -c -o MvaEvaluator.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

TauWP.o: TauWP.cc TauWP.h
	g++ TauWP.cc -c -o TauWP.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

TauCand.o: TauCand.cc TauCand.h
	g++ TauCand.cc -c -o TauCand.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

GenTools.o: GenTools.cc GenTools.h
	g++ GenTools.cc -c -o GenTools.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

ElectronId.o: ElectronId.cc ElectronId.h
	g++ ElectronId.cc -c -o ElectronId.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

TrackIsoReader.o: TrackIsoReader.cc TrackIsoReader.h
	g++ TrackIsoReader.cc -c -o TrackIsoReader.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

LumiFilter.o: LumiFilter.cc LumiFilter.h
	g++ LumiFilter.cc -c -o LumiFilter.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

Cutflow.o: Cutflow.cc Cutflow.h
	g++ Cutflow.cc -c -o Cutflow.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

IrredMVA.o: IrredMVA.cc IrredMVA.h
	g++ IrredMVA.cc -c -o IrredMVA.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

CombinedMVA.o: CombinedMVA.cc CombinedMVA.h
	g++ CombinedMVA.cc -c -o CombinedMVA.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

METWeight.o: METWeight.cc METWeight.h
	g++ METWeight.cc -c -o METWeight.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

LltVariable.o: LltVariable.cc LltVariable.h
	g++ LltVariable.cc -c -o LltVariable.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

MuonId.o: MuonId.cc MuonId.h
	g++ MuonId.cc -c -o MuonId.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

kNNEvaluator.o: kNNEvaluator.cc kNNEvaluator.h
	g++ kNNEvaluator.cc -c -o kNNEvaluator.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3

ScaleFactors.o: ScaleFactors.cc ScaleFactors.h
	g++ ScaleFactors.cc -c -o ScaleFactors.o -I.. -I /nfs/dust/cms/user/aburgm/tautau-fullyhadronic/taucode_new/include `root-config --cflags` -O3
