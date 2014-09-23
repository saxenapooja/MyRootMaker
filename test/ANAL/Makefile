main: main.cc ../Common/Common.a
	g++ main.cc -o main -lTMVA -lXMLIO -lMLP -lTreePlayer -lGenVector `root-config --cflags --libs` -I.. -I /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/include -L /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/lib ../Common/Common.a -lAnalysisTool -O3

main_skim: main_skim.cc	 ../Common/Common.a
	g++ main_skim.cc -o main_skim `root-config --cflags --libs` -I.. -I /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/include -L /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/lib ../Common/Common.a -lAnalysisTool -lsqlite3 -O3
