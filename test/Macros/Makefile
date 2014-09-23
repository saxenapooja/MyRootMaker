BOOST_LIBS = -lboost_regex -lboost_thread -lboost_filesystem
ROOT_LIBS = $(shell root-config --libs) -lMinuit -lMathCore -lrt
libAnalysisTool_la_LIBADD = $(ROOT_LIBS) $(BOOST_LIBS)
merger_LDADD = $(ROOT_LIBS) $(BOOST_LIBS) -lAnalysisTool

VertexAnalysis: VertexAnalysis.C
	g++ VertexAnalysis.C -o script -lTMVA -lXMLIO -lMLP -lTreePlayer $(BOOST_LIBS) `root-config --cflags --libs`
# -I.. -I /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/include -L /nfs/dust/cms/user/pooja/scratch/plot-macro/lib_AnalysisTool/lib ../Common/Common.a -lAnalysisTool -O3

