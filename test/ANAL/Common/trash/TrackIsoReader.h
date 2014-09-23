
#ifndef _WHANALYSIS_TRACKISOREADER_H_
#define _WHANALYSIS_TRACKISOREADER_H_

#include "AnalysisTool/DataForm.h"

#include <sqlite3.h>

class TrackIsoReader
{
public:
	TrackIsoReader(const char* filename);
	~TrackIsoReader();

	std::vector<double> readTrackIso(unsigned int run, unsigned int lumi, unsigned int event, const std::vector<Tau>& taus);
private:
	sqlite3* db;
	sqlite3_stmt* select_stmt;
};

#endif // _WHANALYSIS_TRACKISOREADER_H_

