#include "Common/TrackIsoReader.h"

#include <Math/VectorUtil.h>

#include <stdexcept>

TrackIsoReader::TrackIsoReader(const char* filename)
{
	if(sqlite3_open(filename, &db) != SQLITE_OK)
		throw std::runtime_error(sqlite3_errmsg(db));

	const char query[] = "SELECT pt, eta, phi, iso FROM TrackIso WHERE run=? AND lumi=? AND event=?";
	if(sqlite3_prepare(db, query, sizeof(query), &select_stmt, NULL) != SQLITE_OK)
	{
		std::string msg = sqlite3_errmsg(db);
		sqlite3_close(db);
		throw std::runtime_error(msg);
	}
}

TrackIsoReader::~TrackIsoReader()
{
	sqlite3_finalize(select_stmt);
	sqlite3_close(db);
}

std::vector<double> TrackIsoReader::readTrackIso(unsigned int run, unsigned int lumi, unsigned int event, const std::vector<Tau>& taus)
{
	sqlite3_bind_int(select_stmt, 1, run);
	sqlite3_bind_int(select_stmt, 2, lumi);
	sqlite3_bind_int(select_stmt, 3, event);

	std::vector<double> trackIso(taus.size(), -1.);
	int step_result;

	while((step_result = sqlite3_step(select_stmt)) == SQLITE_ROW)
	{
		const double pt = sqlite3_column_double(select_stmt, 0);
		const double eta = sqlite3_column_double(select_stmt, 1);
		const double phi = sqlite3_column_double(select_stmt, 2);
		const double trackiso = sqlite3_column_double(select_stmt, 3);

		TLorentzVector vec;
		vec.SetPtEtaPhiM(pt, eta, phi, 1.7);

		const Tau* bestTau = NULL;
		for(unsigned int i = 0; i < taus.size(); ++i)
		{
			const double deltaR = ROOT::Math::VectorUtil::DeltaR(vec, taus[i]);
			if(deltaR < 0.3)
			{
				if(!bestTau || deltaR < ROOT::Math::VectorUtil::DeltaR(vec, *bestTau))
				{
					trackIso[i] = trackiso;
					bestTau = &taus[i];
				}
			}
		}
	}

	if(step_result == SQLITE_ERROR)
		throw std::runtime_error(sqlite3_errmsg(db));
	if(sqlite3_reset(select_stmt) != SQLITE_OK)
		throw std::runtime_error(sqlite3_errmsg(db));

	return trackIso;
}
