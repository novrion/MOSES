#pragma once

namespace moses {

class Moses;
class Government;

class Labour {
	friend class Moses;
	friend class Firm;
	friend class Government;
	friend class Household;
	friend class ui;
private:
	    
	// References
	const Moses& moses;


	/* ===== Parameters ===== */
    
	const double entry;  // (csv)
    const double ret;    // (csv)
	const double theta;  // (csv)
    const double ksi;    // (csv)
    const double gamma;  // (csv)
	const double rlu;    // (csv)
    const double skrepa; // (csv)


	/* ===== Variables ===== */
    
	double lf;
    double lu;     // (csv)
   	double ru;     // (csv)
    double qchru;
    double qdwind; // (csv)
   	
public:
	Labour() = default;
	~Labour() = default;

    Labour(const Moses& moses_ref,
			double lu, double entry, double ret,
			double qdwind, double theta, double ksi,
			double gamma, double ru, double rlu,
			double skrepa);

private:
	void print_initialisation();

	void luupdate(Government& g);
};

} // namespace moses
