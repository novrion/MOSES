#include "moses.h"
#include "moses_impl.h"
#include "params.h"
#include "bank.h"
#include "government.h"
#include "household.h"
#include "labour.h"
#include "market.h"
#include "external.h"
#include "firm.h"
#include "csv.h"
#include "ui.h"
#include <iostream>
#include <algorithm>

namespace moses {

Moses::Moses() = default;
Moses::~Moses() = default;

int Moses::start_simulation() {
	return ui::start_simulation();
}

void Moses::print_simulation() {
	ui::print_simulation(*this);
}

void Moses::initialise(const int verbose) {
	std::cout << "[ Initialising MOSES... ]\n";

	try {
		const auto params = csv::parse_params(csv::read_file(csv::ParamsPath).at(0));
		niter = params.niter;
		marketiter = params.marketiter;

		bank = std::move(csv::parse_bank      (*this, csv::read_file(csv::BankPath).at(0)));
		gov  = std::move(csv::parse_government(*this, csv::read_file(csv::GovernmentPath).at(0)));
		hh   = std::move(csv::parse_household (*this, csv::read_file(csv::HouseholdPath).at(0)));
		lab  = std::move(csv::parse_labour    (*this, csv::read_file(csv::LabourPath).at(0)));

		for (const auto& line : csv::read_file(csv::FirmPath)) {
			firms.push_back(csv::parse_firm(*this, params, line));
		}

		for (const auto& line : csv::read_file(csv::MarketPath)) {
			markets.push_back(csv::parse_market(*this, params, line));
		}

		for (const auto& line : csv::read_file(csv::ExternalPath)) {
			externals.push_back(csv::parse_external(*this, line));
		}

		initialise_extra();
	} catch (const std::exception& e) {
		throw std::runtime_error("Failed to initialise MOSES: " + std::string(e.what()));
	}

	std::cout << "[ Successfully initialised MOSES ]\n";

	if (verbose) {
		print_initialisation(verbose);
	}
}

void Moses::initialise_extra() {
	for (auto& firm : firms) {
			firm->initialise_extra(*markets[firm->market_id]);
		}

	bank->initialise_extra(*gov, *hh);
	gov->initialise_extra(*bank);
}

void Moses::print_initialisation(const int verbose) {
	std::cout << "=============== Initialisation Values ===============" << std::endl;
	std::cout << std::endl;

	std::cout << "========== Params ==========" << std::endl;
	std::cout << "niter:      " << niter << std::endl;
	std::cout << "marketiter: " << marketiter << std::endl;
	std::cout << std::endl;
	std::cout << "e1      (0): " << firms[0]->e1	  << "   (" << firms.size() << "): " << firms[firms.size()-1]->e1 << std::endl;
	std::cout << "e2      (0): " << firms[0]->e2      << "   (" << firms.size() << "): " << firms[firms.size()-1]->e2 << std::endl;
	std::cout << "eps     (0): " << firms[0]->eps     << "   (" << firms.size() << "): " << firms[firms.size()-1]->eps << std::endl;
	std::cout << "r       (0): " << firms[0]->r       << "   (" << firms.size() << "): " << firms[firms.size()-1]->r << std::endl;
	std::cout << "smp     (0): " << firms[0]->smp     << "   (" << firms.size() << "): " << firms[firms.size()-1]->smp << std::endl;
	std::cout << "sms     (0): " << firms[0]->sms     << "   (" << firms.size() << "): " << firms[firms.size()-1]->sms << std::endl;
	std::cout << "smt     (0): " << firms[0]->smt     << "   (" << firms.size() << "): " << firms[firms.size()-1]->smt << std::endl;
	std::cout << "smw     (0): " << firms[0]->smw     << "   (" << firms.size() << "): " << firms[firms.size()-1]->smw << std::endl;
	std::cout << "fip     (0): " << firms[0]->fip     << "   (" << firms.size() << "): " << firms[firms.size()-1]->fip << std::endl;
	std::cout << "fis     (0): " << firms[0]->fis     << "   (" << firms.size() << "): " << firms[firms.size()-1]->fis << std::endl;
	std::cout << "fiw     (0): " << firms[0]->fiw     << "   (" << firms.size() << "): " << firms[firms.size()-1]->fiw << std::endl;
	std::cout << "expxdp  (0): " << firms[0]->expxdp  << "   (" << firms.size() << "): " << firms[firms.size()-1]->expxdp << std::endl;
	std::cout << "expxds  (0): " << firms[0]->expxds  << "   (" << firms.size() << "): " << firms[firms.size()-1]->expxds << std::endl;
	std::cout << "expxdw  (0): " << firms[0]->expxdw  << "   (" << firms.size() << "): " << firms[firms.size()-1]->expxdw << std::endl;
	std::cout << "loss    (0): " << firms[0]->loss    << "   (" << firms.size() << "): " << firms[firms.size()-1]->loss << std::endl;
	std::cout << "rho     (0): " << firms[0]->rho     << "   (" << firms.size() << "): " << firms[firms.size()-1]->rho << std::endl;
	std::cout << "resdown (0): " << firms[0]->resdown << "   (" << firms.size() << "): " << firms[firms.size()-1]->resdown << std::endl;
	std::cout << "resmax  (0): " << firms[0]->resmax  << "   (" << firms.size() << "): " << firms[firms.size()-1]->resmax << std::endl;
	std::cout << "tmsto   (0): " << firms[0]->tmsto   << "   (" << firms.size() << "): " << firms[firms.size()-1]->tmsto << std::endl;
	std::cout << "tmimsto (0): " << firms[0]->tmimsto << "   (" << firms.size() << "): " << firms[firms.size()-1]->tmimsto << std::endl;
	std::cout << "beta    (0): " << firms[0]->beta    << "   (" << firms.size() << "): " << firms[firms.size()-1]->beta << std::endl;
	std::cout << "imbeta  (0): " << firms[0]->imbeta  << "   (" << firms.size() << "): " << firms[firms.size()-1]->imbeta << std::endl;
	std::cout << "wtix    (0): " << firms[0]->wtix    << "   (" << firms.size() << "): " << firms[firms.size()-1]->wtix << std::endl;
	std::cout << "iota    (0): " << firms[0]->iota    << "   (" << firms.size() << "): " << firms[firms.size()-1]->iota << std::endl;
	std::cout << "rhobook (0): " << firms[0]->rhobook << "   (" << firms.size() << "): " << firms[firms.size()-1]->rhobook << std::endl;
	std::cout << "rtd     (0): " << firms[0]->rtd     << "   (" << firms.size() << "): " << firms[firms.size()-1]->rtd << std::endl;
	std::cout << "alfabw  (0): " << firms[0]->alfabw  << "   (" << firms.size() << "): " << firms[firms.size()-1]->alfabw << std::endl;
	std::cout << "betabw  (0): " << firms[0]->betabw  << "   (" << firms.size() << "): " << firms[firms.size()-1]->betabw << std::endl;
	std::cout << "redchbw (0): " << firms[0]->redchbw << "   (" << firms.size() << "): " << firms[firms.size()-1]->redchbw << std::endl;
	std::cout << "elinv   (0): " << firms[0]->elinv   << "   (" << firms.size() << "): " << firms[firms.size()-1]->elinv << std::endl;
	std::cout << "utref   (0): " << firms[0]->utref   << "   (" << firms.size() << "): " << firms[firms.size()-1]->utref << std::endl;
	std::cout << std::endl;
	std::cout << "maxdp (0): " << markets[0]->maxdp << "   (" << markets.size() << "): " << markets[markets.size()-1]->maxdp << std::endl;
	std::cout << std::endl;

	bool print_params = false;
	bool print_vectors = false;
	if (verbose > 1) {
		print_params = true;
		print_vectors = true;
	}

	bank->print_initialisation();
	gov->print_initialisation();
	hh->print_initialisation();
	lab->print_initialisation();

	if (verbose > 2) {
		for (const auto& firm : firms) {
			firm->print_initialisation(print_params, print_vectors);
		}
	}

	for (const auto& mkt : markets) {
		mkt->print_initialisation(print_params, print_vectors);
	}

	for (const auto& ext : externals) {
		ext->print_initialisation(print_vectors);
	}
}

double Moses::qdpk() const {
	double sum = 0.0;
	for (const auto& mkt : markets) {
		sum += gov->omega[mkt->id] * (mkt->qdpdom - (gov->qchtxva2() - gov->qchtxva1()));
	}
	for (const auto& ext : externals) {
		sum += gov->omega[ext->id] * (ext->qdpdom - (gov->qchtxva2() - gov->qchtxva1()));
	}
	return sum;
}

// Total quarterly wage sum (firm)
double Moses::qwsf() const {
	double sum = 0.0;
	for (const auto& firm : firms) {
		sum += firm->l * firm->qw/4.0;
	}
	return sum;
}

// Total wage sum (firm) - txw adjusted
double Moses::wstx() const {
	double sum = 0.0;
	for (const auto& firm : firms) {
		sum += firm->l * firm->qw * (1 - gov->txw());
	}
	return sum;
}

// Sum qbuy from a specific sector
double Moses::sum_qbuy(const int mkt_id) const {
	double sum = gov->qbuy[mkt_id] + hh->qbuy[mkt_id] + bank->qbuy[mkt_id];
	for (const auto& mkt : markets) {
		sum += mkt->qbuy[mkt_id];
	}
	for (const auto& ext : externals) { 
		sum += ext->qbuy[mkt_id];
	}// TODO: external should not be included for compute_external_buying()? Should be included for computation of qtbuy though. 
	return sum;
}

// Sum qdeschbw over firms, but only take positive values
double Moses::sum_qdeschbw_positive() const {
    double sum = 0.0;
    for (const auto& firm : firms) {
        sum += std::max(0.0, firm->qdeschbw); 
    }
    return sum;
}	

void Moses::nullify_firm(const size_t index) {
	if (firms.size() <= 1) {
		throw std::runtime_error("Attempting to nullify last firm");
	}

	lab->lu += firms[index]->l;
	firms.erase(firms.begin() + index);
}

void Moses::check_market_health(const bool verbose) {
	std::vector<int> market_counts(markets.size(), 0);
	for (const auto& firm : firms) {
		market_counts[firm->market_id]++;
	}

	if (verbose) {
		std::cout << "\nmkt     ntot\n";
		std::cout << "============\n";
	}
	
	for (size_t i = 0; i < markets.size(); i++) {
		if (market_counts[i] == 0) {
			throw std::runtime_error("Market crash");
		}

		if (verbose) {
			std::cout << i << "       " << market_counts[i] << "\n";
		}
	}
}

size_t Moses::choose_labour_target(const LabourSearchData& data) {
    double scaled_unemployment = data.ll.back() * lab->skrepa;
    double total_labour = std::accumulate(data.ll.begin(), data.ll.end() - 1, 0.0) + scaled_unemployment;
    
    double rand_val = static_cast<double>(std::rand()) / RAND_MAX * total_labour;
    
    size_t target = data.ll.size() - 1;  // Default to unemployment pool
    double cumsum = 0;
    for (size_t j = 0; j < data.ll.size(); j++) {
        // Use scaled value for unemployment in cumsum calculation
        cumsum += (j == data.ll.size() - 1) ? scaled_unemployment : data.ll[j];
        if (cumsum > rand_val) {
            target = j;
            break;
        }
    }
    return target;
}

// Gauss-Jordan Matrix Inversion Algorithm
std::vector<std::vector<double>> Moses::invert_matrix(std::vector<std::vector<double>> m) {
	const int n = (int) m.size();
	std::vector<std::vector<double>> inv = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
	for (int i = 0; i < n; i++) inv[i][i] = 1;

	for (int i = 0; i < n; i++) {
		// find max value pivot element in column
		double max_val = std::abs(m[i][i]);
		double pivot = i;
		for (int j = i + 1; j < n; j++) {
			if (std::abs(m[j][i]) > max_val) {
				max_val = std::abs(m[j][i]);
				pivot = j;
			}
		}

		// check for singular matrix (cannot be inverted)
		if (max_val < 1e-6) {
			throw std::runtime_error("Couldn't invert singular matrix.");
		}

		// swap rows if necessary
		if (pivot != i) {
			std::swap(m[i], m[pivot]);
			std::swap(inv[i], inv[pivot]);		
		}

		// scale row so that the pivot element is 1
		double pivot_scale = m[i][i];
		for (int j = 0; j < n; j++) {
			m[i][j] /= pivot_scale;
			inv[i][j] /= pivot_scale;
		}

		// row elimination
		for (int j = 0; j < n; j++) {
			if (i != j) {
				double factor = m[j][i];
				for (int k = 0; k < n; k++) {
					m[j][k] -= factor * m[i][k];
					inv[j][k] -= factor * inv[i][k];
				}
			}
		}
	}

	return inv;
}

void Moses::luupdate() {
	lab->luupdate(*gov);
	for (auto& firm : firms) {
		firm->luupdate(*lab);
	}
}

void Moses::prodfront() {
	for (auto& mkt : markets) mkt->mtec *= (1 + mkt->qdmtec);
	for (auto& firm : firms) firm->prodfront(*markets[firm->market_id]);
}

void Moses::quarterly_exp() {
	for (auto& firm : firms) {
		firm->quarterly_exp(relative_quarter == 0);
		firm->quarterly_targ();
	}

	for (auto& mkt : markets) {
		mkt->quarterly_exp(*gov);
	}

	for (auto& ext : externals) {
		ext->quarterly_exp(*gov);
	} 
}

void Moses::prodplan() {
	luupdate();
	prodfront();
	for (size_t i = 0; i < firms.size(); i++) {
		firms[i]->init_prodplan();
		if (!firms[i]->target_search()) {
			nullify_firm(i);
			i--;
		}
	}

	check_market_health();

	for (auto& firm : firms) {
		firm->qplanqsave = firm->qplanq;
	}
}

LabourSearchData Moses::labour_search_input() {
	LabourSearchData data;

    data.chl.reserve(firms.size());
	data.ww.reserve(firms.size());
	data.ll.reserve(firms.size() + 1);

	for (const auto& firm : firms) {
		firm->labour_search_input(data);
	}
	data.ll.push_back(lab->lu);

	return data;
}

void Moses::handle_successful_attack(LabourSearchData& data, size_t attacker, size_t target) {
	double now = std::min(lab->theta * data.ll[target], data.chl[attacker]);

	data.ll[attacker] += now;
	data.chl[attacker] -= now;
	data.ll[target] -= now;

	// Is the target a firm or the pool of unemployment?
	if (target < firms.size()) {
		data.chl[target] += now;
	}
}

void Moses::confront(LabourSearchData& data) {
	std::vector<size_t> order(firms.size());
	std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), //should only happen once. Not every nit!
			[&](size_t a, size_t b) {
				return (data.chl[a] / data.ll[a]) > (data.chl[b] / data.ll[b]);
			});

	for (int iter = 0; iter < niter; iter++) {
		for (size_t i : order) {
			if (data.chl[i] <= 0) continue;

			size_t target = choose_labour_target(data);

			if (target == data.ll.size() - 1) {
				handle_successful_attack(data, i, target);
			} else {
				if (data.ww[i] > data.ww[target] * (1 + lab->gamma)) {
					handle_successful_attack(data, i , target);
				} else {
					data.ww[i] += lab->ksi * (data.ww[target] * (1 + lab->gamma) - data.ww[i]);
				}
			}
		}
	}
}

void Moses::labour_search_output(LabourSearchData& data) {
	for (size_t i = 0; i < firms.size(); i++) {
		firms[i]->qchl = data.ll[i] - firms[i]->l;
		firms[i]->qchw = data.ww[i] - firms[i]->qw;
	}

	lab->lu = data.ll.back();

	for (auto& firm : firms) {
		firm->handle_labour_search_layoffs();
	}
}

void Moses::labour_search() {
	LabourSearchData data = labour_search_input();
	confront(data);
	labour_search_output(data);
}

void Moses::labour_update() {
	for (auto& firm : firms) {
		firm->handle_labour_update_layoffs(*lab);
	}

	double oldqw = sum_avg(&Firm::l, &Firm::qw);
	double newqw = 0.0;
	for (const auto& firm : firms) {
		newqw += (firm->l + firm->qchl) * (firm->qw + firm->qchw);
	}
	newqw = newqw / (sum(&Firm::l) + sum(&Firm::qchl));
	lab->qdwind = newqw / oldqw - 1;

	for (auto& firm : firms) {
		firm->update_labour_force_and_wage();
	}

	lab->qchru = lab->lu / (lab->lu + gov->l + sum(&Firm::l)) - lab->ru;
	lab->ru += lab->qchru;
}

void Moses::indalabour() {
	labour_search();
	labour_update();
}

void Moses::labour_market() {
	gov->glabour(*lab);
	indalabour();
	for (auto& firm : firms) {
		firm->planqrevise(*markets[firm->market_id]);
	}
}

void Moses::export_market() {
	for (auto& firm : firms) {
		firm->adjust_foreign(*markets[firm->market_id], *gov);
	}

	for (auto& mkt : markets) {
		mkt->qpfor *= (1 + mkt->qdpfor());
	}

	for (auto& firm : firms) {
		firm->qsfor = firm->qsufor * markets[firm->market_id]->qpfor * (1 + markets[firm->market_id]->rsubs);
	}

	gov->qsubsfor = 0.0;
	qexport = 0.0;
	for (const auto& firm : firms) {
		gov->qsubsfor += firm->qsufor * markets[firm->market_id]->qpfor * markets[firm->market_id]->rsubs;
		qexport += firm->qsfor / (1 + markets[firm->market_id]->rsubs);
	}
}

void Moses::market_entrance() {
	for (auto& firm : firms) {
		firm->qoptsudom = (1 - firm->x) * firm->qoptsu;
	}

	for (auto& mkt : markets) {
		double sum_ch_qoptsudom = 0.0;
		for (const auto& firm : firms) {
			if (firm->market_id == mkt->id) {
				sum_ch_qoptsudom += firm->qoptsudom * (firm->qexpp / firm->qp);
			}
		}

		mkt->qprelpdom = mkt->qpdom * (1 + gov->qchtxva2()) * sum_ch_qoptsudom / sum(&Firm::qoptsudom, mkt->id);
	}
}

std::vector<double> Moses::compute_trial_prices() {
	std::vector<double> pt(NSEC);
	for (auto& mkt : markets) {
		pt[mkt->id] = mkt->compute_trial_price(*gov);
	}
	for (auto& ext : externals) {
		pt[ext->id] = ext->compute_trial_price();
	}

	return pt;
}

void Moses::compute_external_production() {

	// populate io3 matrix with io
	std::vector<std::vector<double>> io3 = std::vector<std::vector<double>>(NEXT, std::vector<double>(NEXT, 0.0));
	for (const auto& ext : externals) {
		for (size_t i = 0; i < NEXT; i++) {
			io3[ext->id - NMKT][i] = ext->io[NMKT + i];
		}
	}

	// modify io3 before matrix inversion
	for (size_t i = 0; i < NEXT; i++) {
		for (size_t j = 0; j < NEXT; j++) {

			// apply ratio
			io3[i][j] *= (1 - externals[j]->imp)/(1 - externals[j]->x);
			
			// subtract from I (identity matrix)
			if (i == j) {
				io3[i][j] = 1 - io3[i][j];
			} else {
				io3[i][j] = -io3[i][j]; // io3[i][j] = 0 - io3[i][j]
			}
		}
	}

	std::vector<std::vector<double>> inverted = invert_matrix(io3);

	// compute buying from each external sector, with a ratio applied
	std::array<double, NEXT> sum_qbuy_modified{};
	for (const auto& ext : externals) {
		const int id = ext->id;
		sum_qbuy_modified[id - NMKT] = sum_qbuy(id) * (1 - ext->imp)/(1 - ext->x);
	}

	// matrix multiplication of total buying from each external sector and the inverted input-output matrix
	for (auto& ext : externals) {
		ext->qq = 0;
		const int relative_id = ext->id - NMKT;
		for (size_t i = 0; i < NEXT; i++) {
			ext->qq += inverted[i][relative_id] * sum_qbuy_modified[i];
		}
	}
}

void Moses::compute_total_buying() {
	for (auto& mkt : markets) {
		mkt->qtbuy = sum_qbuy(mkt->id);
	}
	for (auto& ext : externals) {
		ext->qtbuy = sum_qbuy(ext->id);
	}
}

void Moses::compute_buying(const std::vector<double>& qsp, const std::vector<double>& pt) {
	for (auto& mkt : markets) {
		mkt->compute_buying();
	}
	gov->compute_buying(pt);
	hh->compute_buying(qsp, pt);
	bank->compute_buying(*gov, pt);

	compute_external_production();
	for (auto& ext : externals) {
		ext->compute_buying();
	}

	compute_total_buying();	
}

void Moses::market_confront(std::vector<double>& qsp, std::vector<double>& pt) {
	pt = compute_trial_prices();
	gov->compute_consumption(pt);

	for (int iter = 0; iter < marketiter; iter++) {
		qsp = hh->compute_expenditures(pt);
		compute_buying(qsp, pt);
		if (iter < marketiter - 1) {
			for (auto& mkt : markets) {
				mkt->price_adjust(pt);
			}
		}
	}
}

void Moses::compute_imports(const std::vector<double>& pt) {
	qimport = 0.0;
	for (auto& mkt : markets) {
		mkt->compute_imports(pt);
	}
	for (auto& ext : externals) {
		ext->compute_imports(*gov, pt);
	}
}

void Moses::domestic_result(const std::vector<double>& pt) {
	for (auto& mkt : markets) {
		mkt->domestic_result(pt);
	}
	for (auto& ext : externals) {
		ext->domestic_result(pt);
	}
}

void Moses::external_sectors(const std::vector<double>& pt) {
	hh->qinpay = 0.0;
	for (auto& ext : externals) {
		ext->external_sectors(*gov, *hh, pt);
	}
	hh->qinpay -= (gov->qinvbld + gov->qinvin);
}

void Moses::indirect_taxes(const std::vector<double>& qsp) {
	gov->qvataximp = 0.0;
	for (auto& mkt : markets) {
		mkt->indirect_taxes(*bank, *gov, *hh, qsp);
	}
	for (auto& ext : externals) {
		ext->indirect_taxes(*bank, *gov, *hh, qsp);
	}
}

void Moses::domestic_market() {
	market_entrance();
	hh->init(*bank, *gov, *lab);

	std::vector<double> pt;
	std::vector<double> qsp;
	market_confront(qsp, pt);
	
	compute_imports(pt);
	domestic_result(pt);
	external_sectors(pt);
	hh->household_update(qsp, pt);
	indirect_taxes(qsp);
}

void Moses::firm_sto() {
	std::vector<double> limsto, lower, upper;
	for (auto& firm : firms) {
		double lim = firm->sto + firm->qq - firm->qsufor;
		limsto.push_back(lim); 
		lower.push_back(std::min(lim, firm->minsto));
		upper.push_back(std::min(lim, firm->maxsto));
	}

	gov->qsubsdom = 0;
	for (auto& mkt : markets) {
		mkt->firm_sto(*gov, limsto, upper, lower);
	}
}

void Moses::sto_system() {
	firm_sto();
	for (auto& firm : firms) {
		firm->reference_inventory_levels(*markets[firm->market_id]);
	}
}

void Moses::invfin() {
	for (auto& firm : firms) {
		firm->invfin(*bank, *gov);
	}

	qintf = 0;
	qintk2 = 0;
	gov->qctax = 0;
	hh->qtdiv = 0;
	gov->qsubscash = 0;
	for (const auto& firm : firms) {
		qintf += firm->bw * bank->rif/4.0;
		qintk2 += firm->k2 * bank->rik2/4.0;
		gov->qctax += firm->qtax;
		hh->qtdiv += firm->qdiv;
		gov->qsubscash += firm->qs * gov->rsubscash;
	}
}

void Moses::invfin_adjustments() {
	double oldinv = sum(&Firm::qinvlag);

	for (size_t i = 0; i < firms.size(); i++) {
		if (!firms[i]->invfin_adjustments()) {
			nullify_firm(i);
			i--;
		}
	}

	check_market_health();

	bank->qtchbw = sum(&Firm::qchbw);
	bank->qtchinv = sum(&Firm::qinvlag) - oldinv;
	bank->qtchk2 = sum(&Firm::qchk2);
}

void Moses::monetary_sector() {
	bank->bank_transactions(*gov, *hh);
	bank->credit_market(*gov, *hh);
	invfin_adjustments();
	bank->bank_update(*gov, *hh);
}

void Moses::national_accounting() {
	std::array<double, NGNPCUR> qgnpcur;
	std::array<double, NGNPFIX> qgnpfix;
	qgnpcur.fill(0.0);
	qgnpfix.fill(0.0);

	// Populate temporary vectors for easier calculations
	std::vector<double> pref;
	std::vector<double> qpdom;
	std::vector<double> qchtsto;
	std::vector<double> qtbuyfor;
	for (const auto& mkt : markets) {
		pref.push_back(mkt->pref);
		qpdom.push_back(mkt->qpdom);
		qchtsto.push_back(mkt->qchtsto);
		qtbuyfor.push_back(mkt->qtbuyfor);
	}
	for (const auto& ext : externals) {
		pref.push_back(ext->pref);
		qpdom.push_back(ext->qpdom);
		qtbuyfor.push_back(ext->qtbuyfor);
	}

	// Explicit Sector Production
	for (auto& mkt : markets) {
		mkt->national_accounting(qgnpcur, qgnpfix);
	}

	// External Sector Production
	for (auto& ext : externals) {
		ext->national_accounting(qgnpcur, qgnpfix);
	}

	size_t cur_i = NSEC; // index for qgnpcur (start after explicit and external sector production)
	size_t fix_i = NSEC; // index for qgnpfix (start after explicit and external sector production)

	// Indirect Taxes
	qgnpcur[cur_i++] = sum(gov->qvatax) + sum(&Market::qchtstocurm) - sum(&Market::qchtstocurf) - gov->qvataximp;

	// Subsidies
	qgnpcur[cur_i++] = -gov->qsubs - gov->qsubscash;

	// Government Wages
	qgnpcur[cur_i++] = gov->qws;

	// Government Purchases
	qgnpcur[cur_i++] = sum(gov->qpurch);

	// Private Consumption
	qgnpcur[cur_i++] = hh->nh * sum(hh->qsp) - hh->qsav;

	// Investments made by the explicit model sectors
	qgnpcur[cur_i++] = sum(&Firm::qinvlag) - bank->qtchinv;

	// Investments made by external sectors, exclusive of housing
	qgnpcur[cur_i++] = gov->qinvin;

	// Investments for residential construction
	qgnpcur[cur_i++] = gov->qinvbld;

	// Government investments
	qgnpcur[cur_i++] = gov->qinv;

	// Inventory changes
	qgnpcur[cur_i++] = sum(&Market::qchtstocurm);

	// Exports
	qgnpcur[cur_i++] = qexport;

	// Imports
	qgnpcur[cur_i++] = -(qimport - gov->qvataximp);



	// Government wages, deflated to the government wage level of the reference year
	qgnpfix[fix_i++] = gov->l * gov->wgref/4.0;

	// Government purchases
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * gov->qpurch[i])/qpdom[i];
	}
	fix_i++;
	
	// Private consumption; sum over non-saving categories
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * hh->qsp[i] * hh->nh)/qpdom[i];
	}
	fix_i++;

	// Investments made by the explicit model sectors
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * gov->omega[i] * (sum(&Firm::qinvlag) - bank->qtchinv))/(qpdom[i] * (1 - gov->txva2)/(1 - gov->txva1));
	}
	fix_i++;

	// Investments made by external sectors, exclusive of housing
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * gov->omegain[i] * gov->qinvin)/(qpdom[i] * (1 - gov->txva2)/(1 - gov->txva1));
	}
	fix_i++;

	// Investments for residential construction
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * gov->omegabld[i] * gov->qinvbld)/(qpdom[i] * (1 - gov->txva2)/(1 - gov->txva1));
	}
	fix_i++;

	// Government investments
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] += (pref[i] * gov->omegag[i] * gov->qinv)/(qpdom[i] * (1 - gov->txva2)/(1 - gov->txva1));
	}
	fix_i++;

	// Inventory changes
	for (size_t i = 0; i < NMKT; i++) {
		qgnpfix[fix_i] += pref[i] * qchtsto[i];
	}
	fix_i++;

	// Exports
	for (const auto& mkt : markets) {
		qgnpfix[fix_i] += pref[mkt->id] * sum(&Firm::qsufor, mkt->id);
	}
	for (const auto& ext : externals) {
		qgnpfix[fix_i] += pref[ext->id] * ext->x * ext->qq;
	}
	fix_i++;

	// Imports
	for (size_t i = 0; i < NSEC; i++) {
		qgnpfix[fix_i] -= pref[i] * qtbuyfor[i];
	}



	// Cumulations
	for (size_t i = 0; i < cumgnpcur.size(); i++) {
		cumgnpcur[i] += qgnpcur[i];
	}	

	for (size_t i = 0; i < cumgnpfix.size(); i++) {
		cumgnpfix[i] += qgnpfix[i];
	}

	cumexport += qexport;
	cumimport += qimport;

	if (relative_quarter == 3) {
		gnpcur = cumgnpcur;
		gnpfix = cumgnpfix;
		export_ = cumexport;
		import_ = cumimport;
	}
}

void Moses::simulate_quarter() {
	quarterly_exp();
	prodplan();
	labour_market();
	export_market();
	domestic_market();
	sto_system();

	for (auto& firm : firms) {
		firm->finalqpqsqm();
		firm->quarterly_cum();
	}

	invfin();
	gov->accounting(*bank, *hh);
	monetary_sector();
	national_accounting();
}

void Moses::yearly_init() {
	cumgnpcur.fill(0.0);
	cumgnpfix.fill(0.0);
	cumexport = 0.0;
	cumimport = 0.0;

	gov->yearly_init();
	for (auto& firm : firms) {
		firm->yearly_init();
	}
}

void Moses::simulate_year() {
	yearly_init();

	for (auto& firm : firms) {
		firm->yearly_exp();
		firm->yearly_targ();
	}

	for (relative_quarter = 0; relative_quarter < 4; relative_quarter++) {
		simulate_quarter();
		quarter++;
	}

	for (auto& firm : firms) {
		firm->yearly_update();
	}
}

void Moses::simulate(const int simulation_length, const int verbose) {
	std::cout << "[ Running MOSES simulation... ]\n";
	for (year = 0; year < simulation_length; year++) {
		simulate_year();
		if (verbose > 1) {
			ui::print_year(*this);
		}
	}
	std::cout << "[ Finished MOSES simulation successfully ]\n";

	if (verbose > 0) {
		print_simulation();
	}
}

} // moses namespace
