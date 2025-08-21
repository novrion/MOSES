#include "ui.h"
#include "moses.h"
#include "labour.h"
#include "firm.h"
#include "household.h"
#include "bank.h"
#include "market.h"
#include "external.h"
#include <iostream>
#include <iomanip>

namespace moses {

int ui::start_simulation() {
	std::cout << "========================================================" << std::endl;
    std::cout << "                        M O S E S                       " << std::endl;
    std::cout << "           Model Of the Swedish Economic System         " << std::endl;
    std::cout << "========================================================" << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "                     Elias Lindstenz                    " << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "                          2025                          " << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "         Originally developed in APL by ...             " << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "         Converted to C++ by Elias Lindstenz            " << std::endl;
	std::cout << "                                                        " << std::endl;
	std::cout << "========================================================" << std::endl;
	std::cout << "                                                        " << std::endl;
	
	int simulation_length = get_simulation_length();
	
	
	return simulation_length;
}

int ui::get_simulation_length() {
	while (true) {
		std::string input;
		std::cout << "             Simulation Length (years): ";

		if (!std::getline(std::cin, input)) {
			std::cout << "                    Invalid Input                        " << std::endl;
			continue;
		}

		try {
			size_t pos = 0;
			int simulation_length = std::stoi(input, &pos);

			if (pos != input.length() || simulation_length <= 0) {
				std::cout << "                    Invalid Input                        " << std::endl;
				continue;
			}

			std::cout << "                                                        " << std::endl;
			std::cout << "========================================================" << std::endl;
			std::cout << "                                                        " << std::endl;

			return simulation_length;

		} catch(const std::exception&) {
			std::cout << "                    Invalid Input                        " << std::endl;
			continue;
		}
	}
}

void ui::print_simulation(const Moses& moses) {
    std::cout << "\n================ MOSES Simulation Results ================\n\n";

    // Print GNP components in current prices
    std::cout << "GNP Accounting (Current Prices):\n";
    std::cout << "--------------------------------------------\n";

    // First print explicit sector production (first 4 values)
    std::cout << "Production (Explicit Sectors):\n";
    for(int i = 0; i < NMKT; i++) {
        std::cout << "  Sector " << i << std::setw(34) 
                  << std::right
                  << moses.gnpcur[i] << "\n";
    }

    // Then print external sector production (next 6 values)
    std::cout << "Production (External Sectors):\n";
    for(int i = NMKT; i < NMKT + NEXT; i++) {
        std::cout << "  Sector " << i << std::setw(34)
                  << std::right
                  << moses.gnpcur[i] << "\n";
    }

    // Print remaining current price components
    const char* gnpcur_components[] = {
        "Indirect Taxes",
        "Subsidies",
        "Government Wages",
        "Government Purchases",
		"Private Consumption",
        "Investments (Explicit Sectors)",
        "Investments (External Sectors)",
        "Investments (Housing)",
        "Government Investments",
        "Inventory Changes",
        "Exports",
        "Imports"
    };
    for (int i = NMKT + NEXT; i < NGNPCUR; i++) {
        std::cout << std::left << std::setw(32) << gnpcur_components[i - (NMKT + NEXT)]
                  << std::right << std::setw(12)
                  << moses.gnpcur[i] << "\n";
    }
    std::cout << "\n";

    // Print GNP components in fixed prices
    std::cout << "GNP Accounting (Fixed Prices):\n";
    std::cout << "--------------------------------------------\n";

    // First print explicit sector production (first 4 values)
    std::cout << "Production (Explicit Sectors):\n";
    for(int i = 0; i < NMKT; i++) {
        std::cout << "  Sector " << i << std::setw(34) 
                  << std::right
                  << moses.gnpfix[i] << "\n";
    }

    // Then print external sector production (next 6 values)
    std::cout << "Production (External Sectors):\n";
    for(int i = NMKT; i < NMKT + NEXT; i++) {
        std::cout << "  Sector " << i << std::setw(34)
                  << std::right
                  << moses.gnpfix[i] << "\n";
    }

    // Print remaining fixed price components
    const char* gnpfix_components[] = {
        "Government Wages",
        "Government Purchases",
        "Private Consumption",
        "Investments (Explicit Sectors)",
        "Investments (External Sectors)",
        "Investments (Housing)",
        "Government Investments",
        "Inventory Changes",
        "Exports",
        "Imports"
    };
    for (int i = NMKT + NEXT; i < NGNPFIX; i++) {
        std::cout << std::left << std::setw(32) << gnpfix_components[i - (NMKT + NEXT)]
                  << std::right << std::setw(12)
                  << moses.gnpfix[i] << "\n";
    }
    std::cout << "\n";

    // Print unemployment rate
    std::cout << "Unemployment Rate: "
              << moses.lab->ru * 100 << "%\n";

    // Print number of active firms per market
    std::cout << "\nActive Firms per Market:\n";
    std::cout << "Market    Firms\n";
    std::cout << "---------------\n";
    std::vector<int> market_counts(moses.markets.size(), 0);
    for (const auto& firm : moses.firms) {
        market_counts[firm->market_id]++;
    }
    for (size_t i = 0; i < market_counts.size(); i++) {
        std::cout << std::setw(2) << i << std::setw(12) << market_counts[i] << "\n";
    }
    std::cout << "\n";

    // Print trade balance
    std::cout << "Trade Balance:\n";
    std::cout << "Exports: " << moses.export_ << "\n";
    std::cout << "Imports: " << moses.import_ << "\n";
    std::cout << "Balance: " << moses.export_ - moses.import_ << "\n\n";

    std::cout << "====================================================\n";
}

void ui::print_year(const Moses& moses) {
    // Print year header
    std::cout << "\n============ Year " << moses.year << " ============\n";

    // Print basic economic indicators
    std::cout << "Unemployment Rate: "
              << moses.lab->ru * 100 << "%\n";

    // Calculate and print total GNP (sum all production)
    double total_gnp_current = 0;
    for (int i = 0; i < NMKT + NEXT; i++) {
        total_gnp_current += moses.gnpcur[i];
    }
    std::cout << "Total GNP (Current): "
              << total_gnp_current << "\n";

    // Print total number of firms
    std::cout << "Total Firms: " << moses.firms.size() << "\n";

    // Print firms per market
    std::cout << "Firms per Market: ";
    std::vector<int> market_counts(moses.markets.size(), 0);
    for (const auto& firm : moses.firms) {
        market_counts[firm->market_id]++;
    }
    for (size_t i = 0; i < market_counts.size(); i++) {
        std::cout << market_counts[i];
        if (i < market_counts.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "\n";

    // Print trade metrics
    std::cout << "Trade Balance: "
              << moses.export_ - moses.import_ << "\n";

    // Print inflation rate (CPI)
    std::cout << "Inflation Rate: "
              << moses.hh->qdcpi * 100 << "%\n";

	// Print interest rate
    std::cout << "Interest Rate: "
              << moses.bank->ri * 100 << "%\n";

    // Calculate and print average wage level across firms
    double total_wage = 0;
    int total_workers = 0;
    for (const auto& firm : moses.firms) {
        total_wage += firm->qw * firm->l;
        total_workers += firm->l;
    }
    double avg_wage = total_workers > 0 ? total_wage / total_workers : 0;
    std::cout << "Average Wage Level: "
              << avg_wage << "\n";

    std::cout << "================================\n";
}

} // namespace moses
