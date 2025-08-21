#pragma once
#include <string>
#include <vector>

struct Params;
class Bank;
class Firm;
class Government;
class Household;
class Labour;
class Market;
class External;

namespace moses::csv {

const std::string ParamsPath     = "data/params.csv";
const std::string BankPath       = "data/bank.csv";
const std::string GovernmentPath = "data/government.csv";
const std::string HouseholdPath  = "data/household.csv";
const std::string LabourPath     = "data/labour.csv";
const std::string FirmPath       = "data/firms.csv";
const std::string MarketPath     = "data/markets.csv";
const std::string ExternalPath   = "data/externals.csv";

std::vector<std::string> split_line(const std::string& line);
std::vector<std::string> read_file(const std::string& path);

Params parse_params(const std::string& line);
std::unique_ptr<Bank> parse_bank(const Moses& moses_ref, const std::string& line);
std::unique_ptr<Government> parse_government(const Moses& moses_ref, const std::string& line);
std::unique_ptr<Household> parse_household(const Moses& moses_ref, const std::string& line);
std::unique_ptr<Labour> parse_labour(const Moses& moses_ref, const std::string& line);
std::unique_ptr<Firm> parse_firm(const Moses& moses_ref, const Params& params, const std::string& line);
std::unique_ptr<Market> parse_market(Moses& moses_ref, const Params& params, const std::string& line);
std::unique_ptr<External> parse_external(Moses& moses_ref, const std::string& line);

} // namespace moses::csv
