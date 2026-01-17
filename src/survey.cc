/*
Copyright (C) 2012-2025 tim cotter. All rights reserved.
*/

#include <algorithm>
#include <curl/curl.h>
#include <iomanip>
#include <sstream>
#include <vector>

#include <aggiornamento/aggiornamento.h>
#include <aggiornamento/log.h>

#include "survey.h"

namespace {

static constexpr char kVotingSurveyUrl[] = "https://docs.google.com/spreadsheets/d/11CS8R4pbYFDQcHaXjGxAquIuVWzb8vjciY_fS-Tz_Rk/export?format=tsv&gid=2041703740";

enum Method {
    kApproval,
    kBordaCount,
    kBucklin,
    kCondorcet,
    kGuthrie,
    kInstantRunoff,
    kScore,
    //
    kNumMethods
};
static constexpr int kNumRounds = kNumMethods;
static constexpr int kMaxScore = 10;

class Rankings {
public:
    int rankings_[kNumMethods];
};

class VotingSurveyImpl {
public:
    VotingSurveyImpl() = default;
    ~VotingSurveyImpl() = default;

    std::vector<std::string> method_names_;
    int longest_name_ = 0;
    std::string buffer_;
    std::stringstream sheet_;
    std::stringstream cells_;
    std::string line_;
    std::string s0_;
    std::string s1_;
    std::string s2_;
    std::string s3_;
    std::string s4_;
    std::string s5_;
    std::string s6_;

    int num_votes_ = 0;
    int majority_ = 0;
    int winners_[kNumMethods];
    int approval_counts_[kNumMethods];
    int borda_counts_[kNumMethods];
    int bucklin_counts_[kNumMethods][kNumRounds];
    std::vector<Rankings> condorcet_ballots_;
    int guthrie_counts_[kNumMethods];
    std::vector<Rankings> instant_ballots_;
    int score_totals_[kNumMethods];

    void run() noexcept {
        init();
        readRows();
        analyze();
    }

    void init() noexcept {
        /** table of method names. **/
        method_names_.reserve(kNumMethods);
        method_names_.push_back("Approval");
        method_names_.push_back("Borda Count");
        method_names_.push_back("Bucklin");
        method_names_.push_back("Condorcet");
        method_names_.push_back("Guthrie");
        method_names_.push_back("Instant Runoff");
        method_names_.push_back("Score");

        /** find the longest name for formatting. **/
        for (auto &&name : method_names_) {
            int len = name.size();
            longest_name_ = std::max(longest_name_, len);
        }

        /** load the data from the google sheet. **/
        loadSurveyDataFromGoogleSheet();

        /** stuff it into a string stream. **/
        sheet_ = std::stringstream(std::move(buffer_));
        buffer_.clear();

        /** skip the first line. **/
        std::getline(sheet_, line_);

        /** initialize the tallies. **/
        num_votes_ = 0;
        majority_ = 0;
        for (int method = 0; method < kNumMethods; ++method) {
            winners_[method] = -1;
            approval_counts_[method] = 0;
            borda_counts_[method] = 0;
            guthrie_counts_[method] = 0;
            score_totals_[method] = 0;

            for (int round = 0; round < kNumMethods; ++round) {
                bucklin_counts_[method][round] = 0;
            }
        }
    }

    void loadSurveyDataFromGoogleSheet() noexcept {
        /** init curl. **/
        curl_global_init(CURL_GLOBAL_DEFAULT);
        auto curl = curl_easy_init();
        if (curl) {
            /** set url. **/
            curl_easy_setopt(curl, CURLOPT_URL, kVotingSurveyUrl);
            /** set callback to write to the buffer. **/
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlCallback);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &buffer_);
            /** follow up to 10 redirects. **/
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1LL);
            curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 10LL);
            /** read it, log error. **/
            auto result = curl_easy_perform(curl);
            if (result != CURLE_OK) {
                auto err = curl_easy_strerror(result);
                LOG("curl failed: "<<err);
            }
            /** exit. **/
            curl_easy_cleanup(curl);
        }
        /** exit. **/
        curl_global_cleanup();
    }

    static size_t curlCallback(
        void *contents,
        size_t size,
        size_t nmemb,
        void *userp
    ) noexcept {
        size *= nmemb;
        auto &buffer = * (std::string *) userp;
        buffer.append((char *) contents, size);
        return size;
    }

    void readRows() noexcept {
        /** read the rows. **/
        int row = 0;
        while (std::getline(sheet_, line_)) {
            //LOG("row["<<row<<"]:");
            ++row;
            ++num_votes_;

            /** replace spaces with underscores. **/
            std::replace(line_.begin(), line_.end(), ' ', '_');

            /** replace empty cells with underscores. **/
            for(;;) {
                auto found = line_.find("\t\t");
                if (found == std::string::npos) {
                    break;
                }
                line_.replace(found, 2, "\t_\t");
            }

            /** read cells. **/
            cells_ = std::stringstream(line_);

            /** date/time **/
            cells_ >> s0_;

            /** approval **/
            tallyApproval();

            /** borda count **/
            tallyBordaCount();

            /** bucklin **/
            tallyBucklin();

            /** condorcet **/
            tallyCondorcet();

            /** guthrie **/
            tallyGuthrie();

            /** instant runoff **/
            tallyInstantRunoff();

            /** score **/
            tallyScore();
        }

        /**
        half of the votes (round down) plus one.
        works for even and odd numbers of votes.
        **/
        majority_ = num_votes_ / 2 + 1;
    }

    void tallyApproval() noexcept {
        /** the cell contains the names of the method. **/
        cells_ >> s0_;
        for (int method = 0; method < kNumMethods; ++method) {
            auto found = s0_.find(method_names_[method]);
            if (found != std::string::npos) {
                ++approval_counts_[method];
            }
        }
    }

    void tallyBordaCount() noexcept {
        /**
        the next 7 cells contain the rankings of each method.
        1 is the highest -> 6 points.
        7 is the lowest -> 0 points.
        **/
        for (int method = 0; method < kNumMethods; ++method) {
            cells_ >> s0_;
            int rank = toNumber(s0_);
            if (rank >= 1 && rank <= kNumMethods) {
                borda_counts_[method] += kNumMethods - rank;
            }
        }
    }

    void tallyBucklin() noexcept {
        /**
        the next 7 cells contain the round each method is approved.
        the method is approved for every round thereafter.
        **/
        for (int method = 0; method < kNumMethods; ++method) {
            cells_ >> s0_;
            int first_round = toNumber(s0_);
            if (first_round >= 1 && first_round < kNumRounds) {
                for (int round = first_round - 1; round < kNumMethods; ++round) {
                    ++bucklin_counts_[method][round];
                }
            }
        }
    }

    void tallyCondorcet() noexcept {
        /**
        the next 7 cells contain the rank of each method.
        **/
        Rankings r;
        for (int method = 0; method < kNumMethods; ++method) {
            cells_ >> s0_;
            r.rankings_[method] = toNumber(s0_);
        }
        condorcet_ballots_.push_back(r);
    }

    void tallyGuthrie() noexcept {
        /** the cell contains the name of the method. **/
        cells_ >> s0_;
        for (int method = 0; method < kNumMethods; ++method) {
            auto found = s0_.find(method_names_[method]);
            if (found != std::string::npos) {
                ++guthrie_counts_[method];
                break;
            }
        }
    }

    void tallyInstantRunoff() noexcept {
        /**
        the next 7 cells contain the rank of each method.
        **/
        Rankings r;
        for (int method = 0; method < kNumMethods; ++method) {
            cells_ >> s0_;
            r.rankings_[method] = toNumber(s0_);
        }
        instant_ballots_.push_back(r);
    }

    void tallyScore() noexcept {
        /**
        the next 7 cells contain the score for each method.
        10 is the highest -> 10 points.
        0 is the lowest -> 0 points.
        **/
        for (int method = 0; method < kNumMethods; ++method) {
            cells_ >> s0_;
            int score = toNumber(s0_);
            if (score >= 0 && score <= kMaxScore) {
                score_totals_[method] += score;
            }
        }
    }

    /**
    return the first non-negative integer in a string.
    skip leading non-digits.
    **/
    int toNumber(
        std::string &str
    ) noexcept {
        int len = str.size();
        int i = 0;
        for ( ; i < len; ++i) {
            int ch = str[i];
            if (ch >= '0' && ch <= '9') {
                break;
            }
        }
        int x = -1;
        if (i < len) {
            x = std::stoi(&str[i]);
        }
        return x;
    }

    void analyze() noexcept {
        LOG("Number of ballots: "<<num_votes_);
        analyzeApproval();
        analyzeBordaCount();
        analyzeBucklin();
        analyzeCondorcet();
        analyzeGuthrie();
        analyzeInstantRunoff();
        analyzeScore();
        reportWinners();
    }

    void analyzeApproval() noexcept {
        LOG("Results by Approval Voting (approvals):");
        int winner = -1;
        int max_count = -1;
        for (int method = 0; method < kNumMethods; ++method) {
            int count = approval_counts_[method];
            LOG("  "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<count);
            if (max_count < count) {
                max_count = count;
                winner = method;
            }
        }
        LOG("  Winner: "<<method_names_[winner]);
        winners_[kApproval] = winner;
    }

    void analyzeBordaCount() noexcept {
        LOG("Results by Borda Count (totals):");
        int winner = -1;
        int max_count = -1;
        for (int method = 0; method < kNumMethods; ++method) {
            int count = borda_counts_[method];
            LOG("  "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<count);
            if (max_count < count) {
                max_count = count;
                winner = method;
            }
        }
        LOG("  Winner: "<<method_names_[winner]);
        winners_[kBordaCount] = winner;
    }

    void analyzeBucklin() noexcept {
        LOG("Results by Bucklin Voting (approvals):");
        int winner = -1;
        int leader = -1;
        for (int round = 0; round < kNumMethods; ++round) {
            LOG("  Round["<<round+1<<"]:");
            int max_count = -1;
            for (int method = 0; method < kNumMethods; ++method) {
                int count = bucklin_counts_[method][round];
                LOG("    "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<count);
                if (max_count < count) {
                    max_count = count;
                    leader = method;
                }
            }
            if (max_count >= majority_) {
                winner = leader;
                break;
            }
        }
        if (winner < 0) {
            winner = leader;
        }
        LOG("  Winner: "<<method_names_[winner]);
        winners_[kBucklin] = winner;
    }

    void analyzeCondorcet() noexcept {
        LOG("Results by Condorcet (wins):");

        int wins[kNumMethods];
        for (int method = 0; method < kNumMethods; ++method) {
            wins[method] = 0;
        }

        /** for every possible pairing. **/
        for (int method_a = 0; method_a < kNumMethods; ++method_a) {
            for (int method_b = method_a + 1; method_b < kNumMethods; ++method_b) {
                /** count ballots ranking a higher than b and vice versa. **/
                int count_a = 0;
                int count_b = 0;
                for (auto &&ballot : condorcet_ballots_) {
                    int rank_a = ballot.rankings_[method_a];
                    int rank_b = ballot.rankings_[method_b];
                    if (rank_a < rank_b) {
                        ++count_a;
                    } else if (rank_b < rank_a) {
                        ++count_b;
                    }
                }
                /** update wins. **/
                if (count_a > count_b) {
                    ++wins[method_a];
                } else if (count_b > count_a) {
                    ++wins[method_b];
                }

                LOG("    "<<method_names_[method_a]<<"["<<count_a<<"] vs "
                    <<method_names_[method_b]<<"["<<count_b<<"]");
            }
        }

        /** the the method that wins against all other methods. **/
        static constexpr int kNeededWins = kNumMethods - 1;
        int winner = -1;
        for (int method = 0; method < kNumMethods; ++method) {
            int w = wins[method];
            if (w >= kNeededWins) {
                winner = method;
            }
            LOG("  "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<w);
        }

        if (winner >= 0) {
            LOG("  Winner: "<<method_names_[winner]);
        } else {
            LOG("  Winner: cycle");
        }
        winners_[kCondorcet] = winner;
    }

    void analyzeGuthrie() noexcept {
        LOG("Results by Guthrie Voting (votes):");
        int winner = -1;
        int max_count = -1;
        for (int method = 0; method < kNumMethods; ++method) {
            int count = guthrie_counts_[method];
            LOG("  "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<count);
            if (max_count < count) {
                max_count = count;
                winner = method;
            }
        }
        if (max_count >= majority_) {
            LOG("  Winner: "<<method_names_[winner]);
        } else {
            winner = -1;
            LOG("  Winner: by negotiation");
        }
        winners_[kGuthrie] = winner;
    }

    void analyzeInstantRunoff() noexcept {
        LOG("Results by Instant Runoff Voting (votes):");

        bool elligible[kNumMethods];
        for (int i = 0; i < kNumMethods; ++i) {
            elligible[i] = true;
        }

        int winner = -1;
        for (int round = 1; round <= kNumRounds; ++round) {
            LOG("  Round["<<round<<"]:");

            /** count first place votes. **/
            int counts[kNumMethods];
            for (int method = 0; method < kNumMethods; ++method) {
                counts[method] = 0;
            }
            for (auto &&ballot : instant_ballots_) {
                for (int method = 0; method < kNumMethods; ++method) {
                    int rank = ballot.rankings_[method];
                    if (rank == 1) {
                        ++counts[method];
                        break;
                    }
                }
            }

            /** find the method with the mosts and fewest first place votes. **/
            winner = -1;
            int loser = -1;
            int max_count = -1;
            int min_count = 1000;
            for (int method = 0; method < kNumMethods; ++method) {
                int count = counts[method];
                if (elligible[method]) {
                    LOG("    "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<count);
                }
                if (max_count < count) {
                    max_count = count;
                    winner = method;
                }
                if (elligible[method]) {
                    if (min_count > count) {
                        min_count = count;
                        loser = method;
                    }
                }
            }

            /** the winner must have a majority. **/
            if (max_count >= majority_) {
                break;
            }

            /** remove the loser from the ballots. **/
            LOG("    Loser: "<<method_names_[loser]);
            for (auto &&ballot : instant_ballots_) {
                int loser_rank = ballot.rankings_[loser];

                /** increase the rank of all other lower ranked methods. **/
                for (int method = 0; method < kNumMethods; ++method) {
                    int rank = ballot.rankings_[method];
                    if (rank > loser_rank) {
                        --ballot.rankings_[method];
                    }
                }
                /** eliminate the loser. **/
                ballot.rankings_[loser] = kNumMethods;
                elligible[loser] = false;
            }
        }

        LOG("  Winner: "<<method_names_[winner]);
        winners_[kInstantRunoff] = winner;
    }

    void analyzeScore() noexcept {
        LOG("Results by Score Voting (total):");
        int winner = -1;
        int max_total = -1;
        for (int method = 0; method < kNumMethods; ++method) {
            int total = score_totals_[method];
            LOG("  "<<std::left<<std::setw(longest_name_)<<method_names_[method]<<": "<<total);
            if (max_total < total) {
                max_total = total;
                winner = method;
            }
        }
        LOG("  Winner: "<<method_names_[winner]);
        winners_[kScore] = winner;
    }

    void reportWinners() noexcept {
        LOG("Winners by Method:");
        for (int i = 0; i < kNumMethods; ++i) {
            int winner = winners_[i];
            if (winner >= 0) {
                LOG("  Method["<<std::left<<std::setw(longest_name_)<<method_names_[i]<<"]: "<<method_names_[winner]);
            } else {
                LOG("  Method["<<std::left<<std::setw(longest_name_)<<method_names_[i]<<"]: tbd");
            }
        }
    }
};

} // namespace

void VotingSurvey::run() noexcept {
    auto impl = new(std::nothrow) VotingSurveyImpl;
    impl->run();
    delete impl;
}
