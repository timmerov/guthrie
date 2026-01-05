/*
Copyright (C) 2012-2025 tim cotter. All rights reserved.
*/

#include <algorithm>
#include <curl/curl.h>
#include <iomanip>
#include <sstream>

#include <aggiornamento/aggiornamento.h>
#include <aggiornamento/log.h>

#include "survey.h"

namespace {

static constexpr char kVotingSurveyUrl[] = "https://docs.google.com/spreadsheets/d/11CS8R4pbYFDQcHaXjGxAquIuVWzb8vjciY_fS-Tz_Rk/export?format=tsv&gid=2041703740";

enum class Method {
    kApproval,
    kBordaCount,
    kBucklin,
    kCondorcet,
    kGuthrie,
    kInstantRunoff,
    kScore
};

class VotingSurveyImpl {
public:
    VotingSurveyImpl() = default;
    ~VotingSurveyImpl() = default;

    std::string buffer_;

    void run() noexcept {
        loadSurveyDataFromGoogleSheet();

        std::stringstream ss(std::move(buffer_));
        buffer_.clear();

        /** skip the first line. **/
        std::string line;
        std::getline(ss, line);

        /** read the rows. **/
        int row = 0;
        while (std::getline(ss, line)) {
            LOG("row["<<row<<"]:");
            ++row;

            /** replace spaces with underscores. **/
            std::replace(line.begin(), line.end(), ' ', '_');

            /** replace empty cells with underscores. **/
            for(;;) {
                auto found = line.find("\t\t");
                if (found == std::string::npos) {
                    break;
                }
                line.replace(found, 2, "\t_\t");
            }

            /** read cells. **/
            std::stringstream cells(line);
            std::string s0;
            std::string s1;
            std::string s2;
            std::string s3;
            std::string s4;
            std::string s5;
            std::string s6;

            /** date/time **/
            cells >> s0;

            /** approval **/
            cells >> s0;
            LOG("Approval: "<<s0);

            /** borda count **/
            cells >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;
            LOG("Borda Count: "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s5<<" "<<s6<<" ");

            /** bucklin **/
            cells >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;
            LOG("Bucklin: "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s5<<" "<<s6<<" ");

            /** condorcet **/
            cells >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;
            LOG("Condorcet: "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s5<<" "<<s6<<" ");

            /** guthrie **/
            cells >> s0;
            LOG("Guthrie: "<<s0);

            /** instant runoff **/
            cells >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;
            LOG("Instant Runoff: "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s5<<" "<<s6<<" ");

            /** score **/
            cells >> s0 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;
            LOG("Score: "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s5<<" "<<s6<<" ");
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
};

} // namespace

void VotingSurvey::run() noexcept {
    auto impl = new(std::nothrow) VotingSurveyImpl;
    impl->run();
    delete impl;
}
