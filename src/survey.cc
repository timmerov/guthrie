/*
Copyright (C) 2012-2025 tim cotter. All rights reserved.
*/

#include <curl/curl.h>

#include "survey.h"

#include <aggiornamento/aggiornamento.h>
#include <aggiornamento/log.h>

namespace {

static constexpr char kVotingSurveyUrl[] = "https://docs.google.com/spreadsheets/d/11CS8R4pbYFDQcHaXjGxAquIuVWzb8vjciY_fS-Tz_Rk/export?format=csv&gid=2041703740";

class VotingSurveyImpl {
public:
    VotingSurveyImpl() = default;
    ~VotingSurveyImpl() = default;

    std::string buffer_;

    void run() noexcept {
        loadSurveyDataFromGoogleSheet();

        LOG("csv="<<buffer_);
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
