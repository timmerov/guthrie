/*
Copyright (C) 2012-2026 tim cotter. All rights reserved.
*/

/**
analyze data from a web survey.
**/

#pragma once

class VotingSurvey {
public:
    VotingSurvey() noexcept = default;
    ~VotingSurvey() noexcept = default;

    void run() noexcept;

private:
    void *impl_ = nullptr;
};
