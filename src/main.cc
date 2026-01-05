/*
Copyright (C) 2012-2026 tim cotter. All rights reserved.
*/

/**
guthrie voting.
**/

#include "guthrie.h"
#include "survey.h"

#include <aggiornamento/aggiornamento.h>
#include <aggiornamento/log.h>

int main(
    int argc, char *argv[]
) noexcept {
    (void) argc;
    (void) argv;

    agm::log::init(AGM_TARGET_NAME ".log", false);

    VotingSurvey survey;
    survey.run();

    /*GuthrieVoting guthrie;
    guthrie.run();*/

    return 0;
}
