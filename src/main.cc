/*
Copyright (C) 2012-2025 tim cotter. All rights reserved.
*/

/**
guthrie voting.
**/

#include "guthrie.h"

#include <aggiornamento/aggiornamento.h>
#include <aggiornamento/log.h>

int main(
    int argc, char *argv[]
) noexcept {
    (void) argc;
    (void) argv;

    agm::log::init(AGM_TARGET_NAME ".log", false);

    GuthrieVoting guthrie;
    guthrie.run();

    return 0;
}
