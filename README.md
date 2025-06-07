
Implement Guthrie Voting in a simulated electorate.
Guthrie voting defined and evangelized here:
https://docs.google.com/document/d/1GL__lJMoX5Cku35h4BLXhJHQ_NxuzGaA5tN-OORVdmw/edit?tab=t.0

Loosely speaking, Guthrie voting is any system where:
1) Voters cast an anonymous ballot with a single vote for their favorite candidate.
2) if any candidate has a majority they win.
3) Otherwise, the candidates negotiate a winner by a set of formal rules where:
4) There are one or more rounds of negotiation.
5) Each round, candidates declare how they intend to vote.
6) Candidates may publicly change their votes.
7) A final tally is made and a winner is determined when there are no more changes, i.e. when a Nash equilibrium has been achieved.

For the most part, the set of formal rules for the contingent election doesn’t matter too much.
The candidates will informally agree to unwritten rules that more/less follow Coombs' method.
Coombs' method is used in the code.

In Coombs's method, there are a series of voting rounds alternating between:
A) Voting for a majority winner.
If a candidate has a majority, they win.
B) Voting to eliminate a candidate from the ballot.
Eliminated candidates participate in subsequent rounds to choose the winner.
A candidate has voting power (asset) equal to the number of votes they got from the electorate.


How to build and run:

$ ./gen-unix64.sh
$ ./make-unix64.sh
$ ./bin/guthrie

Configure the parameters by changing these constants at the start of the file ./src/guthrie.cc.

kNTrials
kNVoters
kNCandidates
kElectorateMethod
kNClusters
kNAxes
kAxisWeightDecay
kCandidateMethod
kPrimaryPower
kSeedChoice

Configure the output by change these constants at the start of the file ./src/guthrie.cc.

kFindTheoreticalBestCandidate
kShowElectorateDistribution
kShowVoterBlocs
kShowCoombsRounds
