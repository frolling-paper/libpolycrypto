#include <polycrypto/PolyOps.h>
#include <polycrypto/PolyCrypto.h>
#include <polycrypto/RootsOfUnityEval.h>
#include <polycrypto/Dkg.h>
#include <polycrypto/Lagrange.h>
#include <polycrypto/FFThresh.h> // for random_subset()

#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>

#include <xutils/Log.h>
#include <xutils/Utils.h>
#include <xutils/Timer.h>
#include <xassert/XAssert.h>

#include <polycrypto/Utils.h>

using namespace std;

void testScheme(const std::vector<Dkg::AbstractPlayer*>& players, bool isDkg) {
    if(isDkg) {
        // deal
        for(auto p : players) {
            //logdbg << "Player #" << p->id << " is dealing..." << endl;
            p->deal();
        }

        // incorporate deals from other players
        for(auto p : players) {
            //logdbg << "Player #" << p->id << " is verifying other players..." << endl;
            testAssertTrue(p->verifyOtherPlayers(players, false));
            testAssertTrue(p->verifyOtherPlayers(players, true));
        }

        // verify each player's share of the final agreed-upon secret
        for(auto p : players) {
            //logdbg << "Player #" << p->id << " is verifying his own final share..." << endl;
            testAssertTrue(p->verifyFinalShareProof());
        }
    } else {
        Dkg::AbstractPlayer* dealer = players[0];

        // test if dealing twice works, since our benchmarks need this
        dealer->deal();
        dealer->deal();

        for(auto p : players) {
            //logdbg << "Player #" << p->id << " is verifying other players..." << endl;
            testAssertTrue(p->verifyOtherPlayers({ dealer }, false));
        }
    }

    auto reconstructor = players[0];
    // WARNING: We do not actually test the reconstruction of the final DKG secret, just of the reconstructing
    // player's secret, since performance will be the same but logic is simplified
    Fr secret = reconstructor->getSecret();
    size_t n = reconstructor->params.n;
    size_t t = reconstructor->params.t;
    std::vector<size_t> subset;

    // Worst case: you get n shares and you have to verify all of them before you find a subset of t correct ones 
    subset = random_subset(t, n);
    std::sort(subset.begin(), subset.end());
    testAssertTrue(reconstructor->reconstructionVerify(subset, false));
    Fr r1 = reconstructor->interpolate(subset);

    // Best case: you get all correct shares; here we give the reconstructor exactly t correct shares
    subset = random_subset(t, n);
    std::sort(subset.begin(), subset.end());
    testAssertTrue(reconstructor->reconstructionVerify(subset, true));
    Fr r2 = reconstructor->interpolate(subset);

    testAssertEqual(r1, secret);
    testAssertEqual(r2, secret);
    logdbg << "Secret: " << secret << endl;
}

void testFeldman(const Dkg::DkgParams& params, const Dkg::FeldmanPublicParameters& fpp, bool isDkgPlayer) {
    std::vector<Dkg::AbstractPlayer*> players;
    bool isSimulated = false;

    for(size_t i = 0; i < params.n; i++) {
        //logdbg << "Initializing player #" << i+1 << endl;
        players.push_back(new Dkg::FeldmanPlayer(params, fpp, i, isSimulated, isDkgPlayer));
    }

    testScheme(players, isDkgPlayer);
}

template<class PlayerType>
void testKateSimScheme(const Dkg::DkgParams& params, const Dkg::KatePublicParameters& kpp, bool isDkgPlayer) {
    std::vector<Dkg::AbstractPlayer*> players;
        
    // only VSS dealer is simulated
    players.push_back(new PlayerType(
        params, kpp, 0, 
        true,
        isDkgPlayer));

    for(size_t i = 1; i < params.n; i++) {
        //logdbg << "Initializing player #" << i << endl;
        // NOTE: for Kate, we test simulated players too
        players.push_back(new PlayerType(
            params, kpp, i, 
            false,
            isDkgPlayer));
    }

    testScheme(players, isDkgPlayer);
}

template<class PlayerType>
void testKateBasedScheme(const Dkg::DkgParams& params, const Dkg::KatePublicParameters& kpp, bool isDkgPlayer, bool shouldSimulate) {
    std::vector<Dkg::AbstractPlayer*> players;

    for(size_t i = 0; i < params.n; i++) {
        //logdbg << "Initializing player #" << i << endl;
        // NOTE: for Kate, we test simulated players too
        bool kateIsSimulated = (i % 2 ? false : true);
        players.push_back(new PlayerType(
            params, kpp, i, 
            shouldSimulate ? kateIsSimulated : false,
            isDkgPlayer));
    }

    testScheme(players, isDkgPlayer);
}

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    libpolycrypto::initialize(nullptr, 0);
    srand(static_cast<unsigned int>(time(NULL)));

    //size_t minT = 32;
    size_t maxT = 1<<20;
    size_t maxN = 2*maxT - 1;

    Dkg::KatePublicParameters kpp = Dkg::KatePublicParameters::getRandom(maxT - 1);

    cout << "Computing accumulators..." << endl;
    AccumulatorTree accs(maxN);
    AuthAccumulatorTree authAccs(accs, kpp, maxT);
    size_t t ;
    int i =5;
    for(t=32;t<= maxT;t*=2){
        size_t n = 2*t-1;
        Dkg::DkgParams pp(t,n,1);
        pp.setAuthAccumulators(&authAccs);
        // Dkg::KatePlayer dealer_k(pp,kpp,0,0,0);
        Dkg::KatePlayer dealer_k(pp,kpp,0,0,0);
        // cout<< dealer_k.id<<endl;
        dealer_k.f_id=libpolycrypto::random_field_elems(pp.t);
        cout<<"now commit"<<endl;

        auto startTime1 = std::chrono::high_resolution_clock::now();
        dealer_k.comm=multiExp<G1>(
            kpp.g1si.begin(),
            kpp.g1si.begin() + static_cast<long>(dealer_k.f_id.size()),
            dealer_k.f_id.begin(),
            dealer_k.f_id.end()
        );
        auto endTime1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms1 = endTime1 - startTime1;

        std::cout << "commit time when" << "t = 2^" << i << ": " << fp_ms1.count() << std::endl;
        
        cout<<"now evaluate"<<endl;
        auto startTime2 = std::chrono::high_resolution_clock::now();
        dealer_k.evaluate();
        auto endTime2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms2 = endTime2 - startTime2;
        std::cout << "evaluate time when" << "t=2^" << i << ":  " << fp_ms2.count() << std::endl;
        cout<<"now proof"<<endl;
        auto startTime3 = std::chrono::high_resolution_clock::now();
        dealer_k.computeRealProofs();
        auto endTime3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms3 = endTime3 - startTime3;
        std::cout << "proof time when" << "t = 2^" << i << ":  " << fp_ms3.count() << std::endl;
        cout<<"now verify"<<endl;
        auto startTime4 = std::chrono::high_resolution_clock::now();
        dealer_k.verifyShareFromDealer(dealer_k);
        auto endTime4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms4 = endTime4 - startTime4;
        std::cout << "verify time when" << "t=2^" << i << ":   " << fp_ms4.count() << std::endl;
        i++;
    }


    


    // for(size_t t = minT; t <= maxT; t*=2) {
    //     size_t n = 2*t-1;
    //         Dkg::DkgParams params(t, n, true);
    //         //Dkg::FeldmanPublicParameters fpp(params);
    //         params.setAuthAccumulators(&authAccs);

    //         for(bool isDkg : { false }) {
    //             cout << "Simulating " << t << " out of " << n << " " << (isDkg ? "DKG" : "VSS") << " ..." << endl;

    //             //cout << " * Feldman..." << endl;
    //             //testFeldman(params, fpp, isDkg);

    //             cout << " * Kate..." << endl;
    //             testKateBasedScheme<Dkg::KatePlayer>(params, kpp, isDkg, true);

    //             //cout << " * KateSim..." << endl;
    //             //testKateSimScheme<Dkg::KatePlayer>(params, kpp, isDkg);

    //             cout << " * AMT..." << endl;
    //             testKateBasedScheme<Dkg::MultipointPlayer>(params, kpp, isDkg, false);

    //             //cout << " * AMTSim..." << endl;
    //             //testKateSimScheme<Dkg::MultipointPlayer>(params, kpp, isDkg);

    //             cout << endl;
    //         }
        
    // }

    cout << "All tests succeeded!" << endl;

    return 0;
}
