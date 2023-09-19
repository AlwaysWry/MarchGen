#include "mwvc.h"
#include <sstream>

extern "C" __declspec(dllexport) int MWVC(char graph[], char result_file[], char seed_num[], char cutoff_time_num[], char mode_num[]);


extern "C" __declspec(dllexport) int MWVC(char graph[], char result_file[], char seed_num[], char cutoff_time_num[], char mode_num[]) {

    if (BuildInstance(graph)) {
        cerr << "Open instance graph failed." << endl;
        return 1;
    }

    stringstream ss;
    ss << seed_num;
    ss >> seed;
    ss.clear();
    ss << cutoff_time_num;
    ss >> cutoff_time;
    ss.clear();
    ss << mode_num;
    ss >> mode;
    ss.clear();

    //default seed
    if (seed < 0U || seed > ~0U) {
        seed = 10;
    }

    //default cutoff_time
    if (cutoff_time < 0 || cutoff_time > (int) (~0U >> 1)) {
        cutoff_time = 10;
    }

    if (mode < 0 || mode > 3) {
        mode = 0;
    }

    srand(seed);

    cout << "Parsing graph file " << graph << "..." << endl;

    start = chrono::steady_clock::now();

    ConstructVC();
    LocalSearch();

    if (CheckSolution() == 1) {
        cout << "Solve finished." << endl;
    } else {
        cout << ", the solution is wrong." << endl;
    }

    int weight_sum = 0;
    int original_weight = 0;

    //write results to mwvc_log.txt
    ofstream result(result_file);

    if (!result){
        cerr << "Open output file failed." << endl;
        return 1;
    }

    for (int i = 1; i <= v_num; i++) {
        if (best_v_in_c[i] == 1) {
            //cout << i << ',' << v_weight[i] << endl;
            result << i << endl;
            weight_sum += v_weight[i];
        }

        original_weight += v_weight[i];
    }

    result.close();

    cout << "best size: " << best_c_size << ", " << endl;
    cout << "best weight: " << weight_sum << ", " << endl;
    cout << "reduction: "<< original_weight - weight_sum << endl;

    FreeMemory();

    return 0;
}
