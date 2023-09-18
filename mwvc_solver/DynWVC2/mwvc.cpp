#include "mwvc.h"
#include <sstream>
extern "C" __declspec(dllexport) int DynWVC2(int arg_num, char file[], char seed_num[], char cutoff_time_num[], char mode_num[]);



extern "C" __declspec(dllexport) int DynWVC2(int arg_num, char file[], char seed_num[], char cutoff_time_num[], char mode_num[]) {
    uint seed;

    if (arg_num == 1) {
        cout << "Minimum Weighted Vertex Cover Problem solver." << endl;
        cout << "Usage: ./mwvc_solver [Graph file] [Seed] [Cutoff time] [CC mode]" << endl;
        return 1;
    }

    if (arg_num < 4) {
        cerr << "Missing argument(s)." << endl;
        cout << "Usage: ./mwvc_solver [Graph file] [Seed] [Cutoff time] [CC mode]" << endl;
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

    if (BuildInstance(file) != 0) {
        cerr << "Open instance file failed." << endl;
        return 1;
    }

    if (seed < 0U || seed > ~0U) {
        seed = 10;
    }

    if (cutoff_time < 0 || cutoff_time > (int) (~0U >> 1)) {
        cutoff_time = 1000;
    }

    if (mode < 0 || mode > 3) {
        mode = 0;
    }

    srand(seed);

    cout << file;

    start = chrono::steady_clock::now();

    ConstructVC();
    LocalSearch();

    if (CheckSolution() == 1) {
        cout << "best_w:" << best_weight << " best_time: " << best_comp_time << endl;
    } else {
        cout << ", the solution is wrong." << endl;
    }

    int weight_sum = 0;
    int original_weight = 0;

    for (int i = 1; i <= v_num; i++) {
        if (best_v_in_c[i] == 1) {
            cout << i << ',' << v_weight[i] << endl;
            weight_sum += v_weight[i];
        }

        original_weight += v_weight[i];
    }

    cout << best_c_size << "," << weight_sum << "," << original_weight << endl;

    FreeMemory();

    return 0;
}
