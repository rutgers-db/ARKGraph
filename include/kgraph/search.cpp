/* 
    Copyright (C) 2010,2011 Wei Dong <wdong.pku@gmail.com>. All Rights Reserved.

    DISTRIBUTION OF THIS PROGRAM IN EITHER BINARY OR SOURCE CODE FORM MUST BE
    PERMITTED BY THE AUTHOR.
*/
#ifndef KGRAPH_VALUE_TYPE
#define KGRAPH_VALUE_TYPE float
#endif


#include <cctype>
#include <type_traits>
#include <iostream>
#define timer timer_for_boost_progress_t
#include <boost/progress.hpp>
#undef timer
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include <kgraph.h>
#include <kgraph-data.h>

using namespace std;
using namespace boost;
using namespace kgraph;
namespace po = boost::program_options; 

typedef KGRAPH_VALUE_TYPE value_type;

int main (int argc, char *argv[]) {
    vector<string> params;
    string input_path;
    string index_path;
    string query_path;
    string output_path;
    string init_path;
    string eval_path;
    unsigned K, M, P, S, T;

    po::options_description desc_visible("General options");
    desc_visible.add_options()
    ("help,h", "produce help message.")
    ("version,v", "print version number.")
    ("data", po::value(&input_path), "input path")
    ("index", po::value(&index_path), "index path")
    ("query", po::value(&query_path), "query path")
    ("output", po::value(&output_path), "output path")
    ("init", po::value(&init_path), "init path")
    ("eval", po::value(&eval_path), "eval path")
    (",K", po::value(&K)->default_value(default_K), "")
    (",M", po::value(&M)->default_value(default_M), "")
    (",P", po::value(&P)->default_value(default_P), "")
    (",S", po::value(&S)->default_value(default_S), "")
    (",T", po::value(&T)->default_value(default_T), "")
    ("linear", "")
    ("l2norm", "l2-normalize data, so as to mimic cosine similarity")
    ;

    po::options_description desc("Allowed options");
    desc.add(desc_visible);

    po::positional_options_description p;
    p.add("data", 1);
    p.add("index", 1);
    p.add("query", 1);
    p.add("output", 1);

    po::variables_map vm; 
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm); 

    if (vm.count("version")) {
        cout << "KGraph version " << KGraph::version() << endl;
        return 0;
    }

    if (vm.count("help") || vm.count("data") == 0 || vm.count("index") == 0 || vm.count("query") == 0) {
        cout << "search <data> <index> <query> [output]" << endl;
        cout << desc_visible << endl;
        return 0;
    }

    if (P < K) {
        P = K;
    }

    Matrix<value_type> data;
    Matrix<value_type> query;
    Matrix<unsigned> result; //(query.size(), U);
    unsigned init = 0;

    data.load_lshkit(input_path);
    query.load_lshkit(query_path);
    if (init_path.size()) {
        result.load_lshkit(init_path);
        BOOST_VERIFY(result.size() == query.size());
        init = result.dim();
        BOOST_VERIFY(init >= K);
    }
    if (vm.count("l2norm")) {
        cerr << "L2-normalizing data..." << endl;
        data.normalize2();
        query.normalize2();
    }
    MatrixOracle<value_type, metric::l2sqr> oracle(data);
    float recall = 0;
    float cost = 0;
    float time = 0;
    if (vm.count("linear")) {
        result.resize(query.size(), K);
        boost::progress_display progress(query.size(), cerr);
        {   // to ensure auto_cpu_timer destruct before "time" computation
            boost::timer::auto_cpu_timer timer;
#pragma omp parallel for
            for (unsigned i = 0; i < query.size(); ++i) {
                oracle.query(query[i]).search(K, default_epsilon, result[i]);
#pragma omp critical
                ++progress;
            }
            time = timer.elapsed().wall / 1e9;
        }
        cost = 1.0;
    }
    else {
        result.resize(query.size(), K);
        KGraph::SearchParams params;
        params.K = K;
        params.M = M;
        params.P = P;
        params.S = S;
        params.T = T;
        params.init = init;
        KGraph *kgraph = KGraph::create();
        kgraph->load(index_path.c_str());
        cerr << "Searching..." << endl;
        boost::progress_display progress(query.size(), cerr);
        {
            boost::timer::auto_cpu_timer timer;
#pragma omp parallel for reduction(+:cost)
            for (unsigned i = 0; i < query.size(); ++i) {
                KGraph::SearchInfo info;
                kgraph->search(oracle.query(query[i]), params, result[i], &info);
#pragma omp critical
                ++progress;
                cost += info.cost;
            }
            time = timer.elapsed().wall / 1e9;
        }
        cost /= query.size();
        //cerr << "Cost: " << cost << endl;
        delete kgraph;
    }
    if (output_path.size()) {
        result.save_lshkit(output_path);
    }
    if (eval_path.size()) {
        Matrix<unsigned> gs;
        gs.load_lshkit(eval_path);
        BOOST_VERIFY(gs.dim() >= K);
        BOOST_VERIFY(gs.size() >= query.size());
        kgraph::Matrix<float> gs_dist(query.size(), K);
        kgraph::Matrix<float> result_dist(query.size(), K);
#pragma omp parallel for
        for (unsigned i = 0; i < query.size(); ++i) {
            auto Q = oracle.query(query[i]);
            float *gs_dist_row = gs_dist[i];
            float *result_dist_row = result_dist[i];
            unsigned const *gs_row = gs[i];
            unsigned const *result_row = result[i];
            for (unsigned k = 0; k < K; ++k) {
                gs_dist_row[k] = Q(gs_row[k]);
                result_dist_row[k] = Q(result_row[k]);
            }
            sort(gs_dist_row, gs_dist_row + K);
            sort(result_dist_row, result_dist_row + K);
        }
        recall = AverageRecall(gs_dist, result_dist, K);
    }
    cout << "Time: " << time << " Recall: " << recall << " Cost: " << cost << endl;

    return 0;
}
