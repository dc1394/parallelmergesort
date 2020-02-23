/*! \file parallelstablesort.cpp
	\brief スレッド並列化した安定ソートのパフォーマンスをチェックする

	Copyright © 2018 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include <algorithm>				                // for std::inplace_merge, std::stable_sort
#include <chrono>					                // for std::chrono
#include <cstdint>					                // for std::int32_t
#ifndef __clang__
    #include <execution>                            // for std::execution
#endif
#include <fstream>					                // for std::ifstream, std::ofstream
#include <iostream>					                // for std::cout, std::cerr
#include <iterator>                                 // for std::distance
#include <thread>					                // for std::thread
#include <utility>					                // for std::make_pair, std::pair
#include <vector>					                // for std::vector

#include <boost/archive/text_iarchive.hpp>          // for boost::archive::text_iarchive
#include <boost/format.hpp>                         // for boost::format
#include <boost/process.hpp>                        // for boost::process
// ReSharper disable once CppUnusedIncludeDirective
#include <boost/serialization/utility.hpp>          // for boost::serialization::access
// ReSharper disable once CppUnusedIncludeDirective
#include <boost/serialization/vector.hpp>           // for boost::serialization::access
#include <boost/thread.hpp>                         // for boost::thread

#if defined(__INTEL_COMPILER) || (__GNUC__ >= 5 && __GNUC__ < 8)
    #include <cilk/cilk.h>                          // for cilk_spawn, cilk_sync
#endif

#include <pstl/algorithm>
#include <pstl/execution>			                // for pstl::execution::par_unseq

#include <tbb/parallel_invoke.h>                    // for tbb::parallel_invoke

namespace {
    // 型エイリアス

    using mypair = std::pair<std::int32_t, std::int32_t>;

    //! A enumerated type
    /*!
        パフォーマンスをチェックする際の対象配列の種類を定義した列挙型
    */
    enum class Checktype : std::int32_t {
        // 完全にランダムなデータ
        RANDOM = 0,

        // あらかじめソートされたデータ
        SORT = 1,

        // 最初の1/4だけソートされたデータ
        QUARTERSORT = 2
    };

    //! A global variable (constant expression).
    /*!
        計測する回数
    */
    static auto constexpr CHECKLOOP = 10;

    //! A global variable (constant expression).
    /*!
        ソートする配列の要素数の最初の数
    */
    static auto constexpr N = 100;

    //! A global variable (constant).
    /*!
        再帰するかどうかの閾値
    */
    static auto const THRESHOLD = 3;
        
    //! A function.
    /*!
        並列化されたソート関数のパフォーマンスをチェックする
        \param checktype パフォーマンスをチェックする際の対象配列の種類
        \param ofs 出力用のファイルストリーム
    */
    void check_performance(Checktype checktype, std::ofstream & ofs);

    //! A function.
    /*!
        引数で与えられたstd::functionの実行時間をファイルに出力する
        \param checktype パフォーマンスをチェックする際の対象配列の種類
        \param func 実行するstd::function
        \param n 配列のサイス
        \param ofs 出力用のファイルストリーム
    */
    void elapsed_time(Checktype checktype, std::function<void(std::vector<mypair> &)> const & func, std::int32_t n, std::ofstream & ofs);

#if defined(__INTEL_COMPILER) || (__GNUC__ >= 5 && __GNUC__ < 8)
    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void stable_sort_cilk(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const len = std::distance(first, last);

        if (len <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 現在の再帰の深さが閾値以下のときだけ並列化させる
        if (reci <= THRESHOLD) {
            auto const middle = first + len / 2;

            // 下部をソート（別スレッドで実行）
            cilk_spawn stable_sort_cilk(first, middle, reci);

            // 上部をソート（別スレッドで実行）
            cilk_spawn stable_sort_cilk(middle, last, reci);

            // 二つのスレッドの終了を待機
            cilk_sync;

            // ソートされた下部と上部をマージ
            std::inplace_merge(first, middle, last);
        }
        else {
            // C++標準の安定ソートの関数を呼び出す
            std::stable_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする（Cilkで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void stable_sort_cilk(RandomIter first, RandomIter last)
    {
        // 再帰ありの並列安定ソートの関数を呼び出す
        stable_sort_cilk(first, last, 0);
    }
#endif

#if _OPENMP >= 200805
    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void stable_sort_openmp(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const len = std::distance(first, last);

        if (len <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 現在の再帰の深さが閾値以下のときだけ並列化させる
        if (reci <= THRESHOLD) {
            auto const middle = first + len / 2;

            // 次の関数をタスクとして実行
#pragma omp task
            // 下部をソート
            stable_sort_openmp(first, middle, reci);

            // 次の関数をタスクとして実行
#pragma omp task
            // 上部をソート
            stable_sort_openmp(middle, last, reci);

            // 二つのタスクの終了を待機
#pragma omp taskwait

            // ソートされた下部と上部をマージ
            std::inplace_merge(first, middle, last);
        }
        else {
            // C++標準の安定ソートの関数を呼び出す
            std::stable_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void stable_sort_openmp(RandomIter first, RandomIter last)
    {
#pragma omp parallel    // OpenMP並列領域の始まり
#pragma omp single      // task句はsingle領域で実行
        // 再帰ありの並列安定ソートの関数を呼び出す
        stable_sort_openmp(first, last, 0);
    }
#endif

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void stable_sort_tbb(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const len = std::distance(first, last);

        if (len <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 現在の再帰の深さが閾値以下のときだけ並列化させる
        if (reci <= THRESHOLD) {
            auto const middle = first + len / 2;

            // 二つのラムダ式を別スレッドで実行
            tbb::parallel_invoke(
                // 下部をソート
                [first, middle, reci]() { stable_sort_tbb(first, middle, reci); },
                // 上部をソート
                [middle, last, reci]() { stable_sort_tbb(middle, last, reci); });

            // ソートされた下部と上部をマージ
            std::inplace_merge(first, middle, last);
        }
        else {
            // C++標準の安定ソートの関数を呼び出す
            std::stable_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする（TBBで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void stable_sort_tbb(RandomIter first, RandomIter last)
    {
        // 再帰ありの並列安定ソートの関数を呼び出す
        stable_sort_tbb(first, last, 0);
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void stable_sort_thread(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const len = std::distance(first, last);

        if (len <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 現在の再帰の深さが閾値以下のときだけ並列化させる
        if (reci <= THRESHOLD) {
            auto middle = first + len / 2;

            // 下部をソート（別スレッドで実行）
            auto th1 = std::thread([first, middle, reci]() { stable_sort_thread(first, middle, reci); });

            // 上部をソート（別スレッドで実行）
            auto th2 = std::thread([middle, last, reci]() { stable_sort_thread(middle, last, reci); });

            // 二つのスレッドの終了を待機
            th1.join();
            th2.join();

            // ソートされた下部と上部をマージ
            std::inplace_merge(first, middle, last);
        }
        else {
            // C++標準の安定ソートの関数を呼び出す
            std::stable_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素を安定ソートする（std::threadで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void stable_sort_thread(RandomIter first, RandomIter last)
    {
        // 再帰ありの並列安定ソートの関数を呼び出す
        stable_sort_thread(first, last, 0);
    }

#ifdef DEBUG
    //! A template function.
    /*!
        与えられた二つのstd::vectorのすべての要素が同じかどうかチェックする
        \param v1 一つ目のstd::vector
        \param v2 二つめのstd::vector
        \return 与えられた二つのstd::vectorのすべての要素が同じならtrue、そうでないならfalse
    */
    bool vec_check(std::vector<mypair> const & v1, std::vector<mypair> const & v2);
#endif
}

int main()
{
    std::cout << "物理コア数: " << boost::thread::physical_concurrency();
    std::cout << ", 論理コア数: " << boost::thread::hardware_concurrency() << std::endl;

    std::ofstream ofsrandom("完全にシャッフルされたデータ.csv");
    std::ofstream ofssort("あらかじめソートされたデータ.csv");
    std::ofstream ofsquartersort("最初の1_4だけソートされたデータ.csv");
    
    std::cout << "完全にシャッフルされたデータを計測中...\n";
    check_performance(Checktype::RANDOM, ofsrandom);

    std::cout << "\nあらかじめソートされたデータを計測中...\n";
    check_performance(Checktype::SORT, ofssort);

    std::cout << "\n最初の1/4だけソートされたデータを計測中...\n";
    check_performance(Checktype::QUARTERSORT, ofsquartersort);

    return 0;
}

namespace {
    void check_performance(Checktype checktype, std::ofstream & ofs)
    {
#ifndef _MSC_VER
        std::array< std::uint8_t, 3 > const bom = { 0xEF, 0xBB, 0xBF };
        ofs.write(reinterpret_cast<const char *>(bom.data()), sizeof(bom));
#endif

#if defined(__INTEL_COMPILER) || (__GNUC__ >= 5 && __GNUC__ < 8)
        ofs << u8"配列の要素数,std::stable_sort,std::thread,OpenMP,TBB,CilkPlus,std::stable_sort (Parallelism TS),std::stable_sort (Parallel STLのParallelism TS)\n";
#elif defined(_MSC_VER)
        ofs << "配列の要素数,std::stable_sort,std::thread,TBB,std::stable_sort (MSVC内蔵のParallelism TS),std::stable_sort (Parallel STLのParallelism TS)\n";
#elif _OPENMP < 200805
        ofs << u8"配列の要素数,std::stable_sort,std::thread,TBB,std::stable_sort (Parallel STLのParallelism TS)\n";
#elif __clang__
        ofs << u8"配列の要素数,std::stable_sort,std::thread,OpenMP,TBB,std::stable_sort (Parallel STLのParallelism TS)\n";
#else
        ofs << u8"配列の要素数,std::stable_sort,std::thread,OpenMP,TBB,std::stable_sort (Parallelism TS),std::stable_sort (Parallel STLのParallelism TS)\n";
#endif

        auto n = N;
        for (auto i = 0; i < 7; i++) {
            for (auto j = 0; j < 2; j++) {
                std::cout << n << "個を計測中...\n";

                ofs << n << ',';

                elapsed_time(checktype, [](auto & vec) { std::stable_sort(vec.begin(), vec.end()); }, n, ofs);
                elapsed_time(checktype, [](auto & vec) { stable_sort_thread(vec.begin(), vec.end()); }, n, ofs);

#if _OPENMP >= 200805
                elapsed_time(checktype, [](auto & vec) { stable_sort_openmp(vec.begin(), vec.end()); }, n, ofs);
#endif
                elapsed_time(checktype, [](auto & vec) { stable_sort_tbb(vec.begin(), vec.end()); }, n, ofs);

#if defined(__INTEL_COMPILER) || (__GNUC__ >= 5 && __GNUC__ < 8)
                elapsed_time(checktype, [](auto & vec) { stable_sort_cilk(vec.begin(), vec.end()); }, n, ofs);
#endif

#ifndef __clang__
                elapsed_time(checktype, [](auto & vec) { std::sort(std::execution::par, vec.begin(), vec.end()); }, n, ofs);
#endif

                elapsed_time(checktype, [](auto & vec) { std::stable_sort(pstl::execution::par, vec.begin(), vec.end()); }, n, ofs);

                ofs << std::endl;

                if (i == 6) {
                    break;
                }
                else if (!j) {
                    n *= 5;
                }
            }

            n *= 2;
        }
    }

    void elapsed_time(Checktype checktype, std::function<void(std::vector<mypair> &)> const & func, std::int32_t n, std::ofstream & ofs)
    {
        using namespace std::chrono;

        std::vector< mypair > vec(n);
        auto const program_name = "makestablesortdata";

        switch (checktype) {
        case Checktype::RANDOM:
            {
                auto const filename = (boost::format("sortdata_%d_rand.dat") % n).str();
                std::ifstream ifs(filename);

                if (!ifs.is_open()) {
                    boost::process::child(program_name + (boost::format(" 0 %d") % n).str()).wait();
                    ifs.open(filename);
                }

                boost::archive::text_iarchive ia(ifs);
                ia >> vec;
            }
            break;

        case Checktype::SORT:
            {
                auto const filename = (boost::format("sortdata_%d_already.dat") % n).str();
                std::ifstream ifs(filename);

                if (!ifs.is_open()) {
                    boost::process::child(program_name + (boost::format(" 1 %d") % n).str()).wait();
                    ifs.open(filename);
                }

                boost::archive::text_iarchive ia(ifs);
                ia >> vec;
            }
            break;

        case Checktype::QUARTERSORT:
            {
                auto const filename = (boost::format("sortdata_%d_quartersort.dat") % n).str();
                std::ifstream ifs(filename);

                if (!ifs.is_open()) {
                    boost::process::child(program_name + (boost::format(" 2 %d") % n).str()).wait();
                    ifs.open(filename);
                }

                boost::archive::text_iarchive ia(ifs);
                ia >> vec;
            }
            break;

        default:
            BOOST_ASSERT(!"switchのdefaultに来てしまった！");
            break;
        }

#ifdef DEBUG
        std::vector< mypair > vecback(vec);
#endif

        auto elapsed_time = 0.0;
        for (auto i = 1; i <= CHECKLOOP; i++) {
            auto const beg = high_resolution_clock::now();
            func(vec);
            auto const end = high_resolution_clock::now();

            elapsed_time += (duration_cast<duration<double>>(end - beg)).count();
        }

        ofs << boost::format(u8"%.10f") % (elapsed_time / static_cast<double>(CHECKLOOP)) << ',';

#ifdef DEBUG
        std::stable_sort(pstl::execution::par_unseq, vecback.begin(), vecback.end());

        if (!vec_check(vec, vecback)) {
            std::cerr << "エラー発見！" << std::endl;
        }
#endif
    }
    
#ifdef DEBUG
    bool vec_check(std::vector<mypair> const & v1, std::vector<mypair> const & v2)
    {
        auto const size = v1.size();
        BOOST_ASSERT(size == v2.size());

        for (auto i = 0UL; i < size; i++) {
            if (v1[i] != v2[i]) {
                std::cerr << "Error! i = " << i << std::endl;
                return false;
            }
        }

        return true;
    }
#endif
}
