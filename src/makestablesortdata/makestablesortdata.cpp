#include <cstdint>				                    // for std::int32_t
#include <cstdlib>                                  // for EXIT_FAILURE, EXIT_SUCCESS
#include <fstream>                                  // for std::ofstream
#include <random>				                    // for std::mt19937, std::random_device
#include <vector>				                    // for std::vector

#include <boost/archive/text_oarchive.hpp>          // for boost::archive::text_oarchive
#include <boost/format.hpp>		                    // for boost::format
#include <boost/serialization/utility.hpp>          // for boost::serialization::access
#include <boost/serialization/vector.hpp>           // for boost::serialization::access

#include <pstl/algorithm>
#include <pstl/execution>                           // for pstl::execution::par

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

    //! A function.
    /*!
        ソート対象のデータを生成する
        \param checktype 生成するデータの種類
        \param n 生成するデータの個数
        \return 生成が成功したかどうか
    */
	bool make_sortdata(Checktype checktype, std::int32_t n);
}

int main(int argc, char * argv[])
{
    if (argc != 3) {
        std::exit(EXIT_FAILURE);
    }

    auto const type = static_cast<Checktype>(std::stoi(argv[1]));
    auto const n = std::stoi(argv[2]);

    switch (type) {
    case Checktype::RANDOM:
    {
        if (!make_sortdata(Checktype::RANDOM, n)) {
            std::exit(EXIT_FAILURE);
        }
    }
    break;

    case Checktype::SORT:
    {
        if (!make_sortdata(Checktype::SORT, n)) {
            std::exit(EXIT_FAILURE);
        }
    }
    break;

    case Checktype::QUARTERSORT:
    {
        if (!make_sortdata(Checktype::QUARTERSORT, n)) {
            std::exit(EXIT_FAILURE);
        }
    }
    break;

    default:
        BOOST_ASSERT(!"switchのdefaultに来てしまった！");
        break;
    }

    return EXIT_SUCCESS;
}

namespace {
	bool make_sortdata(Checktype checktype, std::int32_t n)
	{
        // ランダムデバイス
        std::random_device rnd;

        // 乱数エンジン
        auto randengine = std::mt19937(rnd());

        // 乱数の分布（1～n / 10までの一様分布）
        std::uniform_int_distribution<std::int32_t> distribution(1, n / 10);

		std::vector< mypair > vec(n);
        for (auto j = 0; j < n; j++) {
            vec[j] = std::make_pair(distribution(randengine), j);
        }

		switch (checktype) {
		case Checktype::RANDOM:
		    {
			    std::ofstream ofs((boost::format("sortdata_%d_rand.dat") % n).str());
                if (!ofs.is_open()) {
                    return false;
                }

                boost::archive::text_oarchive oa(ofs);
                oa << vec;
		    }
		    break;

		case Checktype::SORT:
    		{
                std::stable_sort(pstl::execution::par_unseq, vec.begin(), vec.end());
            
                std::ofstream ofs((boost::format("sortdata_%d_already.dat") % n).str());
                if (!ofs.is_open()) {
                    return false;
                }
                
                boost::archive::text_oarchive oa(ofs);
                oa << vec;
            }
            break;

		case Checktype::QUARTERSORT:
    		{
                std::stable_sort(pstl::execution::par_unseq, vec.begin(), vec.begin() + n / 4);
                
                std::ofstream ofs((boost::format("sortdata_%d_quartersort.dat") % n).str());
                if (!ofs.is_open()) {
                    return false;
                }

                boost::archive::text_oarchive oa(ofs);
                oa << vec;
            }
            break;

		default:
			BOOST_ASSERT(!"switchのdefaultに来てしまった！");
			break;
		}

        return true;
	}
}
