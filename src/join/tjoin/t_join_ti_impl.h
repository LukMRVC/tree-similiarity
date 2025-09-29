// The MIT License (MIT)
// Copyright (c) 2017 Thomas Huetter
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// Implements the TJoin tree similarity join.

#pragma once
#include "omp.h"

template <typename Label, typename VerificationAlgorithm>
TJoinTI<Label, VerificationAlgorithm>::TJoinTI() {
    ld_ = label::LabelDictionary<Label>();
    pre_candidates_ = 0;
    sum_subproblem_counter_ = 0;
    number_of_labels_ = 0;
    il_lookups_ = 0;
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::execute_join(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    // Convert trees to sets and get the result.
    convert_trees_to_sets(trees_collection, sets_collection);
    std::cout << "Trees converted to sets!\n";

    // Retrieves candidates from the candidate index.
    retrieve_candidates(sets_collection, candidates, distance_threshold);

    // Use the label guided mapping upper bound to send candidates immediately .
    upperbound(trees_collection, candidates, join_result, distance_threshold);

    // Verify all computed join candidates and return the join result.
    verify_candidates(trees_collection, candidates, join_result,
                      distance_threshold);
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::convert_trees_to_sets(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection) {

    // Convert trees to sets and get the result.
    label_set_converter::Converter<Label> lsc;
    lsc.assignFrequencyIdentifiers(trees_collection, sets_collection);
    number_of_labels_ = lsc.get_number_of_labels();
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::retrieve_candidates(
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection,
    std::vector<std::pair<int, int>>& candidates,
    const double distance_threshold) {

    // Initialize candidate index.
    candidate_index::CandidateIndex c_index;

    // Retrieve candidates from the candidate index.
    c_index.lookup(sets_collection, candidates, number_of_labels_,
                   distance_threshold);

    // Copy the number of pre-candidates.
    pre_candidates_ = c_index.get_number_of_pre_candidates();
    // Copy the number of inverted list lookups.
    il_lookups_ = c_index.get_number_of_il_lookups();
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::retrieve_candidates(
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection,
    std::vector<std::pair<int, int>>& candidates,
    const double distance_threshold,
    std::vector<std::chrono::microseconds>& ted_times) {

    // Initialize candidate index.
    candidate_index::CandidateIndex c_index;

    // Retrieve candidates from the candidate index.
    c_index.lookup(sets_collection, candidates, number_of_labels_,
                   distance_threshold, ted_times);

    // Copy the number of pre-candidates.
    pre_candidates_ = c_index.get_number_of_pre_candidates();
    // Copy the number of inverted list lookups.
    il_lookups_ = c_index.get_number_of_il_lookups();
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::upperbound(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    typename VerificationAlgorithm::AlgsCostModel cm(ld_);
    auto total_candidates_to_verify = candidates.size();
    std::atomic<int> i = 0;

    std::vector<std::vector<std::pair<int, int>>>
        per_thread_candidates_to_keep(omp_get_max_threads());
    std::vector<std::vector<join::JoinResultElement>> per_thread_join_results(
        omp_get_max_threads());

#pragma omp parallel for schedule(guided, 1000) shared(ld_, cm, trees_collection, candidates, per_thread_candidates_to_keep, per_thread_join_results)
    for (const auto& pair : candidates) {
        ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
            lgm_algorithm(cm);
        node::TreeIndexLGM ti_1;
        node::TreeIndexLGM ti_2;

        node::index_tree(ti_1, trees_collection[pair.first], ld_, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld_, cm);
        double ub_value = lgm_algorithm.ted_k(ti_1, ti_2, distance_threshold);

        int current_count = i.fetch_add(1, std::memory_order_relaxed);
        if ((current_count + 1) % 10'000 == 0) {
#pragma omp critical
            {
                std::cout << "UB Verified " << (current_count + 1)
                          << " out of " << total_candidates_to_verify << "\n";
            }
        }

        if (ub_value <= distance_threshold) {
            per_thread_join_results[omp_get_thread_num()].emplace_back(
                pair.first, pair.second, ub_value);
        } else {
            per_thread_candidates_to_keep[omp_get_thread_num()].push_back(pair);
        }
    }

    // Clear original candidates and join_result before merging
    candidates.clear();
    // join_result is appended to, so no need to clear.

    // Merge results from all threads
    for (const auto& thread_results : per_thread_join_results) {
        join_result.insert(join_result.end(), thread_results.begin(),
                           thread_results.end());
    }
    for (const auto& thread_candidates : per_thread_candidates_to_keep) {
        candidates.insert(candidates.end(), thread_candidates.begin(),
                          thread_candidates.end());
    }
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::verify_candidates(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    typename VerificationAlgorithm::AlgsCostModel cm(ld_);
    auto i = 0;
    // Verify each pair in the candidate set

    std::vector<std::vector<join::JoinResultElement>> per_each_thread(
        omp_get_max_threads());

#pragma omp parallel for schedule(guided, 5000) shared(ld_, cm, per_each_thread)
    for (const auto& pair : candidates) {
        {
#pragma omp atomic
            i += 1;
        }
        VerificationAlgorithm ted_algorithm(cm);
        typename VerificationAlgorithm::AlgsTreeIndex ti_1;
        typename VerificationAlgorithm::AlgsTreeIndex ti_2;

        if (i % 10'000 == 0) {
            std::cout << "Verified " << i << " out of " << candidates.size()
                      << "\n";
        }
        node::index_tree(ti_1, trees_collection[pair.first], ld_, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld_, cm);
        double ted_value = ted_algorithm.ted_k(ti_1, ti_2, distance_threshold);

        if (ted_value <= distance_threshold) {
            // #pragma omp critical
            //             join_result.emplace_back(pair.first, pair.second,
            //             ted_value);
            per_each_thread[omp_get_thread_num()].emplace_back(
                pair.first, pair.second, ted_value);
        }

        // Sum up all number of subproblems
        sum_subproblem_counter_ += ted_algorithm.get_subproblem_count();
    }

    for (auto& res : per_each_thread) {
        join_result.insert(join_result.end(), res.begin(), res.end());
    }
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_number_of_pre_candidates() const {
    return pre_candidates_;
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_subproblem_count() const {
    return sum_subproblem_counter_;
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_number_of_il_lookups() const {
    return il_lookups_;
}
