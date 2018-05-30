// The MIT License (MIT)
// Copyright (c) 2018 Mateusz Pawlik.
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

/// \file ted_ub/greedy_ub_impl.h
///
/// \details
/// Contains the implementation of the greedy mapping algorithm.

#ifndef TREE_SIMILARITY_TED_UB_GREEDY_UB_IMPL_H
#define TREE_SIMILARITY_TED_UB_GREEDY_UB_IMPL_H

template <typename Label, typename CostModel>
GreedyUB<Label, CostModel>::GreedyUB() : c_() {}

template <typename Label, typename CostModel>
double GreedyUB<Label, CostModel>::verify(const node::Node<Label>& t1,
    const node::Node<Label>& t2, double similarity_threshold) {
  double cost = mapping_cost(lb_mapping_fill_gaps(t1, t2, similarity_threshold));
  if (cost <= static_cast <int> (std::ceil(similarity_threshold))) {
    return cost;
  }
  return std::numeric_limits<double>::infinity();
};

template <typename Label, typename CostModel>
double GreedyUB<Label, CostModel>::mapping_cost(
    const std::vector<std::pair<int, int>>& mapping) const {
  double result = 0.0;
  for (auto m : mapping) {
    result += c_.ren(t1_node_[m.first], t2_node_[m.second]); // renames
  }
  result += t1_input_size_ - mapping.size(); // deletions
  result += t2_input_size_ - mapping.size(); // insertions
  return result;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::lb_mapping(
    const node::Node<Label>& t1, const node::Node<Label>& t2, const int k) {
  init(t1, t2);
  std::vector<std::pair<int, int>> mapping;
  int cand_id = 0;
  int start_pos = 0;
  int end_pos = 0;
  int pos = 0;
  int t1_label_id = 0;
  for (int i = 0; i < t1_input_size_; ++i) { // Loop in postorder.
    std::vector<int>& candidate_ids = t2_label_il_[t1_label_[i]];
    // Use 2k+1 window.
    t1_label_id = t1_label_[i];
    start_pos = t2_label_il_start_pos_[t1_label_id];
    end_pos = std::min(start_pos + 2 * k, (int)candidate_ids.size()-1);
    pos = start_pos;
    while (pos <= end_pos) {
      cand_id = candidate_ids[pos];
      if (cand_id - i > k) {
        // We're certainly outside the window from the right side.
        break;
      }
      if (cand_id == -1 || i - cand_id > k) {
        // The label has been mapped or we're outside the window from the left side.
        ++t2_label_il_start_pos_[t1_label_id]; // For the consecutive i we don't need to start before.
        ++pos;
        continue;
      }
      if (k_relevant(i, cand_id, k)) {
        // The pair can be mapped.
        mapping.push_back({i, cand_id}); // postorder
        candidate_ids[pos] = -1;
        break;
      }
      ++pos;
    }
  }
  mapping = to_ted_mapping(mapping);
  return mapping;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::lb_mapping_deprecated(
    const node::Node<Label>& t1, const node::Node<Label>& t2, const int k) {
  init(t1, t2);
  std::vector<std::pair<int, int>> mapping;
  for (int i = 0; i < t1_input_size_; ++i) { // Loop in postorder.
    std::list<int>& candidate_ids = t2_label_il_hash_[t1_node_[i].get().label().to_string()];
    std::list<int>::iterator cand_it = candidate_ids.begin();
    int cand_id = 0;
    for ( ; cand_it != candidate_ids.end(); ) {
      cand_id = *cand_it;
      if (k_relevant(i, cand_id, k)) {
        mapping.push_back({i, cand_id}); // postorder
        cand_it = candidate_ids.erase(cand_it); // Assignment not needed due to break.
        break;
      }
      ++cand_it;
    }
  }
  mapping = to_ted_mapping(mapping);
  return mapping;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::lb_mapping_fill_gaps(
    const node::Node<Label>& t1, const node::Node<Label>& t2, const int k) {
  std::vector<std::pair<int, int>> mapping = lb_mapping(t1, t2, k);
  mapping = fill_gaps_in_mapping(mapping, k);
  return mapping;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::lb_mapping_fill_gaps_deprecated(
    const node::Node<Label>& t1, const node::Node<Label>& t2, const int k) {
  std::vector<std::pair<int, int>> mapping = lb_mapping_deprecated(t1, t2, k);
  mapping = fill_gaps_in_mapping_deprecated(mapping, k);
  return mapping;
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::update_desc_when_not_mapped(
    const int node, std::vector<int>& count_mapped_desc,
    const std::vector<int>& parent, const int input_size) const {
  if (node < input_size - 1) { // Root has no parent nor the right leaf, and
                               // the nodes of dummy mapping in fill_gaps_in_mapping
                               // do not exist.
    count_mapped_desc[parent[node]] += count_mapped_desc[node];
  }
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::update_desc_when_mapped(
    const int node, std::vector<int>& count_mapped_desc,
    const std::vector<int>& parent, const int input_size) const {
  if (node < input_size - 1) { // Root has no parent nor the right leaf, and
                               // the nodes of dummy mapping in fill_gaps_in_mapping
                               // do not exist.
    count_mapped_desc[parent[node]] += count_mapped_desc[node] + 1;
  }
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::update_desc_and_left_when_not_mapped(
    const int node, std::vector<int>& count_mapped_desc,
    std::vector<int>& count_mapped_left, const std::vector<int>& parent,
    const std::vector<int>& rl, const std::vector<int>& size,
    const std::vector<int>& post_to_pre, const int input_size) const {
  if (node < input_size - 1) { // Root has no parent nor the right leaf, and
                               // the nodes of dummy mapping in fill_gaps_in_mapping
                               // do not exist.
    // Descendants.
    count_mapped_desc[parent[node]] += count_mapped_desc[node];
    // Nodes to the left.
    if (rl[node] > -1) { // if the node has a right leaf
      if (size[node] == 1) { // leaf node
        count_mapped_left[rl[node]] = count_mapped_left[node];
      }
    }
    if (post_to_pre[parent[node]] + 1 == post_to_pre[node]) { // leftmost child of its parent
                                                              // NOTE: Could be precomputed.
      count_mapped_left[parent[node]] = count_mapped_left[node];
    }
  }
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::update_desc_and_left_when_mapped(
    const int node, std::vector<int>& count_mapped_desc,
    std::vector<int>& count_mapped_left, const std::vector<int>& parent,
    const std::vector<int>& rl, const std::vector<int>& size,
    const std::vector<int>& post_to_pre, const int input_size) const {
  if (node < input_size - 1) { // Root has no parent nor the right leaf, and
                               // the nodes of dummy mapping in fill_gaps_in_mapping
                               // do not exist.
    // Descendants.
    count_mapped_desc[parent[node]] += count_mapped_desc[node] + 1;
    // Nodes to the left.
    if (rl[node] > -1) { // if the node has a right leaf
      if (size[node] == 1) { // leaf node
        count_mapped_left[rl[node]] = count_mapped_left[node] + 1;
      } else { // inner node
        count_mapped_left[rl[node]] += 1;
      }
    }
    if (post_to_pre[parent[node]] + 1 == post_to_pre[node]) { // leftmost child of its parent
                                                              // NOTE: Could be precomputed.
      count_mapped_left[parent[node]] = count_mapped_left[node];
    }
  }
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::get_mapped_ancestors_counts(
    std::vector<std::pair<int, int>>& mapping,
    std::vector<int>& t1_count_mapped_anc,
    std::vector<int>& t2_count_mapped_anc) const {
  std::vector<int> t1_mapped_nodes(t1_input_size_); // indexed in postorder
  std::vector<int> t2_mapped_nodes(t2_input_size_); // indexed in postorder
  for (auto m : mapping) {
    t1_mapped_nodes[m.first] = 1;
    t2_mapped_nodes[m.second] = 1;
  }
  int i_in_post = 0;
  // T1
  for (int i = 1; i < t1_input_size_; ++i) { // loop in preorder; root has no parent and 0 mapped ancestors
    i_in_post = t1_pre_to_post_[i];
    t1_count_mapped_anc[i_in_post] = t1_mapped_nodes[t1_parent_[i_in_post]];
    t1_mapped_nodes[i_in_post] += t1_mapped_nodes[t1_parent_[i_in_post]];
  }
  // std::cout << "t1_count_mapped_anc : ";
  // for (auto e : t1_count_mapped_anc) {
  //   std::cout << e << " ";
  // } std::cout << std::endl;
  // T2
  for (int i = 1; i < t2_input_size_; ++i) { // loop in preorder; root has no parent and 0 mapped ancestors
    i_in_post = t2_pre_to_post_[i];
    t2_count_mapped_anc[i_in_post] = t2_mapped_nodes[t2_parent_[i_in_post]];
    t2_mapped_nodes[i_in_post] += t2_mapped_nodes[t2_parent_[i_in_post]];
  }
  // std::cout << "t2_count_mapped_anc : ";
  // for (auto e : t2_count_mapped_anc) {
  //   std::cout << e << " ";
  // } std::cout << std::endl;
};

template <typename Label, typename CostModel>
bool GreedyUB<Label, CostModel>::if_in_corresponding_regions(int t1_begin_gap,
    int i, int t1_end_gap, int t2_begin_gap, int j, int t2_end_gap) const {
  
  // Begin or end pairs of a gap may not exist.
  
  if (t1_begin_gap == -1 && t1_end_gap == t1_input_size_) { // also t2_begin_gap == -1 && t2_end_gap == t2_input_size_
    // Any mapping is fine - no defined regions.
    return true;
  }
  
  int i_pre = t1_post_to_pre_[i];
  int j_pre = t2_post_to_pre_[j];
  
  if (t1_begin_gap == -1) { // also t2_begin_gap == 1
    if (i_pre > t1_post_to_pre_[t1_end_gap] && j_pre > t2_post_to_pre_[t2_end_gap]) {
      return true;
    }
    if (i_pre < t1_post_to_pre_[t1_end_gap] && j_pre < t2_post_to_pre_[t2_end_gap]) {
      return true;
    }
  }
  
  if (t1_end_gap == t1_input_size_) { // also t2_end_gap == t2_input_size_
    if (i_pre > t1_post_to_pre_[t1_begin_gap] && j_pre > t2_post_to_pre_[t2_begin_gap]) {
      return true;
    }
    if (i_pre < t1_post_to_pre_[t1_begin_gap] && j_pre < t2_post_to_pre_[t2_begin_gap]) {
      return true;
    }
  }
  
  int t1_begin_gap_pre = t1_post_to_pre_[t1_begin_gap];
  int t2_begin_gap_pre = t2_post_to_pre_[t2_begin_gap];
  int t1_end_gap_pre = t1_post_to_pre_[t1_end_gap];
  int t2_end_gap_pre = t2_post_to_pre_[t2_end_gap];
  
  // Three different cases!
  if (t1_end_gap_pre < i_pre && i_pre < t1_begin_gap_pre) {
    if (t2_end_gap_pre < j_pre && j_pre < t2_begin_gap_pre) {
      return true;
    }
  }
  if (t1_end_gap_pre < i_pre && t1_begin_gap_pre < i_pre) {
    if (t2_end_gap_pre < j_pre && t2_begin_gap_pre < j_pre) {
      return true;
    }
  }
  if (t1_begin_gap_pre < i_pre && i_pre < t1_end_gap_pre) {
    if (t2_begin_gap_pre < j_pre && j_pre < t2_end_gap_pre) {
      return true;
    }
  }
  if (i_pre < t1_begin_gap_pre && i_pre < t1_end_gap_pre) {
    if (j_pre < t2_begin_gap_pre && j_pre < t2_end_gap_pre) {
      return true;
    }
  }
  
  return false;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::fill_gaps_in_mapping(
    std::vector<std::pair<int, int>>& mapping, const int k) const {
  
  // In result_mapping we store the output of this function.
  std::vector<std::pair<int, int>> result_mapping;
  
  // The counts for mapped descendants and nodes to the left must be maintained.
  std::vector<int> t1_count_mapped_desc(t1_input_size_);
  std::vector<int> t2_count_mapped_desc(t2_input_size_);
    
  // The counts for mapped ancestors.
  std::vector<int> t1_count_mapped_anc(t1_input_size_);
  std::vector<int> t2_count_mapped_anc(t2_input_size_);
  get_mapped_ancestors_counts(mapping, t1_count_mapped_anc, t2_count_mapped_anc);
  
  // Add a dummy mapping pair to cover the gap from the last mapping until the
  // ends of the trees.
  mapping.push_back({t1_input_size_, t2_input_size_});
  
  std::pair<int, int> begin_gap = {-1, -1};
  std::pair<int, int> end_gap;
  std::vector<std::pair<int, int>>::iterator end_gap_it = mapping.begin();
  
  while (end_gap_it != mapping.end()) {
    end_gap = *end_gap_it;
    
    // If there is no gap, continue with the next mapped pair.
    if (end_gap.first - begin_gap.first == 1 || end_gap.second - begin_gap.second == 1) {
      // Update the counts for nodes between begin_gap and end_gap, if any,
      // as not-mapped nodes.
      if (end_gap.first - begin_gap.first > 1) { // There is a gap between first elements.
        for (int not_mapped_first = begin_gap.first + 1; not_mapped_first < end_gap.first; ++not_mapped_first) {
          update_desc_when_not_mapped(not_mapped_first, t1_count_mapped_desc, t1_parent_, t1_input_size_);
        }
      }
      // Update the counts for nodes in end_gap as mapped nodes.
      update_desc_when_mapped(end_gap.first, t1_count_mapped_desc, t1_parent_, t1_input_size_);
      
      // Update the counts for nodes between begin_gap and end_gap, if any,
      // as not-mapped nodes.
      if (end_gap.second - begin_gap.second > 1) { // There is a gap only between second elements.
        for (int not_mapped_second = begin_gap.second + 1; not_mapped_second < end_gap.second; ++not_mapped_second) {
          update_desc_when_not_mapped(not_mapped_second, t2_count_mapped_desc, t2_parent_, t2_input_size_);
        }
      }
      // Update the counts for nodes in end_gap as mapped nodes.
      update_desc_when_mapped(end_gap.second, t2_count_mapped_desc, t2_parent_, t2_input_size_);
      
      // Push the end of the current gap to the result.
      result_mapping.push_back(end_gap);    
      // Update the begin of the next gap.
      begin_gap = end_gap;
      ++end_gap_it;
      continue;
    }
    
    // There is a gap --> try to map nodes --> at the first node pair mapped,
    // change the gap and continue the loop.

    int i = begin_gap.first + 1;
    int i_last = end_gap.first - 1;
    int j = begin_gap.second + 1;
    int j_last = end_gap.second - 1;
    // std::pair<int, int> last_mapped = begin_gap;
    bool mapped_in_gap = false;
    while (i <= i_last && j <= j_last) {
      if (k_relevant(i, j, k) &&
          (t1_count_mapped_desc[i] == t2_count_mapped_desc[j] &&
          t1_count_mapped_anc[i] == t2_count_mapped_anc[j]) &&
          if_in_corresponding_regions(begin_gap.first, i, end_gap.first, begin_gap.second, j, end_gap.second)) {
        // Map (i,j) + push to result mapping + update counts for both nodes
        // incremented as mapped.
        update_desc_when_mapped(i, t1_count_mapped_desc, t1_parent_, t1_input_size_);
        update_desc_when_mapped(j, t2_count_mapped_desc, t2_parent_, t2_input_size_);
        result_mapping.push_back({i, j});
        // last_mapped = {i, j};
        begin_gap = {i, j};
        mapped_in_gap = true;
        break;
      } else {
        // Update counts for the incremented node as non-mapped.
        if (i_last - i > j_last - j) {
          update_desc_when_not_mapped(i, t1_count_mapped_desc, t1_parent_, t1_input_size_);
          ++i;
        } else {
          update_desc_when_not_mapped(j, t2_count_mapped_desc, t2_parent_, t2_input_size_);
          ++j;
        }
      }
    }
    if (mapped_in_gap) {
      continue;
    } else {
      // Check the case when only one of i or j incremented beyond the limit.
      // The gap in one of the trees has been used up.
      // Nothing mapped in gap.
      if (i_last - i < 0 && j_last - j >= 0) {
        for (int j_missing = j; j_missing <= j_last; ++j_missing) {
          update_desc_when_not_mapped(j, t2_count_mapped_desc, t2_parent_, t2_input_size_);
        }
      }
      if (j_last - j < 0 && i_last - i >= 0) {
        for (int i_missing = i; i_missing <= i_last; ++i_missing) {
          update_desc_when_not_mapped(i, t1_count_mapped_desc, t1_parent_, t1_input_size_);
        }
      }
      // Update the counts for nodes in end_gap as mapped nodes.
      update_desc_when_mapped(end_gap.first, t1_count_mapped_desc, t1_parent_, t1_input_size_);
      update_desc_when_mapped(end_gap.second, t2_count_mapped_desc, t2_parent_, t2_input_size_);
      // Push the end of the current gap to the result.
      result_mapping.push_back(end_gap);
      // Update the begin of the next gap.
      begin_gap = end_gap;
      ++end_gap_it;
      continue;
    }
  }
    
  // Remove the last dummy mapping.
  result_mapping.pop_back();
  return result_mapping;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::fill_gaps_in_mapping_deprecated(std::vector<std::pair<int, int>>& mapping,
    const int k) const {
  std::vector<std::pair<int, int>> result_mapping;
  
  // The counts for mapped descendants and nodes to the left must be maintained.
  std::vector<int> t1_count_mapped_desc(t1_input_size_);
  std::vector<int> t2_count_mapped_desc(t2_input_size_);
  std::vector<int> t1_count_mapped_left(t1_input_size_);
  std::vector<int> t2_count_mapped_left(t2_input_size_);
    
  // The counts for mapped ancestors.
  std::vector<int> t1_count_mapped_anc(t1_input_size_);
  std::vector<int> t2_count_mapped_anc(t2_input_size_);
  get_mapped_ancestors_counts(mapping, t1_count_mapped_anc, t2_count_mapped_anc);
  
  std::pair<int, int> begin_gap = {-1, -1};
  
  // Add a dummy mapping pair to cover the gap from the last mapping until the
  // ends of the trees.
  mapping.push_back({t1_input_size_, t2_input_size_});
  
  for (auto end_gap : mapping) {
    // If there is no gap, continue with the next mapped pair.
    if (end_gap.first - begin_gap.first == 1 || end_gap.second - begin_gap.second == 1) {
      // Update the counts for nodes between begin_gap and end_gap, if any,
      // as not-mapped nodes.
      if (end_gap.first - begin_gap.first > 1) { // There is a gap only between first elements.
        for (int not_mapped_first = begin_gap.first + 1; not_mapped_first < end_gap.first; ++not_mapped_first) {
          update_desc_and_left_when_not_mapped(not_mapped_first,
              t1_count_mapped_desc, t1_count_mapped_left, t1_parent_, t1_rl_,
              t1_size_, t1_post_to_pre_, t1_input_size_);
        }
      }
      update_desc_and_left_when_mapped(end_gap.first, t1_count_mapped_desc,
          t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
          t1_input_size_);
      
      // Update the counts for nodes between begin_gap and end_gap, if any,
      // as not-mapped nodes.
      if (end_gap.second - begin_gap.second > 1) { // There is a gap only between second elements.
        for (int not_mapped_second = begin_gap.second + 1; not_mapped_second < end_gap.second; ++not_mapped_second) {
          update_desc_and_left_when_not_mapped(not_mapped_second,
              t2_count_mapped_desc, t2_count_mapped_left, t2_parent_, t2_rl_,
              t2_size_, t2_post_to_pre_, t2_input_size_);
        }
      }
      // Update the counts for nodes in end_gap as mapped nodes.
      update_desc_and_left_when_mapped(end_gap.second, t2_count_mapped_desc,
          t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
          t2_input_size_);
      
      // Push the end of the current gap to the result.
      result_mapping.push_back(end_gap);    
      // Update the begin of the next gap.
      begin_gap = end_gap;
      continue;
    }
    
    // There is a gap --> try to map nodes.
    // While filling the gaps the counts have to be updated.

    int i = begin_gap.first + 1;
    int i_last = end_gap.first - 1;
    int j = begin_gap.second + 1;
    int j_last = end_gap.second - 1;
    std::pair<int, int> last_mapped = begin_gap;
    while (i <= i_last && j <= j_last) {
      if (k_relevant(i, j, k) &&
          (t1_count_mapped_desc[i] == t2_count_mapped_desc[j] &&
          t1_count_mapped_left[i] == t2_count_mapped_left[j] &&
          t1_count_mapped_anc[i] == t2_count_mapped_anc[j]) &&
          if_in_corresponding_regions(last_mapped.first, i, end_gap.first, last_mapped.second, j, end_gap.second)) {
        // Map (i,j) + push to result mapping + update counts for both nodes
        // incremented as mapped.
        update_desc_and_left_when_mapped(i, t1_count_mapped_desc,
            t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
            t1_input_size_);
        update_desc_and_left_when_mapped(j, t2_count_mapped_desc,
            t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
            t2_input_size_);
        result_mapping.push_back({i, j});
        last_mapped = {i, j};
        ++i;
        ++j;
      } else {
        // Update counts for the incremented node as non-mapped.
        if (i_last - i > j_last - j) {
          update_desc_and_left_when_not_mapped(i, t1_count_mapped_desc,
              t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_,
              t1_post_to_pre_, t1_input_size_);
          ++i;
        } else {
          update_desc_and_left_when_not_mapped(j, t2_count_mapped_desc,
              t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_,
              t2_post_to_pre_, t2_input_size_);
          ++j;
        }
      }
    }
    // Check the case when only one of i or j incremented beyond the limit.
    // The gap in one of the trees has been used up.
    if (i_last - i < 0 && j_last - j >= 0) {
      for (int j_missing = j; j_missing <= j_last; ++j_missing) {
        update_desc_and_left_when_not_mapped(j, t2_count_mapped_desc,
            t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
            t2_input_size_);
      }
    }
    if (j_last - j < 0 && i_last - i >= 0) {
      for (int i_missing = i; i_missing <= i_last; ++i_missing) {
        update_desc_and_left_when_not_mapped(i, t1_count_mapped_desc,
            t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
            t1_input_size_);
      }
    }
      
    // Update the counts for nodes in end_gap as mapped nodes.
    update_desc_and_left_when_mapped(end_gap.first, t1_count_mapped_desc,
        t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
        t1_input_size_);
    update_desc_and_left_when_mapped(end_gap.second, t2_count_mapped_desc,
        t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
        t2_input_size_);
    // Push the end of the current gap to the result.
    result_mapping.push_back(end_gap);    
    // Update the begin of the next gap.
    begin_gap = end_gap;
  }
    
  // Remove the last dummy mapping.
  result_mapping.pop_back();
  return result_mapping;
};

template <typename Label, typename CostModel>
std::vector<std::pair<int, int>> GreedyUB<Label, CostModel>::to_ted_mapping(
    const std::vector<std::pair<int, int>>& mapping) const {
  std::vector<std::pair<int, int>> ted_mapping;
  
  std::vector<int> t1_count_mapped_desc(t1_input_size_);
  std::vector<int> t2_count_mapped_desc(t2_input_size_);
  std::vector<int> t1_count_mapped_left(t1_input_size_);
  std::vector<int> t2_count_mapped_left(t2_input_size_);
  
  int t1_i = 0;
  int t2_i = 0;
  int cur_t1 = 0; // Node of T1 in the current mapped pair considered.
  int cur_t2 = 0; // Node of T2 in the current mapped pair considered.
  int prev_t2 = -1; // Previously mapped node in T2.
  bool mapped_t1_node_processed  = false;
  bool mapped_t2_node_processed  = false;
  for (auto m : mapping) {
    cur_t1 = m.first;
    cur_t2 = m.second;
    // Increasing postorder id test.
    // The pairs in the input one-to-one mapping are sorted on the first
    // element. The second elements of two consecutive pairs must not have
    // decreasing postorder ids.
    // NOTE: Non-decreasing ids of second elements could be ensured in lb_mapping.
    //       But, this will reduce the number of initially mapped node which may not be
    //       desired. For example, if given M=<(a,b),(a+1,b+2),(a+2,b+1)>, (a+1,b+2) is
    //       removed due to violating TED conditions, (a+2,b+1) would may still be valid, and
    //       have a lower cost than (a+1,b+1) possible from filling gaps.
    if (cur_t2 < prev_t2) {
      continue;
    }
    // Go back with t1_i and t2_i if the current node ids are smaller.
    // If cur_t1 < t1_i => cur_t1 has been processed as a non-mapped node and
    // it has to be processed as a mapped node, rolling back the changes made
    // to it.
    // If cur_t1 == t1_i => cur_t1 has to be processed as a mapped node,
    // it has not been processed before.
    if (cur_t1 < t1_i) {
      t1_i = cur_t1;
      mapped_t1_node_processed = true;
    }
    if (cur_t2 < t2_i) {
      t2_i = cur_t2;
      mapped_t2_node_processed = true;
    }
    if (cur_t1 == t1_i) {
      mapped_t1_node_processed = false;
    }
    if (cur_t2 == t2_i) {
      mapped_t2_node_processed = false;
    }
        
    // Nodes not in mapping in T1.
    while (t1_i < cur_t1) {
      update_desc_and_left_when_not_mapped(t1_i, t1_count_mapped_desc,
          t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
          t1_input_size_);
      ++t1_i;
    }
    // Nodes not in mapping in T2.
    while (t2_i < cur_t2) {
      update_desc_and_left_when_not_mapped(t2_i, t2_count_mapped_desc,
          t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
          t2_input_size_);
      ++t2_i;
    }
    
    // If the current node pair cannot be mapped, the nodes must not be
    // processed as mapped (the next two if statements).
    
    // Mapped descendants test.
    if (t1_count_mapped_desc[cur_t1] != t2_count_mapped_desc[cur_t2]) {
      continue;
    }
    // TODO: Remove everything related to mapped left.
    // // Mapped nodes to the left test.
    // if (t1_count_mapped_left[cur_t1] != t2_count_mapped_left[cur_t2]) {
    //   continue;
    // }
    
    // Mapped node in T1.
    if (!mapped_t1_node_processed) {
      update_desc_and_left_when_mapped(t1_i, t1_count_mapped_desc,
          t1_count_mapped_left, t1_parent_, t1_rl_, t1_size_, t1_post_to_pre_,
          t1_input_size_);
      ++t1_i;
    } else { // The current mapped node has been processed and it's enough to
             // fix the values it updated.
      if (t1_parent_[t1_i] > -1) { // root has no parent
        // descendants
        t1_count_mapped_desc[t1_parent_[t1_i]] += 1; // THE ONLY DIFF TO ABOVE.
      }
      ++t1_i;
    }
    // Mapped node in T2.
    if (!mapped_t2_node_processed) {
      update_desc_and_left_when_mapped(t2_i, t2_count_mapped_desc,
          t2_count_mapped_left, t2_parent_, t2_rl_, t2_size_, t2_post_to_pre_,
          t2_input_size_);
      ++t2_i;
    } else { // The current mapped node has been processed and it's enough to
             // fix the values it updated.
      if (t2_parent_[t2_i] > -1) { // root has no parent
        // descendants
        t2_count_mapped_desc[t2_parent_[t2_i]] += 1; // THE ONLY DIFF TO ABOVE.
      }
      ++t2_i;
    }
    
    ted_mapping.push_back({cur_t1, cur_t2});
    
    prev_t2 = cur_t2;
  }
  
  // NOTE: The gap after the last mapped node until the last node in the tree
  //       does not matter. Here, we only revise the given mapping. If there
  //       are no pairs to revise left, we can stop. 
  
  return ted_mapping;
};

template <typename Label, typename CostModel>
int GreedyUB<Label, CostModel>::index_nodes_recursion(
    const node::Node<Label>& node,
    std::vector<std::vector<int>>& label_il,
    std::unordered_map<std::string, std::list<int>>& label_il_hash,
    std::vector<std::reference_wrapper<const node::Node<Label>>>& nodes,
    std::vector<int>& post_to_pre,
    std::vector<int>& pre_to_post,
    std::vector<int>& parent,
    std::vector<int>& depth,
    std::vector<int>& size,
    std::vector<int>& label,
    int& start_postorder,
    int& start_preorder,
    unsigned int start_depth) {
  // Here, start_preorder holds this node's preorder id.
  int current_preorder = start_preorder;
  
  // Stores number of descendants of this node. Incrementally computed while
  // traversing the children.
  int desc_sum = 0;

  // Increment start_preorder to hold the correct id of the consecutive node
  // in preorder.
  ++start_preorder;
  
  // To store postorder ids of this node's children.
  std::vector<int> nodes_children;
  
  auto children_start_it = std::begin(node.get_children());
  auto children_end_it = std::end(node.get_children());
  while (children_start_it != children_end_it) {
    desc_sum += index_nodes_recursion(*children_start_it, label_il, label_il_hash, nodes, post_to_pre,
                          pre_to_post,
                          parent, depth, size, label, start_postorder, start_preorder,
                          start_depth + 1);
    // Here, start_postorder-1 is the postorder of the current child.
    // Collect children ids.
    nodes_children.push_back(start_postorder-1);
    // Continue to consecutive children.
    ++children_start_it;
  }
  
  // Here, start_postorder holds this node's postorder id here.

  // Parent - dummy element to grow the vector by this node.
  parent.push_back(-1);
  // Set this node's postorder value as parent for all its children.
  for (auto child_id : nodes_children) {
    parent.at(child_id) = start_postorder;
  }
  
  if (node.is_leaf()) {
    // Leaf has size 1.
    size.push_back(1);
  } else {
    // Inner node has size desc_sum+1.
    size.push_back(desc_sum + 1);
  }
  
  // Depth.
  depth.push_back(start_depth);

  // Add current node to the nodes vector.
  nodes.push_back(std::ref(node));
  
  // Add current node postorder id to label inverted list.
  label_il_hash[node.label().to_string()].push_back(start_postorder);
  // New inverted list.
  unsigned int label_id_in_dict = dict_.insert(node.label());
  label_il[label_id_in_dict].push_back(start_postorder);
  label.push_back(label_id_in_dict);
  
  // Store postorder to preorder translation.
  post_to_pre.push_back(current_preorder);
  
  // Store preorder to postorder translation.
  pre_to_post[current_preorder] = start_postorder;
  
  // Increment start_postorder for the consecutive node in postorder have the
  // correct id.
  ++start_postorder;
  
  // Return the number of nodes in the subtree rooted at this node.
  return desc_sum + 1;
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::index_nodes(
    const node::Node<Label>& root,
    std::vector<std::vector<int>>& label_il,
    std::unordered_map<std::string, std::list<int>>& label_il_hash,
    std::vector<std::reference_wrapper<const node::Node<Label>>>& nodes,
    std::vector<int>& post_to_pre,
    std::vector<int>& pre_to_post,
    std::vector<int>& parent,
    std::vector<int>& depth,
    std::vector<int>& size,
    std::vector<int>& label) {
  // Orders start with '0'.
  int start_postorder = 0;
  // NOTE: Preorder is not used. Remove start_preorder. Or
  //       move the template traversal with postorder and preorder to some notes
  //       of how to traverse trees.
  int start_preorder = 0;
  index_nodes_recursion(root, label_il, label_il_hash, nodes, post_to_pre, pre_to_post, parent, depth, size, label, start_postorder, start_preorder, 0);
  // Here, start_postorder and start_preorder store the size of tree minus 1.
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::post_traversal_indexing(const int input_size,
    const std::vector<std::reference_wrapper<const node::Node<Label>>>& nodes,
    std::vector<int>& rl) {
  int current_leaf = -1;
  int post_id = -1;
  for(int i = 0; i < input_size; i++) { // loop in right-to-left preorder
    post_id = input_size - i - 1; // Translation from right-to-left preorder
                                  // to postorder.
    rl[post_id] = current_leaf;
    if(nodes[post_id].get().is_leaf()) { // NOTE: Could be precomputed.
      current_leaf = post_id;
    }
  }
};

template <typename Label, typename CostModel>
void GreedyUB<Label, CostModel>::init(const node::Node<Label>& t1,
                                    const node::Node<Label>& t2) {
  t1_node_.clear();
  t2_node_.clear();
  t1_post_to_pre_.clear();
  t2_post_to_pre_.clear();
  t1_pre_to_post_.clear();
  t2_pre_to_post_.clear();
  t1_label_il_.clear();
  t2_label_il_.clear();
  t1_label_il_hash_.clear();
  t2_label_il_hash_.clear();
  t1_parent_.clear();
  t2_parent_.clear();
  t1_rl_.clear();
  t2_rl_.clear();
  t1_depth_.clear();
  t2_depth_.clear();
  t1_size_.clear();
  t2_size_.clear();
  dict_.clear();
  t1_label_.clear();
  t2_label_.clear();
  
  // TODO: Do not call get_tree_size() that causes an additional tree traversal.
  //       Index subtree sizes instead - they'll be used anyways.
  // TODO: Getting tree size first makes things easier while indexing nodes.
  //       If the size is known upfront it could be immediately assigned.
  t1_input_size_ = t1.get_tree_size();
  t2_input_size_ = t2.get_tree_size();
  
  t1_pre_to_post_.resize(t1_input_size_);
  t2_pre_to_post_.resize(t2_input_size_);
  
  // The size of the inverted list is at most sum of the input trees sizes.
  t1_label_il_.resize(t1_input_size_ + t2_input_size_);
  t2_label_il_.resize(t1_input_size_ + t2_input_size_);
  
  t2_label_il_start_pos_.resize(t1_input_size_ + t2_input_size_);
  
  index_nodes(t1, t1_label_il_, t1_label_il_hash_, t1_node_, t1_post_to_pre_, t1_pre_to_post_, t1_parent_, t1_depth_, t1_size_, t1_label_);
  index_nodes(t2, t2_label_il_, t2_label_il_hash_, t2_node_, t2_post_to_pre_, t2_pre_to_post_, t2_parent_, t2_depth_, t2_size_, t2_label_);
  
  t1_rl_.resize(t1_input_size_);
  t2_rl_.resize(t2_input_size_);
  post_traversal_indexing(t1_input_size_, t1_node_, t1_rl_);
  post_traversal_indexing(t2_input_size_, t2_node_, t2_rl_);
};

template <typename Label, typename CostModel>
bool GreedyUB<Label, CostModel>::k_relevant(const int x, const int y, const int k) const {
  // The lower bound formula (RA + D + L):
  // |(|T1|-(x+1))-(|T2|-(y+1))| + ||T1_x|-|T2_y|| + |((x+1)-|T1_x|)-((y+1)-|T2_y|)| < k
  // New lower bound formula (R + A + D + L):
  // |(|T1|-(x+1)-depth(x))-(|T2|-(y+1)-depth(y))| + |depth(x)-depth(y)| + ||T1_x|-|T2_y|| + |((x+1)-|T1_x|)-((y+1)-|T2_y|)| < k
  int x_size = t1_size_[x];
  int y_size = t2_size_[y];
  // int lower_bound = std::abs((t1_size_.back() - (x+1)) - (t2_size_.back() - (y+1))) +
  //                   std::abs(x_size - y_size) +
  //                   std::abs(((x+1) - x_size) - ((y+1) - y_size));

  int lower_bound = std::abs((t1_size_.back() - (x+1) - t1_depth_[x]) - (t2_size_.back() - (y+1) - t2_depth_[y])) +
                    std::abs(t1_depth_[x] - t2_depth_[y]) +
                    std::abs(x_size - y_size) +
                    std::abs(((x+1) - x_size) - ((y+1) - y_size));

  // NOTE: The pair (x,y) is k-relevant if lower_bound <= k.
  //       lower_bound < k is not correct because then (x,y) would be
  //       k-irrelevant for lower_bound = k. That would further mean that the
  //       subtrees T1_x and T2_y cannot be mapped with the given budget.
  if (lower_bound <= k) {
    return true;
  }
  return false;
};

template <typename Label, typename CostModel>
const typename GreedyUB<Label, CostModel>::TestItems GreedyUB<Label, CostModel>::get_test_items() const {
  TestItems test_items = {
    t1_label_il_,
    t1_post_to_pre_,
    t1_pre_to_post_,
    t1_parent_,
    t1_rl_,
    t1_depth_,
    t1_size_,
    dict_,
  };
  return test_items;
};

#endif // TREE_SIMILARITY_TED_UB_GREEDY_UB_IMPL_H
