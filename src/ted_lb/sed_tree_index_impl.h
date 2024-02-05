// The MIT License (MIT)
// Copyright (c) 2018 Nikolaus Augsten, Mateusz Pawlik.
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

/// Contains the implementation of the String Edit Distance algorithm.

#pragma once

template <typename CostModel, typename TreeIndex>
double SEDTreeIndex<CostModel, TreeIndex>::ted(const TreeIndex& t1, const TreeIndex& t2) {
      
  // Reset subproblem counter.
  subproblem_counter_ = 0;
  
  // TODO: Swap input trees for memory optimisation, as below.
  // if (x.length() > y.length()) {
  //   String tmp = x;
  //   x = y;
  //   y = tmp;
  // }
  const int kT1Size = t1.tree_size_;
  const int kT2Size = t2.tree_size_;
  
  std::vector<double> cp(kT2Size + 1);
  std::vector<double> cc(kT2Size + 1);
  
  // Edit distance on postorder traversal.
  for (int j = 0; j <= kT2Size; ++j) {
    cp[j] = j;
  }
  for (int i = 1; i <= kT1Size; ++i) {
    cc[0] = i;
    for (int j = 1; j <= kT2Size; ++j) {
      ++subproblem_counter_;
      unsigned int c1 = t1.postl_to_label_id_[i - 1];
      unsigned int c2 = t2.postl_to_label_id_[j - 1];
      cc[j] = std::min({
          cp[j - 1] + ((c1 == c2) ? 0 : 1),
          cc[j - 1] + 1,
          cp[j] + 1
          });
    }
    std::swap(cp, cc);
  }
  double sed_post = cp[kT2Size];
  
  // Edit distance on preorder traversal.
  for (int j = 0; j <= kT2Size; ++j) {
    cp[j] = j;
  }
  for (int i = 1; i <= kT1Size; ++i) {
    cc[0] = i;
    for (int j = 1; j <= kT2Size; ++j) {
      ++subproblem_counter_;
      unsigned int c1 = t1.prel_to_label_id_[i - 1];
      unsigned int c2 = t2.prel_to_label_id_[j - 1];
      cc[j] = std::min({
          cp[j - 1] + ((c1 == c2) ? 0 : 1),
          cc[j - 1] + 1,
          cp[j] + 1
          });
    }
    std::swap(cp, cc);
  }
  double sed_pre = cp[kT2Size];
  
  return std::max(sed_post, sed_pre); // note: cc is asigned to cp, cc is overwritten!
}

template <typename T>
double bounded_sed(const T* a, const T* b, int64_t alen, int64_t blen, int64_t k) {
    // I assume a is always at max the size as b (t2 is always bigger)
    // suffix trimming
    while (alen > 0 && a[alen - 1] == b[blen - 1]) {
        alen--;
        blen--;
    }

    // prefix trimming
    while (alen > 0 && *a == *b) {
        a++;
        b++;
        alen--;
        blen--;
    }
    k = std::min(blen, k);
    // if the shorter string is gone, return b_size or threshold
    if (alen == 0) {
        return k;
    }
    auto size_d = blen - alen;

    if (size_d > k) {
        return k;
    }

    int64_t ZERO_K = std::min(k, alen) / 2 + 2;
    auto array_size = size_d + ZERO_K * 2 + 2;

    std::vector<int64_t> current_row(array_size, -1);
    std::vector<int64_t> next_row(array_size, -1);

    int64_t i = 0;
    int64_t condition_row = size_d + ZERO_K;
    int64_t end_max = condition_row * 2;

    do {
        i++;

        std::swap(current_row, next_row);

        int64_t start;
        int64_t previous_cell;
        int64_t current_cell = -1;
        int64_t next_cell;

        if (i <= ZERO_K) {
            start = -i + 1;
            next_cell = i - 2;
        } else {
            start = i - ZERO_K * 2 + 1;
            next_cell = current_row[ZERO_K + start];
        }

        int64_t end;
        if (i <= condition_row) {
            end = i;
            next_row[ZERO_K + i] = -1;
        } else {
            end = end_max - i;
        }

        for (int64_t q = start, row_index = start + ZERO_K; q < end; q++, row_index++) {
            previous_cell = current_cell;
            current_cell = next_cell;
            next_cell = current_row[row_index + 1];

            int64_t t = std::max({
                                         current_cell + 1,
                                         previous_cell,
                                         next_cell + 1
                                 });

            while (t < alen && t + q < blen && a[t] == b[t + q]) {
                t++;
            }

            next_row[row_index] = t;
        }
    } while (next_row[condition_row] < alen && i <= k);

    return i - 1;
}

template <typename CostModel, typename TreeIndex>
double SEDTreeIndex<CostModel, TreeIndex>::ted_k(const TreeIndex& t1, const TreeIndex& t2, int64_t k) {
    // I assume t1 is always at max the size as t2 (t2 is always bigger)
    auto alen = t1.prel_to_label_id_.size();
    auto blen = t2.prel_to_label_id_.size();

    auto sed_post = bounded_sed(t1.postl_to_label_id_.data(), t2.postl_to_label_id_.data(), alen, blen, k);
    auto sed_pre = bounded_sed(t1.prel_to_label_id_.data(), t2.prel_to_label_id_.data(), alen, blen, k);
    return std::max(sed_post, sed_pre);
}
