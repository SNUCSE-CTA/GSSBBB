#pragma once
#include "Base.h"

void VectorIntersection(std::vector<int> &A, std::vector<int> &B,
                        std::vector<int> &results) {
  int a_idx = 0, b_idx = 0;
  while (a_idx < A.size() && b_idx < B.size()) {
    if (A[a_idx] < B[b_idx]) {
      a_idx++;
    } else if (A[a_idx] > B[b_idx]) {
      b_idx++;
    } else {
      results.push_back(A[a_idx]);
      a_idx++;
      b_idx++;
    }
  }
}

int VectorIntersectionSize(const std::vector<int> &A,
                           const std::vector<int> &B) {
  int intersection_size = 0;
  int a_idx = 0, b_idx = 0;
  while (a_idx < A.size() && b_idx < B.size()) {
    if (A[a_idx] < B[b_idx]) {
      a_idx++;
    } else if (A[a_idx] > B[b_idx]) {
      b_idx++;
    } else {
      intersection_size++;
      a_idx++;
      b_idx++;
    }
  }
  return intersection_size;
}
// note: research log marker 6
