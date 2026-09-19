// ---------------------------------------------------------------------------
// Starter tests.
//
// Three tests are here to show the pattern and to give you a green baseline.
// They are properties of Shannon entropy that hold regardless of how the code
// is implemented, so they stay valid through every refactor you do. That is the
// point: these are the tests you can refactor *behind*.
//
// See TASK.md for what you should be adding here as you go.
// ---------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "entropy/entropy.h"

namespace {

// The dataset layout the code assumes: column-major, so variable `i` of sample
// `s` lives at d[s + i * nsamples]. Building test data by hand is the fastest
// way to find out whether you actually understand the layout.
std::vector<int> makeColumnMajor(const std::vector<std::vector<int>>& columns,
                                 int nsamples) {
  std::vector<int> d;
  d.reserve(columns.size() * static_cast<std::size_t>(nsamples));
  for (const auto& col : columns) {
    for (int s = 0; s < nsamples; ++s) {
      d.push_back(col[static_cast<std::size_t>(s)]);
    }
  }
  return d;
}

constexpr double kTol = 1e-9;

}  // namespace

// A column that never varies carries no information.
TEST(EntropyProperties, ConstantColumnHasZeroEntropy) {
  constexpr int nsamples = 8;
  const std::vector<int> col(nsamples, 1);
  const std::vector<int> d = makeColumnMajor({col}, nsamples);

  int sel[4] = {0, -1, -1, -1};
  const double h = entropyFast(d.data(), nsamples, 1, 0, sel);

  EXPECT_NEAR(h, 0.0, kTol);
}

// A balanced binary column carries exactly one nat (natural log, not base 2).
TEST(EntropyProperties, BalancedBinaryColumnIsLogTwo) {
  constexpr int nsamples = 8;
  const std::vector<int> col = {0, 1, 0, 1, 0, 1, 0, 1};
  const std::vector<int> d = makeColumnMajor({col}, nsamples);

  int sel[4] = {0, -1, -1, -1};
  const double h = entropyFast(d.data(), nsamples, 1, 0, sel);

  EXPECT_NEAR(h, std::log(2.0), kTol);
}

// Joining a variable with itself adds nothing: H(X, X) == H(X).
TEST(EntropyProperties, JointWithSelfEqualsMarginal) {
  constexpr int nsamples = 9;
  const std::vector<int> col = {0, 1, 2, 0, 1, 2, 0, 0, 2};
  const std::vector<int> d = makeColumnMajor({col, col}, nsamples);

  int marginal[4] = {0, -1, -1, -1};
  int joint[4] = {0, 1, -1, -1};

  const double hX = entropyFast(d.data(), nsamples, 2, 0, marginal);
  const double hXX = entropyFast(d.data(), nsamples, 2, 0, joint);

  EXPECT_NEAR(hX, hXX, kTol);
}

// ---------------------------------------------------------------------------
// Your turn. Things worth pinning down before you change anything (TASK.md
// Stage 1 onwards). No hints here on purpose - work out what the right answer
// is from the maths, not from what the current code happens to return.
//
//  - H(X, Y) >= max(H(X), H(Y)) for any two columns
//  - H(X, Y) <= H(X) + H(Y), with equality iff X and Y are independent
//  - a column with k equally likely values has entropy log(k)
//  - the estimator argument `c` selects between four estimators, but only one
//    of them is reachable from the CLI. Pin down what each should return, from
//    the maths - not from what the current code happens to return.
//  - reading the same file twice in one process gives the same data both times
//  - a search over a fixed small dataset produces the same triples, in the
//    same order, whether run on 1 thread or 8
// ---------------------------------------------------------------------------
