#include "astap/quads.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "astap/astro_math.h"

namespace astap {
  namespace {
    // Sorts the six distances in place, largest first. The unrolled comparison
    // network of the original is kept: it is faster than a generic sort for six
    // elements and reproduces the exact same ordering.
    inline void sort6_descending(double &d1, double &d2, double &d3, double &d4, double &d5,
                                 double &d6) {
      double t;
      if (d2 > d1) {
        t = d1;
        d1 = d2;
        d2 = t;
      }
      if (d3 > d2) {
        t = d2;
        d2 = d3;
        d3 = t;
      }
      if (d4 > d3) {
        t = d3;
        d3 = d4;
        d4 = t;
      }
      if (d5 > d4) {
        t = d4;
        d4 = d5;
        d5 = t;
      }
      if (d6 > d5) {
        t = d5;
        d5 = d6;
        d6 = t;
      }

      if (d2 > d1) {
        t = d1;
        d1 = d2;
        d2 = t;
      }
      if (d3 > d2) {
        t = d2;
        d2 = d3;
        d3 = t;
      }
      if (d4 > d3) {
        t = d3;
        d3 = d4;
        d4 = t;
      }
      if (d5 > d4) {
        t = d4;
        d4 = d5;
        d5 = t;
      }

      if (d2 > d1) {
        t = d1;
        d1 = d2;
        d2 = t;
      }
      if (d3 > d2) {
        t = d2;
        d2 = d3;
        d3 = t;
      }
      if (d4 > d3) {
        t = d3;
        d3 = d4;
        d4 = t;
      }

      if (d2 > d1) {
        t = d1;
        d1 = d2;
        d2 = t;
      }
      if (d3 > d2) {
        t = d2;
        d2 = d3;
        d3 = t;
      }

      if (d2 > d1) {
        t = d1;
        d1 = d2;
        d2 = t;
      }
    }

    // All C(n,4) index combinations used by find_many_quads, in the same order as
    // the case statement of the original.
    const int kQuadIndices5[5][4] = {
      {0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}, {1, 2, 3, 4},
    };
    const int kQuadIndices6[15][4] = {
      {0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 2, 5}, {0, 1, 3, 4}, {0, 1, 3, 5},
      {0, 1, 4, 5}, {0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5},
      {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}, {2, 3, 4, 5},
    };
    const int kQuadIndices7[35][4] = {
      {0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 2, 5}, {0, 1, 2, 6}, {0, 1, 3, 4}, {0, 1, 3, 5},
      {0, 1, 3, 6}, {0, 1, 4, 5}, {0, 1, 4, 6}, {0, 1, 5, 6}, {0, 2, 3, 4}, {0, 2, 3, 5},
      {0, 2, 3, 6}, {0, 2, 4, 5}, {0, 2, 4, 6}, {0, 2, 5, 6}, {0, 3, 4, 5}, {0, 3, 4, 6},
      {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 3, 6}, {1, 2, 4, 5},
      {1, 2, 4, 6}, {1, 2, 5, 6}, {1, 3, 4, 5}, {1, 3, 4, 6}, {1, 3, 5, 6}, {1, 4, 5, 6},
      {2, 3, 4, 5}, {2, 3, 4, 6}, {2, 3, 5, 6}, {2, 4, 5, 6}, {3, 4, 5, 6},
    };
  } // namespace

  void quicksort_starlist(RowList &a, long lo, long hi) {
    if (lo >= hi) return;
    long Lo = lo;
    long Hi = hi;
    double pivot = a(0, static_cast<size_t>((Lo + Hi) / 2));
    do {
      while (a(0, static_cast<size_t>(Lo)) < pivot) Lo++; // sort in X only
      while (a(0, static_cast<size_t>(Hi)) > pivot) Hi--;
      if (Lo <= Hi) {
        std::swap(a(0, static_cast<size_t>(Lo)), a(0, static_cast<size_t>(Hi)));
        std::swap(a(1, static_cast<size_t>(Lo)), a(1, static_cast<size_t>(Hi)));
        Lo++;
        Hi--;
      }
    } while (Lo <= Hi);
    if (Hi > lo) quicksort_starlist(a, lo, Hi);
    if (Lo < hi) quicksort_starlist(a, Lo, hi);
  }

  void find_many_quads(const RowList &starlist, RowList &quads, int mode) {
    const size_t nrstars = starlist.count();

    int num_closest;
    int num_quads_per_group;
    const int (*quad_indices)[4];
    switch (mode) {
      case 6:
        num_closest = 6;
        num_quads_per_group = 15; // C(6,4)
        quad_indices = kQuadIndices6;
        break;
      case 7:
        num_closest = 7;
        num_quads_per_group = 35; // C(7,4)
        quad_indices = kQuadIndices7;
        break;
      case 5:
      default:
        num_closest = 5;
        num_quads_per_group = 5;
        quad_indices = kQuadIndices5;
        break;
    }

    size_t nrquads = 0;
    quads.resize(8, nrstars * static_cast<size_t>(num_quads_per_group));

    std::vector<long> closest_indices(static_cast<size_t>(num_closest));
    std::vector<double> closest_distances(static_cast<size_t>(num_closest));

    const double *StarsX = starlist.data(0);
    const double *StarsY = starlist.data(1);

    for (size_t i = 0; i < nrstars; i++) {
      closest_distances[0] = 0; // reference star distance is zero
      closest_indices[0] = static_cast<long>(i);
      for (int j = 1; j < num_closest; j++) {
        closest_indices[static_cast<size_t>(j)] = -1;
        closest_distances[static_cast<size_t>(j)] = 1E99;
      }

      double x1 = StarsX[i]; // reference star
      double y1 = StarsY[i];

      for (size_t j = 0; j < nrstars; j++) {
        if (i == j) continue;
        double dx = StarsX[j] - x1;
        double dy = StarsY[j] - y1;
        double distance = dx * dx + dy * dy;
        if (distance <= 1) continue; // not an identical star

        // Insertion sort: find the position, then shift.
        int insert_pos = -1;
        for (int k = num_closest - 1; k >= 1; k--) {
          if (distance < closest_distances[static_cast<size_t>(k)])
            insert_pos = k;
          else
            break;
        }
        if (insert_pos >= 0) {
          for (int k = num_closest - 1; k >= insert_pos + 1; k--) {
            closest_distances[static_cast<size_t>(k)] = closest_distances[static_cast<size_t>(k - 1)];
            closest_indices[static_cast<size_t>(k)] = closest_indices[static_cast<size_t>(k - 1)];
          }
          closest_distances[static_cast<size_t>(insert_pos)] = distance;
          closest_indices[static_cast<size_t>(insert_pos)] = static_cast<long>(j);
        }
      }

      if (closest_indices[static_cast<size_t>(num_closest - 1)] == -1) continue;

      for (int q = 0; q < num_quads_per_group; q++) {
        double x1q = StarsX[closest_indices[static_cast<size_t>(quad_indices[q][0])]];
        double y1q = StarsY[closest_indices[static_cast<size_t>(quad_indices[q][0])]];
        double x2 = StarsX[closest_indices[static_cast<size_t>(quad_indices[q][1])]];
        double y2 = StarsY[closest_indices[static_cast<size_t>(quad_indices[q][1])]];
        double x3 = StarsX[closest_indices[static_cast<size_t>(quad_indices[q][2])]];
        double y3 = StarsY[closest_indices[static_cast<size_t>(quad_indices[q][2])]];
        double x4 = StarsX[closest_indices[static_cast<size_t>(quad_indices[q][3])]];
        double y4 = StarsY[closest_indices[static_cast<size_t>(quad_indices[q][3])]];

        double xt = (x1q + x2 + x3 + x4) * 0.25; // quad centre
        double yt = (y1q + y2 + y3 + y4) * 0.25;

        bool identical_quad = false;
        for (size_t k = 0; k < nrquads; k++) {
          if (std::fabs(xt - quads(6, k)) < 1 && std::fabs(yt - quads(7, k)) < 1) {
            identical_quad = true;
            break;
          }
        }
        if (identical_quad) continue;

        double dx, dy;
        dx = x1q - x2;
        dy = y1q - y2;
        double dist1 = std::sqrt(dx * dx + dy * dy);
        dx = x1q - x3;
        dy = y1q - y3;
        double dist2 = std::sqrt(dx * dx + dy * dy);
        dx = x1q - x4;
        dy = y1q - y4;
        double dist3 = std::sqrt(dx * dx + dy * dy);
        dx = x2 - x3;
        dy = y2 - y3;
        double dist4 = std::sqrt(dx * dx + dy * dy);
        dx = x2 - x4;
        dy = y2 - y4;
        double dist5 = std::sqrt(dx * dx + dy * dy);
        dx = x3 - x4;
        dy = y3 - y4;
        double dist6 = std::sqrt(dx * dx + dy * dy);

        sort6_descending(dist1, dist2, dist3, dist4, dist5, dist6);

        quads(0, nrquads) = dist1;
        quads(1, nrquads) = dist2 / dist1;
        quads(2, nrquads) = dist3 / dist1;
        quads(3, nrquads) = dist4 / dist1;
        quads(4, nrquads) = dist5 / dist1;
        quads(5, nrquads) = dist6 / dist1;
        quads(6, nrquads) = xt;
        quads(7, nrquads) = yt;
        nrquads++;
      }
    }

    quads.resize(8, nrquads);
  }

  void find_quads(int nrstars_image, RowList &starlist, RowList &quads) {
    constexpr int kBucketCapacity = 5; // max quads per bucket, grows when needed
    constexpr double kGridInv = 0.2; // pre-calculated inverse of grid_size (1.0 / 5.0)

    size_t nrstars = starlist.count();

    // Base the quad group size selection on the number of stars in the *image*
    // and not on the number of database stars, since the database field could be
    // larger.
    if (nrstars_image < 15 && nrstars > 6) {
      find_many_quads(starlist, quads, 7); // 35 times more quads
      return;
    }
    if (nrstars_image < 30 && nrstars > 5) {
      find_many_quads(starlist, quads, 6); // 15 times more quads
      return;
    }
    if (nrstars_image < 60 && nrstars > 4) {
      find_many_quads(starlist, quads, 5); // 5 times more quads
      return;
    }

    if (nrstars < 4) {
      // not enough stars for quads
      quads.resize(8, 0);
      return;
    }

    size_t bandw;
    if (nrstars >= 150) {
      quicksort_starlist(starlist, 0, static_cast<long>(nrstars) - 1); // sort in X only
      // The resulting tolerance band is about twice the average star distance,
      // assuming the stars are equally distributed.
      bandw = static_cast<size_t>(pround(2 * std::sqrt(static_cast<double>(nrstars))));
    } else {
      bandw = nrstars; // switch off the pre-filtering in X
    }

    const double *StarsX = starlist.data(0);
    const double *StarsY = starlist.data(1);

    // In hash table design the number of buckets is often set to 1-2 times the
    // expected number of entries, giving a load factor of 0.5-1.0.
    const size_t hash_table_len = nrstars * 2;
    std::vector<std::vector<size_t> > hash_table(hash_table_len);
    for (auto &b: hash_table) b.reserve(kBucketCapacity);

    size_t nrquads = 0;
    quads.resize(8, nrstars);

    for (size_t i = 0; i < nrstars; i++) {
      double distance1 = 1E99; // distance closest star
      double distance2 = 1E99; // distance second closest star
      double distance3 = 1E99; // distance third closest star
      size_t j_index1 = 0, j_index2 = 0, j_index3 = 0;

      size_t Sstart = i > bandw ? i - bandw : 0;
      // Search in a limited X band only. The star list is sorted in X, which
      // increases the search speed by about 30%.
      size_t Send = std::min(nrstars - 1, i + bandw);

      double x1 = StarsX[i]; // first star position of the quad
      double y1 = StarsY[i];

      for (size_t j = Sstart; j <= Send; j++) {
        if (j == i) continue; // do not check the star with itself
        double disty = sqr(StarsY[j] - y1);
        if (disty >= distance3) continue; // pre-check to increase the processing speed
        double distance = sqr(StarsX[j] - x1) + disty; // square distances are used
        if (distance <= 1) continue; // not an identical star

        if (distance < distance1) {
          distance3 = distance2;
          j_index3 = j_index2;
          distance2 = distance1;
          j_index2 = j_index1;
          distance1 = distance;
          j_index1 = j;
        } else if (distance < distance2) {
          distance3 = distance2;
          j_index3 = j_index2;
          distance2 = distance;
          j_index2 = j;
        } else if (distance < distance3) {
          distance3 = distance;
          j_index3 = j;
        }
      }

      if (distance3 >= 1E99) continue; // did not find 4 stars in the restricted area

      double x2 = StarsX[j_index1]; // second star position of the quad
      double y2 = StarsY[j_index1];
      double x3 = StarsX[j_index2];
      double y3 = StarsY[j_index2];
      double x4 = StarsX[j_index3];
      double y4 = StarsY[j_index3];

      double xt = (x1 + x2 + x3 + x4) * 0.25; // mean x position of the quad
      double yt = (y1 + y2 + y3 + y4) * 0.25; // mean y position of the quad

      // Check for a duplicate quad using the hash table. This makes the routine
      // about 25% faster than comparing against all quads found so far.
      long hash_x = ptrunc(xt * kGridInv);
      long hash_y = ptrunc(yt * kGridInv);
      size_t idx = static_cast<size_t>(std::labs(hash_x * 31 + hash_y)) % hash_table_len;

      bool identical_quad = false;
      for (size_t k: hash_table[idx]) {
        if (std::fabs(xt - quads(6, k)) < 1 && std::fabs(yt - quads(7, k)) < 1) {
          identical_quad = true;
          break;
        }
      }
      if (identical_quad) continue;

      double dist1 = std::sqrt(distance1); // star1-star2, reuse the value already calculated
      double dist2 = std::sqrt(distance2); // star1-star3
      double dist3 = std::sqrt(distance3); // star1-star4

      double dx, disty;
      dx = x2 - x3;
      disty = y2 - y3;
      double dist4 = std::sqrt(dx * dx + disty * disty);
      dx = x2 - x4;
      disty = y2 - y4;
      double dist5 = std::sqrt(dx * dx + disty * disty);
      dx = x3 - x4;
      disty = y3 - y4;
      double dist6 = std::sqrt(dx * dx + disty * disty);

      sort6_descending(dist1, dist2, dist3, dist4, dist5, dist6);

      quads(0, nrquads) = dist1; // largest distance
      quads(1, nrquads) = dist2 / dist1; // scale relative to the largest distance
      quads(2, nrquads) = dist3 / dist1;
      quads(3, nrquads) = dist4 / dist1;
      quads(4, nrquads) = dist5 / dist1;
      quads(5, nrquads) = dist6 / dist1;
      quads(6, nrquads) = xt;
      quads(7, nrquads) = yt;

      hash_table[idx].push_back(nrquads);
      nrquads++;
    }

    quads.resize(8, nrquads); // adapt to the number found
  }
} // namespace astap
