#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void fill_array(NumericVector A_array,
                IntegerVector dims,
                IntegerVector from,
                IntegerVector to,
                IntegerVector time,
                NumericVector value,
                int net_idx,
                int trial_idx,
                bool symmetric) {

    int N = from.size();
    int D0 = dims[0]; // networks
    int D1 = dims[1]; // trials
    int D2 = dims[2]; // timesteps
    int D3 = dims[3]; // from
    int D4 = dims[4]; // to

    for (int i = 0; i < N; ++i) {
        int f = from[i] - 1;
        int t_ = to[i] - 1;
        int ts = time[i] - 1;
        int k = trial_idx;
        int n = net_idx;

        if (f < 0 || t_ < 0 || ts < 0 || f >= D3 || t_ >= D4 || ts >= D2) continue;

        int index = t_ * D3 * D2 * D1 * D0 +
            f * D2 * D1 * D0 +
            ts * D1 * D0 +
            k * D0 +
            n;

        A_array[index] = value[i];

        if (symmetric && f != t_) {
            int index_sym = f * D3 * D2 * D1 * D0 +
                t_ * D2 * D1 * D0 +
                ts * D1 * D0 +
                k * D0 +
                n;
            A_array[index_sym] = value[i];
        }
    }
}
