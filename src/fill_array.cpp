#include <Rcpp.h>
using namespace Rcpp;

//' @useDynLib STbayes, .registration = TRUE
//' @export
// [[Rcpp::export]]
void fill_array(NumericVector A_array,
                IntegerVector dims,
                IntegerVector focal,
                IntegerVector other,
                IntegerVector time,
                NumericVector value,
                int net_idx,
                int trial_idx,
                bool symmetric) {

    int N = focal.size();
    int D0 = dims[0]; // networks
    int D1 = dims[1]; // trials
    int D2 = dims[2]; // timesteps
    int D3 = dims[3]; // focal
    int D4 = dims[4]; // other

    for (int i = 0; i < N; ++i) {
        int f = focal[i] - 1;
        int t_ = other[i] - 1;
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
