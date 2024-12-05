#include "fft_distributed.h"
#include "core/get_time.h"

#include <unordered_map>

// Distributed version of the iterative fast fourier transform from ICP
// This function could perform FFT or IFFT, determined by isInverse
// Precondition: size of X and number of processes are power of 2, size of X >= number of processes
// The input X may contain partial or full input array
std::vector<std::complex<double>> icpFftDistributed(const std::vector<std::complex<double>> &X, bool isInverse) {
    int processes, currProcId;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);  // Total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcId);  // Current process ID

	const int n = X.size();
	const int r = std::log2(n);
    const int logP = std::log2(processes);
    const int elePerProc = n / processes;
    const int startIdx = currProcId * elePerProc; // Start index of elements owned by this process
    const int endIdx = startIdx + elePerProc; // Exclusive end index

	const std::complex<double> img(0.0, 1.0);
	double exponentSign = isInverse ? 1.0 : -1.0;
	std::complex<double> omegaBase = std::exp(exponentSign * 2.0 * img * std::numbers::pi / (double)n);

	std::vector<std::complex<double>> R(n); // Result array
    std::copy(X.begin() + startIdx, X.begin() + endIdx, R.begin() + startIdx);

    std::vector<std::complex<double>> S(n); // Auxillary array to hold previous value of R
    std::copy(R.begin() + startIdx, R.begin() + endIdx, S.begin() + startIdx);
    
    std::unordered_map<int, std::complex<double>> omegaCache; // Store previously calculated omega

	// Outer loop O(log n), m represent stage
    for (int m = 0; m < r; ++m) {
        // Check if data exchange is needed at this stage
        if (m < (logP)) {
            std::bitset<32> partnerProcIdBits(currProcId);
            partnerProcIdBits.flip(logP - m - 1);
            int partnerProcId = partnerProcIdBits.to_ulong();
            // Exchange data with partner
            MPI_Sendrecv(
                R.data() + startIdx, elePerProc, MPI_DOUBLE_COMPLEX, partnerProcId, m,
                S.data() + partnerProcId * elePerProc, elePerProc, MPI_DOUBLE_COMPLEX, partnerProcId, m, MPI_COMM_WORLD, MPI_STATUS_IGNORE 
            );
        }

        // Inner loop O(n)
		for (int i = startIdx; i < endIdx; ++i) {
			// // Let (b_0 b_1 ... b_r-1) be the binary representation of i 
            unsigned int j = i & ~(1 << (r - 1 - m)); // Clear the bit, j:= (b_0 ... b_m−1 0 b_m+1 ... b_r−1)
            unsigned int k = i | (1 << (r - 1 - m));  // Set the bit, k:= (b_0 ... b_m−1 1 b_m+1 ... b_r−1)

			unsigned int omegaExp = getFirstNBits(reverseBits(i, r), m + 1) << (r - 1 - m); // (b_m b_m-1 ... b-0 ... 0 ... 0)
            std::complex<double> omega;
            if (omegaCache.find(omegaExp) != omegaCache.end()) {
                omega = omegaCache[omegaExp];
            } else {
                omega = std::pow(omegaBase, omegaExp);
                omegaCache[omegaExp] = omega;
            }
            R[i] = S[j] + S[k] * omega;
		}
        // Update S with new values
        std::copy(R.begin() + startIdx, R.begin() + endIdx, S.begin() + startIdx);
    }

    // All processes gather result from other processes
    MPI_Allgather(R.data() + startIdx, elePerProc, MPI_DOUBLE_COMPLEX, R.data(), elePerProc, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

    // In-place bit-reversal reordering
    for (int i = 0; i < n; ++i) {
        int reversedIndex = reverseBits(i, r);
        if (i < reversedIndex) {
            std::swap(R[i], R[reversedIndex]);
        }
    }

	// Normalize result if performing IFFT
    if (isInverse) {
        for (int i = 0; i < n; ++i) {
            R[i] /= n;
        }
    }

    return R;
}
