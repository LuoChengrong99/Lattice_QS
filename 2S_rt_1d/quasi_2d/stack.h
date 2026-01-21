void hstack(
    const unsigned int* A, const unsigned int* B, unsigned int* C,
    std::size_t M, std::size_t N1, std::size_t N2
) {
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N1; ++j)
            C[i*(N1 + N2) + j] = A[i*N1 + j]; 
        for (std::size_t j = 0; j < N2; ++j)
            C[i*(N1 + N2) + (N1 + j)] = B[i*N2 + j]; 
    }
    return; 
}
void vstack(
    const unsigned int* A, const unsigned int* B, unsigned int* C,
    std::size_t M1, std::size_t M2, std::size_t N
) {
    for (std::size_t i = 0; i < M1; ++i)
        for (std::size_t j = 0; j < N; ++j)
            C[i*N + j] = A[i*N + j]; 

    for (std::size_t i = 0; i < M2; ++i)
        for (std::size_t j = 0; j < N; ++j)
            C[(M1 + i)*N + j] = B[i*N + j]; 
    return; 
}