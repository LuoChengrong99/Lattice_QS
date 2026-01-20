struct Ran {
    unsigned long long u, v, w; 
    Ran(unsigned long long j) : v(4101842887655102017LL), w(1) {
        u = j ^ v; int64(); 
        v = u; int64(); 
        w = v; int64(); 
    }
    inline unsigned long long int64() {
        u = u * 2862933555777941757LL + 7046029254386353087LL; 
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8; 
        w = 4294957665U * (w & 0xffffffff) + (w >> 32); 
        unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4; 
        return (x + v) ^ w; 
    }
    inline double doub() { return 5.42101086242752217E-20 * int64(); } 
    inline unsigned int int32() { return (unsigned int)int64(); } 
}; 
void shuffle_idx(std::vector<unsigned int> &vec_idx, Ran &myran) 
{
    std::size_t size_idx = vec_idx.size(); 
    for(std::size_t i = size_idx - 1; i > 0; i--) 
    { 
        std::size_t rand_j = static_cast<std::size_t>(myran.doub() * (i + 1)); 
        unsigned int temp_idx = vec_idx[i]; 
        vec_idx[i] = vec_idx[rand_j]; 
        vec_idx[rand_j] = temp_idx; 
    }
    return; 
}