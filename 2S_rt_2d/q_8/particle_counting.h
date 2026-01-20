#include "particle_state.h" 
void initial_setting(
    const unsigned int Lx, const unsigned int Ly, const unsigned int N, const unsigned int NA, 
    std::vector<ParticleStates>& sigma, std::vector<int>& x, std::vector<int>& y, 
    std::vector<unsigned int>& R1A, std::vector<unsigned int>& R2A, std::vector<unsigned int>& R3A, std::vector<unsigned int>& R4A, 
    std::vector<unsigned int>& R5A, std::vector<unsigned int>& R6A, std::vector<unsigned int>& R7A, std::vector<unsigned int>& R8A, std::vector<unsigned int>& TA, 
    std::vector<unsigned int>& R1B, std::vector<unsigned int>& R2B, std::vector<unsigned int>& R3B, std::vector<unsigned int>& R4B, 
    std::vector<unsigned int>& R5B, std::vector<unsigned int>& R6B, std::vector<unsigned int>& R7B, std::vector<unsigned int>& R8B, std::vector<unsigned int>& TB, 
    Ran& myran 
) 
{
    for (std::size_t i = 0; i < N; ++i) 
    {
        double r = myran.doub(); 
        x[i] = static_cast<int> (myran.doub() * Lx); 
        y[i] = static_cast<int> (myran.doub() * Ly); 
        unsigned int ii = static_cast<unsigned int>(y[i] * Lx + x[i]); 
        if (i < NA) 
        {
            if (r < 0.5 * 1/8) 
            {
                R1A[ii] += 1; 
                sigma[i] = ParticleStates::Running_1_state; 
            }
            else if (r < 0.5 * 2/8) 
            {
                R2A[ii] += 1; 
                sigma[i] = ParticleStates::Running_2_state; 
            }
            else if (r < 0.5 * 3/8) 
            {
                R3A[ii] += 1; 
                sigma[i] = ParticleStates::Running_3_state; 
            }
            else if (r < 0.5 * 4/8) 
            {
                R4A[ii] += 1; 
                sigma[i] = ParticleStates::Running_4_state; 
            }
            else if (r < 0.5 * 5/8) 
            {
                R5A[ii] += 1; 
                sigma[i] = ParticleStates::Running_5_state; 
            }
            else if (r < 0.5 * 6/8) 
            {
                R6A[ii] += 1; 
                sigma[i] = ParticleStates::Running_6_state; 
            }
            else if (r < 0.5 * 7/8) 
            {
                R7A[ii] += 1; 
                sigma[i] = ParticleStates::Running_7_state; 
            }
            else if (r < 0.5 * 8/8) 
            {
                R8A[ii] += 1; 
                sigma[i] = ParticleStates::Running_8_state; 
            }
            else 
            {
                TA[ii] += 1; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
        }
        else 
        {
            if (r < 0.5 * 1/8) 
            {
                R1B[ii] += 1; 
                sigma[i] = ParticleStates::Running_1_state; 
            }
            else if (r < 0.5 * 2/8) 
            {
                R2B[ii] += 1; 
                sigma[i] = ParticleStates::Running_2_state; 
            }
            else if (r < 0.5 * 3/8) 
            {
                R3B[ii] += 1; 
                sigma[i] = ParticleStates::Running_3_state; 
            }
            else if (r < 0.5 * 4/8) 
            {
                R4B[ii] += 1; 
                sigma[i] = ParticleStates::Running_4_state; 
            }
            else if (r < 0.5 * 5/8) 
            {
                R5B[ii] += 1; 
                sigma[i] = ParticleStates::Running_5_state; 
            }
            else if (r < 0.5 * 6/8) 
            {
                R6B[ii] += 1; 
                sigma[i] = ParticleStates::Running_6_state; 
            }
            else if (r < 0.5 * 7/8) 
            {
                R7B[ii] += 1; 
                sigma[i] = ParticleStates::Running_7_state; 
            }
            else if (r < 0.5 * 8/8) 
            {
                R8B[ii] += 1; 
                sigma[i] = ParticleStates::Running_8_state; 
            }
            else 
            {
                TB[ii] += 1; 
                sigma[i] = ParticleStates::Tumbling_0_state; 
            }
        }
    }

    return; 
}
void counting_running_A(
    const unsigned int Lx, const unsigned int Ly, 
    const int jj, const int x, const int y, 
    const ParticleStates s, const ParticleStates s_next, 
    std::vector<unsigned int>& R1A, std::vector<unsigned int>& R2A, std::vector<unsigned int>& R3A, std::vector<unsigned int>& R4A, 
    std::vector<unsigned int>& R5A, std::vector<unsigned int>& R6A, std::vector<unsigned int>& R7A, std::vector<unsigned int>& R8A, std::vector<unsigned int>& TA 
) 
{
    if (s == ParticleStates::Running_1_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R1A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R1A[jj] -= 1; 
            R1A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_2_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R2A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R2A[jj] -= 1; 
            R2A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_3_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R3A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R3A[jj] -= 1; 
            R3A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_4_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R4A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R4A[jj] -= 1; 
            R4A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_5_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R5A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R5A[jj] -= 1; 
            R5A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_6_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R6A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R6A[jj] -= 1; 
            R6A[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_7_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R7A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R7A[jj] -= 1; 
            R7A[y * Lx + x] += 1; 
        }
    }
    else 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R8A[jj] -= 1; 
            TA[y * Lx + x] += 1; 
        }
        else 
        {
            R8A[jj] -= 1; 
            R8A[y * Lx + x] += 1; 
        }
    }

    return; 
}
void counting_tumbling_A(
    const unsigned int Lx, const unsigned int Ly, 
    const int jj, const int x, const int y, 
    const ParticleStates s, const ParticleStates s_next, 
    std::vector<unsigned int>& R1A, std::vector<unsigned int>& R2A, std::vector<unsigned int>& R3A, std::vector<unsigned int>& R4A, 
    std::vector<unsigned int>& R5A, std::vector<unsigned int>& R6A, std::vector<unsigned int>& R7A, std::vector<unsigned int>& R8A, std::vector<unsigned int>& TA 
)
{
    if (s_next == ParticleStates::Running_1_state) 
    {
        TA[jj] -= 1; 
        R1A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_2_state) 
    {
        TA[jj] -= 1; 
        R2A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_3_state) 
    {
        TA[jj] -= 1; 
        R3A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_4_state) 
    {
        TA[jj] -= 1; 
        R4A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_5_state) 
    {
        TA[jj] -= 1; 
        R5A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_6_state) 
    {
        TA[jj] -= 1; 
        R6A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_7_state) 
    {
        TA[jj] -= 1; 
        R7A[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_8_state) 
    {
        TA[jj] -= 1; 
        R8A[y * Lx + x] += 1; 
    }
    else 
    {
        TA[jj] -= 1; 
        TA[y * Lx + x] += 1; 
    }

    return; 
} 
void counting_running_B(
    const unsigned int Lx, const unsigned int Ly, 
    const int jj, const int x, const int y, 
    const ParticleStates s, const ParticleStates s_next, 
    std::vector<unsigned int>& R1B, std::vector<unsigned int>& R2B, std::vector<unsigned int>& R3B, std::vector<unsigned int>& R4B, 
    std::vector<unsigned int>& R5B, std::vector<unsigned int>& R6B, std::vector<unsigned int>& R7B, std::vector<unsigned int>& R8B, std::vector<unsigned int>& TB 
)
{
    if (s == ParticleStates::Running_1_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R1B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R1B[jj] -= 1; 
            R1B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_2_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R2B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R2B[jj] -= 1; 
            R2B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_3_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R3B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R3B[jj] -= 1; 
            R3B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_4_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R4B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R4B[jj] -= 1; 
            R4B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_5_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R5B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R5B[jj] -= 1; 
            R5B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_6_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R6B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R6B[jj] -= 1; 
            R6B[y * Lx + x] += 1; 
        }
    }
    else if (s == ParticleStates::Running_7_state) 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R7B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R7B[jj] -= 1; 
            R7B[y * Lx + x] += 1; 
        }
    }
    else 
    {
        if (s_next == ParticleStates::Tumbling_0_state) 
        {
            R8B[jj] -= 1; 
            TB[y * Lx + x] += 1; 
        }
        else 
        {
            R8B[jj] -= 1; 
            R8B[y * Lx + x] += 1; 
        }
    }

    return; 
}
void counting_tumbling_B(
    const unsigned int Lx, const unsigned int Ly, 
    const int jj, const int x, const int y, 
    const ParticleStates s, const ParticleStates s_next, 
    std::vector<unsigned int>& R1B, std::vector<unsigned int>& R2B, std::vector<unsigned int>& R3B, std::vector<unsigned int>& R4B, 
    std::vector<unsigned int>& R5B, std::vector<unsigned int>& R6B, std::vector<unsigned int>& R7B, std::vector<unsigned int>& R8B, std::vector<unsigned int>& TB 
)
{
    if (s_next == ParticleStates::Running_1_state) 
    {
        TB[jj] -= 1; 
        R1B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_2_state) 
    {
        TB[jj] -= 1; 
        R2B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_3_state) 
    {
        TB[jj] -= 1; 
        R3B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_4_state) 
    {
        TB[jj] -= 1; 
        R4B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_5_state) 
    {
        TB[jj] -= 1; 
        R5B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_6_state) 
    {
        TB[jj] -= 1; 
        R6B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_7_state) 
    {
        TB[jj] -= 1; 
        R7B[y * Lx + x] += 1; 
    }
    else if (s_next == ParticleStates::Running_8_state) 
    {
        TB[jj] -= 1; 
        R8B[y * Lx + x] += 1; 
    }
    else 
    {
        TB[jj] -= 1; 
        TB[y * Lx + x] += 1; 
    } 

    return; 
} 