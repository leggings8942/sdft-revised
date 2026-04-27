# Knowledge Priority（v1.4.2）

1. Treat `SDFT_Unified_Complete_v1.4.2.html` as the canonical theoretical source.
2. Treat `sdft_unified_v1_4_A_reference.py` as the canonical implementation reference (v1.4.2 affordance-extended version).
3. Treat directly uploaded campaign / facility / regional / customer data as higher priority than generic assumptions.
4. When multiple data files conflict:
   - prefer explicitly observed values over derived summaries
   - prefer recent operational data over stale data
   - preserve uncertainty instead of forcing certainty
5. Separate all outputs into four evidence tiers:
   - **Observed**：directly measured from uploaded data
   - **Derived**：deterministically computed from observed values
   - **Assumed**：analyst-set constants or scenario parameters
   - **Hypothesis**：interpretive or behavioral inference
6. Keep notation stable:
   - S: entropy
   - D: fractal dimension
   - H: Hurst exponent
   - Θ: temperature-like parameter（used in F = U - ΘS; do NOT confuse with 𝓣）
   - 𝓣: regime tension（structural diff + operator diff + gauge distortion + normative_barrier + cognitive_barrier）
   - F = U - ΘS: free energy
   - q: quaternion state
   - W_ij: link operator
   - A: affordance score（0–1; do NOT confuse with A_vec = vector potential in Hamiltonian）
   - A_vec: vector potential in Hamiltonian equations
   - P_behavior: behavior probability = sigmoid(α_A·A + α_F·F + α_H·H - α_T·T_reg)
   - φ_A: affordance potential = ν_A · A
   - E_total: combined electric field = E_struct + γ_A · E_aff（v1.4.2）
7. Version conflict resolution: v1.4.2 always takes precedence over v1.4.1 and earlier versions.
