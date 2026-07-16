"""
KAIKU -- is the whisper arithmetic, or just a magnitude-distribution echo?
==========================================================================
(kaiku: Finnish for echo. A whisper that is only a reflection of the
weight histogram is an echo, not a voice.)

Jyrkanne found a reproducible ~3.5-3.8 sigma lift at log5/log7 in the
prime-weighted spectrum vs a flat control. The registered suspicion: it is
a Kac-Rice effect -- a beat between the HEAVY-TAILED magnitude DISTRIBUTION
of {log p / sqrt(p)} and the level density -- NOT the arithmetic of the
primes. This file discriminates, with a control ladder built to be the
decisive test.

The logic. The prime ladder is TWO things at once:
  (1) a MULTISET of magnitudes {log p / sqrt(p) : p prime} -- a histogram.
  (2) an ORDERING of those magnitudes onto edges by prime index.
Kac-Rice cares only about (1) the histogram. Arithmetic cares about (2)
the specific values/ordering being the actual primes.

The discriminator -- four weight sets, ALL with the SAME MULTISET of
magnitudes (so identical histogram, identical Kac-Rice content), differing
only in what the magnitudes MEAN:

  PRIME     : w_e = log p_e / sqrt(p_e), p_e = e-th prime (Jyrkanne's set)
  SHUFFLE   : the SAME multiset, randomly permuted onto edges
              -> identical histogram, arithmetic ordering destroyed.
  SORTED    : the same multiset, sorted descending along a fixed edge order
              -> identical histogram, monotone (non-arithmetic) ordering.
  SURROGATE : a DIFFERENT multiset with the SAME shape -- draw w ~ log x /
              sqrt(x) for x uniform, matched mean/var -- primes replaced by
              a smooth density of the same tail.
              -> nearly identical histogram, NO primes at all.

Prediction table (registered):
  If the whisper is KAC-RICE (magnitude distribution):
     SHUFFLE, SORTED, SURROGATE all reproduce the log5/log7 lift within
     ~1 sigma of PRIME. The lift tracks the histogram, not the arithmetic.
  If the whisper is ARITHMETIC (the primes qua primes):
     PRIME shows the lift; SHUFFLE / SORTED / SURROGATE do NOT (lift drops
     by >= 2 sigma). Only the correctly-identified primes sing.

=== REGISTERED PREDICTION (before running) ===
K1  KAC-RICE (the bet). Mean log5/log7 sigma for SHUFFLE, SORTED, SURROGATE
    each within 1.5 sigma of PRIME's. -> the whisper is a magnitude-
    distribution echo; it is NOT about primes; Jyrkanne's hint is
    demoted to an artifact and the honest ledger records it as such.
  ARITHMETIC (the long shot). PRIME's log5/log7 lift exceeds each of
    SHUFFLE/SORTED/SURROGATE by >= 2 sigma. -> the ordering matters, the
    primes qua primes are doing something. If THIS fires, re-run 5 seeds
    and, if it holds, it is a genuine (small) arithmetic signal and a
    much longer look is warranted.

6 seeds per weight set, N=400, full Dynamo dynamics (current-driven A from
zero), alkuluku prime-trace on each final spectrum. Report per-set mean
sigma at log2/3/5/7 with seed scatter.

Do not hype. Do not lie. Just show. Betting on the echo -- and measuring
it rather than assuming it.
"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
# reuse Jyrkanne machinery
exec(open(os.path.join(os.path.dirname(__file__), 'jyrkanne.py')).read().split('if __name__')[0])

def prime_multiset(ne):
    ps = primes_upto_count(ne)
    w = np.log(ps) / np.sqrt(ps)
    return w / w.mean()

def edge_list(A):
    N = A.shape[0]
    iu = np.array(np.triu_indices(N, 1))
    mask = A[iu[0], iu[1]] != 0
    return iu[0][mask], iu[1][mask]

def weights_from_vector(A, wvec, ei, ej):
    N = A.shape[0]; W = np.zeros((N, N))
    W[ei, ej] = wvec; W[ej, ei] = wvec
    return W

def build_weightset(kind, A, ei, ej, rng):
    ne = len(ei)
    base = prime_multiset(ne)                    # the PRIME multiset (index order)
    if kind == 'prime':
        w = base
    elif kind == 'shuffle':
        w = rng.permutation(base)                # same multiset, random order
    elif kind == 'sorted':
        w = np.sort(base)[::-1]                  # same multiset, monotone order
    elif kind == 'surrogate':
        x = rng.uniform(2, ne + 2, ne)           # smooth density, no primes
        w = np.log(x) / np.sqrt(x); w = w / w.mean()
        # match mean/var to base for a fair histogram
        w = (w - w.mean()) / w.std() * base.std() + base.mean()
        w = np.clip(w, 1e-3, None)
    else:
        raise ValueError(kind)
    return weights_from_vector(A, w, ei, ej)

def run_set(kind, seed, N=400, K=0.1, mu=0.05, equil=1000, drive=2500):
    G, Ag, tris = make(N=N, seed=seed)
    ei, ej = edge_list(Ag)
    rng = np.random.default_rng(7000 + seed)
    C = build_weightset(kind, Ag, ei, ej, rng)
    r2 = np.random.default_rng(200 + seed)
    th = r2.uniform(0, 2*np.pi, N); tho = th.copy()
    def phase_step(th, tho, Aconn):
        d = th[:,None]-th[None,:]-Aconn
        dd = np.arctan2(np.sin(d), np.cos(d))
        beta = np.sum(C*dd**2, axis=1); gamma = 1/(1+TAU*beta+1e-12)**2
        force = np.sum(C*np.sin(-dd), axis=1)
        v = (1-0.05*DT)*(th-tho)
        return np.mod(th+v+DT**2*gamma**2*force, 2*np.pi), th.copy()
    Az = np.zeros((N,N))
    for _ in range(equil):
        th, tho = phase_step(th, tho, Az)
    Aconn = np.zeros((N,N)); edge = (Ag != 0)
    for _ in range(drive):
        th, tho = phase_step(th, tho, Aconn)
        d = th[:,None]-th[None,:]-Aconn
        Jc = C*np.sin(d); Jc = 0.5*(Jc-Jc.T)
        cc = curlcurl(Aconn, tris, N)
        Aconn = (Aconn + DT*(Jc - K*cc - mu*Aconn))*edge
        Aconn = 0.5*(Aconn-Aconn.T)*edge
    d = th[:,None]-th[None,:]-Aconn
    M = -C*np.cos(d)*np.exp(1j*Aconn); M = 0.5*(M+M.conj().T)
    np.fill_diagonal(M, np.sum(C*np.cos(d), axis=1))
    return prime_trace_sigmas(eigvalsh(M))

if __name__ == '__main__':
    kinds = ['prime', 'shuffle', 'sorted', 'surrogate']
    seeds = [1, 2, 3, 4, 5, 6]
    freqs = ['log2', 'log3', 'log5', 'log7']
    data = {k: {f: [] for f in freqs} for k in kinds}
    for seed in seeds:
        for kind in kinds:
            sig = run_set(kind, seed)
            for f in freqs:
                data[kind][f].append(sig[f])
        print(f"seed {seed}: " + "  ".join(
            f"{k}[log5={np.mean(data[k]['log5'][-1:]):+.1f},"
            f"log7={np.mean(data[k]['log7'][-1:]):+.1f}]" for k in kinds))

    summary = {}
    for k in kinds:
        summary[k] = {f: dict(mean=round(float(np.mean(data[k][f])), 2),
                              std=round(float(np.std(data[k][f])), 2))
                      for f in freqs}
    print("\n=== per-set mean sigma (6 seeds) ===")
    for k in kinds:
        print(f"{k:10s}  " + "  ".join(
            f"{f}={summary[k][f]['mean']:+.2f}+/-{summary[k][f]['std']:.2f}" for f in freqs))

    # discriminator: compare PRIME vs each control at log5, log7
    def lift_gap(f):
        pm = np.mean(data['prime'][f])
        gaps = {c: pm - np.mean(data[c][f]) for c in ['shuffle', 'sorted', 'surrogate']}
        return gaps
    g5, g7 = lift_gap('log5'), lift_gap('log7')
    # pooled seed std for a rough significance of the gap
    def pooled_std(f, c):
        return np.sqrt(np.var(data['prime'][f]) + np.var(data[c][f]))
    arithmetic = all(
        (np.mean(data['prime'][f]) - np.mean(data[c][f])) >= 2.0
        for f in ['log5', 'log7'] for c in ['shuffle', 'sorted', 'surrogate'])
    kacrice = all(
        abs(np.mean(data['prime'][f]) - np.mean(data[c][f])) <= 1.5
        for f in ['log5', 'log7'] for c in ['shuffle', 'sorted', 'surrogate'])
    verdict = dict(
        prime_log5=summary['prime']['log5']['mean'],
        prime_log7=summary['prime']['log7']['mean'],
        gap_log5_vs_controls={c: round(g5[c], 2) for c in g5},
        gap_log7_vs_controls={c: round(g7[c], 2) for c in g7},
        K1_KACRICE_echo=bool(kacrice),
        ARITHMETIC_signal=bool(arithmetic),
        HEADLINE=(
            "ARITHMETIC: the correctly-ordered primes lift log5/log7 above "
            "same-histogram controls by >=2 sigma -- re-run 5 seeds, this is "
            "a genuine (small) arithmetic signal" if arithmetic else
            "KAC-RICE ECHO: same-histogram controls (shuffle/sorted/surrogate) "
            "reproduce the log5/log7 lift within scatter -- the whisper is a "
            "magnitude-DISTRIBUTION artifact, NOT arithmetic. Jyrkanne's hint "
            "is demoted to an echo. The cliff is unscaled and the near side "
            "holds no secret arithmetic." if kacrice else
            "MIXED -- see per-set numbers; neither clean threshold met"))
    print("\n" + json.dumps(verdict, indent=2))
    json.dump(dict(summary=summary, verdict=verdict,
                   raw={k: data[k] for k in kinds}),
              open(os.path.join(os.path.dirname(__file__), '..', 'results',
                                'kaiku_results.json'), 'w'), indent=1)
