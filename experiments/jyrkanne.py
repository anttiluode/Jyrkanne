"""
JYRKANNE -- the amplitude wall (jyrkanne: Finnish for cliff)
============================================================
The last experiment firable from the near side of the RH cliff, with the
verdict conditions written to measure what it CAN decide, not what one
would like it to.

What is already done (do not re-run):
  Dynamo: the phase current sin(th_i - th_j - A_ij) already drives A; the
  glass already self-generates flux and self-lifts to GUE. "Close the loop"
  is finished. This file does the ONE new thing: force the couplings toward
  the prime ladder log p / sqrt(p) and ask whether the self-generated GUE
  survives contact with arithmetic.

The honest question (Berry-Keating's amplitude-hermiticity tension):
  When the bare coupling magnitudes are pushed from flat -> prime-ladder,
  does the self-generated flux keep the spectrum REAL and GUE, or do
  eigenvalues peel off the real line?

What this CAN show (registered):
  - whether the amplitude-hermiticity tension bites in this medium, and if
    so at what interpolation strength (where reality breaks).
  - whether prime-shaped weights preserve, degrade, or destroy the GUE
    class the dynamo produces.
What this CANNOT show, stated so it is never claimed:
  - it does NOT realize Lemma 5.2 (exact prime amplitudes on primitive
    loops with composite-loop cancellation). A taper toward the ladder is
    not the explicit formula.
  - GUE + real + prime-magnitudes is NECESSARY, not sufficient. The
    ARBITER of "knows the primes" is the prime-trace test (alkuluku), NOT
    the r-statistic. So alkuluku runs on the final spectrum, and a GUE-but-
    no-prime-peaks outcome is reported as generic chaos in arithmetic
    clothing -- the null this ecosystem expects.

Setup: BirthOfClockfield graph. Map each edge to a "prime" by ladder order;
give it target magnitude w_p = log p / sqrt(p). Interpolate the coupling
    C_ij(s) = (1-s) * flat + s * prime_ladder,   s in [0,1]
Run Dynamo dynamics (current-driven A from zero) at each s. Measure:
  r-statistic of the covariant spectrum; max|Im(eig)| (reality);
  self-generated flux; and at s=1, the alkuluku prime-trace sigmas.

=== REGISTERED PREDICTIONS (before any run) ===
J1  REALITY HOLDS BY CONSTRUCTION (check, not discovery). The covariant
    operator is built Hermitian at every s, so max|Im| < 1e-9 throughout.
    If it FAILS the construction is broken. (The interesting reality
    question -- does a NON-Hermitian realization of the prime amplitudes
    leave the line -- is a DIFFERENT experiment; here we hold Hermiticity
    and ask whether GUE survives, which is the cleanly separable half.)
J2  DOES GUE SURVIVE THE LADDER? Registered outcomes, all clean:
    (a) r stays >= GOE + 0.5(GUE-GOE) at s=1  -> prime weights preserve
        the class; the amplitude taper is GUE-compatible.
    (b) r falls below that midpoint monotonically in s -> prime weights
        DEGRADE the class; arithmetic disorder pulls it back toward GOE/
        Poisson. Report the s where it crosses.
    No KILL -- both (a) and (b) are results. This measures a fact.
J3  THE ARBITER (alkuluku on the s=1 spectrum). Registered NULL: none of
    the six log-prime frequencies exceeds 4 sigma; prime-ladder MAGNITUDES
    do not manufacture prime PEAKS, because peaks need the loop PHASES
    (log p transit times), not just amplitudes. If the null holds, the
    honest headline is "prime magnitudes shape the weights but the spectrum
    still does not know the primes -- the arbiter reads generic."
    MIRACLE CLAUSE: if 4+ log-prime freqs exceed 4 sigma at s=1 and not in
    a flat-weight control, STOP, re-run 3 seeds, and if it survives, that
    is the phone call. (Registered so the bar cannot move afterward.)

Interpretation locked in advance:
  Most likely: J1 pass, J2 = degrade (b), J3 null. Headline: "prime-ladder
  magnitudes are GUE-compatible-to-degrading and carry NO arithmetic
  fingerprint; the amplitude wall, in its Hermitian half, is about class
  robustness, not about primes -- and the primes remain absent." That is
  the near side of the cliff, honestly mapped.

Do not hype. Do not lie. Just show.
"""
import numpy as np, json
import networkx as nx
from numpy.linalg import eigvalsh, eigvals

TAU, DT = 2.737, 0.02

def make(N=400, k=8, p=0.15, seed=1):
    G = nx.watts_strogatz_graph(N, k, p, seed=seed)
    A = nx.to_numpy_array(G)
    tris = [tuple(c) for c in nx.enumerate_all_cliques(G) if len(c) == 3]
    return G, A, tris

def primes_upto_count(n):
    """First n primes."""
    out, c = [], 2
    while len(out) < n:
        if all(c % p for p in out if p * p <= c):
            out.append(c)
        c += 1
    return np.array(out)

def prime_ladder_weights(A):
    """Assign each edge a target magnitude log p / sqrt(p), p indexed by a
    fixed ladder ordering of the edges."""
    N = A.shape[0]
    iu = np.array(np.triu_indices(N, 1))
    mask = A[iu[0], iu[1]] != 0
    ei, ej = iu[0][mask], iu[1][mask]
    ne = len(ei)
    ps = primes_upto_count(ne)
    w = np.log(ps) / np.sqrt(ps)                 # log p / sqrt(p), decreasing
    # normalize to unit mean so s-interpolation is magnitude-fair vs flat=1
    w = w / w.mean()
    W = np.zeros((N, N))
    W[ei, ej] = w; W[ej, ei] = w
    return W

def flux_tri(Aconn, tris):
    return np.array([Aconn[i, j] + Aconn[j, k] + Aconn[k, i] for i, j, k in tris])

def curlcurl(Aconn, tris, N):
    F = flux_tri(Aconn, tris)
    out = np.zeros((N, N))
    for t, (i, j, k) in enumerate(tris):
        f = F[t]
        out[i, j] += f; out[j, i] -= f
        out[j, k] += f; out[k, j] -= f
        out[k, i] += f; out[i, k] -= f
    return out

def r_stat(ev):
    ev = np.sort(ev.real); n = len(ev); ev = ev[int(.05*n):int(.95*n)]
    s = np.diff(ev); s = s[s > 1e-12]
    return float((np.minimum(s[1:], s[:-1])/np.maximum(s[1:], s[:-1])).mean())

def calibrate(N=400, n=8, seed0=555):
    rg, ru = [], []
    for m in range(n):
        rng = np.random.default_rng(seed0+m)
        a = rng.standard_normal((N,N)); H=(a+a.T)/np.sqrt(2); rg.append(r_stat(eigvalsh(H)))
        b = rng.standard_normal((N,N))+1j*rng.standard_normal((N,N)); ru.append(r_stat(eigvalsh((b+b.conj().T)/2)))
    return float(np.mean(rg)), float(np.mean(ru))

def prime_trace_sigmas(ev):
    ev = np.sort(ev.real); n = len(ev); lo, hi = int(.05*n), int(.95*n)
    lam, idx = ev[lo:hi], np.arange(lo, hi, dtype=float)
    coef = np.polyfit(lam, idx, 9)
    grid = np.linspace(lam[0], lam[-1], 4096)
    dN = np.interp(grid, lam, idx) - np.polyval(coef, grid)
    dN -= dN.mean(); dN *= np.hanning(len(dN))
    F = np.abs(np.fft.rfft(dN))
    om = 2*np.pi*np.fft.rfftfreq(len(grid), grid[1]-grid[0])
    sig = {}
    for name, w in {'log2':np.log(2),'log3':np.log(3),'log5':np.log(5),
                    'log7':np.log(7),'log11':np.log(11),'log13':np.log(13)}.items():
        i = np.argmin(np.abs(om - w)); band = F[max(1,i-40):i+40]
        sig[name] = round(float((F[i]-np.median(band))/max(band.std(),1e-12)), 2)
    return sig

def run_s(s, Wflat, Wprime, G, Ag_graph, tris, N, seed=1,
          K=0.1, mu=0.05, equil=1000, drive=2500):
    C = (1-s)*Wflat + s*Wprime            # interpolated coupling magnitudes
    rng = np.random.default_rng(200+seed)
    th = rng.uniform(0, 2*np.pi, N); tho = th.copy()
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
    Aconn = np.zeros((N,N))
    edge = (Ag_graph != 0)
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
    im = float(np.abs(eigvals(M).imag).max())
    ev = eigvalsh(M)
    phi = float(np.abs(flux_tri(Aconn, tris)).sum())
    return dict(s=round(s,2), r=round(r_stat(ev),4), max_im=im,
                flux=round(phi,1), ev=ev)

if __name__ == '__main__':
    GOE, GUE = calibrate()
    mid = GOE + 0.5*(GUE-GOE)
    print(f"calibration GOE={GOE:.4f} GUE={GUE:.4f} midpoint={mid:.4f}")
    G, Ag, tris = make(N=400, seed=1)
    N = 400
    Wflat = (Ag != 0).astype(float)
    Wprime = prime_ladder_weights(Ag)
    print(f"prime-ladder weights: mean={Wprime[Ag!=0].mean():.3f} "
          f"range [{Wprime[Ag!=0].min():.3f}, {Wprime[Ag!=0].max():.3f}]")

    rows = []
    for s in [0.0, 0.25, 0.5, 0.75, 1.0]:
        r = run_s(s, Wflat, Wprime, G, Ag, tris, N)
        rows.append(r)
        print(f"s={s:.2f}: r={r['r']:.4f}  max|Im|={r['max_im']:.1e}  flux={r['flux']}")

    # arbiter: alkuluku prime-trace on s=1 vs flat control
    sig_prime = prime_trace_sigmas(rows[-1]['ev'])
    sig_flat = prime_trace_sigmas(rows[0]['ev'])
    n_peaks_prime = sum(v > 4 for v in sig_prime.values())
    n_peaks_flat = sum(v > 4 for v in sig_flat.values())
    print(f"alkuluku s=1 prime sigmas: {sig_prime}")
    print(f"alkuluku s=0 (flat) sigmas: {sig_flat}")

    r_end = rows[-1]['r']; r_start = rows[0]['r']
    verdict = dict(
        GOE=round(GOE,4), GUE=round(GUE,4), midpoint=round(mid,4),
        J1_reality_holds=bool(max(x['max_im'] for x in rows) < 1e-9),
        J2_outcome=("preserves_GUE" if r_end >= mid else "degrades_toward_GOE"),
        J2_r_start=r_start, J2_r_end=r_end,
        J2_r_trajectory=[x['r'] for x in rows],
        J3_prime_peaks_s1=n_peaks_prime,
        J3_prime_peaks_flat=n_peaks_flat,
        J3_arbiter_null=bool(n_peaks_prime < 4),
        MIRACLE=bool(n_peaks_prime >= 4 and n_peaks_flat < 4),
        HEADLINE=(
            "PHONE CALL: prime-ladder weights produced log-prime spectral "
            "peaks -- re-run 3 seeds immediately" if (n_peaks_prime>=4 and n_peaks_flat<4)
            else "Prime-ladder MAGNITUDES are GUE-{} and carry NO arithmetic "
                 "fingerprint (arbiter null): the Hermitian half of the amplitude "
                 "wall is about class robustness, not primes. The primes remain "
                 "absent; the cliff is unscaled.".format(
                     "compatible" if r_end>=mid else "degrading")))
    print(json.dumps(verdict, indent=2))
    save = dict(calibration=dict(GOE=GOE,GUE=GUE),
                rows=[{k:v for k,v in r.items() if k!='ev'} for r in rows],
                alkuluku=dict(s1=sig_prime, flat=sig_flat), verdict=verdict)
    json.dump(save, open('../results/jyrkanne_results.json','w'), indent=1)
