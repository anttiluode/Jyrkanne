"""
ULAM -- can the indexing MAP steer the graph, and is anything arithmetic?
=========================================================================
Kaiku proved the Jyrkanne whisper was MONOTONE ORDERING: the graph feels
the macroscopic GRADIENT of the weight ladder, not the primes. That was
with a LINEAR edge index (primes laid 2,3,5,7,... in order).

The Ulam spiral works because a QUADRATIC/spiral indexing of the integers
aligns with prime-rich polynomials (n^2+n+41, etc). The registered
question: if we index the couplings by a quadratic/spiral map instead of
linearly, does the graph feel structure the linear order could not -- and
crucially, does ANY of it survive de-arithmetization (the Kaiku SORTED
control), i.e. is it primes or is it just a different geometric ramp?

Honest prior (stated before running): almost certainly NO arithmetic. The
primes are labels on edges, not participants in the dynamics; the graph
can only feel the SHAPE of the weight field across its geometry, whatever
map produced that shape. Kaiku's SORTED control is the discriminator that
keeps us honest, ported here.

===================== THE INDEXING MAPS ==========================
Assign edge e (in a FIXED canonical edge order) a coupling magnitude
w_e = g(idx(e)) where g(m) = log(m-th prime)/sqrt(m-th prime), and idx()
is the map from edge position to ladder index:

  LINEAR    : idx(e) = e                      (Jyrkanne/Kaiku baseline)
  SPIRAL    : lay edges on a sqrt(E) x sqrt(E) grid in an Ulam spiral;
              idx follows spiral order -> the ladder is written along the
              spiral, so diagonals carry quadratic subsequences.
  QUADRATIC : idx(e) = (e^2 mod E) reindexed  (a quadratic scramble of the
              ladder position -- Ulam's polynomial structure without the 2D)

For each map, FOUR Kaiku-style arms (all sharing the histogram):
  PRIME     : real prime ladder under the map
  SORTED    : same multiset, monotone along the map's index (de-arithmetized
              but SAME geometric gradient) -> isolates ordering/gradient.
  SHUFFLE   : same multiset, random -> isolates histogram.
  DShuffle  : "diagonal shuffle" -- permute only WITHIN spiral diagonals,
              preserving the spiral's coarse structure but scrambling the
              fine arithmetic. (spiral map only)

===================== REGISTERED PREDICTIONS =====================
U1  MAP CHANGES THE GEOMETRIC SIGNATURE. The spiral/quadratic maps produce
    a DIFFERENT set of elevated FFT frequencies than linear's log5 ramp
    (because a different index map = a different weight-field shape). This
    is expected and NOT arithmetic: it just shows the graph is steerable by
    the indexing. Report which frequencies light up per map.
U2  STILL NOT PRIMES (the Kaiku test, ported). For EVERY map, PRIME ~ SORTED
    within 1.5 sigma at every lit frequency: the signal is the map's
    geometric gradient, de-arithmetized SORTED reproduces it, so it is not
    the primes. KILL of the null (= the interesting outcome): PRIME exceeds
    SORTED by >= 2 sigma at some frequency under the spiral map AND that
    excess dies under DShuffle -> the spiral's quadratic arithmetic is
    doing something the monotone gradient cannot. That would be an Ulam-
    style signal and warrants many seeds. Expected: does NOT fire.
U3  ARBITER. Run the alkuluku log-prime scan on every PRIME spectrum. The
    log-prime peaks (log2..log13) should remain sub-4-sigma for all maps
    (the spiral aligns with quadratics, not with the log-p transit-time
    comb the explicit formula needs). Report; miracle clause as always.

Interpretation locked in advance:
  Most likely: U1 yes (maps steer the geometry), U2 null (PRIME=SORTED
  everywhere -> geometric not arithmetic), U3 null (no log-p comb).
  Headline: "the indexing map steers WHICH geometric resonance the graph
  feels, but no map makes it feel the PRIMES -- Ulam's spiral reveals
  arithmetic in the plane because primes SIT on the quadratics; here the
  primes are only labels, so the spiral reveals only its own geometry."

This is designed so YOUR heavy run = same code, bigger knobs: raise N,
n_seeds, and add maps. Nothing here is hardcoded to N=400.

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, sys, os
import networkx as nx
from numpy.linalg import eigvalsh

TAU, DT = 2.737, 0.02

# ------------------------------------------------------------ graph + primes
def make(N=400, k=8, p=0.15, seed=1):
    G = nx.watts_strogatz_graph(N, k, p, seed=seed)
    A = nx.to_numpy_array(G)
    tris = [tuple(c) for c in nx.enumerate_all_cliques(G) if len(c) == 3]
    return G, A, tris

def primes_upto_count(n):
    out, c = [], 2
    while len(out) < n:
        if all(c % p for p in out if p*p <= c): out.append(c)
        c += 1
    return np.array(out)

def prime_ladder(ne):
    ps = primes_upto_count(ne); w = np.log(ps)/np.sqrt(ps); return w/w.mean()

def edge_list(A):
    N = A.shape[0]; iu = np.array(np.triu_indices(N, 1)); m = A[iu[0], iu[1]] != 0
    return iu[0][m], iu[1][m]

# ------------------------------------------------------------ indexing maps
def spiral_order(E):
    """Ulam spiral visiting order on a square grid covering E cells.
    Returns a permutation: spiral_pos -> linear grid index, truncated to E."""
    L = int(np.ceil(np.sqrt(E)))
    grid = -np.ones((L, L), int)
    x = y = L // 2
    grid[y, x] = 0
    steps, num, d = 1, 1, 0
    dirs = [(1, 0), (0, -1), (-1, 0), (0, 1)]   # R, U, L, D (spiral out)
    while num < L*L:
        for _ in range(2):
            dx, dy = dirs[d % 4]
            for _ in range(steps):
                x += dx; y += dy
                if 0 <= x < L and 0 <= y < L and num < L*L:
                    grid[y, x] = num; num += 1
            d += 1
        steps += 1
    order = np.argsort(grid.ravel())            # cell -> spiral index
    # map each of the first E linear cells to its spiral rank
    spiral_rank = np.empty(L*L, int)
    spiral_rank[grid.ravel()[grid.ravel() >= 0]] = np.arange((grid >= 0).sum())
    r = spiral_rank[:E]
    # remap to a dense 0..E-1 ranking so indices are always in-range
    return np.argsort(np.argsort(r)), L

def build_index(kind, E):
    if kind == 'linear':
        return np.arange(E), None
    if kind == 'spiral':
        rank, L = spiral_order(E)
        return rank, L
    if kind == 'quadratic':
        e = np.arange(E)
        q = (e*e) % max(E, 1)
        return np.argsort(np.argsort(q)), None    # rank of quadratic scramble
    raise ValueError(kind)

def diagonal_labels(E, L):
    """For the spiral grid, label each edge by its diagonal (i+j) for DShuffle."""
    rank, _ = spiral_order(E)
    # invert: spiral index -> (row,col)
    coords = np.zeros((E, 2), int)
    # rebuild grid to get coords
    grid = -np.ones((L, L), int); x = y = L//2; grid[y, x] = 0
    steps, num, d = 1, 1, 0; dirs = [(1,0),(0,-1),(-1,0),(0,1)]
    pos = {0: (y, x)}
    while num < L*L:
        for _ in range(2):
            dx, dy = dirs[d % 4]
            for _ in range(steps):
                x += dx; y += dy
                if 0 <= x < L and 0 <= y < L and num < L*L:
                    grid[y, x] = num; pos[num] = (y, x); num += 1
            d += 1
        steps += 1
    diag_by_spidx = np.array([pos[i][0] + pos[i][1] for i in range(E)])
    # rank[edge_pos] = spiral index of that edge; diag per edge position:
    return diag_by_spidx[rank]

# ------------------------------------------------------------ weight arms
def build_weights(map_kind, arm, A, ei, ej, rng):
    ne = len(ei)
    idx, L = build_index(map_kind, ne)
    base = prime_ladder(ne)                        # ladder magnitudes by index
    wmap = base[np.argsort(np.argsort(idx))]       # place ladder along the map
    if arm == 'prime':
        w = wmap
    elif arm == 'sorted':
        # same multiset, monotone along the SAME map index (de-arithmetized,
        # preserves the geometric gradient the map induces)
        order = np.argsort(idx)
        w = np.empty(ne); w[order] = np.sort(base)[::-1]
    elif arm == 'shuffle':
        w = rng.permutation(base)
    elif arm == 'dshuffle' and map_kind == 'spiral':
        diag = diagonal_labels(ne, L)
        w = wmap.copy()
        for dlab in np.unique(diag):
            mask = diag == dlab
            w[mask] = rng.permutation(w[mask])     # scramble within diagonal
    else:
        w = wmap
    W = np.zeros((A.shape[0], A.shape[0]))
    W[ei, ej] = w; W[ej, ei] = w
    return W

# ------------------------------------------------------------ dynamics (Dynamo)
def flux_tri(Ac, tris): return np.array([Ac[i,j]+Ac[j,k]+Ac[k,i] for i,j,k in tris])
def curlcurl(Ac, tris, N):
    F = flux_tri(Ac, tris); out = np.zeros((N, N))
    for t, (i, j, k) in enumerate(tris):
        f = F[t]; out[i,j]+=f; out[j,i]-=f; out[j,k]+=f; out[k,j]-=f; out[k,i]+=f; out[i,k]-=f
    return out

def prime_trace_sigmas(ev, freqs=None):
    if freqs is None:
        freqs = {'log2':np.log(2),'log3':np.log(3),'log5':np.log(5),
                 'log7':np.log(7),'log11':np.log(11),'log13':np.log(13)}
    ev = np.sort(ev.real); n = len(ev); lo, hi = int(.05*n), int(.95*n)
    lam, idx = ev[lo:hi], np.arange(lo, hi, dtype=float)
    coef = np.polyfit(lam, idx, 9)
    grid = np.linspace(lam[0], lam[-1], 4096)
    dN = np.interp(grid, lam, idx) - np.polyval(coef, grid)
    dN -= dN.mean(); dN *= np.hanning(len(dN)); F = np.abs(np.fft.rfft(dN))
    om = 2*np.pi*np.fft.rfftfreq(len(grid), grid[1]-grid[0]); sig = {}
    for name, w in freqs.items():
        i = np.argmin(np.abs(om - w)); band = F[max(1, i-40):i+40]
        sig[name] = round(float((F[i]-np.median(band))/max(band.std(), 1e-12)), 2)
    return sig, F, om

def dominant_freqs(F, om, k=4, om_max=2.5):
    """Report the top-k FFT frequencies below om_max (map's geometric signature)."""
    m = (om > 0.2) & (om < om_max)
    idx = np.argsort(F[m])[::-1][:k]
    return [(round(float(om[m][i]), 3), round(float(F[m][i]), 1)) for i in idx]

def run(map_kind, arm, seed, N=400, K=0.1, mu=0.05, equil=1000, drive=2500):
    G, Ag, tris = make(N=N, seed=seed); ei, ej = edge_list(Ag)
    rng = np.random.default_rng(7000 + seed)
    C = build_weights(map_kind, arm, Ag, ei, ej, rng)
    r2 = np.random.default_rng(200 + seed); th = r2.uniform(0, 2*np.pi, N); tho = th.copy()
    def step(th, tho, Ac):
        d = th[:,None]-th[None,:]-Ac; dd = np.arctan2(np.sin(d), np.cos(d))
        beta = np.sum(C*dd**2, axis=1); gamma = 1/(1+TAU*beta+1e-12)**2
        force = np.sum(C*np.sin(-dd), axis=1); v = (1-0.05*DT)*(th-tho)
        return np.mod(th+v+DT**2*gamma**2*force, 2*np.pi), th.copy()
    Az = np.zeros((N, N))
    for _ in range(equil): th, tho = step(th, tho, Az)
    Ac = np.zeros((N, N)); edge = (Ag != 0)
    for _ in range(drive):
        th, tho = step(th, tho, Ac); d = th[:,None]-th[None,:]-Ac
        Jc = C*np.sin(d); Jc = 0.5*(Jc-Jc.T); cc = curlcurl(Ac, tris, N)
        Ac = (Ac + DT*(Jc - K*cc - mu*Ac))*edge; Ac = 0.5*(Ac-Ac.T)*edge
    d = th[:,None]-th[None,:]-Ac; M = -C*np.cos(d)*np.exp(1j*Ac); M = 0.5*(M+M.conj().T)
    np.fill_diagonal(M, np.sum(C*np.cos(d), axis=1))
    sig, F, om = prime_trace_sigmas(eigvalsh(M))
    return sig, dominant_freqs(F, om)

# ------------------------------------------------------------ driver
def sweep(maps, arms, seeds, N=400):
    out = {}
    for mk in maps:
        out[mk] = {}
        arms_here = arms + (['dshuffle'] if mk == 'spiral' else [])
        for arm in arms_here:
            sigs, doms = [], []
            for s in seeds:
                sig, dom = run(mk, arm, s, N=N)
                sigs.append(sig); doms.append(dom)
            keys = sigs[0].keys()
            out[mk][arm] = {f: dict(mean=round(float(np.mean([x[f] for x in sigs])), 2),
                                    std=round(float(np.std([x[f] for x in sigs])), 2))
                            for f in keys}
            out[mk][arm]['_dominant'] = doms[0]
            print(f"[{mk}/{arm}] " + " ".join(
                f"{f}={out[mk][arm][f]['mean']:+.1f}" for f in ['log3','log5','log7']))
    return out

if __name__ == '__main__':
    # LIGHT default (container). Heavy run: raise N and seeds on real hardware.
    N = int(os.environ.get('ULAM_N', 300))
    nseed = int(os.environ.get('ULAM_SEEDS', 3))
    maps = ['linear', 'spiral', 'quadratic']
    arms = ['prime', 'sorted', 'shuffle']
    seeds = list(range(1, nseed + 1))
    print(f"N={N}, seeds={seeds}, maps={maps}")
    out = sweep(maps, arms, seeds, N=N)

    # verdict: for each map, is PRIME ~ SORTED (geometric) or PRIME >> SORTED (arithmetic)?
    verdict = {}
    for mk in maps:
        rows = out[mk]
        checks = {}
        for f in ['log5', 'log7']:
            p = rows['prime'][f]['mean']; so = rows['sorted'][f]['mean']; sh = rows['shuffle'][f]['mean']
            checks[f] = dict(prime=p, sorted=so, shuffle=sh,
                             prime_vs_sorted=round(p - so, 2),
                             prime_vs_shuffle=round(p - sh, 2))
        arithmetic = any(abs(checks[f]['prime_vs_sorted']) >= 2.0 for f in checks)
        if mk == 'spiral' and 'dshuffle' in rows:
            for f in ['log5', 'log7']:
                checks[f]['dshuffle'] = rows['dshuffle'][f]['mean']
        verdict[mk] = dict(checks=checks, arithmetic_signal=bool(arithmetic),
                           geometric_only=bool(not arithmetic))
    headline = ("Some map shows PRIME >> SORTED (arithmetic beyond gradient) -- "
                "investigate with many seeds"
                if any(verdict[m]['arithmetic_signal'] for m in maps) else
                "Every map: PRIME ~ SORTED. The indexing steers WHICH geometric "
                "resonance the graph feels, but no map makes it feel the PRIMES. "
                "Ulam's spiral reveals arithmetic because primes SIT on the plane's "
                "quadratics; here primes are only labels, so the spiral reveals only "
                "its own geometry. Cliff unscaled.")
    verdict['HEADLINE'] = headline
    print("\n" + json.dumps(verdict, indent=2, default=str))
    os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'results'), exist_ok=True)
    json.dump(dict(sweep=out, verdict=verdict),
              open(os.path.join(os.path.dirname(__file__), '..', 'results',
                                'ulam_results.json'), 'w'), indent=1, default=str)
