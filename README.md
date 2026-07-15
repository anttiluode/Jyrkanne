# Jyrkänne — the amplitude wall, and a whisper

*jyrkänne (Finnish): cliff, precipice. The last experiment firable from the near side. It does not scale the cliff. It does hear something faint coming off it.*

**Force the Clockfield's couplings onto the prime ladder `log p / √p` and let the dynamo run. The self-generated GUE survives contact with arithmetic — spectrum real to 10⁻¹⁴, class preserved. And the arbiter (alkuluku) hears a whisper it did not hear from flat weights: a reproducible, four-seed lift at log 5 and log 7, ~3.5–3.8σ above a control that sits at zero. Not the 4-of-6 phone-call bar. Not primes-realized. But prime *magnitudes alone*, with no loop-phase wiring, carry a partial arithmetic fingerprint — and that is more than the null this was registered to expect.**

*PerceptionLab / Antti Luode with Claude (Opus 4.8). Helsinki, July 2026.*
*The frontier where wall one (symmetry class) meets wall two (amplitudes). Downstream of `Dynamo` (the self-generating gauge field) and `alkuluku` (the arbiter). Firmly on the near side of the RH cliff — see the guard, which matters more here than anywhere.*

> Do not hype. Do not lie. Just show.

---

## What was actually new here (and what Gemini got wrong)

Gemini proposed "close the loop and let the phase current drive the connection." That is **Dynamo**, already done — the current already drives A, the glass already self-generates flux and self-lifts to GUE. Its points 1 and 2 described a finished result as future work. The one genuinely new thing it pointed at — buried under an over-claim — was point 3: **push the couplings to the prime ladder and see if GUE survives arithmetic.** This repo does exactly and only that.

And it does it with the verdict conditions written against Gemini's seduction. Gemini's closing line — *"if the spectrum stays real and GUE under prime weights, you have a physical mechanism for the Riemann operator"* — is **false**, and it is the precise trap this ecosystem exists to refuse. GUE + real + prime-*magnitudes* is necessary, not sufficient. So the arbiter of "knows the primes" is not the r-statistic; it is the **prime-trace test** (`alkuluku`), and it ran on the final spectrum with a miracle clause fixed in advance.

## The experiment

BirthOfClockfield graph, N=400. Interpolate the coupling magnitudes from flat to the prime ladder:

```
C_ij(s) = (1−s)·flat  +  s·(log p / √p),   s ∈ [0,1]
```

Run the full **Dynamo** dynamics at each s — phase current `sin(θ_i − θ_j − A_ij)` drives the connection A from zero, Maxwell cost pushes back, glass boils and settles. Measure the covariant spectrum's r-statistic, its reality `max|Im|`, the self-generated flux, and — at s=1 vs a flat control — the alkuluku log-prime FFT peaks.

## Results

Self-calibrated: GOE 0.509, GUE 0.606, midpoint 0.558.

| s (flat → prime) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 |
|---|--:|--:|--:|--:|--:|
| r-statistic | 0.605 | 0.592 | 0.599 | 0.569 | **0.600** |
| max\|Im(eig)\| | \multicolumn{5}{c}{< 7×10⁻¹⁵ throughout} |
| self-gen flux | 1097 | 1037 | 998 | 950 | 842 |

- **[V] J1 — reality holds.** max|Im| < 7×10⁻¹⁵ at every s. The covariant operator is Hermitian by construction, so pushing to prime weights does not, in this Hermitian realization, take the spectrum off the line. (The *non*-Hermitian realization of the amplitudes — where the tension could genuinely break reality — is a separate experiment, flagged below.)
- **[V] J2 — GUE survives the ladder.** r ends at 0.600, above the midpoint, GUE-class preserved across the whole interpolation (a shallow dip at s=0.75, back up at s=1). Prime-ladder magnitudes are GUE-compatible. The dynamo keeps self-generating flux throughout (842 at s=1, from zero).
- **[~→!] J3 — the arbiter hears a whisper.** Registered expectation was a clean null (magnitudes can't make peaks; peaks need loop phases = log-p transit times). The null *technically holds* — 2 of 6 frequencies over 4σ at s=1, below the 4-of-6 miracle bar. But the pooled four-seed picture is not null:

| frequency | prime weights (s=1) | flat control (s=0) | lift |
|---|--:|--:|--:|
| log 2 | +0.38σ | −0.94σ | +1.32 |
| log 3 | +1.99σ | −0.79σ | +2.78 |
| **log 5** | **+3.90σ** | +0.11σ | **+3.79** |
| **log 7** | **+3.61σ** | +0.11σ | **+3.50** |

  The lift **grows through the low primes** and lands at ~3.5–3.8σ at log 5 and log 7, over a control pinned at zero — reproducible across all four seeds. The prime *magnitudes alone*, carrying no transit-time phase structure, produce a partial arithmetic fingerprint. It is a whisper, not a voice: it does not clear the phone-call bar, log 11/log 13 are unreliable (those primes barely fit the edge-ladder), and the individual peaks wander seed to seed. But it is real, it beats its control, and it was not expected.

## What this means — stated at exactly its strength

**Licensed:** In this medium, the symmetry-class wall and the reality wall are *both* passable under prime-shaped amplitudes — GUE and a real spectrum coexist with `log p / √p` magnitudes, self-generated flux and all. That is the Hermitian half of Berry–Keating's amplitude–hermiticity tension, and in this half the tension does **not** bite: you can have the arrow, the real spectrum, and the prime magnitudes at once. And — the genuinely new thing — prime magnitudes are not spectrally inert: they push power toward the log-prime frequencies, faintly but reproducibly.

**Forbidden (and this is where Gemini's framing had to be refused):** this is **not** a Riemann operator, not a physical mechanism for one, and the whisper is **not** the primes realized. Three hard reasons: (1) The arbiter's miracle clause did not fire — 2 of 6, not 4 of 6. A whisper at log 5/log 7 is a hint, not the explicit formula's comb. (2) Prime *magnitudes* `log p/√p` are not the Prime Orbit Condition. Lemma 5.2 needs those amplitudes on **primitive prime loops** with **transit times** log p (the phases, which is where the peaks actually come from) and **all composite loops cancelling** — an infinite tower of constraints, of which this taper realizes essentially none. (3) The whisper most likely comes from the *magnitude distribution* of the ladder (heavy-tailed, `log p/√p` has structure) beating against the spectral density, not from arithmetic closure. That is a Kac–Rice-flavored effect, not zeta — and distinguishing the two is the honest next question, not a phone call.

The cliff is unscaled. What moved is that we now know the near-side walls (class, reality) are both passable together under prime weights, and that prime magnitudes leave a faint trace the arbiter can detect — which sharpens exactly where the real work is: the **phases** (transit times) and the **cancellation** (composite loops), neither touched here.

---

## Ledger

**[V] Verified.** Under prime-ladder coupling magnitudes, the Dynamo spectrum stays real to 7×10⁻¹⁵ (J1) and GUE-class (r=0.600, J2) across the full flat→prime interpolation, with self-generated flux throughout. Four-seed reproducible lift at log 5 (+3.79) and log 7 (+3.50) over a zero control (J3) — prime magnitudes are not spectrally inert.

**[K] Killed.** Gemini's claim that GUE+real+prime-weights would constitute "a physical mechanism for the Riemann operator" — false; explicitly refused, and the arbiter installed precisely to refuse it. The registered clean-null expectation for J3 — partially overturned: it is a whisper, not silence, and the honest report says so rather than rounding to "null."

**[~] Gray / open.** The whisper does not clear the miracle bar (2/6, not 4/6) and is most plausibly a magnitude-distribution (Kac–Rice) effect rather than arithmetic — untested which. log 11/log 13 unreliable (edge-ladder truncation). Single graph family, N=400, 4 seeds for the arbiter. The interpolation is a smooth taper, not the exact ladder on ordered primitive loops. Hermitian realization only.

**[B] Bet — the two untouched walls.** The **phase** wall: wire transit times log p as actual loop delays (not just magnitudes) and re-run the arbiter — that is where real peaks, if any, must come from, and it is the Prime Orbit Condition's actual content. The **non-Hermitian** wall: realize the amplitudes without forcing Hermiticity and watch whether reality breaks — that is the true amplitude–hermiticity tension, the half this repo deliberately held fixed. And, always: whether the whisper is zeta or Kac–Rice. None of these is scaled here; all are now sharply posed.

---

## Thoughts

I went in expecting a clean null and registered for one, and the honest thing is that a clean null is not quite what came back. The spectrum did not learn the primes — the arbiter is right, and the phone did not ring. But prime magnitudes turned out not to be inert either: they lean on the log-prime frequencies, faintly, reproducibly, growing through log 5 and log 7, over a control that sits flat at zero. The correct word for that is *whisper*, and the discipline is to report a whisper as a whisper — neither burying it under "null" because it missed the bar, nor inflating it toward the miracle because it landed at the right frequencies. It is a hint that the amplitude distribution of the primes has spectral consequences even without their arithmetic closure, and the most likely mundane explanation (a Kac–Rice beat between a heavy-tailed weight distribution and the level density) is checkable and probably correct and would make the whisper *not* about zeta at all. That check is the next experiment, and I would bet on the mundane explanation — which is exactly why it has to be run rather than assumed.

The larger shape of where this leaves the program: the near-side walls are mapped and both passable. Symmetry class — passable (Nuoli). Self-generation — done (Dynamo). Ghost-freedom — clean (Aave). Reality under prime magnitudes — holds (here). What is *not* touched, and what the whisper points at by its very incompleteness, is the far side: the phases (transit times, the Prime Orbit Condition's real content) and the cancellation (Lemma 5.2's infinite tower). Those are not a night's work and they are not passable by a taper. The value of Jyrkänne is that it stands at the edge, confirms the ground behind it is solid, hears something faint on the wind, and correctly declines to call a faint sound an arrival. Gemini wanted to call it an arrival. The arbiter, and the four seeds, and the honest word "whisper," are why we don't.

And the guard, load-bearing as never before: a whisper at log 5 is not the Riemann Hypothesis, is not a Riemann operator, and is quite possibly not even about primes. `alkuluku` heard it and `alkuluku` also says it is below the bar. Both of those are the same instrument being honest. The morgue holds Gemini's "physical mechanism for the Riemann operator." The trophy case holds a real spectrum, a preserved class, and a faint, well-controlled, probably-mundane whisper at the right frequencies. Knowing the difference is the whole program.

---

## Reproduce

```bash
pip install numpy scipy networkx matplotlib
python experiments/jyrkanne.py          # interpolation sweep + arbiter (1 seed)
# arbiter across seeds: rerun run_s at s=0,1 for seeds 2-4 and pool the sigmas
```

Registered predictions J1–J3 with the miracle clause are in the docstring above the code, fixed before the run. `results/jyrkanne_results.json` and `results/jyr_arbiter_pooled.json` hold every number; `figs/jyrkanne.png` shows GUE surviving the ladder (left) and the four-seed whisper at log 5/log 7 (right).

## References

`Dynamo` — the self-generating gauge field this drives. · `alkuluku` — the arbiter and its miracle clause. · `Nuoli`/`Mittari`/`Aave` — the near-side walls (class, gauge structure, ghost). · The `HilbertPolyaReintepretation` paper, Lemma 5.2 and the Prime Orbit Condition — the far-side walls (phases, cancellation) this does not touch. · Kac–Rice — the mundane explanation the whisper must be tested against. · Berry & Keating (1999) — the amplitude–hermiticity tension, whose Hermitian half is mapped here.

---

*Gemini said: reach the prime weights and you have the Riemann operator. We reached the prime weights. The spectrum stayed real, the class held, the dynamo kept spinning — and the arbiter heard a whisper at log 5 and log 7 that the control did not. It is not the operator. It is not even certainly the primes. It is a faint, honest, reproducible sound at the edge of a cliff we did not climb, reported at exactly the volume it was heard. That volume — no louder — is the whole discipline.*
