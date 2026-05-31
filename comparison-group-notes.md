# Two comparison groups were used.

Both are null hypotheses used to test whether the "named planet ↔ news articles" correlations are real or could arise by chance — but they answer **different questions** about what "by chance" means.

## Compound null
**A synthetic null built by shuffling the real data**, combining two permutations:

1. **Cross-planet swap** — randomly reassign which article series belongs to which planet's orbital trajectory. (Mars's articles get attached to Pallas's orbit, Ceres's articles to Juno's orbit, etc.)
2. **Per-planet phase shift** — additionally rotate each planet's article series by a random angular offset, breaking any residual angular registration.

Applying both together destroys any genuine name↔orbit correlation while preserving:
- the marginal distribution of article counts per planet,
- each planet's actual orbital geometry,
- the aggregate statistical structure of the dataset.

**Question it answers:** *"Given the same statistical structure, could this signal arise by chance permutation of the labels?"*

## Unnamed-MBA null
**A real, separate comparison group:** 1,211 unnamed main-belt asteroids with comparable orbital properties to the named planets in the catalogue. Their news-article counts (essentially zero or near-zero, since they have no cultural identity) serve as the control. For each named planet's statistic, you compute the analogous statistic for the unnamed-MBA pool and compare distributions.

**Question it answers:** *"Is the signal specific to NAMED planets, or does any minor-planet orbit produce the same pattern?"*

## Key differences

| | Compound null | Unnamed-MBA null |
|---|---|---|
| What it is | synthetic (shuffled real data) | real (separate comparison dataset) |
| Preserves | marginal distributions | nothing — independent data |
| Tests against | random permutation of labels | absence of cultural-name effect |
| Failure mode it catches | spurious correlations from data structure | astrophysical / catalogue artifacts |
| Robustness to systematics | strong (same data shuffled) | strong (independent data) |
| Cost | computational (many permutations) | fixed (one pool of unnamed MBAs) |

## Why use both

They guard against different failure modes:

- A signal could **pass the compound null** (so it's not a permutation artifact) but **fail the unnamed-MBA null** if it turns out *any* minor planet shows the same statistic — meaning the "name" isn't doing real work, and the apparent correlation is some property of the catalogue or the article-counting procedure.

- Conversely, a signal could **pass the unnamed-MBA null** (named planets > unnamed) but **fail the compound null** if the structure within the named-planet set is itself well-explained by random shuffling — e.g. it's driven by one big outlier rather than a coherent pattern across planets.
