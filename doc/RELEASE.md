<!-- TEMPLATE
# v0.0.0 (yyyy-mm-dd)
## 💥 Breaking
## 📝 Documentation
## 🚀 Performance
## 🌟 New Features
## 🐛 Bug Fixes
## 🔧 Maintenance
## ⛔️ Deprecated
+ commitMessage. Issue #XX [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/xxxxx)]
-->

<!-- ############################################################# # -->

# v0.12.1 (2025-12-05)

## 🐛 Bug Fixes

- Pass keep through transform to preserve requested fields such as FLAG


-------------------------------------------------------------

# Past Releases

<!-- <details>
<summary> v0.X.X (2024-MM-DD) </summary>

</details> -->



<details>
<summary> v0.12.0 (2025-12-05) </summary>

## 💥 Breaking

+ `main.transform` now requires `path_sam: Path | str`; the old API that accepted a SAM iterator is removed.

+ Output dictionaries from `main.transform` now expose `MIDSV` (formerly `CSSPLIT`), so callers must read `MIDSV`.

+ Reference Ns are always emitted as `=N`/`=n`.

## 🔧 Maintenance

+ Documentation refreshed to outline the breaking changes and differences from v0.11.x.

</details>


<details>
<summary> v0.11.1 (2024-12-20) </summary>

## 🔧 Maintenance

+ Replace setup.py to pyproject.toml. Issue #4 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/16b861ac3851f70e02e7c6565280caff7b1f63b7)]

+ Resolve the errors of testing. Issue #5 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c85694ebeb093c2d2d9febb373d745420310edb7)]

</details>
