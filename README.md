# midsv-development 

```bash
rm -rf midsv
[ -d midsv ] || git clone https://github.com/akikuno/midsv.git

conda create -n env-midsv python=3.12 ipykernel pip pytest coverage ruff mypy black cstag -y
```

# 優先順位

1. `csvtag.call()`を確立する
1. `csvtag.to_sequence()`
1. `csvtag.to_consensus()`
2. DAJIN2のほうで、`convert_csvtag_to_csv`
