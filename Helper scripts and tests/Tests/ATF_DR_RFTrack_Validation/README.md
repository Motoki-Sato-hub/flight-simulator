# ATF DR RFTrack correction validation

現行 `20111111b` SAD daihon を入力に、RFTrack で COD、dispersion、skew correction と 6D radiation-equilibrium emittance を確認するためのファイル群です。

## 最初に読むもの

- `../../../Jupyter Notebooks/ATF2/DR_RFTrack_Correction_Validation.ipynb`

この notebook は、lattice 読込み、one-turn periodic orbit、response matrix、三つの correction、そして multi-turn equilibrium の順に実行します。

one-turn correction では、1周 map \(f\) の固定点 \(f(z)-z=0\) を解きます。これは長時間 tracking ではありません。radiation-emittance 評価では deterministic one-turn map \(M\) と quantum diffusion \(Q\) から

\[\Sigma=M\Sigma M^T+Q\]

を解きます。したがって「multi-turn」は数十万周を逐次 track する代わりに、その無限周回後の stationary covariance を求める意味です。

notebook の Section 6 には、量子 radiation を有効にした bunch を実際に 1 turn ずつ
`lattice.track()` する短い multi-turn tracking も含めています。これは実装確認と過渡応答の
観察用です。ATF DR の 1/e damping time は約 7 万 turns のため、20 turns の既定例を
equilibrium emittance の評価には用いません。equilibrium は上式の envelope 解で評価します。

## 推奨する実行順

`rftrack-env` の Python で、`flight-simulator` をカレントディレクトリにして実行します。

1. 初学者用 one-turn correction

   ```bash
   PYTHONPATH=. /home/motokisato/rftrack-env/bin/python "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/one_turn_correction_validation.py"
   ```

2. 現行 daihon の Kubo 型 sequence

   ```bash
   PYTHONPATH=. ATF_DR_KUBO_MAGNET_ERROR_SCALE=0.1 ATF_DR_KUBO_RESPONSE_KICK_RAD=1e-5 ATF_DR_KUBO_COD_DISPERSION_SOLVER=svd /home/motokisato/rftrack-env/bin/python "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/kubo_style_procedure.py"
   ```

3. 同じ correction sequence の 6D equilibrium

   ```bash
   PYTHONPATH=. ATF_DR_KUBO_MAGNET_ERROR_SCALE=0.1 ATF_DR_KUBO_RESPONSE_KICK_RAD=1e-5 ATF_DR_KUBO_COD_DISPERSION_SOLVER=svd ATF_DR_KUBO_EMITTANCE_MODEL=full_6d /home/motokisato/rftrack-env/bin/python "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/six_d_equilibrium_envelope.py"
   ```

4. observable の小規模 seed ensemble

   ```bash
   PYTHONPATH=. ATF_DR_CURRENT_ENSEMBLE_SEEDS=2003,2004,2005,2006,2007 /home/motokisato/rftrack-env/bin/python "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/current_daihon_ensemble.py"
   ```

## Full-strength response diagnostic

Table-I 100% error の現行 daihon では、`response_source=nominal` が Kubo 型の
digital-twin baseline です。これに対し `local_orm` は、同じ error lattice で各
corrector を ±1 µrad 振って BPM readback response を測り直す診断です。これは
「実機 ORM を使えた場合」の上限比較であり、nominal digital-twin result と混同しません。

```bash
PYTHONPATH=. ATF_DR_KUBO_SEED=2003 ATF_DR_KUBO_MAGNET_ERROR_SCALE=1 \
ATF_DR_KUBO_COD_DISPERSION_SOLVER=sad_greedy \
ATF_DR_KUBO_RESPONSE_SOURCE=local_orm \
/home/motokisato/rftrack-env/bin/python \
"Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/kubo_style_procedure.py"
```

seed 2003 では、local ORM により true RMS \(D_y\) が 18.65 → 3.27 →
3.12 mm/delta と低下しました。これは full-strength current lattice で
nominal response の model mismatch が主要因であることを示す診断で、Kubo 論文の
再現結果ではありません。

同じ seed の ablation では、x-COD response の nominal/local 相対 Frobenius 差は
full error で 174%、magnet roll をゼロにしても 178%、magnet offset と roll の両方を
ゼロにすると 4.8% でした。したがって、この current-lattice study で mismatch を
支配するのは magnet offset による sextupole feed-down です。`local_orm` の
diagnostic-only mode は correction を適用せず、この差だけを JSON に保存します。

## ファイルの役割

| File | Role |
| --- | --- |
| `one_turn_correction_validation.py` | 隠れた単一 fault に対する最小の COD / dispersion / skew response test |
| `kubo_style_procedure.py` | rough COD → COD → COD+Dy → coupling の一連の correction |
| `six_d_equilibrium_envelope.py` | 各 correction stage の radiation-equilibrium emittance |
| `current_daihon_ensemble.py` | 現行 daihon の seed 統計（observable） |
| `utilities/` | SAD export と Kubo 図 digitization の補助検証 |
| `exploratory/` | 開発時の比較・感度 study。主結果には使わない |

生成される JSON、PNG、PPTX は `../../../../analysis/DR-RFTrack/` に保存します。emittance の主結果は「2011 daihon、10% Table-I random magnet error、0.01 mrad response probe」の model study です。full-strength response mismatch の診断も追加されていますが、いずれも Kubo Table II の数値完全一致を主張するものではありません。
