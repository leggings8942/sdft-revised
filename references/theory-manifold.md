# SDFT v0.2 多様体理論の前提

## 1. 背景理論：リーマン多様体の離散近似

SDFT v0.2 は「軌道が描く一般のリーマン多様体」を背景理論とする。
連続版のリーマン計量 $g_{ij}$ は有限データからの推定がノイジーなため、
$k$-NN グラフ上の組み合わせ的手法で幾何量を定義する。

### 1.1 理論的位置づけ

```
v0.1 → 位相幾何学（連結性・ベッチ数・永続ホモロジー）
v0.2 → 多様体（Takens 埋め込み・離散曲率・スペクトルギャップ）← 本版
v0.3 → 情報幾何学（Fisher 計量・α-接続・自然勾配）
```

v0.2 は確率分布を経由せず、時系列データを直接多様体に持ち上げて解析する。
情報幾何学（Fisher 行列 $G$, 温度 $T$, 自由エネルギー $F = U - TS$）は
v0.3 まで使用しない。

### 1.2 v0.0.6 との対応関係

| v0.0.6（情報幾何学版） | v0.2（多様体版） | 変更理由 |
|---|---|---|
| ヒストグラム → 経験分布 | **Takens 遅延座標埋め込み** | 分布仮定不要・力学系構造を保存 |
| Fisher 情報行列 $G$ | **k-NN グラフ + Ollivier-Ricci** | ノイズ耐性・非パラメトリック |
| Bhattacharyya 距離 | **測地距離（Dijkstra on k-NN）** | 多様体上の真の距離 |
| $T$（情報的温度） | **廃止** | $k_B$ アナロジー不要 |
| $F = U - TS$ | **$F_{\text{geom}} = \lambda_2$（スペクトルギャップ）** | 系の結合強度の直接指標 |
| $H$（Hurst 指数） | **$\lambda_1$（最大リアプノフ指数）** | カオス/安定の直接測定 |
| $D$（Higuchi フラクタル次元） | **$D$（Grassberger-Procaccia 相関次元）** | 多様体上の真の次元推定 |
| $S = k_B \ln 2 \cdot H_{bit}$ | **$S_g$（軌道体積エントロピー）** | 分布仮定不要 |

## 2. Takens 遅延座標埋め込み

### 2.1 Takens の埋め込み定理（1981）

$d \geq 2n + 1$（$n$ はアトラクタの次元）のとき、
遅延座標マップ $\Phi_\tau: x(t) \mapsto (x(t), x(t+\tau), \ldots, x(t+(d-1)\tau))$
は generic な条件下で埋め込み（injective immersion）となる。

### 2.2 パラメータ推定

- **遅延 $\tau$**: 相互情報量の最初の極小値（Fraser & Swinney 1986）
- **次元 $d$**: False Nearest Neighbors（Kennel, Brown & Abarbanel 1992）

## 3. 離散リーマン幾何

### 3.1 Ollivier-Ricci 曲率

$$\kappa(x, y) = 1 - \frac{W_1(m_x, m_y)}{d(x, y)}$$

$m_x$: 点 $x$ の $k$-近傍上の一様分布。$W_1$: Wasserstein-1 距離。

- **References**: Ollivier (2009), Lin-Lu-Yau (2011)
- **計算**: 最適割当問題（Hungarian algorithm）で $W_1$ を計算

### 3.2 スペクトルギャップ $\lambda_2$

正規化グラフ Laplacian $L = I - D^{-1/2}WD^{-1/2}$ の
第 2 最小固有値（Fiedler 値）。

- **Cheeger 不等式**: $\lambda_2/2 \leq h(G) \leq \sqrt{2\lambda_2}$
- $h(G)$: Cheeger 定数（等周定数の離散版）
- **References**: Chung (1997), Belkin & Niyogi (2003)

## 4. グランドポテンシャル Φ の正当化

### 4.1 方向 A：分解構造の保存

$$\Phi = F_{\text{geom}} - \mu \cdot N$$

- $F_{\text{geom}} = \lambda_2$: 系の構造エネルギー
- $\mu = T_2$（測地距離）: 1人の移動コスト
- $N$: 移動人数（Observed）

### 4.2 「放置時に減少」の保証

| 項 | 放置時の変化 | Φ への影響 |
|---|---|---|
| $\lambda_2$ 低下 | 構造が緩む → 結合が弱まる | Φ 減少 ✓ |
| $\mu N$ 増大 | レジーム乖離 → 距離増大 | Φ 減少 ✓ |

### 4.3 v0.3 での完成

$F_{\text{geom}} - \mu N$ の「引き算」は、v0.3 で情報幾何を導入すると
Legendre 変換として正当化される。v0.2 時点では empirical な正当化のみ。

## 5. T₃, T₄ の簡易版と精密化予定

### T₃（平行移動の不整合）

- **v0.2（簡易版）**: 局所 PCA フレーム間の回転角
- **v0.2.1（精密版）**: Singer-Wu vector diffusion maps

### T₄（ホロノミー）

- **v0.2（簡易版）**: 曲率差 $|\kappa_A - \kappa_B| / \max(|\kappa_A|, |\kappa_B|)$
- **v0.2.1（精密版）**: 離散ホロノミー（閉ループに沿った平行移動の蓄積差）

## 6. 参考文献

- Takens, F. (1981). "Detecting strange attractors in turbulence."
- Fraser, A.M., Swinney, H.L. (1986). "Independent coordinates for strange attractors from mutual information."
- Kennel, M.B., Brown, R., Abarbanel, H.D.I. (1992). "Determining embedding dimension for phase-space reconstruction..."
- Grassberger, P., Procaccia, I. (1983). "Measuring the strangeness of strange attractors."
- Rosenstein, M.T., Collins, J.J., De Luca, C.J. (1993). "A practical method for calculating largest Lyapunov exponents..."
- Ollivier, Y. (2009). "Ricci curvature of Markov chains on metric spaces."
- Lin, Y., Lu, L., Yau, S.-T. (2011). "Ricci curvature of graphs."
- Chung, F.R.K. (1997). "Spectral Graph Theory."
- Belkin, M., Niyogi, P. (2003). "Laplacian eigenmaps for dimensionality reduction and data representation."
